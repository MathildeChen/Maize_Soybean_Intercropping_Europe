# -------------------------------------------------------------------------
#
#       03-1 - Simulations for soybean allocation in Europe (main analysis)
#       Author: M. Chen, Inrae, 2024
#         
# -------------------------------------------------------------------------

# ----------------------------------------
# Packages & tools
library(tidyverse)
library(stringr)
library(lubridate)
library(CCMHr)
library(parallel) ; library(doParallel); library(foreach)

# > Function for allocation
#source("E:/POSTDOC INRAE/ANALYSES/B_OPTIMISATION/00_0_Function_for_allocation_3.R")
source("E:/POSTDOC INRAE/ANALYSES/B_OPTIMISATION/scripts_R2/00_0_Function_for_allocation_3.R")

# ----------------------------------------
# Data

# Train data for soybean yields 
# Yields in Argentina, Brazil, Canada, China, India, Italy, and USA
# 1981-2016
load("E:/POSTDOC INRAE/ANALYSES/B_OPTIMISATION/00_Data/00_tab_soybean.rda")

# Coordinates of pixels in EU
load("E:/POSTDOC INRAE/ANALYSES/B_OPTIMISATION/00_Data/00_dat_coords_EU27.rda")
dat_coords_eu27$gridcode = paste0(dat_coords_eu27$x, "_", dat_coords_eu27$y)
dim(dat_coords_eu27) # N=2699 

load("E:/POSTDOC INRAE/ANALYSES/B_OPTIMISATION/00_Data/00_dat_coords_EU42.rda")
dat_coords_EU <- dat_coords_eu42
dat_coords_eu42$gridcode = paste0(dat_coords_eu42$x, "_", dat_coords_eu42$y)
dim(dat_coords_eu42) # N=4192

# Projected yields in Europe
load("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_preds/04_EU/Ya_pred_eu_2000_2023.rda")

# > EU27
Ya_pred_eu %>%
  filter(id_eu27==1) %>% 
  distinct(x,y) %>% dim(.) # N=2699

# > EU42
Ya_pred_eu %>% 
  distinct(x,y) %>% dim(.) # N=4192    

# > total surface cropland in the EU27 estimated 
cropland_surf <- Ya_pred_eu %>% 
  filter(id_eu27==1) %>%
  distinct(x,y,cropland_area_ha) %>% 
  pull(cropland_area_ha) %>%   
  sum(.) ; (cropland_surf)*10^-6 # ~105.629 Mha

# ----------------------------------------
# Minimum observed yield in Italy and USA
# > what is the current minimum "accepted" yield? 
min_Ya <- tab_soybean %>% 
  filter(country_name %in% c("Italy", "United States of America")) %>% 
  group_by(x,y) %>% 
  summarise(mean_Ya = mean(Ya)) %>% 
  ungroup() %>%
  pull(mean_Ya) %>% 
  min(.) ; min_Ya # min_Ya = 1.018906 t/ha

# Sets of pixels with yields > minimum accepted yield (min_Ya)
pixels_with_enough_soybean_1 <- Ya_pred_eu %>% 
  filter(model=="pca.m.2", crop=="soybean") %>% 
  group_by(x, y) %>% 
  summarise(mean_Ya_pred = mean(Ya_pred)) %>%
  filter(mean_Ya_pred >= min_Ya) ; dim(pixels_with_enough_soybean_1) #  3217 (among 4192) pixels with yields >= min_Ya t/ha

# Mean average soybean grain yield in the EU27 (FAOSTATS 2018-2022)
# > https://www.fao.org/faostat/fr/#data/QCL
yields_FAO_EU27 <- c(30481, 30986, 28444, 28872, 22895) # yields in the 2018-2022 period expressed in 100 g / ha 
mean_Ya <- mean((yields_FAO_EU27*100)*(10^-6)) ; mean_Ya # mean_Ya = 2.83356 t/ha

# > pixels with yields >= mean yield in average in EU27 (FAOSTATS 2018-2022)
pixels_with_enough_soybean_2 <- Ya_pred_eu %>% 
  filter(model=="pca.m.2", crop=="soybean") %>% 
  group_by(x, y) %>% 
  summarise(mean_Ya_pred = mean(Ya_pred)) %>%
  filter(mean_Ya_pred >= mean_Ya) ; dim(pixels_with_enough_soybean_2) #  417 (among 4192) pixels with yields >= mean_Ya t/ha

# ----------------------------------------
# Prepare data use for allocation

# MAIN ANALYSIS: Pixels with yields >= min_Ya t/ha in the EU27
data_for_allocation_1 <- Ya_pred_eu %>% 
  filter(id_eu27==1, model=="pca.m.2") %>% 
  spread(key=crop, value=Ya_pred) %>% 
  left_join(pixels_with_enough_soybean_1, by=c("x", "y")) %>% 
  # > remove the pixels with < min_Ya t/ha in the EU27
  filter(is.na(mean_Ya_pred)==F) %>% 
  dplyr::select(x, y, year, 
                "area"="cropland_area_ha", 
                "Ya_pred"="soybean", 
                "Ya_pred_2"="maize")

# Selected pixels with yields >= min_Ya t/ha used for the analyses 
# > total number of pixels
data_for_allocation_1 %>% 
  distinct(x,y) %>% dim(.) # 2034 (among 2699) pixels with yields >= min_Ya t/ha in the EU27

# ----------------
# Pixels with yields >= mean_Ya t/ha in the EU27
data_for_allocation_2 <- Ya_pred_eu %>% 
  filter(model=="pca.m.2") %>% 
  spread(key=crop, value=Ya_pred) %>% 
  left_join(pixels_with_enough_soybean_2, by=c("x", "y")) %>% 
  filter(id_eu27==1, is.na(mean_Ya_pred)==F) %>% 
  dplyr::select(x, y, year, 
                "area"="cropland_area_ha", 
                "Ya_pred"="soybean", 
                "Ya_pred_2"="maize")

# Selected pixels with yields >= mean_Ya t/ha used for the analyses 
# > total number of pixels
data_for_allocation_2 %>% 
  distinct(x,y) %>% dim(.) # 317 (among 2699) pixels with yields >= min_Ya t/ha in the EU27

# ----------------------------------------
# Simulation plan regarding pLER values 

# > 

tab_pLERs_ref <- data_for_allocation_1 %>% 
  distinct(x,y) %>% 
  mutate(pLER_crop_xy=0.504716,
         pLER_crop_1_xy=0.504716,
         pLER_crop_2_xy=0.78328)




# > 75th percentile of mean nitrogen use for maize
#   in the EU in 2018
#   source of data: 
#      Ludemann, C.I., Gruere, A., Heffer, P. et al. 
#      Global data on fertilizer use by crop and by country. 
#      Sci Data 9, 501 (2022). https://doi.org/10.1038/s41597-022-01592-z
N_maize_eu_data <- readxl::read_xlsx("C:/Users/benni/Desktop/Sauvegarde 05 2025/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/EU27_N_FERTI_IFA.xlsx", col_types = c("text", "text", rep("numeric", 10)))
N_maize_eu_perc50 <- as.numeric(as.character(quantile(N_maize_eu_data$N_rate_kg_per_ha, na.rm=T, probs = 0.5))) ; N_maize_eu_perc50   # Median quantity of N applied/ha in maize in EU: 125 kg/ha 
N_maize_eu_perc75 <- as.numeric(as.character(quantile(N_maize_eu_data$N_rate_kg_per_ha, na.rm=T, probs = 0.75))) ; N_maize_eu_perc75  # 75th percentile of N applied/ha in maize in EU: 148 kg/ha 
  
# > Partial Land Equivalent Ratios (pLERs) of soybean and maize
# based on the meta-analysis of Xu et al., 2020
# variation according to nitrogen fertilisation in the EU 
num_pLER_s_sensi <- 0.591 - 0.000583 * c(N_maize_eu_perc50, N_maize_eu_perc75) ; num_pLER_s_sensi  
# in the median situation: 0.518125; in the 75th percentile situation: 0.504716
# -> for comparison: num_pLER_s = 0.56

num_pLER_m_sensi <- 0.767 + 0.00011  * c(N_maize_eu_perc50, N_maize_eu_perc75) ; num_pLER_m_sensi  
# in the median situation: 0.78075 ; in the 75th percentile situation: 0.78328
# -> for comparison: num_pLER_m = 0.79

# ----------------------------------------
# Other features of the simulation plan 

# > European Union 27 soybean supply
# Based on FAOSTATS - 2018-2022
num_eu_supply <- 36.3*10^6

# > Surface max allocated in Mha (25% of the total EU croplands)
# EUROSTATS DATA (based on CORINE land cover): ~100 Mha (see: https://ec.europa.eu/eurostat/databrowser/view/lan_lcv_ovw__custom_13006546/default/table?lang=en)
# --> 1/4 of total cropland area = 100/4 ~ 25 Mha 
num_eu_cropland25 <- 25*10^6

# > Combinaisons of pLERs, pLERm, crop frequencies and soybean target production
#   median N use (125 kg/ha): pLER_s = 0.518125 and pLER_m = 0.78075
simul_plan_restricted_perc50 <- expand.grid(pLER_s_i   = num_pLER_s_sensi[1],
                                            pLER_m_i   = num_pLER_m_sensi[1],          
                                            freqs_i    = c(0.14, 0.16, 0.20, 0.25, 0.33, 0.5),
                                            target_i   = c(num_eu_supply*0.25, num_eu_supply*0.5, num_eu_supply*0.75, num_eu_supply),
                                            max_surf_i = num_eu_cropland25) %>% 
  unite(col = "scenario", c(pLER_s_i, pLER_m_i, freqs_i, target_i, max_surf_i), sep = "_", remove = F) %>% 
  mutate(pLER_s_i = as.numeric(as.character(pLER_s_i)),
         pLER_m_i = as.numeric(as.character(pLER_m_i)))

dim(simul_plan_restricted_perc50)
#24 scenarios
unique(simul_plan_restricted_perc50$pLER_s_i) # 0.518125
unique(simul_plan_restricted_perc50$pLER_m_i) # 0.78075

#   75th percentile of N use (148 kg/ha): pLER_s = 0.504716 and pLER_m = 0.78328
simul_plan_restricted_perc75 <- expand.grid(pLER_s_i   = num_pLER_s_sensi[2],
                                            pLER_m_i   = num_pLER_m_sensi[2],          
                                            freqs_i    = c(0.14, 0.16, 0.20, 0.25, 0.33, 0.5),
                                            target_i   = c(num_eu_supply*0.25, num_eu_supply*0.5, num_eu_supply*0.75, num_eu_supply),
                                            max_surf_i = num_eu_cropland25) %>% 
  unite(col = "scenario", c(pLER_s_i, pLER_m_i, freqs_i, target_i, max_surf_i), sep = "_", remove = F) %>% 
  mutate(pLER_s_i = as.numeric(as.character(pLER_s_i)),
         pLER_m_i = as.numeric(as.character(pLER_m_i)))

dim(simul_plan_restricted_perc75)
#24 scenarios
unique(simul_plan_restricted_perc75$pLER_s_i) # 0.504716
unique(simul_plan_restricted_perc75$pLER_m_i) # 0.78328

#   TND = 0: pLER_s = 0.525 and pLER_m = 0.755
simul_plan_restricted_tnd0 <- expand.grid(pLER_s_i   = 0.525,
                                          pLER_m_i   = 0.755,          
                                          freqs_i    = c(0.14, 0.16, 0.20, 0.25, 0.33, 0.5),
                                          target_i   = c(num_eu_supply*0.25, num_eu_supply*0.5, num_eu_supply*0.75, num_eu_supply),
                                          max_surf_i = num_eu_cropland25) %>% 
  unite(col = "scenario", c(pLER_s_i, pLER_m_i, freqs_i, target_i, max_surf_i), sep = "_", remove = F) %>% 
  mutate(pLER_s_i = as.numeric(as.character(pLER_s_i)),
         pLER_m_i = as.numeric(as.character(pLER_m_i)))

dim(simul_plan_restricted_tnd0)
#24 scenarios
unique(simul_plan_restricted_tnd0$pLER_s_i) # 0.525
unique(simul_plan_restricted_tnd0$pLER_m_i) # 0.755

#############

# N fertilizer use for maize in the EU
load("E:/POSTDOC INRAE/ANALYSES/B_OPTIMISATION/00_Data/03_Nrate_maize.rda")

simul_pLERs_ferti <- data_for_allocation_1 %>% 
  distinct(x,y) %>% 
  #add Nrate for each cell
  left_join(fertilization_m_tab, by=c("x", "y")) %>% 
  rename(Nrate_m_i = Nrate_m) %>% 
  mutate(Nrate_m_i = if_else(Nrate_m_i < 0, 0, Nrate_m_i)) %>% 
  #compute pLER for each cell
  mutate(pLER_crop_xy = 0.591 - 0.000583 * Nrate_m_i,
         pLER_crop_1_xy=pLER_crop_xy,
         pLER_crop_2_xy = 0.767 + 0.00011  * Nrate_m_i) %>% 
  dplyr::select(x, y, pLER_crop_xy, pLER_crop_1_xy, pLER_crop_2_xy)


simul_plan_restricted_ferti <- expand.grid(freqs_i    = c(0.14, 0.16, 0.20, 0.25, 0.33, 0.5),
                                          target_i   = c(num_eu_supply*0.25, num_eu_supply*0.5, num_eu_supply*0.75, num_eu_supply),
                                          max_surf_i = num_eu_cropland25) %>% 
  mutate(split_id = 1:n()) %>% 
  split(.$split_id) %>% 
  map_dfr(., ~{
    
    simul_pLERs_ferti$freqs_i <- .x$freqs_i
    simul_pLERs_ferti$target_i <- .x$target_i
    simul_pLERs_ferti$max_surf_i <- .x$max_surf_i
    
    simul_pLERs_ferti$pLER_s_scenario <- .x$split_id
    simul_pLERs_ferti$pLER_m_scenario <- .x$split_id
    
    simul_pLERs_ferti
    
  }, .id = "split_id") %>%
  unite(col = "scenario", c(pLER_s_scenario, pLER_m_scenario, freqs_i, target_i, max_surf_i), sep = "_", remove = F) %>% 
  mutate(pLER_crop_1_xy = as.numeric(as.character(pLER_crop_1_xy)),
         pLER_crop_2_xy = as.numeric(as.character(pLER_crop_2_xy)))



# > Path to save the results 
path_alloc <- "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/08_allocations/allocations_with_max_surf/"

# ----------------------------------------
# ----------------------------------------
#      Allocation based on each scenario 
# -----{restricted set of simulation}-----
  

list_simul <- list(
  "N_ferti" = list(data = data_for_allocation_1, tested_scenarios = simul_plan_restricted_ferti, save_name = "sensi/sensi_5_3")
  #"N_perc50" = list(data = data_for_allocation_1, tested_scenarios = simul_plan_restricted_perc50, save_name = "sensi/sensi_5_1"),
  #"N_perc75" = list(data = data_for_allocation_1, tested_scenarios = simul_plan_restricted_perc75, save_name = "sensi/sensi_5_2"),
  #"TND_0"    = list(data = data_for_allocation_1, tested_scenarios = simul_plan_restricted_tnd0,   save_name = "sensi/sensi_6_1")
)

# > Performing 3 first sets of sensitivity analyses
n.cores <- parallel::detectCores() - 3 ; n.cores # 7 on the computer, 11 on the laptop

# >>> create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

# register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

# > all sensitivity analyses (except the one testing all pLERs values*crop frequency)
list_simul %>% 
  map(., ~{
      
      # > data for allocation
      sensi_data <- .x$data
      # > list of scenarios to be tested
      sensi_scenarios <- .x$tested_scenarios
      
      # > allocate soybean and maize according in each scenario
      sensi_allocations <- sensi_scenarios %>%
        split(.$scenario) %>% 
        map(., ~{
          
          # > Retrieve values from the selected scenario
          target_i_x   <- unique(.x$target_i)
          #pLER_s_i_x   <- unique(.x$pLER_s_i)
          #pLER_m_i_x   <- unique(.x$pLER_m_i)
          freqs_i_x    <- unique(.x$freqs_i)
          max_surf_i_x <- unique(.x$max_surf_i)
          
          # > Table with pLERs values
          tab_pLERs_i <- data_for_allocation_1 %>% 
            distinct(x,y) %>% 
            mutate(pLER_crop_xy   = .x$pLER_crop_1_xy,
                   pLER_crop_1_xy = .x$pLER_crop_1_xy,
                   pLER_crop_2_xy = .x$pLER_crop_2_xy)
          
          
          # > Allocation
          alloc_i  <- suppressMessages(allocation_intercropping_solecropping(data              = sensi_data,
                                                                             target_production = target_i_x,  
                                                                             max_surf          = max_surf_i_x,
                                                                             tab_pLER_crop     = tab_pLERs_i, 
                                                                             freq_crop         = freqs_i_x))
          
          # > Results of allocation
          res_i <- suppressMessages(results_allocation(data_allocation = alloc_i, 
                                                       tab_pLER_crops  = tab_pLERs_i, 
                                                       freq_crops_i    = freqs_i_x, 
                                                       freq_crops_m    = c(freqs_i_x, freqs_i_x)))
          
          # > Format data prior exportation
          # > Full dataset allocated
          res_i$data_res <- res_i$data_res %>% 
            mutate(id_pixel = case_when(pixel_solecropping       == 1 ~ "Soybean sole cropping",
                                        pixel_intercropping_free == 1 ~ "Maize sole cropping",
                                        pixel_landsaving         == 1 ~ "Land saving",
                                        TRUE~NA))
          # > Results 
          # Labels
          res_i$res1 <- res_i$res1 %>% 
            mutate(crop = recode(crop, "crop 1"="Soybean", "crop 2"="Maize")) %>%
            mutate(strategy = recode(strategy, "intercrop"="Intercropping", "sole crop"="Sole crop"))
          # > Output
          res_i   })
      
      # > Save 
      save(sensi_allocations, 
           file = paste0(path_alloc, .x$save_name, ".rda"))
      
      cat("=")
      
    })
  
# >>> stop cluster//
stopCluster(my.cluster)


