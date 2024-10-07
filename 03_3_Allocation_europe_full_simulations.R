# -------------------------------------------------------------------------
#
#       03-2 - Full set of simulations for soybean allocation in Europe
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
source("E:/POSTDOC INRAE/ANALYSES/B_OPTIMISATION/00_0_Function_for_allocation_3.R")

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

# Mean average soybean grain yield in the EU27 (FAOSTATS 2018-2022)
# > https://www.fao.org/faostat/fr/#data/QCL
yields_FAO_EU27 <- c(30481, 30986, 28444, 28872, 22895) # yields in the 2018-2022 period expressed in 100 g / ha 
mean_Ya <- mean((yields_FAO_EU27*100)*(10^-6)) ; mean_Ya # mean_Ya = 2.83356 t/ha

# ----------------------------------------
# Sets of pixels with productivity higher than a certain threshold

# Sensitivity analysis 1 
# > pixels with yields >= mean_Ya t/ha in the EU27
pixels_sensi1 <- Ya_pred_eu %>% 
  filter(model=="pca.m.2", crop=="soybean") %>% 
  group_by(x, y) %>% 
  summarise(mean_Ya_pred = mean(Ya_pred)) %>%
  filter(mean_Ya_pred >= mean_Ya) ; dim(pixels_sensi1) #  417 (among 4192) pixels with yields >= mean_Ya t/ha

# > select pixels with yields >= mean_Ya t/ha used for the analyses 
data_for_allocation_sensi1 <- Ya_pred_eu %>% 
  filter(id_eu27==1, model=="pca.m.2") %>% 
  spread(key=crop, value=Ya_pred) %>% 
  left_join(pixels_sensi1, by=c("x", "y")) %>% 
  filter(is.na(mean_Ya_pred)==F) %>% 
  dplyr::select(x, y, year, 
                "area"="cropland_area_ha", 
                "Ya_pred"="soybean", 
                "Ya_pred_2"="maize")

# > total number of pixels
data_for_allocation_sensi1 %>% 
  distinct(x,y) %>% dim(.) # 317 (among 2699) pixels with yields >= min_Ya t/ha in the EU27

# Sensitivity analysis 2.1 - EU42
# > pixels with yields > minimum accepted yield (min_Ya)
pixels_sensi2_1 <- Ya_pred_eu %>% 
  filter(model=="pca.m.2", crop=="soybean") %>% 
  group_by(x, y) %>% 
  summarise(mean_Ya_pred = mean(Ya_pred)) %>%
  filter(mean_Ya_pred >= min_Ya) ; dim(pixels_sensi2_1) #  3217 (among 4192) pixels with yields >= min_Ya t/ha in the EU42

# > select
data_for_allocation_sensi2_1 <- Ya_pred_eu %>% 
  filter(model=="pca.m.2") %>% 
  spread(key=crop, value=Ya_pred) %>% 
  left_join(pixels_sensi2_1, by=c("x", "y")) %>% 
  # > remove the pixels with < min_Ya t/ha in the EU27
  filter(is.na(mean_Ya_pred)==F) %>% 
  dplyr::select(x, y, year, 
                "area"="cropland_area_ha", 
                "Ya_pred"="soybean", 
                "Ya_pred_2"="maize")

# > total number of pixels
data_for_allocation_sensi2_1 %>% 
  distinct(x,y) %>% dim(.) # 3217 (among 4192) pixels with yields >= min_Ya t/ha in the EU42

# Sensitivity analysis 2.2 - EU42
# > select pixels with yields >= mean_Ya t/ha used for the analyses 
data_for_allocation_sensi2_2 <- Ya_pred_eu %>% 
  filter(model=="pca.m.2") %>% 
  spread(key=crop, value=Ya_pred) %>% 
  left_join(pixels_sensi1, by=c("x", "y")) %>% 
  # > remove the pixels with < min_Ya t/ha in the EU27
  filter(is.na(mean_Ya_pred)==F) %>% 
  dplyr::select(x, y, year, 
                "area"="cropland_area_ha", 
                "Ya_pred"="soybean", 
                "Ya_pred_2"="maize")

# > total number of pixels
data_for_allocation_sensi2_2 %>% 
  distinct(x,y) %>% dim(.) # 417 (among 4192) pixels with yields >= min_Ya t/ha in the EU42

# Sensitivity analysis 3 - predictions from other models 
pixels_sensi3 <- Ya_pred_eu %>% 
  filter(crop=="soybean") %>% 
  group_by(model, x, y) %>% 
  summarise(mean_Ya_pred = mean(Ya_pred)) %>%
  filter(mean_Ya_pred >= min_Ya) 

pixels_sensi3 %>% 
  group_by(model) %>%
  count()

#  model       n
#1 avg.m    3097
#2 avg.s    2628
#3 pca.m.2  3217
#4 pca.m.3  3243

# > avg.m
data_for_allocation_avg_m <- Ya_pred_eu %>% 
  filter(model=="avg.m") %>% 
  spread(key=crop, value=Ya_pred) %>% 
  left_join(pixels_sensi3[which(pixels_sensi3$model=="avg.m"),], by=c("x", "y")) %>% 
  filter(id_eu27==1, is.na(mean_Ya_pred)==F) %>% 
  dplyr::select(x, y, year, 
                "area"="cropland_area_ha", 
                "Ya_pred"="soybean", 
                "Ya_pred_2"="maize")

# Number of pixels with yields > 1 t/ha used for the analyses 
data_for_allocation_avg_m %>% 
  distinct(x,y) %>% dim(.) # N=1977

# > avg.s
data_for_allocation_avg_s <- Ya_pred_eu %>% 
  filter(model=="avg.s") %>% 
  spread(key=crop, value=Ya_pred) %>% 
  left_join(pixels_sensi3[which(pixels_sensi3$model=="avg.s"),], by=c("x", "y")) %>% 
  filter(id_eu27==1, is.na(mean_Ya_pred)==F) %>% 
  dplyr::select(x, y, year, 
                "area"="cropland_area_ha", 
                "Ya_pred"="soybean", 
                "Ya_pred_2"="maize")

# Number of pixels with yields > 1 t/ha used for the analyses 
data_for_allocation_avg_s %>% 
  distinct(x,y) %>% dim(.) # N=1770

# > pca.m.3
data_for_allocation_pca_m_3 <- Ya_pred_eu %>% 
  filter(model=="pca.m.3") %>% 
  spread(key=crop, value=Ya_pred) %>% 
  left_join(pixels_sensi3[which(pixels_sensi3$model=="pca.m.3"),], by=c("x", "y")) %>% 
  filter(id_eu27==1, is.na(mean_Ya_pred)==F) %>% 
  dplyr::select(x, y, year, 
                "area"="cropland_area_ha", 
                "Ya_pred"="soybean", 
                "Ya_pred_2"="maize")

# Number of pixels with yields > 1 t/ha used for the analyses 
data_for_allocation_pca_m_3 %>% 
  distinct(x,y) %>% dim(.) # N=2037

# ----------------------------------------
# Simulation plan 

# > Partial Land Equivalent Ratios (pLERs) of soybean and maize
# based on the meta-analysis of Xu et al., 2020
num_pLER_s <- 0.56
num_pLER_m <- 0.79

# > European Union 27 soybean supply
# Based on FAOSTATS - 2018-2022
num_eu_supply <- 36.3*10^6

# > Surface max allocated in Mha (25% of the total EU croplands)
# EUROSTATS DATA (based on CORINE land cover): ~100 Mha (see: https://ec.europa.eu/eurostat/databrowser/view/lan_lcv_ovw__custom_13006546/default/table?lang=en)
# --> 1/4 of total cropland area = 100/4 ~ 25 Mha 
num_eu_cropland25 <- 25*10^6

# > Combinaisons of pLERs, pLERm, crop frequencies and soybean target production
simul_plan_restricted <- expand.grid(pLER_s_i   = c(seq(0.3, 0.7, by=0.1), 0.56),
                                     pLER_m_i   = 0.79,          
                                     freqs_i    = c(0.14, 0.16, 0.20, 0.25, 0.33, 0.5),
                                     target_i   = c(num_eu_supply*0.25, num_eu_supply*0.5, num_eu_supply*0.75, num_eu_supply),
                                     max_surf_i = num_eu_cropland25) %>% 
  unite(col = "scenario", c(pLER_s_i, pLER_m_i, freqs_i, target_i, max_surf_i), sep = "_", remove = F) %>% 
  mutate(pLER_s_i = as.numeric(as.character(pLER_s_i)),
         pLER_m_i = as.numeric(as.character(pLER_m_i)))

dim(simul_plan_restricted)
#144 scenarios

# > Full set of simulations
simul_plan_full <- expand.grid(pLER_s_i = seq(0.3, 0.7, by=0.01),
                               pLER_m_i = seq(0.5, 0.9, by=0.01),          
                               freqs_i = c(0.14, 0.16, 0.20, 0.25, 0.33, 0.5),
                               target_i = c(num_eu_supply*0.25, 
                                            num_eu_supply, 
                                            num_eu_supply*0.75, 
                                            num_eu_supply*0.5),
                               max_surf_i = num_eu_cropland25) %>% 
  unite(col = "scenario", c(pLER_s_i, pLER_m_i, freqs_i, target_i, max_surf_i), sep = "_", remove = F) %>% 
  mutate(pLER_s_i = as.numeric(as.character(pLER_s_i)),
         pLER_m_i = as.numeric(as.character(pLER_m_i)))

dim(simul_plan_full)
#40344 scenarios

# > Path to save the results 
path_alloc <- "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/08_allocations/allocations_with_max_surf/"

# ----------------------------------------
# -----   {Sensitivity analyses}   -----
# ----------------------------------------

# > sensitivity analysis with different yield threshold
list_sensi_1 <- list(
  "sensi_1" = list(data             = data_for_allocation_sensi1, 
                   tested_scenarios = simul_plan_restricted, 
                   save_name        = "sensi/sensi_1")
)

# > sensitivity analysis including countries outside of EU
list_sensi_2 <- list(
  "sensi_2_1" = list(data             = data_for_allocation_sensi2_1,  
                     tested_scenarios = simul_plan_restricted, 
                     save_name        = "sensi/sensi_2_1"),
  "sensi_2_2" = list(data             = data_for_allocation_sensi2_2,  
                     tested_scenarios = simul_plan_restricted, 
                     save_name        = "sensi/sensi_2_2")
)

# > sensitivity analyses with different yield predictive models 
list_sensi_3 <- list(
  "sensi_3_1" = list(data             = data_for_allocation_avg_m,   
                     tested_scenarios = simul_plan_restricted, 
                        save_name     = "sensi/sensi_3_1_avg_m"),
  "sensi_3_2" = list(data             = data_for_allocation_avg_s,   
                     tested_scenarios = simul_plan_restricted, 
                     save_name        = "sensi/sensi_3_2_avg_s"),
  "sensi_3_3" = list(data             = data_for_allocation_pca_m_3, 
                     tested_scenarios = simul_plan_restricted, 
                     save_name        = "sensi/sensi_3_3_pca_m_3")
)

# > sensitivity analyses examining the performance of intercropping vs sole cropping 
# according to pLERs values (takes a lot of time because 40344 scenarios tested)

split_simul <- simul_plan_full %>% split(list(.$target_i, .$freqs_i))
list_sensi_4 <- list(
  "full_simulation_25_14"   = list(data = data_for_allocation,  tested_scenarios = split_simul$'14600000.0.14', save_name = "sensi_pLERs/sensi_25_14"),
  "full_simulation_25_16"   = list(data = data_for_allocation,  tested_scenarios = split_simul$'14600000.0.16', save_name = "sensi_pLERs/sensi_25_16"),
  "full_simulation_25_20"   = list(data = data_for_allocation,  tested_scenarios = split_simul$'14600000.0.2', save_name = "sensi_pLERs/sensi_25_20"),
  "full_simulation_25_25"   = list(data = data_for_allocation,  tested_scenarios = split_simul$'14600000.0.25', save_name = "sensi_pLERs/sensi_25_25"),
  "full_simulation_25_33"   = list(data = data_for_allocation,  tested_scenarios = split_simul$'14600000.0.33', save_name = "sensi_pLERs/sensi_25_33"),
  "full_simulation_25_50"   = list(data = data_for_allocation,  tested_scenarios = split_simul$'14600000.0.5', save_name = "sensi_pLERs/sensi_25_50"),
  
  "full_simulation_50_14"   = list(data = data_for_allocation,  tested_scenarios = split_simul$'29200000.0.14', save_name = "sensi_pLERs/sensi_50_14"),
  "full_simulation_50_16"   = list(data = data_for_allocation,  tested_scenarios = split_simul$'29200000.0.16', save_name = "sensi_pLERs/sensi_50_16"),
  "full_simulation_50_20"   = list(data = data_for_allocation,  tested_scenarios = split_simul$'29200000.0.2', save_name = "sensi_pLERs/sensi_50_20"),
  "full_simulation_50_25"   = list(data = data_for_allocation,  tested_scenarios = split_simul$'29200000.0.25', save_name = "sensi_pLERs/sensi_50_25"),
  "full_simulation_50_33"   = list(data = data_for_allocation,  tested_scenarios = split_simul$'29200000.0.33', save_name = "sensi_pLERs/sensi_50_33"),
  "full_simulation_50_50"   = list(data = data_for_allocation,  tested_scenarios = split_simul$'29200000.0.5', save_name = "sensi_pLERs/sensi_50_50"),
  
  "full_simulation_75_14"   = list(data = data_for_allocation,  tested_scenarios = split_simul$'43800000.0.14', save_name = "sensi_pLERs/sensi_75_14"),
  "full_simulation_75_16"   = list(data = data_for_allocation,  tested_scenarios = split_simul$'43800000.0.16', save_name = "sensi_pLERs/sensi_75_16"),
  "full_simulation_75_20"   = list(data = data_for_allocation,  tested_scenarios = split_simul$'43800000.0.2', save_name = "sensi_pLERs/sensi_75_20"),
  "full_simulation_75_25"   = list(data = data_for_allocation,  tested_scenarios = split_simul$'43800000.0.25', save_name = "sensi_pLERs/sensi_75_25"),
  "full_simulation_75_33"   = list(data = data_for_allocation,  tested_scenarios = split_simul$'43800000.0.33', save_name = "sensi_pLERs/sensi_75_33"),
  "full_simulation_75_50"   = list(data = data_for_allocation,  tested_scenarios = split_simul$'43800000.0.5', save_name = "sensi_pLERs/sensi_75_50"),
  
  "full_simulation_100_14"   = list(data = data_for_allocation,  tested_scenarios = split_simul$'58400000.0.14', save_name = "sensi_pLERs/sensi_100_14"),
  "full_simulation_100_16"   = list(data = data_for_allocation,  tested_scenarios = split_simul$'58400000.0.16', save_name = "sensi_pLERs/sensi_100_16"),
  "full_simulation_100_20"   = list(data = data_for_allocation,  tested_scenarios = split_simul$'58400000.0.2', save_name = "sensi_pLERs/sensi_100_20"),
  "full_simulation_100_25"   = list(data = data_for_allocation,  tested_scenarios = split_simul$'58400000.0.25', save_name = "sensi_pLERs/sensi_100_25"),
  "full_simulation_100_33"   = list(data = data_for_allocation,  tested_scenarios = split_simul$'58400000.0.33', save_name = "sensi_pLERs/sensi_100_33"),
  "full_simulation_100_50"   = list(data = data_for_allocation,  tested_scenarios = split_simul$'58400000.0.5', save_name = "sensi_pLERs/sensi_100_50")
  
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
list(list_sensi_1,
     list_sensi_2,
     list_sensi_3) %>%
     #list_sensi_4) %>% 
  map(., ~{
    
    map(.x, ~{
      
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
          pLER_s_i_x   <- unique(.x$pLER_s_i)
          pLER_m_i_x   <- unique(.x$pLER_m_i)
          freqs_i_x    <- unique(.x$freqs_i)
          max_surf_i_x <- unique(.x$max_surf_i)
          
          # > Allocation
          alloc_i  <- suppressMessages(allocation_intercropping_solecropping(data              = sensi_data,
                                                                             target_production = target_i_x,  
                                                                             max_surf          = max_surf_i_x,
                                                                             pLER_crop         = pLER_s_i_x, 
                                                                             freq_crop         = freqs_i_x))
          
          # > Results of allocation
          res_i <- suppressMessages(results_allocation(data_allocation = alloc_i, 
                                                       pLERs           = c(pLER_s_i_x, pLER_m_i_x), 
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
    
    cat("\n")
    
  })
  
# >>> stop cluster//
stopCluster(my.cluster)


