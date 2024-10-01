# -------------------------------------------------------------------------
#
#       03 - Full set of simulations for soybean allocation in Europe
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

# Coordinates of pixels in EU
load("E:/POSTDOC INRAE/ANALYSES/B_OPTIMISATION/00_Data/00_dat_coords_EU27.rda")
dat_coords_eu27$gridcode = paste0(dat_coords_eu27$x, "_", dat_coords_eu27$y)
dim(dat_coords_eu27) # N=2699 

load("E:/POSTDOC INRAE/ANALYSES/B_OPTIMISATION/00_Data/00_dat_coords_EU42.rda")
dat_coords_EU <- dat_coords_eu42
dat_coords_eu42$gridcode = paste0(dat_coords_eu42$x, "_", dat_coords_eu42$y)
dim(dat_coords_eu42) # N=4192


# Predicted yields - europe
load("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_preds/04_EU/Ya_pred_eu_2000_2023.rda")

# > EU27
Ya_pred_eu %>% 
  distinct(x,y) %>% dim(.) # N=4192    

# > EU42
Ya_pred_eu %>%
  filter(id_eu27==1) %>% 
  distinct(x,y) %>% dim(.) # N=2699

# ----------------------------------------
# Sets of pixels with productivity higher than a certain threshold
# > pixels with yields >1 t/ha
pixels_with_enough_soybean <- Ya_pred_eu %>% 
  filter(model=="pca.m.2", crop=="soybean") %>% 
  group_by(x, y) %>% 
  summarise(mean_Ya_pred = mean(Ya_pred)) %>%
  filter(mean_Ya_pred >= 2.677) ; dim(pixels_with_enough_soybean) #  3240 (among 4192) pixels with yields > 1 t/ha

# > Sensitivity analysis : pixels with yields >2 t/ha
pixels_with_enough_soybean2 <- Ya_pred_eu %>% 
  filter(model=="pca.m.2", crop=="soybean") %>% 
  group_by(x, y) %>% 
  summarise(mean_Ya_pred = mean(Ya_pred)) %>%
  filter(mean_Ya_pred >= 2) ; dim(pixels_with_enough_soybean2) #  1986 (among 4192) pixels with yields > 2 t/ha

# > Sensitivity analysis : pixels with yields >2.8 t/ha (mean yield of soybean in Eu27 according to FAO)
pixels_with_enough_soybean3 <- Ya_pred_eu %>% 
  filter(model=="pca.m.2", crop=="soybean") %>% 
  group_by(x, y) %>% 
  summarise(mean_Ya_pred = mean(Ya_pred)) %>%
  filter(mean_Ya_pred >= 2.8) ; dim(pixels_with_enough_soybean2) #  461 (among 4192) pixels with yields > 2.8 t/ha

# > Sensitivity analysis : predictions from other models 
pixels_with_enough_soybean3 <- Ya_pred_eu %>% 
  filter(crop=="soybean") %>% 
  group_by(model, x, y) %>% 
  summarise(mean_Ya_pred = mean(Ya_pred)) %>%
  filter(mean_Ya_pred >= 1) 

pixels_with_enough_soybean3 %>% 
  group_by(model) %>%
  count()

#  model       n
#1 avg.m    3108
#2 avg.s    2640
#3 pca.m.2  3240
#4 pca.m.3  3254

# ----------------------------------------
# Prepare data use for allocation

# MAIN ANALYSIS: Pixels with yields > 1 t/ha in the EU27
data_for_allocation <- Ya_pred_eu %>% 
  filter(model=="pca.m.2") %>% 
  spread(key=crop, value=Ya_pred) %>% 
  left_join(pixels_with_enough_soybean, by=c("x", "y")) %>% 
  filter(id_eu27==1, is.na(mean_Ya_pred)==F) %>% 
  dplyr::select(x, y, year, 
                "area"="cropland_area_ha", 
                "Ya_pred"="soybean", 
                "Ya_pred_2"="maize")

# Number of pixels with yields > 1 t/ha used for the analyses 
data_for_allocation %>% 
  distinct(x,y) %>% dim(.) # N=2047

# ----------------
# Pixels with yields > 1 t/ha in the EU42
data_for_allocation_eu42 <- Ya_pred_eu %>% 
  filter(model=="pca.m.2") %>% 
  spread(key=crop, value=Ya_pred) %>% 
  left_join(pixels_with_enough_soybean, by=c("x", "y")) %>% 
  filter(is.na(mean_Ya_pred)==F) %>% 
  dplyr::select(x, y, year, 
                "area"="cropland_area_ha", 
                "Ya_pred"="soybean", 
                "Ya_pred_2"="maize")

# Number of pixels with yields > 1 t/ha used for the analyses 
data_for_allocation_eu42 %>% 
  distinct(x,y) %>% dim(.) # N=3240

# ----------------
# Pixels with yields > 2 t/ha in the EU27
data_for_allocation2 <- Ya_pred_eu %>% 
  filter(model=="pca.m.2") %>% 
  spread(key=crop, value=Ya_pred) %>% 
  left_join(pixels_with_enough_soybean2, by=c("x", "y")) %>% 
  filter(id_eu27==1, is.na(mean_Ya_pred)==F) %>% 
  dplyr::select(x, y, year, 
                "area"="cropland_area_ha", 
                "Ya_pred"="soybean", 
                "Ya_pred_2"="maize")

# Number of pixels with yields > 1 t/ha used for the analyses 
data_for_allocation2 %>% 
  distinct(x,y) %>% dim(.) # N=1387

# ----------------
# Pixels with yields > 2 t/ha in the EU42
data_for_allocation2_eu42 <- Ya_pred_eu %>% 
  filter(model=="pca.m.2") %>% 
  spread(key=crop, value=Ya_pred) %>% 
  left_join(pixels_with_enough_soybean2, by=c("x", "y")) %>% 
  filter(is.na(mean_Ya_pred)==F) %>% 
  dplyr::select(x, y, year, 
                "area"="cropland_area_ha", 
                "Ya_pred"="soybean", 
                "Ya_pred_2"="maize")

# Number of pixels with yields > 1 t/ha used for the analyses 
data_for_allocation2_eu42 %>% 
  distinct(x,y) %>% dim(.) # N=1986

# ----------------
# Pixels with yields > 1 t/ha in the EU27

# > avg.m
data_for_allocation_avg_m <- Ya_pred_eu %>% 
  filter(model=="avg.m") %>% 
  spread(key=crop, value=Ya_pred) %>% 
  left_join(pixels_with_enough_soybean3[which(pixels_with_enough_soybean3$model=="avg.m"),], by=c("x", "y")) %>% 
  filter(id_eu27==1, is.na(mean_Ya_pred)==F) %>% 
  dplyr::select(x, y, year, 
                "area"="cropland_area_ha", 
                "Ya_pred"="soybean", 
                "Ya_pred_2"="maize")

# Number of pixels with yields > 1 t/ha used for the analyses 
data_for_allocation_avg_m %>% 
  distinct(x,y) %>% dim(.) # N=1984

# > avg.s
data_for_allocation_avg_s <- Ya_pred_eu %>% 
  filter(model=="avg.s") %>% 
  spread(key=crop, value=Ya_pred) %>% 
  left_join(pixels_with_enough_soybean3[which(pixels_with_enough_soybean3$model=="avg.s"),], by=c("x", "y")) %>% 
  filter(id_eu27==1, is.na(mean_Ya_pred)==F) %>% 
  dplyr::select(x, y, year, 
                "area"="cropland_area_ha", 
                "Ya_pred"="soybean", 
                "Ya_pred_2"="maize")

# Number of pixels with yields > 1 t/ha used for the analyses 
data_for_allocation_avg_s %>% 
  distinct(x,y) %>% dim(.) # N=1778

# > pca.m.3
data_for_allocation_pca_m_3 <- Ya_pred_eu %>% 
  filter(model=="pca.m.3") %>% 
  spread(key=crop, value=Ya_pred) %>% 
  left_join(pixels_with_enough_soybean3[which(pixels_with_enough_soybean3$model=="pca.m.3"),], by=c("x", "y")) %>% 
  filter(id_eu27==1, is.na(mean_Ya_pred)==F) %>% 
  dplyr::select(x, y, year, 
                "area"="cropland_area_ha", 
                "Ya_pred"="soybean", 
                "Ya_pred_2"="maize")

# Number of pixels with yields > 1 t/ha used for the analyses 
data_for_allocation_pca_m_3 %>% 
  distinct(x,y) %>% dim(.) # N=2042

# ----------------------------------------
# Simulation plan 

# > Partial Land Equivalent Ratios (pLERs) of soybean and maize
# based on the meta-analysis of Xu et al., 2020
num_pLER_s <- 0.56
num_pLER_m <- 0.79

# > European Union 27 soybean supply
# Based on FAOSTATS - 2018-2022
num_eu_supply <- 36.3*10^6

# > Surface max allocated (25% of the total EU croplands)
# ~ 31.8 Mha
num_eu_cropland25 <- data_for_allocation %>% 
  distinct(x,y,area) %>% 
  pull(area) %>% 
  sum(.)*0.25 ; num_eu_cropland25/10^6

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
# ----------------------------------------
#      Allocation based on each scenario 
# -----{restricted set of simulation}-----

# restricted list of scenarios
n.cores <- parallel::detectCores() - 3 ; n.cores # 7 on the computer, 11 on the laptop

# >>> create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

# register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

allocations <- simul_plan_restricted %>% 
  split(.$scenario) %>% 
  map(., ~{
    
    # > Retrieve values from the selected scenario
    target_i_x   <- unique(.x$target_i)
    pLER_s_i_x   <- unique(.x$pLER_s_i)
    pLER_m_i_x   <- unique(.x$pLER_m_i)
    freqs_i_x    <- unique(.x$freqs_i)
    max_surf_i_x <- unique(.x$max_surf_i)
    
    # > Allocation
    alloc_i  <- suppressMessages(allocation_intercropping_solecropping(data              = data_for_allocation,
                                                                       target_production = target_i_x,  
                                                                       max_surf          = max_surf_i_x,
                                                                       pLER_crop         = pLER_s_i_x, 
                                                                       freq_crop         = freqs_i_x))
    
    # > Results of allocation
    res_i <- suppressMessages(results_allocation(data_allocation=alloc_i,
                                                 pLERs      = c(pLER_s_i_x, pLER_m_i_x),
                                                 freq_crops_i = freqs_i_x,
                                                 freq_crops_m = c(freqs_i_x, freqs_i_x)))
    
    # > Format data prior exportation
    # > Full dataset allocated
    res_i$data_res <- res_i$data_res %>% 
      mutate(id_pixel = case_when(
        pixel_solecropping       == 1 ~ "Soybean sole cropping",
        pixel_intercropping_free == 1 ~ "Maize sole cropping",
        pixel_landsaving         == 1 ~ "Land saving",
        TRUE~NA
      ))
    
    # > Results 
    # Labels
    res_i$res1 <- res_i$res1 %>% 
      mutate(crop = recode(crop, "crop 1"="Soybean", "crop 2"="Maize")) %>%
      mutate(strategy = recode(strategy, "intercrop"="Intercropping", "sole crop"="Sole crop"))
    
    # > Output
    res_i
    
  })

# >>> stop cluster//
stopCluster(my.cluster)
beepr::beep(1)

# Save 
save(allocations, file = paste0(path_alloc, "allocations_soybean_maize_eu_restricted_min_rdt.rda"))

# ----------------------------------------
# -----   {Sensitivity analyses}   -----
# ----------------------------------------

# > sensitivity analyses with different yield predictive models 
list_sensi_1 <- list(
  "avg_m"   = list(data = data_for_allocation_avg_m,   tested_scenarios = simul_plan_restricted, save_name = "sensi/sensi_1_avg_m"),
  "avg_s"   = list(data = data_for_allocation_avg_s,   tested_scenarios = simul_plan_restricted, save_name = "sensi/sensi_1_avg_s"),
  "pca_m_3" = list(data = data_for_allocation_pca_m_3, tested_scenarios = simul_plan_restricted, save_name = "sensi/sensi_1_pca_m_3")
)

# > sensitivity analysis with different yield threshold
list_sensi_2 <- list(
  "yield_2" = list(data = data_for_allocation2, tested_scenarios = simul_plan_restricted, save_name = "sensi/sensi_2_yield_2")
)

# > sensitivity analysis including countries outside of EU
list_sensi_3 <- list(
  "eu_42"   = list(data = data_for_allocation_eu42,  tested_scenarios = simul_plan_restricted, save_name = "sensi/sensi_3_eu42")
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
list(#list_sensi_1,
     #list_sensi_2,
     #list_sensi_3,
     list_sensi_4) %>% 
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


