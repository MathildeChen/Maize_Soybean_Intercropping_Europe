# ALLOCATION OF 1 CROP SPECIES
# The objective is to produce a target quantity of a given crop species
# in a given region (e.g., country, continent).  
# - data: a data.frame containing for each site-year
#        > x,y:     coordinates of the pixels
#        > year:    year (here, 2000-2023)
#        > Ya_pred: sole crop yield predictions (t.ha) 
#        > area:    available area to grow the crop 
# - target_production: production target of the species (in t)
# - max_surf: maximum surface allocated
# - freq_crop: the frequence of crop in the rotation (range between 0 and 1)
allocation <- function(data, 
                       target_production,
                       max_surf)
{
  
  # -------------------------
  # 1. Compute the mean production by pixel
  #    and order the pixels by productivity
  #    (the most productive first)
  data_prod <- data %>% 
    group_by(x,y,area) %>% 
    # > mean yield over time
    summarise(mean_Ya = mean(Ya_pred)) %>% 
    # > mean productivity (yield * area)
    mutate(mean_prod = mean_Ya*area) %>% 
    arrange(desc(mean_prod)) %>% 
    ungroup()
  
  # 2. Identify the required area to reach
  #    target production specified 
  data_area_req <- data_prod %>% 
    mutate(cum_prod      = cumsum(mean_prod),
           cum_surf      = cumsum(area)) %>% 
    mutate(pixel_to_keep = if_else(lag(cum_prod)>target_production | lag(cum_surf)>max_surf, 0, 1),          # surface needed to reach target_supply
           pixel_to_keep = if_else(is.na(pixel_to_keep)==T, 1, pixel_to_keep)) %>%  # correct the 1rst pixel
    dplyr::select(x, y, pixel_to_keep)
  
  # 3. Merge with the initial data
  data_out <- data %>% 
    left_join(data_area_req, by = c("x", "y"))
  
  # 4. Output
  return(data_out)
  
}

#----------------------
# ALLOCATION OF 1 CROP SPECIES IN INTERCROPPING AND IN SOLE CROPPING 
# The objective is to produce a target quantity of a given crop species
# in a given region (e.g., country, continent) 
# using 2 different cropping strategies: intercropping or sole cropping
# Arguments: 
# - data: a data.frame containing for each site-year
#        > x,y:     coordinates of the pixels
#        > year:    year (here, 2000-2023)
#        > Ya_pred: sole crop yield predictions (t.ha)  
#        > area:    available area to grow the crop 
# - target_production: production target of the species (in t)
# - max_surf: maximum surface allocated
# - freq_crop: the frequence of crop in the rotation (range between 0 and 1)
allocation_intercropping_solecropping <- function(data,
                                                  target_production, 
                                                  max_surf,
                                                  freq_crop, 
                                                  pLER_crop)
{
  
  # 0. Pre-checks
  # > Enter crop frequency as a proportion
  testthat::expect_false(freq_crop >1, info = "Enter crop frequency as proportion (0 <= freq_crop <= 1)")
  testthat::expect_false(freq_crop <0, info = "Enter crop frequency as proportion (0 <= freq_crop <= 1)")
  
  # > Enter pLER as a proportion
  testthat::expect_false(pLER_crop >1, info = "Enter pLER as proportion (0 <= pLER_crop <= 1)")
  testthat::expect_false(pLER_crop <0, info = "Enter pLER as proportion (0 <= pLER_crop <= 1)")
  
  # -------------------------
  # Allocation in sole cropping 
  # 1. Data format
  #    - area: depending on the frequency of crop in the rotation
  data_solecropping <- data %>% 
    mutate(area = area*freq_crop)
  
  # 2. Allocation 
  allocation_solecropping <- allocation(data = data_solecropping, 
                                        target_production = target_production,
                                        max_surf = max_surf) %>% 
    dplyr::select(x, y, pixel_to_keep) %>% 
    distinct(.) %>% 
    rename("pixel_solecropping"="pixel_to_keep")
  
  # -------------------------
  # Allocation in intercropping 
  # 1. Data format
  #    - area:       depending on the frequency of crop in the rotation
  #    - production: multiply by the pLER of the crop 
  data_intercropping <- data %>% 
    mutate(area    = area*freq_crop,
           Ya_pred = Ya_pred*pLER_crop) 
  
  # 2. Allocation 
  allocation_intercropping <- allocation(data = data_intercropping, 
                                         target_production = target_production,
                                         max_surf = max_surf) %>%
    dplyr::select(x, y, pixel_to_keep) %>% 
    distinct(.) %>% 
    rename("pixel_intercropping"="pixel_to_keep")
  
  # -------------------------
  # Merge both results and identify the pixels
  # not used (pixel_free) and not used in solecropping (pixel_intercropping_free)
  allocation_total <- data %>% 
    left_join(allocation_solecropping,  by = c("x", "y")) %>% 
    left_join(allocation_intercropping, by = c("x", "y")) %>%
    mutate(pixel_intercropping_free = if_else(pixel_intercropping == 1 & pixel_solecropping == 0, 1, 0),
           pixel_free = if_else(pixel_intercropping == 0 & pixel_solecropping == 0 & pixel_intercropping_free == 0, 1, 0))
  
  # Check on pixels allocation
  test <- allocation_total %>% 
    distinct(x, y, pixel_solecropping, pixel_intercropping, pixel_intercropping_free, pixel_free) %>% 
    mutate(check1 = pixel_solecropping + pixel_intercropping_free + pixel_free, 
           check2 = pixel_intercropping + pixel_free) %>% 
    filter(check1 != 1 | check2 != 1)
  testthat::expect_equal(nrow(test), 0, info = "Problem in the allocation; some pixels are allocated twice")
  
  # Output 
  return(allocation_total)
  
}

#----------------------
# EXAMINE RESULTS OF ALLOCATIONS
# Compute the total production and surface resulting of allocation
# in intercropping and in sole cropping 
# - data_allocation: a data.frame containing for each site-year
#        > x,y:     coordinates of the pixels
#        > year:    year (here, 2000-2023)
#        > Ya_pred: sole crop yield predictions (t.ha) of crop 1
#        > Ya_pred_2: sol crop yield predictions of crop 2
#        > area:    available area to grow the crop 
#        > pixel_...: binary variables indicating which pixel is required to 
#          grown the crop: pixel_solecropping pixel_intercropping, pixel_intercropping_free, pixel_free
# - pLER: a vector of partial land equivalent ratios associated with the 2 crops

results_allocation <- function(data_allocation, 
                               pLERs      = c(0.5, 0.5),
                               freq_crops_i = 0.5,
                               freq_crops_m = c(0.5, 0.5))
{
  
  # -------------------
  # 0. Pre-checks
  # > Crop frequency 
  if(freq_crops_m[1] != freq_crops_i | freq_crops_m[2] != freq_crops_i){
    warning("The crop frequencies provided for intercrop and sole crop differ; please check")
  }
  
  # -------------------
  # 1. Based on previous allocation, 
  # compute the mean production per pixel 
  # and the area required to reach the production target
  data_allocation_productions <- data_allocation %>% 
    # > compute mean yields over time
    group_by(x, y, area, pixel_intercropping, pixel_solecropping, pixel_intercropping_free, pixel_free) %>% 
    summarise(mean_Ya_1 = mean(Ya_pred),
              mean_Ya_2 = mean(Ya_pred_2)) %>% 
    # > mean productivity (yield * area)
    mutate(mean_prod_1_i = mean_Ya_1*pLERs[1]*area*freq_crops_i*pixel_intercropping,
           mean_prod_2_i = mean_Ya_2*pLERs[2]*area*freq_crops_i*pixel_intercropping,
           mean_prod_1_m = mean_Ya_1*area*freq_crops_m[1]*pixel_solecropping,
           mean_prod_2_m = mean_Ya_2*area*freq_crops_m[2]*(pixel_intercropping_free+pixel_free)) %>% 
    # > sort pixels according to soybean productivity in intercropping and sole cropping first, and then on maize sole crops
    arrange(desc(mean_prod_1_i), desc(mean_prod_1_m), desc(mean_prod_2_m)) %>% 
    # > total productivity in both strategies
    mutate(total_prod_i = sum(c(mean_prod_1_i, mean_prod_2_i)),
           total_prod_m = sum(c(mean_prod_1_m, mean_prod_2_m))) %>%
    ungroup() %>% 
    mutate(cum_total_prod_i = cumsum(total_prod_i),
           cum_total_prod_m = cumsum(total_prod_m)) %>%
    # > area required to reach the same level of productivity in the 2 strategies
    mutate(pixel_reach_equal_prod = if_else(lag(cum_total_prod_m) > cum_total_prod_i, 0, 1),
           pixel_reach_equal_prod = if_else(is.na(pixel_reach_equal_prod)==T, 1, pixel_reach_equal_prod)) %>% 
    # > land saving from intercropping
    mutate(pixel_landsaving = if_else(pixel_free == 1 & pixel_reach_equal_prod == 1, 1, 0)) %>%
    # Format table 
    dplyr::select(x, y, area, mean_Ya_1, mean_Ya_2, mean_prod_1_i, mean_prod_2_i, mean_prod_1_m, mean_prod_2_m,
                  total_prod_i, total_prod_m, cum_total_prod_i, cum_total_prod_m, 
                  pixel_intercropping, pixel_solecropping, pixel_intercropping_free, pixel_free, pixel_reach_equal_prod, pixel_landsaving)
  
  # -------------------
  # 2. Compute metrics
  # 2.1. For a given surface, how much do we produce for crop 1 and 2?
  # > Productions 
  production_1_i <- sum(data_allocation_productions$mean_prod_1_i)
  production_2_i <- sum(data_allocation_productions$mean_prod_2_i)
  
  production_1_m <- sum(data_allocation_productions$mean_prod_1_m)
  production_2_m <- sum(data_allocation_productions[which(data_allocation_productions$pixel_intercropping_free==1),]$mean_prod_2_m)
  
  # > Surfaces 
  surface_1_i <- sum(data_allocation_productions[which(data_allocation_productions$pixel_intercropping==1),]$area*freq_crops_i)
  surface_2_i <- sum(data_allocation_productions[which(data_allocation_productions$pixel_intercropping==1),]$area*freq_crops_i)
  
  surface_1_m <- sum(data_allocation_productions[which(data_allocation_productions$pixel_solecropping==1),]$area*freq_crops_m[1])
  surface_2_m <- sum(data_allocation_productions[which(data_allocation_productions$pixel_intercropping_free==1),]$area*freq_crops_m[2])
  
  # > Yields
  yield_1_i_pler1 <- (production_1_i/surface_1_i)
  yield_2_i_pler2 <- (production_2_i/surface_2_i)
  
  teta_1_yield_1_m <- (production_1_m/surface_1_m)
  teta_2_yield_2_m <- (production_2_m/surface_2_m)
  
  # 2.2. Comparing efficiency of both strategies
  # > Tetas
  teta_1 <- teta_1_yield_1_m/(yield_1_i_pler1/pLERs[1])
  teta_2 <- teta_2_yield_2_m/(yield_2_i_pler2/pLERs[2])
  
  # > Ratio 1 
  #   if < 1, intercropping produces more maize than monocrop,
  #   if > 1, reverse
  Ratio_1 <- (teta_2 / pLERs[2]) * ( 1 - (pLERs[1] / teta_1))
  
  # > Ratio 2
  #   if > 1, intercropping produces more maize than monocrop,
  #   if < 1, reverse
  Ratio_2 <- (pLERs[1] / teta_1) + (pLERs[2] / teta_2)
  
  # 2.3. Which surface required to produce the same production of crop 1+crop 2?
  # > Total productions per crop
  tot_production_1_i <- sum(data_allocation_productions[which(data_allocation_productions$pixel_intercropping==1),]$mean_prod_1_i)
  tot_production_2_i <- sum(data_allocation_productions[which(data_allocation_productions$pixel_intercropping==1),]$mean_prod_2_i)
  
  tot_production_1_m <- sum(data_allocation_productions[which(data_allocation_productions$pixel_solecropping==1),]$mean_prod_1_m)
  tot_production_2_m <- sum(data_allocation_productions[which(data_allocation_productions$pixel_intercropping_free==1),]$mean_prod_2_m)
  # > production on additional areas
  landsaving_production_2_m <- sum(data_allocation_productions[which(data_allocation_productions$pixel_landsaving==1),]$mean_prod_2_m)
  
  # Checks 
  #tot_production_m <- sum(data_allocation_productions[which(data_allocation_productions$pixel_reach_equal_prod==1),]$total_prod_m)
  #testthat::expect_equal(tot_production_1_m+tot_production_2_m, tot_production_m)
  
  # > Total surfaces per crop
  tot_surface_1_i <- sum(data_allocation_productions[which(data_allocation_productions$pixel_intercropping==1),]$area*freq_crops_i)
  tot_surface_2_i <- sum(data_allocation_productions[which(data_allocation_productions$pixel_intercropping==1),]$area*freq_crops_i)
  
  tot_surface_1_m <- sum(data_allocation_productions[which(data_allocation_productions$pixel_solecropping==1),]$area*freq_crops_m[1])
  tot_surface_2_m <- sum(data_allocation_productions[which(data_allocation_productions$pixel_intercropping_free==1),]$area*freq_crops_m[2])
  # > Land saving 
  landsaving_surface_2_m <- sum(data_allocation_productions[which(data_allocation_productions$pixel_landsaving==1),]$area*freq_crops_m[2])
  
  # Checks (but depends on the crop frequency)
  #tot_surface_m <- sum(data_allocation_productions[which(data_allocation_productions$pixel_reach_equal_prod==1),]$area*freq_crops_m[2])
  #testthat::expect_equal(tot_surface_1_m+tot_surface_2_m+land_saving_2_m, tot_surface_m)
  
  
  # -------------------
  # 3. Output
  
  result_1 <- data.frame(crop       = c("crop 1",       "crop 2",       "crop 1",         "crop 2"),
                         strategy   = c("intercrop",    "intercrop",    "sole crop",      "sole crop"),
                         production = c(production_1_i, production_2_i, production_1_m,   production_2_m),
                         surface    = c(surface_1_i,    surface_2_i,    surface_1_m,      surface_2_m))
  
  result_2 <- data.frame(Teta_1  = teta_1,
                         Teta_2  = teta_2,
                         Ratio_1 = Ratio_1, 
                         Ratio_2 = Ratio_2)
  
  result_3 <- data.frame(crop             = c("crop 1",           "crop 2",           "crop 1",             "crop 2",           "landsaving"),
                         strategy         = c("intercrop",        "intercrop",        "sole crop",          "sole crop",        "sole crop"),
                         total_production = c(tot_production_1_i, tot_production_2_i, tot_production_1_m,   tot_production_2_m, landsaving_production_2_m),
                         total_surface    = c(tot_surface_1_i,    tot_surface_2_i,    tot_surface_1_m,      tot_surface_2_m,    landsaving_surface_2_m))
  
  # -------------------
  # Out
  return(list("data_res" = data_allocation_productions,
              "res1"     = result_1,
              "res2"     = result_2,
              "res3"     = result_3))
  
  
}




