# -------------------------------------------------------------------------
# 
#       02-1. Prepare data containing yield, climate, and irrigation data in EU
#       Author: M. Chen, Inrae, 2024
#         
# -------------------------------------------------------------------------

# ----------------------------------
# > Packages
library(tidyverse) ; library(stringr) ; library(lubridate) ; library(CCMHr)
library(cowplot)
library(terra) ; library(raster)
library(parallel) ; library(doParallel); library(foreach)
# > machine learning
library(caret) ; library(ranger) ; library(fastshap)

# Home-made functions performing the dimension reductions
source(".../00_Functions_dimension_reduction.R")

# > path to climatic data
path_day <- "..."
path_month <- "..."

# ----------------------------------------
# Compute monthly averages from daily averages

monthly_average <- function(var_i, 
                            load = T, 
                            data_var_i = NULL, 
                            crop = NULL){ # crop = "maize" or "soybean"
  
  # > load data 
  if(load==TRUE)
  {
    if(is.null(crop) == T){ message("No crop is provided in the arguments") } 
    data <- loadRDa(paste0(path_day, "/", crop, "/era5daily_data_", var_i, ".rda"))
  }
  
  # > or use the data provided
  if(load==FALSE)
  {
    data <- data_var_i
    
  }
  
  # > compute means for raw and accumulated variables
  data_m <- data %>% 
    # > retrieve month from date 
    mutate(month = month(date),
           year = year(date)) %>% 
    # > compute monthly averages for the variable (raw and accumulated)
    group_by(site_year, x, y, gridcode, country_name, year, month, clim.var) %>% 
    summarise(mean_clim.value = mean(clim.value),
              mean_cum_clim.value = mean(cum_clim.value)) %>% 
    # > rename monthly
    arrange(month)
  
  return(data_m)
  
}

# ----------------------------------
# Data
# - irrigation (SPAM2010)
# - climatic variables (from ERA5 database)

# -------------------
# Irrigation
# retrieved from the SPAM dataset (accessible at: https://doi.org/10.7910/DVN/PRFF8V)

# Soybean 
raster_irrigation_s <- raster::raster(".../SPAM/spam2010v2r0_global_harv_area.geotiff/spam2010V2r0_global_H_SOYB_I.tif"); raster_irrigation_s
#class      : RasterLayer 
#dimensions : 2160, 4320, 9331200  (nrow, ncol, ncell)
#resolution : 0.083333, 0.083333  (x, y)
#extent     : -180, 179.9986, -89.99928, 90  (xmin, xmax, ymin, ymax)
#crs        : +proj=longlat +datum=WGS84 +no_defs 
#source     : spam2010V2r0_global_H_SOYB_I.tif 
#names      : spam2010V2r0_global_H_SOYB_I 

# > resolution reduction (from 0.083333 -> 0.5) 
# check the new resolution
agg_irrigation_s <- terra::aggregate(raster_irrigation_s, fact=6, fun="mean") ; agg_irrigation_s
#class      : RasterLayer 
#dimensions : 360, 720, 259200  (nrow, ncol, ncell)
#resolution : 0.499998, 0.499998  (x, y)
#extent     : -180, 179.9986, -89.99928, 90  (xmin, xmax, ymin, ymax)
#crs        : +proj=longlat +datum=WGS84 +no_defs 
#source     : memory
#names      : spam2010V2r0_global_H_SOYB_I 
#values     : 0, 2215.939  (min, max)

# > table irrigated portion
irrigation_s_tab <- data.frame(irrigated_portion = values(agg_irrigation_s),
                               coordinates(agg_irrigation_s)) %>%  
  # > round x and y to 2 digits to be consistent with yield data 
  mutate(x=round(x,2),
         y=round(y,2)) %>% 
  # transform NA into 0
  mutate(irrigated_portion =if_else(is.na(irrigated_portion)==T, 0, irrigated_portion)) %>% 
  # proportion of area in %., we need to /100 to have %
  mutate(irrigated_portion_perc = irrigated_portion/100) %>% 
  unite("gridcode", x, y, remove=F)

# check values of irrigation 
summary(irrigation_s_tab$irrigated_portion_perc)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.000000  0.000000  0.000000  0.005443  0.000000 22.159389

# > Maize
raster_irrigation_m <- raster::raster(".../SPAM/spam2010v2r0_global_harv_area.geotiff/spam2010V2r0_global_H_MAIZ_I.tif"); raster_irrigation_m
#class      : RasterLayer 
#dimensions : 2160, 4320, 9331200  (nrow, ncol, ncell)
#resolution : 0.083333, 0.083333  (x, y)
#extent     : -180, 179.9986, -89.99928, 90  (xmin, xmax, ymin, ymax)
#crs        : +proj=longlat +datum=WGS84 +no_defs 
#source     : spam2010V2r0_global_H_MAIZ_I.tif 
#names      : spam2010V2r0_global_H_MAIZ_I 

# > resolution reduction (from 0.083333 -> 0.5) 
# check the new resolution
agg_irrigation_m <- terra::aggregate(raster_irrigation_m, fact=6, fun="mean") ; agg_irrigation_m
#class      : RasterLayer 
#dimensions : 360, 720, 259200  (nrow, ncol, ncell)
#resolution : 0.499998, 0.499998  (x, y)
#extent     : -180, 179.9986, -89.99928, 90  (xmin, xmax, ymin, ymax)
#crs        : +proj=longlat +datum=WGS84 +no_defs 
#source     : memory
#names      : spam2010V2r0_global_H_MAIZ_I 
#values     : 0, 5624.472  (min, max)

# > table irrigated portion
irrigation_m_tab <- data.frame(irrigated_portion = values(agg_irrigation_m),
                               coordinates(agg_irrigation_m)) %>%  
  # > round x and y to 2 digits to be consistent with yield data 
  mutate(x=round(x,2),
         y=round(y,2)) %>% 
  # transform NA into 0
  mutate(irrigated_portion =if_else(is.na(irrigated_portion)==T, 0, irrigated_portion)) %>% 
  # proportion of area in %., we need to /100 to have %
  mutate(irrigated_portion_perc = irrigated_portion/100) %>% 
  unite("gridcode", x, y, remove=F)

# check values of irrigation 
summary(irrigation_m_tab$irrigated_portion_perc)
# Min.     1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.00000  0.00000  0.00000  0.03623  0.00000 56.24472

# -------------------
# Climate in Europe 

# > Store data in a list  
list_climate <- list()

# > Climate variables names and abrevations
vars_names <- data.frame(clim.var = c("max_2m_temperature", "min_2m_temperature", 
                                      "et0", 
                                      "surface_net_solar_radiation", 
                                      "total_precipitation", "vapor_pressure_deficit_1")) %>% 
  mutate(clim.var_abb = recode(clim.var, 
                               "min_2m_temperature"         ="min_temp",
                               "max_2m_temperature"         ="max_temp",
                               "et0"                        ="et0",
                               "surface_net_solar_radiation"="rad",
                               "total_precipitation"        ="prec",
                               "vapor_pressure_deficit_1"   ="vpd_1"))  %>% 
  mutate(clim.var_lab = recode(clim.var, 
                               "min_2m_temperature"         ="Minimum temperature",
                               "max_2m_temperature"         ="Maximum temperature",
                               "et0"                        ="Evapotranspiration ref",
                               "surface_net_solar_radiation"="Solar radiations",
                               "total_precipitation"        ="Precipitation",
                               "vapor_pressure_deficit_1"   ="Vapor pressure deficit"))

# > Load data and compute monthly and annual averages 
for(var_i in c('max_2m_temperature', 'min_2m_temperature', 'surface_net_solar_radiation', 
               'et0', 'total_precipitation', 
               'vapor_pressure_deficit_1'))
{
  
  var_i_abb <- vars_names[which(vars_names$clim.var==var_i), "clim.var_abb"]
  
  # > Load climatic data
  era5_var_i <- loadRDa(paste0(path_day, "/era5daily_", var_i_abb, "_EU.rda"))
  
  # -------------------------
  # > Compute monthly averages based on daily data 
  # all variables except precipitations
  if(var_i %in% c('max_2m_temperature', 'min_2m_temperature', 
                  'surface_net_solar_radiation', 'et0', 
                  'vapor_pressure_deficit_1'))
  {
    
    # > Compute monthly averages 
    tab_month_annual <- monthly_average(var_i = var_i, 
                                        load = F, 
                                        data_var_i = era5_var_i) %>% 
      group_by(site_year, x, y, gridcode, country_name, clim.var) %>%
      # > identify month id
      arrange(year, month) %>%
      mutate(month_nb = row_number()) %>% 
      # > wider format 
      ungroup() %>% 
      dplyr::select(-mean_cum_clim.value, -month, -year) %>%
      pivot_wider(names_from = c("clim.var", "month_nb"), 
                  values_from = "mean_clim.value", 
                  names_prefix = "monthly_", 
                  names_sep = "_")
    
  }
  
  # -------------------------
  # > In the case of precipitation, which are already 
  # computed at a monthly time step
  if(var_i == "total_precipitation")
  {
    
    # > Change format of monthly data 
    tab_month_annual <- era5_var_i %>% 
      group_by(x, y, year, gridcode, country_name, clim.var) %>%
      # > identify month id
      arrange(year, month) %>%
      mutate(month_nb = row_number()) %>% 
      # > remove useless column
      ungroup() %>% 
      dplyr::select(-month) %>%
      pivot_wider(names_from = c("clim.var", "month_nb"), 
                  values_from = "clim.value", 
                  names_prefix = "monthly_", 
                  names_sep = "_") %>% 
      mutate(site_year = paste0(gridcode, "_", year))
    
  }
  
  # -------------------------
  # > Compute annual averages 
  # soybean (only 7 months)
  tab_month_annual$year_average_soybean <- tab_month_annual %>% 
    dplyr::select(starts_with("monthly")) %>% 
    dplyr::select(-ends_with("_8")) %>% 
    apply(., MARGIN = 1, mean)
  
  # maize (8 months)
  tab_month_annual$year_average_maize <- tab_month_annual %>% 
    dplyr::select(starts_with("monthly")) %>% 
    apply(., MARGIN = 1, mean)
  
  # -------------------------
  # > Store data
  # > month & years 
  list_climate[[paste0(var_i_abb)]] <- tab_month_annual
  
  # -------------------------
  # > Remove the data 
  rm(var_i_abb, era5_var_i, tab_month_annual)
  
} 

# > check
list_climate %>% map(., ~{ names(.x)})

# > rename year averages 
list_climate$max_temp$year_max_2m_temperature_soybean     <- list_climate$max_temp$year_average_soybean
list_climate$min_temp$year_min_2m_temperature_soybean     <- list_climate$min_temp$year_average_soybean
list_climate$rad$year_surface_net_solar_radiation_soybean <- list_climate$rad$year_average_soybean
list_climate$et0$year_et0_soybean                         <- list_climate$et0$year_average_soybean
list_climate$prec$year_total_precipitation_soybean        <- list_climate$prec$year_average_soybean
list_climate$vpd_1$year_vapor_pressure_deficit_soybean    <- list_climate$vpd_1$year_average_soybean

list_climate$max_temp$year_max_2m_temperature_maize     <- list_climate$max_temp$year_average_maize
list_climate$min_temp$year_min_2m_temperature_maize     <- list_climate$min_temp$year_average_maize
list_climate$rad$year_surface_net_solar_radiation_maize <- list_climate$rad$year_average_maize
list_climate$et0$year_et0_maize                         <- list_climate$et0$year_average_maize
list_climate$prec$year_total_precipitation_maize        <- list_climate$prec$year_average_maize
list_climate$vpd_1$year_vapor_pressure_deficit_maize    <- list_climate$vpd_1$year_average_maize

list_climate %>% map(., ~{ names(.x)})

list_climate %>% map_dfr(., ~{ nrow(.x)}, .id="var")
list_climate %>% map_dfr(., ~{ length(unique(.x$gridcode))}, .id="var")

# > Merge all climatic variables in 1 tab
tab_climate_EU <- list_climate$min_temp %>% dplyr::select(-year_average_soybean, -year_average_maize) %>% 
     left_join(., list_climate$max_temp %>% dplyr::select(-year_average_soybean, -year_average_maize), by = c("site_year", "x", "y", "gridcode", "country_name")) %>% 
     left_join(., list_climate$rad      %>% dplyr::select(-year_average_soybean, -year_average_maize), by = c("site_year", "x", "y", "gridcode", "country_name")) %>% 
     left_join(., list_climate$et0      %>% dplyr::select(-year_average_soybean, -year_average_maize), by = c("site_year", "x", "y", "gridcode", "country_name")) %>% 
     left_join(., list_climate$prec     %>% dplyr::select(-year_average_soybean, -year_average_maize), by = c("site_year", "x", "y", "gridcode", "country_name")) %>% 
     left_join(., list_climate$vpd_1    %>% dplyr::select(-year_average_soybean, -year_average_maize), by = c("site_year", "x", "y", "gridcode", "country_name")) 

summary(tab_climate_EU)

# -------------------------------------------------------------------------
# Transform data based on PCA scores derived from full dataset
# Data 

# Soybean
tab_climate_EU_soybean <- tab_climate_EU %>% 
  dplyr::select(-ends_with("_8"), -starts_with("year_"))

names(tab_climate_EU_soybean)

# Load PCA loads and scores derived from PCA on monthly averages at global scale
load(".../data/00_tab_soybean.rda")
pca_soybean <- loadRDa(paste0(path_day, "/soybean/pca_soybean.rda"))

# change name for vpd_1
pca_soybean$list_pca_per_variable$vpd_1 <- pca_soybean$list_pca_per_variable$vapor_pressure_deficit

# > Check
tab_soybean %>% filter(site_year == "-53.25_-18.75_1981") %>% pull(PC2_month_max_temp) # -0.5896709 ou -0.5768624?
pca_soybean$tab_PCA_scores[rownames(pca_soybean$tab_PCA_scores)== "-53.25_-18.75_1981",]$PC2_month_max_temp # -0.5896709 ou -0.5768624?

# > Apply PCA loads on new data 
list_scores_eu <- list()
for(var_j in unique(vars_names$clim.var)) 
{
  
  clim.var_abb_j <- vars_names[which(vars_names$clim.var==var_j),]$clim.var_abb
  
  # > PCA loads
  loads_world <- pca_soybean$list_pca_per_variable[[paste0(clim.var_abb_j)]]$pca$rotation
  
  # > Mean and sd of the original data 
  mu_world <- tab_soybean %>% ungroup(.) %>%
    dplyr::select(starts_with(paste0("monthly_", var_j))) %>% 
    colMeans(.)
  
  sd_world <- tab_soybean %>% ungroup(.) %>%
    dplyr::select(starts_with(paste0("monthly_", var_j))) %>% 
    apply(., 2, sd)
  
  # > New data (standardized based on mu and sd in the full dataset)
  tab_climate_EU_j <- tab_climate_EU_soybean %>% ungroup(.) %>%
    dplyr::select(site_year, starts_with(paste0("monthly_", var_j)))
  
  # > Scale new data based on mean and sd of the original data 
  z_tab_climate_EU_j <- scale(tab_climate_EU_j[,-1], center = mu_world, scale=sd_world)
  
  list_scores_j <- list()
  # > Compute scores for the var_j
  for(k in 1:ncol(loads_world))
  {
    
    # > compute scores based on loads and scaled data 
    score_j_k <- sapply(1:ncol(z_tab_climate_EU_j), function(x) z_tab_climate_EU_j[,x] * loads_world[x,k] ) %>%
      apply(., 1, sum)
    # > store in list 
    list_scores_j[[paste0(k)]] <- data.frame(site_year = tab_climate_EU_j$site_year,
                                             score    = score_j_k, 
                                             clim.var = clim.var_abb_j)
    
  }
  
  # > Table formatting
  tab_scores_EU_j <- map_dfr(list_scores_j, data.frame, .id = "id_score") %>% 
    mutate(PC = paste0("PC", id_score, "_month_", clim.var_abb_j)) %>% 
    dplyr::select(-id_score, -clim.var) %>% 
    spread(key = PC, value = score)
  
  # > Store table of all scores for var_j
  list_scores_eu[[paste0(clim.var_abb_j)]] <- tab_scores_EU_j
  
  # > Remove unused objects 
  rm(loads_world, mu_world, sd_world, tab_climate_EU_j, z_tab_climate_EU_j, list_scores_j, tab_scores_EU_j)
  
}

# > PC scores in table
tab_scores_EU_soybean <- list_scores_eu %>% map_dfc(., data.frame) %>% 
  dplyr::select(1, starts_with("PC")) %>% 
  rename("site_year" = 1)

summary(tab_scores_EU_soybean)

rm(list_scores_eu, pca_soybean, tab_soybean)

# -------------------
# Maize
tab_climate_EU_maize <- tab_climate_EU %>% 
  dplyr::select(-starts_with("year_"))

names(tab_climate_EU_maize)

# Load PCA loads and scores derived from PCA on monthly averages at global scale
load(".../data/00_tab_maize.rda")
pca_soybean <- loadRDa(paste0(path_day, "/maize/pca_maize.rda"))

pca_maize$list_pca_per_variable$vpd_1 <- pca_maize$list_pca_per_variable$vapor_pressure_deficit

# > Check
tab_maize %>% filter(site_year == "-71.75_45.25_1981") %>% pull(PC2_month_max_temp) #  0.01785573 ou  0.02131632?
pca_maize$tab_PCA_scores[rownames(pca_maize$tab_PCA_scores)== "-71.75_45.25_1981",]$PC2_month_max_temp #  0.01785573 ou  0.02131632?

# > Apply PCA loads on new data 
list_scores_eu <- list()

for(var_j in unique(vars_names$clim.var)) 
{
  
  clim.var_abb_j <- vars_names[which(vars_names$clim.var==var_j),]$clim.var_abb
  
  # > PCA loads
  loads_world <- pca_maize$list_pca_per_variable[[paste0(clim.var_abb_j)]]$pca$rotation
  
  # > Mean and sd of the original data 
  mu_world <- tab_maize %>% ungroup(.) %>%
    dplyr::select(starts_with(paste0("monthly_", var_j))) %>% 
    colMeans(.)
  
  sd_world <- tab_maize %>% ungroup(.) %>%
    dplyr::select(starts_with(paste0("monthly_", var_j))) %>% 
    apply(., 2, sd)
  
  # > New data (standardized based on mu and sd in the full dataset)
  tab_climate_EU_j <- tab_climate_EU_maize %>% ungroup(.) %>%
    dplyr::select(site_year, starts_with(paste0("monthly_", var_j)))
  
  # > Scale new data based on mean and sd of the original data 
  z_tab_climate_EU_j <- scale(tab_climate_EU_j[,-1], center = mu_world, scale=sd_world)
  
  list_scores_j <- list()
  # > Compute scores for the var_j
  for(k in 1:ncol(loads_world))
  {
    
    # > compute scores based on loads and scaled data 
    score_j_k <- sapply(1:ncol(z_tab_climate_EU_j), function(x) z_tab_climate_EU_j[,x] * loads_world[x,k] ) %>%
      apply(., 1, sum)
    # > store in list 
    list_scores_j[[paste0(k)]] <- data.frame(site_year = tab_climate_EU_j$site_year,
                                             score    = score_j_k, 
                                             clim.var = clim.var_abb_j)
    
  }
  
  # > Table formatting
  tab_scores_EU_j <- map_dfr(list_scores_j, data.frame, .id = "id_score") %>% 
    mutate(PC = paste0("PC", id_score, "_month_", clim.var_abb_j)) %>% 
    dplyr::select(-id_score, -clim.var) %>% 
    spread(key = PC, value = score)
  
  # > Store table of all scores for var_j
  list_scores_eu[[paste0(clim.var_abb_j)]] <- tab_scores_EU_j
  
  # > Remove unused objects 
  rm(loads_world, mu_world, sd_world, tab_climate_EU_j, z_tab_climate_EU_j, list_scores_j, tab_scores_EU_j)
  
}

# > PC scores in table
tab_scores_EU_maize <- list_scores_eu %>% map_dfc(., data.frame) %>% 
  dplyr::select(1, starts_with("PC")) %>% 
  rename("site_year" = 1)

summary(tab_scores_EU_maize)

# -------------------
# Merge all data

# > Soybean
tab_soybean_EU <- tab_climate_EU %>% 
  dplyr::select(-ends_with("_8"), -ends_with("_maize")) %>%
  # > rename annual averages
  rename(
    "year_max_2m_temperature"="year_max_2m_temperature_soybean",
    "year_min_2m_temperature"="year_min_2m_temperature_soybean",
    "year_surface_net_solar_radiation"="year_surface_net_solar_radiation_soybean",
    "year_et0"="year_et0_soybean",
    "year_total_precipitation"="year_total_precipitation_soybean",
    "year_vapor_pressure_deficit"="year_vapor_pressure_deficit_soybean") %>% 
  # > add PCA scores
  left_join(., tab_scores_EU_soybean, by = "site_year") %>% 
  # > add irrigation 
  left_join(., irrigation_s_tab, by = c("gridcode", "x", "y"))

dim(tab_soybean_EU) # 100608
length(unique(tab_soybean_EU$gridcode)) # 4192
summary(tab_soybean_EU)

# > Maize
tab_maize_EU <- tab_climate_EU %>% 
  dplyr::select(-ends_with("_soybean")) %>%
  # > rename annual averages
  rename(
    "year_max_2m_temperature"="year_max_2m_temperature_maize",
    "year_min_2m_temperature"="year_min_2m_temperature_maize",
    "year_surface_net_solar_radiation"="year_surface_net_solar_radiation_maize",
    "year_et0"="year_et0_maize",
    "year_total_precipitation"="year_total_precipitation_maize",
    "year_vapor_pressure_deficit"="year_vapor_pressure_deficit_maize") %>% 
  # > add PCA scores
  left_join(., tab_scores_EU_maize, by = "site_year") %>% 
  # > add irrigation 
  left_join(., irrigation_m_tab, by = c("gridcode", "x", "y"))

dim(tab_maize_EU) # 100608    
length(unique(tab_maize_EU$gridcode)) # 4192
summary(tab_maize_EU)

# ---------------------
# Save
save(tab_soybean_EU, file = ".../data/02_tab_eu_soybean.rda")
save(tab_maize_EU, file = ".../data/02_tab_eu_maize.rda")

# ---------------------
stop()

  



