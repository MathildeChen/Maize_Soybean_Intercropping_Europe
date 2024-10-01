# -------------------------------------------------------------------------
# 
#       01-1. Prepare data containing yield, climate, and irrigation data  
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

# Homemade function
source("E:/POSTDOC INRAE/DATA/01_CLIMATE/ERA5/functions_to_read_era5.R")

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
    data <- loadRDa(paste0("C:/Users/benni/Documents/Post doc/ERA5_daily/", crop, "/era5daily_data_", var_i, ".rda"))
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
# - climatic variables (from ERA5 database)
# - yield data (from GHDY)

# -------------------
# Yield
load("E:/POSTDOC INRAE/DATA/02_YIELDS/GDHY_v1.3/yield_no.trend_full.rda")

# > check if there is some duplicate (normally no)
data_crop %>% 
  distinct(gridcode, x, y, country_name) %>% 
  group_by(gridcode) %>% 
  summarise(n = n()) %>% 
  filter(n>1) # 0 line: OK

# > yield data for soybean
yield_s <- data_crop %>% 
  filter(crop == "soybean") %>% 
  #filter(country_name == "China") %>% 
  mutate(site_year = paste0(gridcode, "_", year))

length(unique(yield_s$gridcode)) # in total 3444
#length(unique(yield_s[which(yield_s$country_name != "Desert"),]$gridcode)) # 2784
#length(unique(yield_s[which(yield_s$country_name == "Desert"),]$gridcode)) # 660

# > yield data for maize
yield_m <- data_crop %>% 
  filter(crop == "maize") %>% 
  #filter(country_name == "China") %>% 
  mutate(site_year = paste0(gridcode, "_", year))

length(unique(yield_m$gridcode)) # 2139
#length(unique(yield_m[which(yield_m$country_name != "Desert"),]$gridcode)) # 1734
#length(unique(yield_m[which(yield_m$country_name == "Desert"),]$gridcode)) # 405

# > load 1 initial yield file to resample era5 data 
yield_ref <- rast("E:/POSTDOC INRAE/DATA/02_YIELDS/GDHY_v1.3/gdhy_v1.2_v1.3_20190128/maize/yield_1981.nc4")

# -------------------------------------------------------------------------
# Outliers in the yield datasets

# Some site-years with yield > 20 t/ha 
out_1 <- yield_m %>% 
  filter(Ya > 20)
to_remove_1 <- unique(out_1$site_year)

# Some site-years with constant yields over time
# > soybean
out_2 <- yield_s %>% 
  filter(country_name != "Desert") %>% 
  group_by(gridcode) %>% 
  summarise(sd_Ya = sd(Ya)) %>% 
  filter(sd_Ya < 10^-6) # 158
to_remove_2 <- unique(out_2$gridcode)

# > maize
out_3 <- yield_m %>% 
  filter(country_name != "Desert") %>% 
  group_by(gridcode) %>% 
  summarise(sd_Ya = sd(Ya)) %>% 
  filter(sd_Ya < 10^-6) # 74
to_remove_3 <- unique(out_3$gridcode)

# Remove the outliers 
# > soybean 
yield_s_corrected <- yield_s %>% 
  mutate(to_remove = if_else(gridcode %in% to_remove_2, 1, 0)) %>% 
  filter(to_remove != 1) 

dim(yield_s_corrected) # 116516

# > maize
yield_m_corrected <- yield_m %>% 
  mutate(to_remove = if_else(site_year %in% to_remove_1 | gridcode %in% to_remove_3, 1, 0)) %>% 
  filter(to_remove != 1) %>% 
  # site with no climate data in the ERA5-land daily estimates
  filter(gridcode != "179.75_65.25")

dim(yield_m_corrected) # 72964

# How much sites and site-years lost?
# > soybean
length(unique(yield_s$site_year))-length(unique(yield_s_corrected$site_year)) # -5000 site-years
length(unique(yield_s$gridcode))-length(unique(yield_s_corrected$gridcode)) # -158 site

# > maize
length(unique(yield_m$site_year))-length(unique(yield_m_corrected$site_year)) # -2598 site-years
length(unique(yield_m$gridcode))-length(unique(yield_m_corrected$gridcode)) # -75 site

# What is the proportion of Desert in the new data? (should be ~20%)
# > soybean 
yield_s_corrected %>% 
  mutate(is_desert=if_else(country_name=="Desert", 1, 0)) %>% 
  group_by(is_desert) %>% 
  count() %>% 
  ungroup() %>%
  mutate(total=sum(n)) %>% 
  mutate(freq=n/total)
#   is_desert     n  total  freq
# 1         0 92756 116516 0.796
# 2         1 23760 116516 0.204

# > maize
yield_m_corrected %>% 
  mutate(is_desert=if_else(country_name=="Desert", 1, 0)) %>% 
  group_by(is_desert) %>% 
  count() %>% 
  ungroup() %>%
  mutate(total=sum(n)) %>% 
  mutate(freq=n/total)
#  is_desert     n  total  freq
# 1         0 58420 72964 0.801
# 2         1 14544 72964 0.199

# -------------------
# Irrigation
# Soybean 
raster_irrigation_s <- raster::raster("E:/POSTDOC INRAE/DATA/02_YIELDS/SPAM/spam2010v2r0_global_harv_area.geotiff/spam2010V2r0_global_H_SOYB_I.tif"); raster_irrigation_s
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
raster_irrigation_m <- raster::raster("E:/POSTDOC INRAE/DATA/02_YIELDS/SPAM/spam2010v2r0_global_harv_area.geotiff/spam2010V2r0_global_H_MAIZ_I.tif"); raster_irrigation_m
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
# Climate

# > Select the crop (soybean or maize)
#crop <- "soybean"
crop <- "maize"

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
  
  # -------------------------
  # > Compute monthly averages based on daily data 
  # all variables except precipitations
  if(var_i %in% c('max_2m_temperature', 'min_2m_temperature', 
                  'surface_net_solar_radiation', 'et0', 
                  'vapor_pressure_deficit_1'))
  {
    
    # > Load climatic data
    era5_var_i <- loadRDa(paste0("C:/Users/benni/Documents/Post doc/ERA5_daily/", crop ,"/era5daily_data_", var_i_abb, ".rda")) # %>% filter(site_year %in% unique(yield_s$site_year))
    
    # > Compute monthly averages 
    tab_month_annual <- monthly_average(var_i = var_i, 
                                         load = F, 
                                         data_var_i = era5_var_i) %>% 
      group_by(site_year, x, y, gridcode, country_name, clim.var) %>%
      # > correct year 
      #mutate(year = if_else(country_name == "Brazil" & month %in% c(11,12), year+1, year)) %>% 
      # > identify month id
      arrange(year, month) %>%
      mutate(month_nb = row_number()) %>% 
      # to remove 1980 and 2017 site years
      filter(substr(site_year, nchar(site_year)-3, nchar(site_year)) %in% as.character(1981:2016)) %>% 
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
    
    # > Load climatic data
    era5_var_i <- loadRDa(paste0("C:/Users/benni/Documents/Post doc/ERA5_daily/", crop, "/era5monthly_data_prec.rda"))
    
    # > Change format of monthly data 
    tab_month_annual <- era5_var_i %>% 
      group_by(site_year, x, y, gridcode, country_name, clim.var) %>%
      # > identify month id
      arrange(year, month) %>%
      mutate(month_nb = row_number()) %>% 
      # > remove useless column
      ungroup() %>% 
      dplyr::select(-month, -year) %>%
      pivot_wider(names_from = c("clim.var", "month_nb"), 
                  values_from = "clim.value", 
                  names_prefix = "monthly_", 
                  names_sep = "_")
    
  }
  
  # -------------------------
  # > Compute annual averages 
  tab_month_annual$year_average <- tab_month_annual %>% 
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

list_climate$max_temp$year_max_2m_temperature     <- list_climate$max_temp$year_average
list_climate$min_temp$year_min_2m_temperature     <- list_climate$min_temp$year_average
list_climate$rad$year_surface_net_solar_radiation <- list_climate$rad$year_average
list_climate$et0$year_et0                         <- list_climate$et0$year_average
list_climate$prec$year_total_precipitation        <- list_climate$prec$year_average
list_climate$vpd_1$year_vapor_pressure_deficit    <- list_climate$vpd_1$year_average

list_climate %>% map(., ~{ names(.x)})

list_climate %>% map_dfr(., ~{ nrow(.x)}, .id="var")
list_climate %>% map_dfr(., ~{ length(unique(.x$site_year))}, .id="var")

# > Merge all climatic variables in 1 tab
tab_climate <- list_climate$min_temp %>% 
  left_join(., list_climate$max_temp %>% dplyr::select(-year_average), by = c("site_year", "x", "y", "gridcode", "country_name")) %>% 
  left_join(., list_climate$rad      %>% dplyr::select(-year_average), by = c("site_year", "x", "y", "gridcode", "country_name")) %>% 
  left_join(., list_climate$et0      %>% dplyr::select(-year_average), by = c("site_year", "x", "y", "gridcode", "country_name")) %>% 
  left_join(., list_climate$prec     %>% dplyr::select(-year_average), by = c("site_year", "x", "y", "gridcode", "country_name")) %>% 
  left_join(., list_climate$vpd_1    %>% dplyr::select(-year_average), by = c("site_year", "x", "y", "gridcode", "country_name")) 

summary(tab_climate)

#tab_climate_s <- tab_climate %>% 
#  filter(site_year %in% unique(yield_s_corrected$site_year))
#dim(tab_climate_s) # 116516     

tab_climate_m <- tab_climate %>% 
  filter(site_year %in% unique(yield_m_corrected$site_year))
dim(tab_climate_m) # 72964       

# -------------------------------------------------------------------------
# PCA on monthly averages 
# Home-made functions performing the dimension reductions
source("E:/POSTDOC INRAE/ANALYSES/A_MODEL_COMP/00_Functions_dimension_reduction.R")

# Only select the site_years that are in the Yield data 
# > soybean 
tab_climate_pca <- tab_climate_s

dim(tab_climate_pca) # 116516          
summary(tab_climate_pca) # no NA

message("SOYBEAN - PRINCIPAL COMPONENT ANALYSIS")
pca_scores <- function_pca(type_data = "M", 
                           vars_names = vars_names,
                           data = tab_climate_pca, 
                           scale = T, 
                           cum_clim.var = F)

# > save pca object to retrieve the scores
save(pca_scores, file = paste0("C:/Users/benni/Documents/Post doc/ERA5_daily/soybean/pca_soybean.rda"))

# add site-year label
pca_scores$tab_PCA_scores$site_year <- tab_climate_pca$site_year

summary(pca_scores$tab_PCA_scores) # mean 0

pca_scores_s <- pca_scores

rm(tab_climate_pca, pca_scores)

# > maize
tab_climate_pca <- tab_climate_m 

dim(tab_climate_pca) # 72964
summary(tab_climate_pca) # no NA

message("MAIZE - PRINCIPAL COMPONENT ANALYSIS")
pca_scores <- function_pca(type_data = "M", 
                           vars_names = vars_names,
                           data = tab_climate_pca, 
                           scale = T, 
                           cum_clim.var = F)

summary(pca_scores$tab_PCA_scores) # mean = 0 
names(pca_scores$tab_PCA_scores) # only 7 scores for min temp
pca_scores$list_pca_per_variable$min_temp$pca$rotation # 7 PC 

# add site-year label
pca_scores$tab_PCA_scores$site_year <- tab_climate_pca$site_year

# > save pca object to retrieve the scores
save(pca_scores, file = paste0("C:/Users/benni/Documents/Post doc/ERA5_daily/maize/pca_maize.rda"))

summary(pca_scores$tab_PCA_scores) # mean 0

pca_scores_m <- pca_scores

rm(tab_climate_pca, pca_scores)

# -------------------------------------------------------------------------
# Merge to yield and irrigation data 

# > soybean
tab_soybean <- yield_s_corrected %>% 
  mutate(year = as.numeric(as.character(year))) %>% 
  left_join(tab_climate_s, by=c("site_year", "x", "y", "gridcode", "country_name")) %>% 
  left_join(pca_scores_s$tab_PCA_scores, by=c("site_year")) %>% 
  left_join(irrigation_s_tab, by=c("x", "y", "gridcode")) %>% 
  dplyr::select(-to_remove)

summary(tab_soybean)

# > maize
tab_maize <- yield_m_corrected %>% 
  mutate(year = as.numeric(as.character(year))) %>% 
  left_join(tab_climate_m, by=c("site_year", "x", "y", "gridcode", "country_name")) %>% 
  left_join(pca_scores_m$tab_PCA_scores, by=c("site_year")) %>% 
  left_join(irrigation_m_tab, by=c("x", "y", "gridcode")) %>% 
  dplyr::select(-to_remove)

summary(tab_maize)

# -------------------------------------------------------------------------
# Save 
save(tab_soybean, file = "E:/POSTDOC INRAE/ANALYSES/B_OPTIMISATION/00_Data/00_tab_soybean.rda")
save(tab_maize, file = "E:/POSTDOC INRAE/ANALYSES/B_OPTIMISATION/00_Data/00_tab_maize.rda")

