# -------------------------------------------------------------------------
#
#       02-0. ERA5-Land daily daily climatic data 2000-2023 for EU
#       Merge data from individual .nc to global dataset 
#       Author: M. Chen, Inrae, 2024 
#         
# -------------------------------------------------------------------------

# ----------------------------------------
# Packages & tools
library(tidyverse)
library(stringr)
library(lubridate)
library(terra) ; library(rnaturalearth)
library(parallel) ; library(doParallel); library(foreach)
library(CCMHr) 

# Homemade function to read daily climate data from ERA5-land dataset. 
source(".../functions_to_read_era5.R")

# ----------------------------------------
# Data 

# > path to climatic data
path_day <- "..."
path_month <- "..."

# > load 1 initial yield file to resample era5 data 
#   yield file is from the GDHY dataset (accessible here: https://doi.pangaea.de/10.1594/PANGAEA.909132)
yield_ref <- rast(".../GDHY_v1.3/gdhy_v1.2_v1.3_20190128/maize/yield_1981.nc4")

# > coordinates
load(".../data/00_dat_coords_EU42.rda")
dat_coords_EU <- dat_coords_eu42

dat_coords_EU <- dat_coords_EU %>%
  unite(col = "gridcode", c("x", "y"), remove = F, sep = "_") %>%
  dplyr::select(-continent, -region)

# > check if there is some duplicate (normally no)
dat_coords_EU %>% 
  distinct(gridcode, x, y, country_name) %>% 
  group_by(gridcode) %>% 
  summarise(n = n()) %>% 
  filter(n>1) # 0 line: OK

# > set for data loading
dat_coords_EU_temp <- dat_coords_EU
dat_coords_EU <- dat_coords_EU_temp %>% 
  dplyr::select(gridcode, x, y, country_name) %>%
  mutate(country_name = "Italy")

# > count nb of cells per crop (excluding desert)
dim(dat_coords_EU) # 4192 sites

# ----------------------------------------
# Individual .nc files with climatic ERA5 data 

# > extract all the names of the files
filenames <- list.files(path, pattern="*.nc", full.names = TRUE)

# > split the files among the different variables 
filetable <- data.frame(filename = filenames) %>% 
  # > add variable 
  mutate(var = case_when(
    str_detect(filename, "10m_u_component_of_wind")         == T ~ "10m_u_component_of_wind",
    str_detect(filename, "10m_v_component_of_wind")         == T ~ "10m_v_component_of_wind",
    str_detect(filename, "mean_total_precipitation")        == T ~ "mean_precipitation",
    str_detect(filename, "maximum_total_precipitation")     == T ~ "total_precipitation",
    str_detect(filename, "mean_2m_temperature")             == T ~ "2m_temperature",
    str_detect(filename, "mean_2m_dewpoint_temperature")    == T ~ "2m_dewpoint_temperature",
    str_detect(filename, "minimum_2m_temperature")          == T ~ "min_2m_temperature",
    str_detect(filename, "maximum_2m_temperature")          == T ~ "max_2m_temperature",
    str_detect(filename, "minimum_2m_dewpoint_temperature") == T ~ "min_2m_dewpoint_temperature",
    str_detect(filename, "maximum_2m_dewpoint_temperature") == T ~ "max_2m_dewpoint_temperature",
    str_detect(filename, "surface_pressure")                == T ~ "surface_pressure",
    str_detect(filename, "surface_net_solar_radiation")     == T ~ "surface_net_solar_radiation"
  )) %>% 
  # > add month and year 
  mutate(year  = substr(substr(filename, nchar(filename)-10, nchar(filename)), 2, 5),
         month = substr(substr(filename, nchar(filename)-10, nchar(filename)), 7, 8))

# > examine data 
filetable %>% 
  group_by(var, year) %>%
  summarise(n_files=n()) %>% 
  ggplot(., aes(x=as.numeric(as.character(year)), y=var, fill=as.factor(n_files))) +
  geom_tile(colour="white") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        panel.grid = element_blank()) +
  scale_fill_viridis_d(direction = -1, name="Number of months available") + 
  labs(x="Years")

# For the moment, daily data are split between months
# 1 file = daily data for each pixel / month / year
# Temporal range: 1980-2017
# Spatial coverage: global, 0.5° resolution

# ----------------------------------------
# PREPARATION FOR DATA LOADING
# > name variables and abrevations
vars_names <- data.frame(clim.var = c("max_2m_temperature", "min_2m_temperature", 
                                      "et0", "surface_net_solar_radiation", 
                                      "total_precipitation", "vapor_pressure_deficit")) %>% 
  mutate(clim.var_abb = recode(clim.var, 
                               "min_2m_temperature"         ="min_temp",
                               "max_2m_temperature"         ="max_temp",
                               "et0"                        ="et0",
                               "surface_net_solar_radiation"="rad",
                               "total_precipitation"        ="prec",
                               "vapor_pressure_deficit"     ="vpd_1"))  %>% 
  mutate(clim.var_lab = recode(clim.var, 
                               "min_2m_temperature"         ="Minimum temperature",
                               "max_2m_temperature"         ="Maximum temperature",
                               "et0"                        ="Evapotranspiration ref",
                               "surface_net_solar_radiation"="Solar radiations",
                               "total_precipitation"        ="Precipitation",
                               "vapor_pressure_deficit"   ="Vapor pressure deficit"))

# > ERA5 variables to compute VPD, ET0
var_vpd_1 <- c("min_2m_temperature", "max_2m_temperature", "min_2m_dewpoint_temperature", "max_2m_dewpoint_temperature")

var_et0 <- c("10m_u_component_of_wind", "10m_v_component_of_wind", "min_2m_temperature", "max_2m_temperature", 
             "2m_dewpoint_temperature", "surface_net_solar_radiation", "surface_pressure")

# > select files to merge for each variable 
files_to_merge <- list()
files_to_merge[[paste0("max_temp")]] <- filetable %>% mutate(to_keep = case_when(#var == "max_2m_temperature"          & year == 2000 ~ 1, 
                                                                                 var == "max_2m_temperature"          & year %in% 2000:2023 ~ 1, 
                                                                                 TRUE ~ 0)) %>% filter(to_keep == 1) 
files_to_merge[[paste0("min_temp")]] <- filetable %>% mutate(to_keep = case_when(#var == "min_2m_temperature"          & year == 2000 ~ 1, 
                                                                                 var == "min_2m_temperature"          & year %in% 2000:2023 ~ 1, 
                                                                                 TRUE ~ 0)) %>% filter(to_keep == 1) 
files_to_merge[[paste0("rad")]]      <- filetable %>% mutate(to_keep = case_when(#var == "surface_net_solar_radiation" & year == 2000 ~ 1, 
                                                                                 var == "surface_net_solar_radiation" & year %in% 2000:2023 ~ 1, 
                                                                                 TRUE ~ 0)) %>% filter(to_keep == 1) 
files_to_merge[[paste0("et0")]]      <- filetable %>% mutate(to_keep = case_when(#var %in% var_et0                     & year == 2000 ~ 1, 
                                                                                 var %in% var_et0                     & year %in% 2000:2023 ~ 1, 
                                                                                 TRUE ~ 0)) %>% filter(to_keep == 1) 
files_to_merge[[paste0("vpd_1")]]    <- filetable %>% mutate(to_keep = case_when(#var %in% var_vpd_1                   & year == 2000 ~ 1, 
                                                                                 var %in% var_vpd_1                   & year %in% 2000:2023 ~ 1, 
                                                                                 TRUE ~ 0)) %>% filter(to_keep == 1) 


# Load, merge and recompute variable for the set of gridcells used for test
for(v in c("max_temp", "min_temp", "rad"))
{
  
  # > Load data
  era5daily_init <- merge_era5_data(var = v,
                                    crop = "Maize",
                                    files_to_merge = files_to_merge[[paste0(v)]]$filename, 
                                    dat.coords = dat_coords_EU,
                                    yield_ref = yield_ref, 
                                    save_output = F)
  
  # > (Re)Compute the variables
  era5daily_correct<- correct_era5_data(clim.var = v, 
                                        data.clim.var = era5daily_init, 
                                        cum.value = T)
  
  # > save
  save(era5daily_correct, 
       file = paste0(path_daily, "/era5daily_", v, "_EU.rda"))
  
  # > remove
  rm(era5daily_init, era5daily_correct)
  
}

rm(era5daily_init, era5daily_correct)

# -----------
# VPD
# vpd_1 takes too much time, need to split loading

v <- "vpd_1"
for(y in c(2000:2023))
{
  
  # > files to load
  files_y <- files_to_merge[[paste0(v)]][which(files_to_merge[[paste0(v)]]$year==y),]$filename
  # > Load data
  era5daily_init <- merge_era5_data(var = v,
                                    crop = "Maize",
                                    files_to_merge = files_y, 
                                    dat.coords = dat_coords_EU,
                                    yield_ref = yield_ref, 
                                    save_output = F)
  
  # > (Re)Compute the variables
  era5daily_correct<- correct_era5_data(clim.var = v, 
                                        data.clim.var = era5daily_init, 
                                        cum.value = T)
  
  
  # > save
  save(era5daily_correct, 
       file = paste0(path_daily, "/era5daily_", v, "_", y, "_EU.rda"))
  
  # > remove unused files
  rm(files_y, era5daily_init, era5daily_correct)
  
}

# Merge vpd files together
list_vpd_1_temp <- list()

# load & merge 
for(y in 2000:2023)
{
  
  # > Load data
  era5daily_correct_y <- loadRDa(paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/01_days/temp_eu/vpd_1_temp/era5daily_vpd_1_", y, "_EU.rda"))
  
  # > Store in a list
  list_vpd_1_temp[[paste0(y)]] <- era5daily_correct_y
  
  # > remove unused files
  rm(era5daily_correct_y)
  
}

era5daily_correct <- map_dfr(list_vpd_1_temp, data.frame)

# checks
unique(era5daily_correct$clim.var) # "vapor_pressure_deficit"

length(unique(era5daily_correct$site_year))

era5daily_correct %>% 
  mutate(year=year(date)) %>% 
  pull(year) %>% unique(.) # check if all year are here

era5daily_correct %>% 
  mutate(year=year(date)) %>% 
  group_by(year, site_year) %>% 
  count() %>% 
  summary() # 244 lines per site-year per year 

save(era5daily_correct, 
     file = paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/01_days/temp_eu/era5daily_vpd_1_EU.rda"))

rm(era5daily_init, era5daily_correct)

# -----------
# ET0
# et0 takes too much time, need to split loading

v <- "et0"

for(y in c(2000:2023))
{
  
  # > files to load
  files_y <- files_to_merge[[paste0(v)]][which(files_to_merge[[paste0(v)]]$year==y),]$filename
  # > Load data
  era5daily_init <- merge_era5_data(var = v,
                                    crop = "Maize",
                                    files_to_merge = files_y, 
                                    dat.coords = dat_coords_EU,
                                    yield_ref = yield_ref, 
                                    save_output = F)
  
  # > (Re)Compute the variables
  era5daily_correct<- correct_era5_data(clim.var = v, 
                                        data.clim.var = era5daily_init, 
                                        cum.value = T)
  
  
  # > save
  save(era5daily_correct, 
       file = paste0(path_daily, "/", v, "_temp/era5daily_", v, "_", y, "_EU.rda"))
  
  # > remove unused files
  rm(files_y, era5daily_init, era5daily_correct)
  
}

# Merge et0 files together
list_et0_temp <- list()

# load & merge 
for(y in 2000:2023)
{
  
  # > Load data
  era5daily_correct_y <- loadRDa(paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/01_days/temp_eu/et0_temp/era5daily_et0_", y, "_EU.rda"))
  
  # > Store in a list
  list_et0_temp[[paste0(y)]] <- era5daily_correct_y
  
  # > remove unused files
  rm(era5daily_correct_y)
  
}

era5daily_correct <- map_dfr(list_et0_temp, data.frame)

era5daily_correct %>% 
  group_by(site_year) %>% 
  count() %>% 
  summary() # 244 lines per site-year

save(era5daily_correct, 
     file = paste0(path_daily, "/era5daily_et0_EU.rda"))

# -----------
# TOTAL PRECIPITATIONS 
# (directly downloaded from ERA5 Land portal)

# > extract all the names of the files
filenames_month <- list.files(path_month, pattern="*.nc", full.names = TRUE)
# > split the files among the different variables 
filetable_month <- data.frame(filename = filenames_month) %>% 
  # > add variable 
  mutate(var = case_when(
    str_detect(filename, "10m_u_component_of_wind")         == T ~ "10m_u_component_of_wind",
    str_detect(filename, "10m_v_component_of_wind")         == T ~ "10m_v_component_of_wind",
    str_detect(filename, "mean_total_precipitation")        == T ~ "mean_precipitation",
    str_detect(filename, "total_precipitation")             == T ~ "total_precipitation",
    str_detect(filename, "mean_2m_temperature")             == T ~ "2m_temperature",
    str_detect(filename, "mean_2m_dewpoint_temperature")    == T ~ "2m_dewpoint_temperature",
    str_detect(filename, "minimum_2m_temperature")          == T ~ "min_2m_temperature",
    str_detect(filename, "maximum_2m_temperature")          == T ~ "max_2m_temperature",
    str_detect(filename, "minimum_2m_dewpoint_temperature") == T ~ "min_2m_dewpoint_temperature",
    str_detect(filename, "maximum_2m_dewpoint_temperature") == T ~ "max_2m_dewpoint_temperature",
    str_detect(filename, "surface_pressure")                == T ~ "surface_pressure",
    str_detect(filename, "surface_net_solar_radiation")     == T ~ "surface_net_solar_radiation"
  )) %>% 
  # > add month and year 
  mutate(year  = substr(filename, nchar(filename)-6, nchar(filename)-3))


# > 2000 - 2022 
files_to_merge[[paste0("prec")]]     <- filetable_month %>% mutate(to_keep = case_when(var == "total_precipitation" & year %in% 2000:2022 ~ 1, 
                                                                                       TRUE ~ 0)) %>% filter(to_keep == 1) 

# > Load data
era5monthly_init_2000_2022 <- merge_era5_data(var = "prec",
                                              crop = "Maize",
                                              files_to_merge = files_to_merge[[paste0("prec")]], 
                                              dat.coords = dat_coords_EU,
                                              yield_ref = yield_ref, 
                                              save_output = F, 
                                              monthly = TRUE)

# > 2023 (not clean data)
# > export the data from .nc file
raster_i_init <- rast(paste0(path_month, "/download_monthly_total_precipitation_2023.nc")) ; raster_i_init
raster_i <- terra::aggregate(raster_i_init, fact=2, fun="mean") ; raster_i

# > add projection
crs(raster_i) <- "epsg:4326"

# > realign era5 data on yield data
raster_i_resample <- resample(raster_i, yield_ref)

# > transform it into a data.frame
tab_raster_i <- as.data.frame(raster_i_resample, xy=T) %>% 
  mutate(x = if_else(x>180, x-360, x)) %>% 
  right_join(., dat_coords_EU, by=c("x", "y")) %>% 
  # rename columns
  rename("tp_01"="tp_expver=1_1", 
         "tp_02"="tp_expver=1_2", 
         "tp_03"="tp_expver=1_3", 
         "tp_04"="tp_expver=1_4", 
         "tp_05"="tp_expver=1_5", 
         "tp_06"="tp_expver=1_6", 
         "tp_07"="tp_expver=1_7",
         "tp_08"="tp_expver=1_8", 
         "tp_09"="tp_expver=1_9", 
         "tp_10"="tp_expver=1_10",
         "tp_11"="tp_expver=5_11") %>% 
  # remove weird columns data
  dplyr::select(-starts_with("tp_expver"))

# > long format data.frame 
crop = "Maize"
era5monthly_init_2023 <- tab_raster_i %>% 
  # > set in long format
  gather(key = variable, value = clim.value, -x, -y, -gridcode, -country_name) %>% 
  # > add real variable and month number
  separate(variable, c("clim.var", "month"), remove = F) %>% 
  # > rename clim.var
  mutate(clim.var = "total_precipitation") %>% 
  # > add year
  mutate(year  = 2023) %>% 
  mutate(year = if_else(crop == "Soybean" & country_name %in% c("Argentina", "Brazil") & month %in% c("11", "12"), year+1, year)) %>% 
  # > separate by year
  split(.$year) %>% 
  # > identify months to keep for each region (soybean/maize growing season)
  map_dfr(., ~{ 
    
    .x %>% 
      mutate(to_keep = 0) %>%
      mutate(to_keep = case_when(
        # Soybean 
        crop == "Soybean" & country_name %in% c("Argentina", "Brazil")                                  & month %in% c("11", "12", "01", "02", "03", "04", "05")  ~ 1,  
        crop == "Soybean" & country_name %in% c("Desert", "China", "United States of America", "Italy") & month %in% c("04", "05", "06", "07", "08", "09", "10")  ~ 1, 
        crop == "Soybean" & country_name %in% c("Canada")                                               & month %in% c("05", "06", "07", "08", "09", "10", "11")  ~ 1, 
        crop == "Soybean" & country_name %in% c("India")                                                & month %in% c("06", "07", "08", "09", "10", "11", "12")  ~ 1,
        # Maize 
        crop == "Maize"   & month %in% c("04", "05", "06", "07", "08", "09", "10", "11") ~ 1,
        TRUE ~ 0
      )) %>% 
      filter(to_keep == 1)
    
  }, .id = "year") %>%
  # > create site*year combination as unique ID
  mutate(site_year = paste0(gridcode, "_", year)) %>% 
  # > remove useless variables
  dplyr::select(site_year, x, y, gridcode, country_name, month, year, clim.var, clim.value)

# > Merge 2000-2022 and 2023 data 
era5monthly_init <- rbind(era5monthly_init_2000_2022,era5monthly_init_2023) 

dim(era5monthly_init) # 673728                 

era5monthly_init %>% 
  group_by(year, month) %>% 
  count() %>% 
  summary(.) # 3509 lines : ok

# > (Re)Compute the variables
era5daily_correct <- era5monthly_init %>% 
  mutate(clim.value = clim.value*1e3)

# > save
save(era5daily_correct, 
     file = paste0(path_daily, "/era5daily_prec_EU.rda"))

rm(era5monthly_init, era5daily_correct)

