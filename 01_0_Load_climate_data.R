# -------------------------------------------------------------------------
#
#       01-0. ERA5-Land daily daily climatic data 1980-2017
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

# Homemade function
source("E:/POSTDOC INRAE/DATA/01_CLIMATE/ERA5/functions_to_read_era5.R")

# ----------------------------------------
# Data 

# > path to climatic data
path <- "C:/Users/benni/Documents/Post doc/Test"

# > load 1 initial yield file to resample era5 data 
yield_ref <- rast("E:/POSTDOC INRAE/DATA/02_YIELDS/GDHY_v1.3/gdhy_v1.2_v1.3_20190128/maize/yield_1981.nc4")

# > yield data to retrieve coordinates
load("E:/POSTDOC INRAE/DATA/02_YIELDS/GDHY_v1.3/yield_no.trend_full.rda")

# > check if there is some duplicate (normally no)
data_crop %>% 
  distinct(gridcode, x, y, country_name) %>% 
  group_by(gridcode) %>% 
  summarise(n = n()) %>% 
  filter(n>1) # 0 line: OK

# > count nb of cells per crop (excluding desert)
data_crop %>% 
  filter(country_name !="Desert") %>% 
  distinct(gridcode, x, y, crop) %>%
  group_by(crop) %>% 
  summarise(n_gridcells = n())
#    crop    n_gridcells
# 1 Maize          730
# 2 Soybean        2784

# > count nb of cells per country (+ desert) per crop
# > Soybean
data_crop %>% 
  filter(crop != "maize") %>% 
  distinct(gridcode, x, y, country_name) %>%
  group_by(country_name) %>% 
  summarise(n_gridcells = n()) %>% 
  mutate(tot_gridcells = sum(n_gridcells)) %>% 
  mutate(freq_gridcells = (n_gridcells/tot_gridcells)*100) %>%
  mutate(labs=paste0(n_gridcells, " (", round(freq_gridcells, 1), "%)")) %>% 
  dplyr::select(country_name, labs)

#  country_name             Soybean    
#1 Argentina                283 (8.2%) 
#2 Brazil                   424 (12.3%)
#3 Canada                   28 (0.8%)  
#4 China                    994 (28.8%)
#5 India                    199 (5.8%) 
#6 Italy                    26 (0.8%)  
#7 United States of America 830 (24.1%)
#8 Desert                   663 (19.2%)

# > Maize
data_crop %>% 
  filter(crop == "maize") %>% 
  distinct(gridcode, x, y, country_name) %>%
  group_by(country_name) %>% 
  summarise(n_gridcells = n()) %>% 
  mutate(tot_gridcells = sum(n_gridcells)) %>% 
  mutate(freq_gridcells = (n_gridcells/tot_gridcells)*100) %>%
  arrange(desc(freq_gridcells)) %>% 
  mutate(labs=paste0(n_gridcells, " (", round(freq_gridcells, 1), "%)")) %>% 
  dplyr::select(country_name, labs)

# Top 10
#   country_name              labs       
#  1 United States of America 984 (46%)  
#  2 Desert                   405 (18.9%)
#  3 China                    192 (9%)   
#  4 France                   148 (6.9%) 
#  5 Spain                    65 (3%)    
#  6 Italy                    52 (2.4%)  
#  7 Hungary                  44 (2.1%)  
#  8 Canada                   32 (1.5%)  
#  9 Republic of Serbia       31 (1.4%)  
# 10 Turkey                   31 (1.4%)  

# > retrieve coordinates from grid-cells
dat.coords_soybean <- data_crop %>% 
  filter(crop != "maize") %>% 
  distinct(gridcode, x, y, country_name)

dat.coords_maize <- data_crop %>% 
  filter(crop == "maize") %>% 
  distinct(gridcode, x, y, country_name)

# > number of site-years
dat.site_year_soybean <- data_crop %>% 
  filter(crop != "maize") %>% 
  distinct(gridcode, x, y, year, country_name, Ya) %>% 
  unite("site_year", gridcode, year, remove = F)

length(unique(dat.site_year_soybean$site_year)) # n=122121
dat.site_year_soybean %>% group_by(site_year) %>% count() %>% filter(n!=1) # 0 line: OK

dat.site_year_maize <- data_crop %>% 
  filter(crop == "maize") %>% 
  distinct(gridcode, x, y, year, country_name, Ya) %>% 
  unite("site_year", gridcode, year, remove = F)

length(unique(dat.site_year_maize$gridcode)) # n=2139
length(unique(dat.site_year_maize$site_year)) # n=75562
dat.site_year_maize %>% group_by(site_year) %>% count() %>% filter(n!=1) # 0 line: OK

#dat.coords_soybean_test <- dat.coords_soybean[1:10,]

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
# FILES TO LOAD
# > ERA5 variables to compute VPD, ET0
var_vpd_1 <- c("min_2m_temperature", "max_2m_temperature", "min_2m_dewpoint_temperature", "max_2m_dewpoint_temperature")

var_et0 <- c("10m_u_component_of_wind", "10m_v_component_of_wind", "min_2m_temperature", "max_2m_temperature", 
             "2m_dewpoint_temperature", "surface_net_solar_radiation", "surface_pressure")

# > select files to merge for each variable 
files_to_merge <- list()

files_to_merge[[paste0("max_temp")]] <- filetable %>% mutate(to_keep = case_when(var %in% "max_2m_temperature" & year == 1980 & month %in% c("11", "12") ~ 1, 
                                                                                 var %in% "max_2m_temperature" & year %in% 1981:2016 ~ 1, 
                                                                                 TRUE ~ 0)) %>%   filter(to_keep == 1)  
files_to_merge[[paste0("min_temp")]] <- filetable %>% mutate(to_keep = case_when(var %in% "min_2m_temperature" & year == 1980 & month %in% c("11", "12") ~ 1, 
                                                                                 var %in% "min_2m_temperature" & year %in% 1981:2016 ~ 1, 
                                                                                 TRUE ~ 0)) %>%   filter(to_keep == 1) 
files_to_merge[[paste0("rad")]]      <- filetable %>% mutate(to_keep = case_when(var %in% "surface_net_solar_radiation" & year == 1980 & month %in% c("11", "12") ~ 1, 
                                                                                 var %in% "surface_net_solar_radiation" & year %in% 1981:2016 ~ 1, 
                                                                                 TRUE ~ 0)) %>%   filter(to_keep == 1) 
files_to_merge[[paste0("et0")]]      <- filetable %>% mutate(to_keep = case_when(var %in% var_et0 & year == 1980 & month %in% c("11", "12") ~ 1, 
                                                                                 var %in% var_et0 & year %in% 1981:2016 ~ 1, 
                                                                                 TRUE ~ 0)) %>%   filter(to_keep == 1)  
files_to_merge[[paste0("vpd_1")]]    <- filetable %>% mutate(to_keep = case_when(var %in% var_vpd_1 & year == 1980 & month %in% c("11", "12") ~ 1, 
                                                                                 var %in% var_vpd_1 & year %in% 1981:2016 ~ 1, 
                                                                                 TRUE ~ 0)) %>%   filter(to_keep == 1) 

# ----------------------------------------
# SOLAR RADIATION, MINIMUM AND MAXIMUM TEMPERATURE (°C)

# > check the selected files 
plot_grid(# > Temperature max
          files_to_merge[[paste0("max_temp")]] %>% group_by(var, year) %>% count() %>% 
            ggplot(., aes(x=as.numeric(as.character(year)), y=n, color=var)) + geom_line() + geom_point() + theme_bw() + theme(legend.position = "none") + facet_wrap(.~var),
          # > Temperature min
          files_to_merge[[paste0("min_temp")]] %>% group_by(var, year) %>% count() %>% 
            ggplot(., aes(x=as.numeric(as.character(year)), y=n, color=var)) + geom_line() + geom_point() + theme_bw() + theme(legend.position = "none") + facet_wrap(.~var),
          # > Radiation
          files_to_merge[[paste0("rad")]] %>% group_by(var, year) %>% count() %>% 
            ggplot(., aes(x=as.numeric(as.character(year)), y=n, color=var)) + geom_line() + geom_point() + theme_bw() + theme(legend.position = "none") + facet_wrap(.~var), 
          ncol=3)

# > Soybean
# Load, merge and recompute variable for the set of gridcells used for test
for(v in c("max_temp", "min_temp", "rad"))
{
  
  # > Load data
  era5daily_init <- merge_era5_data(var = v,
                                    crop = "Soybean",
                                    files_to_merge = files_to_merge[[paste0(v)]]$filename, 
                                    dat.coords = dat.coords_soybean,
                                    yield_ref = yield_ref, 
                                    save_output = F)
  
  # > (Re)Compute the variables
  era5daily_correct<- correct_era5_data(clim.var = v, 
                                        data.clim.var = era5daily_init, 
                                        cum.value = T)
  
  # > save
  save(era5daily_correct, 
       file = paste0("C:/Users/benni/Documents/Post doc/ERA5_daily/soybean/era5daily_", v, ".rda"))
  
  # > remove
  rm(era5daily_init, era5daily_correct)
  
  
}

# > Maize
# Load, merge and recompute variable for the set of gridcells used for test
for(v in c("max_temp", "min_temp", "rad"))
{
  
  # > Load data
  era5daily_init <- merge_era5_data(var = v,
                                    crop = "Maize",
                                    files_to_merge = files_to_merge[[paste0(v)]][which(files_to_merge[[paste0(v)]]$year!=1980),]$filename, 
                                    dat.coords = dat.coords_maize,
                                    yield_ref = yield_ref, 
                                    save_output = F)
  
  # > (Re)Compute the variables
  era5daily_correct<- correct_era5_data(clim.var = v, 
                                        data.clim.var = era5daily_init, 
                                        cum.value = T)
  
  # > save
  save(era5daily_correct, 
       file = paste0("C:/Users/benni/Documents/Post doc/ERA5_daily/maize/era5daily_", v, ".rda"))
  
  # > remove
  rm(era5daily_init, era5daily_correct)
  
}

rm(v)

# ----------------------------------------
# VAPOUR PRESSURE DEFICIT (vpd)

# variable name
v <- "vpd_1"

# > check selected files 
files_to_merge[[paste0(v)]] %>% group_by(var, year) %>% count() %>% 
  ggplot(., aes(x=as.numeric(as.character(year)), y=n, color=var)) + geom_line() + geom_point() + theme_bw() + theme(legend.position = "none") + facet_wrap(.~var)

# > Soybean
for(y in c(1981:2016))
{
  
  # > files to load
  files_y <- files_to_merge[[paste0(v)]] %>% 
    # > identify the data needed
    mutate(to_keep = case_when(
      year == y-1 & month %in% c("11", "12") ~ 1,
      year %in% y ~ 1,
      TRUE ~ 0)) %>% 
    filter(to_keep == 1)
  
  # > Load data
  era5daily_init <- merge_era5_data(var = v,
                                    crop = "Soybean",
                                    files_to_merge = files_y$filename, 
                                    dat.coords = dat.coords_soybean,
                                    yield_ref = yield_ref, 
                                    save_output = F)
  
  # > (Re)Compute the variables
  era5daily_correct<- correct_era5_data(clim.var = v, 
                                        data.clim.var = era5daily_init, 
                                        cum.value = T)
  
  
  # > save
  save(era5daily_correct, 
       file = paste0("C:/Users/benni/Documents/Post doc/ERA5_daily/soybean/vpd_1/era5daily_data_vpd_1_", y, ".rda"))
  
  # > remove unused files
  rm(files_y, era5daily_init, era5daily_correct)
  
}

# > merge all files
list_vpd_1_temp_s <- list()

# > load & merge 
for(y in 1981:2016)
{
  
  # > Load data
  era5daily_correct_y <- loadRDa(paste0("C:/Users/benni/Documents/Post doc/ERA5_daily/soybean/vpd_1/era5daily_data_vpd_1_", y, ".rda"))
  
  # > Store in a list
  list_vpd_1_temp_s[[paste0(y)]] <- era5daily_correct_y
  
  # > remove unused files
  rm(era5daily_correct_y)
  
}

# checks if the number of lines is similar across years (different number of lines between bissextiles and normal years)
# also check if the data are not similar between years
list_vpd_1_temp_s %>% 
  map_dfr(., ~{
    data.frame(n_lines = nrow(.x), 
               mean_clim.value = mean(.x$clim.value)) }, .id="year")

# merge
era5daily_data_vpd_1_s <- map_dfr(list_vpd_1_temp_s, data.frame)

# > save 
save(era5daily_data_vpd_1_s, file = "C:/Users/benni/Documents/Post doc/ERA5_daily/soybean/era5daily_data_vpd_1.rda")

rm(list_vpd_1_temp_s, era5daily_data_vpd_1_s)

# ----
# > Maize 

for(y in c(1981:2016))
{
  
  # > files to load
  files_y <- files_to_merge[[paste0(v)]] %>% 
    # > identify the data needed
    mutate(to_keep = case_when(
      year %in% y ~ 1,
      TRUE ~ 0)) %>% 
    filter(to_keep == 1)
  
  # > Load data
  era5daily_init <- merge_era5_data(var = v,
                                    crop = "Maize",
                                    files_to_merge = files_y$filename, 
                                    dat.coords = dat.coords_maize %>% filter(gridcode=="179.75_65.25"),
                                    yield_ref = yield_ref, 
                                    save_output = F)
  
  # > (Re)Compute the variables
  era5daily_correct<- correct_era5_data(clim.var = v, 
                                        data.clim.var = era5daily_init, 
                                        cum.value = T)
  
  
  # > save
  save(era5daily_correct, 
       file = paste0("C:/Users/benni/Documents/Post doc/ERA5_daily/maize/vpd_1/era5daily_data_vpd_1_", y, ".rda"))
  
  # > remove unused files
  rm(files_y, era5daily_init, era5daily_correct)
  
}

# > merge all files
list_vpd_1_temp_m <- list()

# > load & merge 
for(y in 1981:2016)
{
  
  # > Load data
  era5daily_correct_y <- loadRDa(paste0("C:/Users/benni/Documents/Post doc/ERA5_daily/maize/vpd_1/era5daily_data_vpd_1_", y, ".rda"))
  
  # > Store in a list
  list_vpd_1_temp_m[[paste0(y)]] <- era5daily_correct_y
  
  # > remove unused files
  rm(era5daily_correct_y)
  
}

# checks if the number of lines is similar across years (different number of lines between bissextiles and normal years)
list_vpd_1_temp_m %>% 
  map_dfr(., ~{
    data.frame(n_lines = nrow(.x), 
               mean_clim.value = mean(.x$clim.value)) }, .id="year")

# merge
era5daily_data_vpd_1_m <- map_dfr(list_vpd_1_temp_m, data.frame)

# > save 
save(era5daily_data_vpd_1_m, file = "C:/Users/benni/Documents/Post doc/ERA5_daily/maize/era5daily_data_vpd_1.rda")

rm(list_vpd_1_temp_m, era5daily_data_vpd_1_m)

rm(v)

# ----------------------------------------
# REFERENCE EVAPOTRANSPIRATION (ET0, in mm.day-1) 

# more information on: https://www.fao.org/3/x0490e/x0490e06.htm

# variable name
v <- "et0"

# > check the selected files 
files_to_merge[[paste0(v)]] %>% group_by(var, year) %>% count() %>% 
  ggplot(., aes(x=as.numeric(as.character(year)), y=n, color=var)) + geom_line() + geom_point() + theme_bw() + theme(legend.position = "none") + facet_wrap(.~var)
# -> some combinations of years and variables have 14 files because we load files for months 11 and 12 twice
# -> these files are not counted twice in the final output

# ----
# > Soybean
for(y in c(1981:2016))
{
  
  # > files to load
  files_y <- files_to_merge[[paste0(v)]] %>% 
    # > identify the data needed
    mutate(to_keep = case_when(
      year == y-1 & month %in% c("11", "12") ~ 1,
      year %in% y ~ 1,
      TRUE ~ 0)) %>% 
    filter(to_keep == 1)
  
  # > Load data
  era5daily_init <- merge_era5_data(var = v,
                                    crop = "Soybean",
                                    files_to_merge = files_y$filename, 
                                    dat.coords = dat.coords_soybean,
                                    yield_ref = yield_ref, 
                                    save_output = F)
  
  # > (Re)Compute the variables
  era5daily_correct<- correct_era5_data(clim.var = v, 
                                        data.clim.var = era5daily_init, 
                                        cum.value = T)
  
  
  # > save
  save(era5daily_correct, 
       file = paste0("C:/Users/benni/Documents/Post doc/ERA5_daily/soybean/et0/era5daily_data_", v, "_", y, ".rda"))
  
  # > remove unused files
  rm(files_y, era5daily_init, era5daily_correct)
  
}

# > merge all files
list_et0_temp_s <- list()

# > load & merge 
for(y in 1981:2016)
{
  
  # > Load data
  era5daily_correct_y <- loadRDa(paste0("C:/Users/benni/Documents/Post doc/ERA5_daily/soybean/et0/era5daily_data_et0_", y, ".rda"))
  
  # > Store in a list
  list_et0_temp_s[[paste0(y)]] <- era5daily_correct_y
  
  # > remove unused files
  rm(era5daily_correct_y)
  
}

# checks if the number of lines is similar across years (different number of lines between bissextiles and normal years)
# also check if the data are not similar between years
list_et0_temp_s %>% 
  map_dfr(., ~{
    data.frame(n_lines = nrow(.x), 
               mean_clim.value = mean(.x$clim.value)) }, .id="year")

# merge
era5daily_data_et0_s <- map_dfr(list_et0_temp_s, data.frame)

# > save
save(era5daily_data_et0_s, file = "C:/Users/benni/Documents/Post doc/ERA5_daily/soybean/era5daily_data_et0.rda")

rm(list_et0_temp_s, era5daily_data_et0_s)


# ----
# > Maize

for(y in c(1981:2016))
{
  
  # > files to load
  files_y <- files_to_merge[[paste0(v)]] %>% 
    # > identify the data needed
    mutate(to_keep = case_when(
      year %in% y ~ 1,
      TRUE ~ 0)) %>% 
    filter(to_keep == 1)
  
  # > Load data
  era5daily_init <- merge_era5_data(var = v,
                                    crop = "Maize",
                                    files_to_merge = files_y$filename, 
                                    dat.coords = dat.coords_maize,
                                    yield_ref = yield_ref, 
                                    save_output = F)
  
  # > (Re)Compute the variables
  era5daily_correct<- correct_era5_data(clim.var = v, 
                                        data.clim.var = era5daily_init, 
                                        cum.value = T)
  
  
  # > save
  save(era5daily_correct, 
       file = paste0("C:/Users/benni/Documents/Post doc/ERA5_daily/maize/et0/era5daily_data_", v, "_", y, ".rda"))
  
  # > remove unused files
  rm(files_y, era5daily_init, era5daily_correct)
  
}

# > merge all files
list_et0_temp_m <- list()

# load & merge 
for(y in 1981:2016)
{
  
  # > Load data
  era5daily_correct_y <- loadRDa(paste0("C:/Users/benni/Documents/Post doc/ERA5_daily/maize/et0/era5daily_data_et0_", y, ".rda"))
  
  # > Store in a list
  list_et0_temp_m[[paste0(y)]] <- era5daily_correct_y
  
  # > remove unused files
  rm(era5daily_correct_y)
  
}

# checks if the number of lines is similar across years (different number of lines between bissextiles and normal years)
# also check if the data are not similar between years
list_et0_temp_m %>% 
  map_dfr(., ~{
    data.frame(n_lines = nrow(.x), 
               mean_clim.value = mean(.x$clim.value, na.rm=T)) }, .id="year")

# merge
era5daily_data_et0_m <- map_dfr(list_et0_temp_m, data.frame)

era5daily_data_et0_m %>%
  filter(is.na(clim.value)==T) %>% 
  distinct(gridcode)

# > save 
save(era5daily_data_et0_m, file = "C:/Users/benni/Documents/Post doc/ERA5_daily/maize/era5daily_data_et0.rda")

rm(list_et0_temp_m, era5daily_data_et0_m)

rm(v)

# ----------------------------------------
# DAILY TOTAL PRECIPITATION averaged by MONTH 

# > extract all the names of the files
filenames_month <- list.files("C:/Users/benni/Documents/Post doc/Test_month", pattern="*.nc", full.names = TRUE)
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

# > 1980-2016
files_to_merge_total_precipitation <- filetable_month %>% mutate(to_keep = case_when(var == "total_precipitation" & year %in% 1980:2016 ~ 1, 
                                                                                     TRUE ~ 0)) %>% filter(to_keep == 1) 

# > Soybean
era5monthly_prec_s <- merge_era5_data(var = "prec",
                                      crop = "Soybean",
                                      files_to_merge = files_to_merge_total_precipitation, 
                                      dat.coords = dat.coords_soybean,
                                      yield_ref = yield_ref, 
                                      save_output = F, 
                                      monthly = TRUE)  %>% 
  # to remove 1980 and 2017 site years
  filter(substr(site_year, nchar(site_year)-3, nchar(site_year)) %in% as.character(1981:2016)) %>% 
  # recompute the variable in mm (instead of m)
  mutate(clim.value = clim.value*1e3)

# > check
era5monthly_prec_s %>% 
  group_by(site_year) %>% 
  count() %>% 
  summary() # only 7 lines per site-year

# > save
save(era5monthly_prec_s, file = paste0("C:/Users/benni/Documents/Post doc/ERA5_daily/soybean/era5monthly_data_prec.rda"))

rm(era5monthly_prec_s)

# ----
# > Maize
era5monthly_prec_m <- merge_era5_data(var = "prec",
                                      crop = "Maize",
                                      files_to_merge = files_to_merge_total_precipitation, 
                                      dat.coords = dat.coords_maize,
                                      yield_ref = yield_ref, 
                                      save_output = F, 
                                      monthly = TRUE) %>% 
  # to remove 1980 and 2017 site years
  filter(substr(site_year, nchar(site_year)-3, nchar(site_year)) %in% as.character(1981:2016)) %>% 
  # recompute the variable in mm (instead of m)
  mutate(clim.value = clim.value*1e3)

# > check
era5monthly_prec_m %>% 
  group_by(site_year) %>% 
  count() %>% 
  summary() # only 8 lines per site-year

# > save
save(era5monthly_prec_m, file = paste0("C:/Users/benni/Documents/Post doc/ERA5_daily/maize/era5monthly_data_prec.rda"))

rm(era5monthly_prec_m)

