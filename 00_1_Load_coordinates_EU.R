# -------------------------------------------------------------------------
#           
#     SET OF PIXELS USED FOR MAKE PROJECTIONS OF YIELD IN EUROPE    
#                       (EU27 and EU27 extended)
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Loading packages 
library(tidyverse) ; library(stringr)
library(cowplot)
library(terra) ; library(raster)
library("rnaturalearth") ; library("rnaturalearthdata") ; library(sf) ; library(sp) ; library(rworldmap) ; library(osmextract)

# -------------------------------------------------------------------------
# Home-made functions
# > add specific path to the project
source(".../00_Functions.R")

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

# Reference rasters (with yield and irrigation data)
# > add specific path to the data from Iizumi, T., Sakai, T. The global dataset of historical yields for major crops 1981–2016. Sci Data 7, 97 (2020). https://doi.org/10.1038/s41597-020-0433-7 ; the dataset is accessible here: https://doi.pangaea.de/10.1594/PANGAEA.909132
ref_rast <- rast(".../GDHY_v1.3/gdhy_v1.2_v1.3_20190128/soybean/yield_2016.nc4")

# > Reference grid with coordinates of all points in lands
ref_grid0 <- terra::rast(res=0.5)

# List of countries
# > EU 27
eu27_country_names <- c('Austria','Belgium','Bulgaria','Croatia','Cyprus','Northern Cyprus',
                        'Czech Republic','Denmark','Estonia','Finland','France',
                        'Germany','Greece','Hungary','Ireland','Italy',
                        'Latvia','Lithuania','Luxembourg', 'Malta', 'Netherlands',
                        'Poland','Portugal','Romania','Slovakia','Slovenia',
                        'Spain','Sweden')
length(eu27_country_names) # 27+1 (Northern cyprus)

# > EU 42 (continental EU+Switzerland)
eu42_country_names <- c(eu27_country_names, 
                        'Albania', 'Armenia', 'Azerbaijan', 'Belarus', 
                        'Bosnia and Herzegovina', 'Georgia', 
                        'Kosovo', 'Macedonia', 'Moldova', 'Montenegro',  
                        'Norway', 'Republic of Serbia', 'Switzerland', 'Turkey', 
                        'Ukraine', 'United Kingdom')
length(eu42_country_names) # 42+1 (Northern Cyprus)+1 (Switzerland)

list_coords <- list()
# > Shape files and coordinates
for(country_i in eu42_country_names)
{
  
  # > shape file for the countries 
  sph_i <- rnaturalearth::ne_countries(country = country_i, 
                                      scale = "medium", 
                                      returnclass = "sf")
  
  # > split into 0.5° pixels, including those in the sea
  ref_grid_i <- rasterize(sph_i, ref_grid0, touches = T)
  
  # > Format shapefile as data.frame (extract coordinates)
  dat_coords_i <- ref_grid_i %>% 
    as.data.frame(., xy=T) %>% 
    distinct(x, y) 
  
  # > Attribute country and world region to each combination of long and lat
  dat_coords_country_i <- coords2continent(points = dat_coords_i)
  
  # > Store
  list_coords[[paste0(country_i)]] <- dat_coords_country_i
}

# Merge all data and identify sets of countries
dat_coords_merged <- plyr::ldply(list_coords, .id="country_name_check") %>%
  mutate(id_eu27   = if_else(country_name_check %in% eu27_country_names, 1, 0),
         id_ext    = if_else(country_name_check %in% c('Albania', 'Armenia', 'Azerbaijan', 'Belarus', 
                                                       'Bosnia and Herzegovina', 'Georgia', 
                                                       'Kosovo', 'Macedonia', 'Moldova', 'Montenegro',  
                                                       'Norway', 'Republic of Serbia', 'Switzerland', 'Turkey', 
                                                       'Ukraine', 'United Kingdom'), 1, 0),
         id_eu_ext = if_else(country_code %in% eu42_country_names, 1, 0)) 

# Coordinates in the EU27
dat_coords_eu27 <- dat_coords_merged %>%
  # > choose only countries in the EU27
  filter(id_eu27 == 1) %>% 
  # > remove the points outside of Europe (French Guiana, Brazil etc.)
  filter(x >-14, x < 51,
         y > 34, y < 71) %>% 
  distinct(x,y, continent, region, country_name, country_code)

dim(dat_coords_eu27) # N=2699 pixels
unique(dat_coords_eu27$country_name) # NAs are in the sea

# Coordinates in the EU42
dat_coords_eu42 <- dat_coords_merged %>%
  # > remove the points outside of Europe (French Guiana, Brazil etc.)
  filter(x >-14, x < 51,
         y > 34, y < 71) %>% 
  distinct(x,y, continent, region, country_name, country_code)

dim(dat_coords_eu42) # N=4192 pixels
unique(dat_coords_eu42$country_name) # NAs are in the sea


# Check
plot_grid(ggplot() +
            geom_sf(data=europe) +
            geom_point(data=dat_coords_eu27, 
                       aes(x=x, y=y, color=country_name), 
                       size=0.25) + 
            theme_bw() +
            theme(legend.position = "none") +
            lims(x = c(-14,51), y=c(34,71)) +
            ggtitle("EU27 set of data"),
          ggplot() +
            geom_sf(data=europe) +
            geom_point(data=dat_coords_eu42,
                       aes(x=x, y=y, color=country_name), 
                       size=0.25) + 
            theme_bw() +
            theme(legend.position = "none") +
            lims(x = c(-14,51), y=c(34,71)) +
            ggtitle("EU42 set of data"))

# Check in Ireland
plot_grid(ggplot() +
            geom_sf(data=europe) +
            geom_point(data=dat_coords_eu27, aes(x=x, y=y, color=country_name), size=1) + 
            theme_bw() +
            theme(legend.position = "none") +
            lims(x = c(-11,0), y=c(51,57)) +
            ggtitle("EU27 set of data"),
          ggplot() +
            geom_sf(data=europe) +
            geom_point(data=dat_coords_eu42, aes(x=x, y=y, color=country_name), size=1) + 
            theme_bw() +
            theme(legend.position = "none") +
            lims(x = c(-11,0), y=c(51,57)) +
            ggtitle("EU42 set of data"))

# Save
save(dat_coords_eu27, file = ".../data/00_dat_coords_EU27.rda")
save(dat_coords_eu42, file = ".../data/00_dat_coords_EU42.rda")
