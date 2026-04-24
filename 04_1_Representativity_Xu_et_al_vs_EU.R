# Packages & tools
library(tidyverse)
library(stringr)
library(lubridate)
library(terra) ; library(raster) ; library(sf) ; library(rnaturalearth)
library(cowplot)
library(wesanderson)
library(metR)

# > Color palettes
# > gradient blue - yellow - red
pal <- wes_palette("Zissou1", 6, type="continuous")

# > countries maps
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
europe <- rnaturalearth::ne_countries(continent = "europe", scale = "medium", returnclass = "sf")
eu27 <- rnaturalearth::ne_countries(country = c('Finland', 'Sweden', 'Estonia', 'Latvia', 'Denmark', 'Lithuania', 
                                                'Ireland', 'Germany', 'Poland', 'Netherlands', 'Belgium', 'France', 
                                                #'Czech Republic', 
                                                'Czechia',
                                                'Luxembourg', 'Slovakia', 'Austria', 'Hungary', 
                                                'Romania', 'Italy', 'Slovenia', 'Croatia', 'Bulgaria', 'Spain', 'Portugal', 
                                                'Greece', 'Northern Cyprus', 'Cyprus'), scale = "medium", returnclass = "sf")
not_eu27 <- rnaturalearth::ne_countries(country = c('Norway', 'United Kingdom', 'Belarus', 'Ukraine', 'Moldova', 'Switzerland', 
                                                    'Republic of Serbia', 'Bosnia and Herzegovina', 'Montenegro', 'Kosovo', 'Georgia', 
                                                    'Albania', 'Macedonia', 'Turkey', 'Azerbaijan', 'Armenia'), scale = "medium", returnclass = "sf")
eu27_ext <- rnaturalearth::ne_countries(country = c('Finland', 'Sweden', 'Estonia', 'Latvia', 'Denmark', 'Lithuania', 
                                                    'Ireland', 'Germany', 'Poland', 'Netherlands', 'Belgium', 'France', 
                                                    #'Czech Republic', 
                                                    'Czechia',
                                                    'Luxembourg', 'Slovakia', 'Austria', 'Hungary', 
                                                    'Romania', 'Italy', 'Slovenia', 'Croatia', 'Bulgaria', 'Spain', 'Portugal', 
                                                    'Greece', 'Northern Cyprus', 'Cyprus', 
                                                    'Norway', 'United Kingdom', 'Belarus', 'Ukraine', 'Moldova', 'Switzerland', 
                                                    'Republic of Serbia', 'Bosnia and Herzegovina', 'Montenegro', 'Kosovo', 'Georgia', 
                                                    'Albania', 'Macedonia', 'Turkey', 'Azerbaijan', 'Armenia'), scale = "medium", returnclass = "sf")
russia <- rnaturalearth::ne_countries(country = 'Russia', scale = "medium", returnclass = "sf")
# > download ocean outlines
ocean <- ne_download(
  scale = 50,
  type = "ocean",
  category = "physical",
  returnclass = "sf")

# Data europe
load("E:/POSTDOC INRAE/ANALYSES/B_OPTIMISATION/00_Data/02_tab_eu_maize.rda")

# Data meta-analysis
library(readxl)
FCR_maize_and_soybean_intercropping <- read_excel("E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/05_NATURE_COMS/01_R2/data xu et al., 2023/FCR-maize and soybean intercropping.xlsx", 
                                                  col_types = c("numeric", "numeric", "numeric", 
                                                                "text", "text", "numeric", "text", 
                                                                "numeric", "numeric", "text", "text", "text",
                                                                "numeric", "numeric")) %>%
  mutate(Continent = factor(Continent, levels = c("Europe", "Africa", "Asia", "Oceania", "North America", "South America")))

# Distribution of studies
table(FCR_maize_and_soybean_intercropping$Continent)
# Europe   Africa    Asia   Oceania  North America  South America 
# 11       106       464    10       63             4

# Identify the unique study-location paires
FCR_maize_and_soybean_intercropping_unique_loc <- FCR_maize_and_soybean_intercropping %>% 
  distinct(LiteratureID, Lon, Lat, Continent, Country) %>% 
  rename("x"="Lon", "y"="Lat", "Continent_init"="Continent", "Country_init"="Country") %>%
  mutate(Continent_init = factor(Continent_init, levels = c("Europe", "Africa", "Asia", "Oceania", "North America", "South America")))

# ----------------------------------------
# > Attribute country and world region based on longitude and latitude 
# The single argument to this function, points, is a data.frame in which:
#   - column 1 contains the longitude in degrees
#   - column 2 contains the latitude in degrees
coords2continent = function(points)
{  
  require(rworldmap)
  countriesSP <- getMap(resolution='li')
  #countriesSP <- getMap(resolution='high') #you could use high res map from rworldxtra if you were concerned about detail
  
  # converting points to a SpatialPoints object
  # setting CRS directly to that from rworldmap
  pointsSP = SpatialPoints(points, proj4string=CRS(proj4string(countriesSP)))  
  
  
  # use 'over' to get indices of the Polygons object containing each point 
  indices = over(pointsSP, countriesSP)
  
  return(data.frame(
    x=points$x,
    y=points$y,
    continent=indices$continent,   # returns the continent (6 continent model)
    region=indices$REGION,   # returns the continent (7 continent model)
    country_name=indices$ADMIN,  #returns country name
    country_code=indices$ISO3 # returns the ISO3 code 
  ))
  
}


# Check the country associated with the points
country_coords_expe_ma <- coords2continent(points = FCR_maize_and_soybean_intercropping_unique_loc %>% 
                                             dplyr::select(x,y)) %>%
  # two NA in the sea close to china
  mutate(region = if_else(is.na(region), "Asia", region))

# Check by continent
# > Africa
FCR_maize_and_soybean_intercropping_unique_loc %>% 
  filter(Continent_init == "Africa") %>% 
  left_join(., country_coords_expe_ma,
            by=c("x", "y"))

# LiteratureID     x       y Continent_init Country_init continent region country_name country_code
#           52  11.6    4.08  Africa         Cameroon     Africa    Africa Cameroon     CMR         
#           58  37.5    0.518 Africa         Kenya        Africa    Africa Kenya        KEN         
#           58  37.9    0.1   Africa         Kenya        Africa    Africa Kenya        KEN         
#           60   8.28   7.41  Africa         Nigeria      Africa    Africa Nigeria      NGA         
#           64  39.3  -15.3   Africa         Mozambique   Africa    Africa Mozambique   MOZ         
#           64  36.7  -15.3   Africa         Mozambique   Africa    Africa Mozambique   MOZ         
#           64  35.2  -13.3   Africa         Mozambique   Africa    Africa Mozambique   MOZ         
#           69  30.9   31.1   Africa         Egypt        Africa    Africa Egypt        EGY         
#           71   6.12   9.75  Africa         Nigeria      Africa    Africa Nigeria      NGA         
#           75   4.00   8.00  Africa         Nigeria      Africa    Africa Nigeria      NGA         
#           79  -0.66   9.96  Africa         Ghana        Africa    Africa Ghana        GHA         
#           79  -1.01  10.8   Africa         Ghana        Africa    Africa Ghana        GHA         
#           82   9.33   4.23  Africa         Cameroon     Africa    Africa Cameroon     CMR         
#           84   7.38   9.02  Africa         Nigeria      Africa    Africa Nigeria      NGA         
#           88   7.75   8.6   Africa         Nigeria      Africa    Africa Nigeria      NGA         
#           89   5.25   8.36  Africa         Nigeria      Africa    Africa Nigeria      NGA  

# > Asia 
FCR_maize_and_soybean_intercropping_unique_loc %>% 
  filter(Continent_init == "Asia", Country_init!="China") %>% 
  left_join(., country_coords_expe_ma,
            by=c("x", "y"))

# LiteratureID     x     y Continent_init Country_init continent region country_name country_code
#           54  84.4 27.6  Asia           Nepal        Eurasia   Asia   Nepal        NPL         
#           57  79.9 21.2  Asia           India        Eurasia   Asia   India        IND         
#           57 121   15.8  Asia           Philippines  Eurasia   Asia   Philippines  PHL         
#           57  80.7  7.28 Asia           Sri Lanka    Eurasia   Asia   Sri Lanka    LKA         
#           57 100.  13.8  Asia           Thailand     Eurasia   Asia   Thailand     THA         
#           59  48.3 36.7  Asia           Iran         Eurasia   Asia   Iran         IRN         
#           61  78.4 17.3  Asia           India        Eurasia   Asia   India        IND         
#           62  76.1 11.5  Asia           India        Eurasia   Asia   India        IND         
#           66  75.1 15.4  Asia           India        Eurasia   Asia   India        IND         
#           70  77.2 28.6  Asia           India        Eurasia   Asia   India        IND         
#           72  53.1 36.6  Asia           Iran         Eurasia   Asia   Iran         IRN         
#           73  77.4 23.3  Asia           India        Eurasia   Asia   India        IND         
#           76  83.4 27.8  Asia           Nepal        Eurasia   Asia   Nepal        NPL         
#           78 103.  16.4  Asia           Thailand     Eurasia   Asia   Thailand     THA         
#           85  70.7 23.2  Asia           India        Eurasia   Asia   India        IND         
#           87  77.2 28.7  Asia           India        Eurasia   Asia   India        IND   
 
FCR_maize_and_soybean_intercropping_unique_loc %>% 
  filter(Continent_init == "Asia", Country_init == "China") %>% 
  left_join(., country_coords_expe_ma,
            by=c("x", "y")) %>% View(.)

# One point is showing the ocean (ID 19)
# One point is not found by the function, but correctly points in China
FCR_maize_and_soybean_intercropping_unique_loc %>% 
  filter(Continent_init == "Asia", Country_init == "China") %>% 
  left_join(., country_coords_expe_ma,
            by=c("x", "y")) %>% 
  filter(is.na(country_name))
 
# > Europe
FCR_maize_and_soybean_intercropping_unique_loc %>% 
  filter(Continent_init=="Europe") %>% 
  left_join(., country_coords_expe_ma,
            by=c("x", "y"))

# LiteratureID     x     y Continent_init Country_init continent region country_name      
#           63  26.3  44.5 Europe         Romania      Eurasia   Europe Romania           
#           83  20.4  44.8 Europe         Serbia       Eurasia   Europe Republic of Serbia
#           90  22.2  55.7 Europe         Lithuania    Eurasia   Europe Lithuania         
#           90  24.2  56.0 Europe         Lithuania    Eurasia   Europe Lithuania         
#           90  24.4  54.2 Europe         Lithuania    Eurasia   Europe Lithuania   

# > North and South America
FCR_maize_and_soybean_intercropping_unique_loc %>% 
  filter(Continent_init %in% c("North America", "South America")) %>% 
  left_join(., country_coords_expe_ma,
            by=c("x", "y"))

# LiteratureID      x     y Continent_init Country_init             continent     region country_name country_code
#           51  -73.4  45.8 North America  Canada                   North America North… Canada       CAN         
#           51  -73.6  45.5 North America  Canada                   North America North… Canada       CAN         
#           51  -73.6  45.5 North America  Canada                   North America North… Canada       CAN         
#           53  -75.7  45.4 North America  Canada                   North America North… Canada       CAN         
#           55  -85.7  32.4 North America  United States of America North America North… United Stat… USA         
#           57 -158.   21.4 North America  United States of America North America North… United Stat… USA         
#           65  -58.3 -37.8 South America  Argentina                South America South… Argentina    ARG         
#           67  -75.7  45.4 North America  Canada                   North America North… Canada       CAN         
#           77  -76.7  39.5 North America  United States of America North America North… United Stat… USA  

# > Oceania
FCR_maize_and_soybean_intercropping_unique_loc %>% 
  filter(Continent_init == "Oceania") %>% 
  left_join(., country_coords_expe_ma,
            by=c("x", "y"))

# LiteratureID     x     y Continent_init Country_init continent region    country_name country_code
#           57  151. -33.9 Oceania        Australia    Australia Australia Australia    AUS         
#           68  145  -36   Oceania        Australia    Australia Australia Australia    AUS      

# ----------------------------------------
# Koppen Geiger zones

# > read Koppen-Geiger classification in raster format
# data are from http://koeppen-geiger.vu-wien.ac.at/present.htm
# August 2022 version
KGmap <- raster("C:/Users/benni/Documents/Post doc/Map_KG-Global/Map_KG-Global/KG_1986-2010.grd") ; KGmap

#class      : RasterLayer 
#dimensions : 2160, 4320, 9331200  (nrow, ncol, ncell)
#resolution : 0.08333333, 0.08333333  (x, y)
#extent     : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#crs        : +proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 
#source     : KG_1986-2010.grd 
#names      : layer 
#values     : 1, 32  (min, max)

# > labels
climate_types <- data.frame(layer=1:32, 
                            layer_name=c('Af', 'Am', 'As', 'Aw', 'BSh', 'BSk', 'BWh', 'BWk', 'Cfa', 'Cfb','Cfc', 'Csa', 'Csb', 'Csc', 'Cwa','Cwb', 'Cwc', 'Dfa', 'Dfb', 'Dfc','Dfd', 'Dsa', 'Dsb', 'Dsc', 'Dsd','Dwa', 'Dwb', 'Dwc', 'Dwd', 'EF','ET', 'Ocean'),
                            layer_1   =c('A', 'A', 'A', 'A', 'B', 'B', 'B', 'B', 'C', 'C','C', 'C', 'C', 'C', 'C','C', 'C', 'D', 'D', 'D','D', 'D', 'D', 'D', 'D','D', 'D', 'D', 'D', 'E','E', 'Ocean'))

# > resolution reduction (from 0.083333 -> 0.5)
#KGmap <- terra::aggregate(KGmap, fact=6, fun="modal") ; KGmap
#KGmap

# > convert as sf object
KGmap_sf=st_as_sf(x = as.data.frame(KGmap, xy=T), coords = c("x", "y"), crs = 4326)

# > Points with the MA experiments
points <- tibble(longitude = FCR_maize_and_soybean_intercropping_unique_loc$x,
                 latitude  = FCR_maize_and_soybean_intercropping_unique_loc$y)

# > convert to spatial
points <- points %>% 
  st_as_sf(coords = c("longitude", "latitude"),
           crs = 4326)

# > check that the two objects have the same CRS
st_crs(points) == st_crs(KGmap_sf)

# > locate for each point the nearest polygon
nearest_polygon <- st_join(points, KGmap_sf,
                           join = st_nearest_feature)

nearest_MA <- data.frame(set="MA", 
                         st_coordinates(nearest_polygon$geometry),
                         layer=nearest_polygon$layer) %>% 
  left_join(climate_types)

# > Points with data in Europe
eu_coords <- tab_maize_EU %>% 
  distinct(gridcode,x,y)

summary(eu_coords$x)

# > coordinates
points_eu <- eu_coords %>% 
  st_as_sf(coords = c("x", "y"),
           crs = 4326)

# > locate for each point the nearest polygon
nearest_polygon_eu <- st_join(points_eu, KGmap_sf,
                              join = st_nearest_feature)

nearest_EU <- data.frame(set="EU", 
                         st_coordinates(nearest_polygon_eu$geometry),
                         layer=nearest_polygon_eu$layer) %>% 
  left_join(climate_types)


# > merge the results 
nearest <- rbind(nearest_MA, nearest_EU)

nearest %>% 
  filter(layer_name!="Ocean") %>% 
  group_by(set, layer_1) %>% 
  count() %>% 
  group_by(set) %>% 
  mutate(tot = sum(n),
         freq=(n/tot)*100) %>%
  ggplot(.) +
  geom_raster(aes(x=layer_1, y=set, fill=freq)) +
  geom_text(aes(x=layer_1, y=set, label=round(freq, 1))) + 
  scale_fill_viridis_c(option="B", direction=-1, begin=0.3, name="Frequence (%)") +
  theme_cowplot() +
  theme(legend.position = "bottom") +
  labs(x="Köppen Geiger classification", y ="")


ggplot() +
  geom_sf(data=world %>% 
            filter(continent !="Antarctica"), color="black", fill="grey95")  + 
  geom_point(data = nearest %>% 
               filter(layer_name!="Ocean"), 
             aes(x=X, y=Y, color = layer_name, size=set))+
  theme_map() +
  theme(legend.position = "bottom") +
  scale_color_viridis_d(option="B", direction=-1) +
  scale_size_manual(values=c(0.1,2))+
  facet_grid(set~.)

# ----------------------------------------
# > Coordinates of experimental studies in the meta-analysis 
expe_coords <- FCR_maize_and_soybean_intercropping %>% 
  distinct(StudyID, Lon, Lat, Planting_time, Harvesting_time) %>% 
  rename("x"="Lon", "y"="Lat", "gridcode"="StudyID") %>% 
  mutate(x=if_else(x < 0, x + 360, x)) %>% 
  mutate(set="expe_ma")

summary(expe_coords$x)
names(expe_coords)

# > Coordinates of the grid-cells with soybean and maize yield projections in the EU 
eu_coords <- tab_maize_EU %>% 
  distinct(gridcode,x,y) %>% 
  mutate(x=if_else(x < 0, x + 360, x)) %>% 
  mutate(Planting_time = 4, Harvesting_time = 11) %>% 
  mutate(set="eu")

summary(eu_coords$x)
names(eu_coords)

# > coordinates
all_coords <- rbind(expe_coords, eu_coords)
dim(all_coords) # 4365

# > load 1 initial yield file to resample era5 data 
yield_ref <- rast("E:/POSTDOC INRAE/DATA/02_YIELDS/GDHY_v1.3/gdhy_v1.2_v1.3_20190128/maize/yield_1981.nc4")

# Load by year the data 

data_all_years <- NULL

for(year_i in c("2000_2010", "2011_2023"))
{
  
  # > export the data from .nc file
  # temperature
  t2m_raster_i_init <- rast(paste0("C:/Users/benni/Documents/Post doc/Test_month/monthly_2m_temperature_", year_i , ".nc")) 
  t2m_raster_i <- terra::aggregate(t2m_raster_i_init, fact=2, fun="mean")
  crs(t2m_raster_i) <- "epsg:4326"# > add projection
  t2m_raster_i_resample <- resample(t2m_raster_i, yield_ref)# > realign era5 data on yield data
  
  # total precipitations
  tp_raster_i_init <- rast(paste0("C:/Users/benni/Documents/Post doc/Test_month/monthly_total_precipitations_", year_i , ".nc")) 
  tp_raster_i <- terra::aggregate(tp_raster_i_init, fact=2, fun="mean")
  crs(tp_raster_i) <- "epsg:4326"# > add projection
  tp_raster_i_resample <- resample(tp_raster_i, yield_ref) # > realign era5 data on yield data
  
  # > extract the data for these points in particular 
  # temperature
  t2m_year_i <- cbind(all_coords, extract(t2m_raster_i_resample, all_coords[,2:3])) %>% 
    gather(key = "variable", value="t2m", starts_with("t2m")) %>% 
    separate(variable, c("clim_var", "date_epoch"), sep="_valid_time=", remove = T) %>%
    mutate(date = as.POSIXct(as.numeric(date_epoch), origin="1970-01-01", tz="GMT")) %>% 
    mutate(month = lubridate::month(date),
           year = lubridate::year(date)) %>% 
    dplyr::select(-date, -date_epoch, -clim_var) %>%
    # > correct value (convert K in °C)
    mutate(t2m = t2m-273.15)
  
  # total precipitations
  tp_year_i <- cbind(all_coords, extract(tp_raster_i_resample, all_coords[,2:3])) %>% 
    gather(key = "variable", value="tp", starts_with("tp")) %>% 
    separate(variable, c("clim_var", "date_epoch"), sep="_valid_time=", remove = T) %>%
    mutate(date = as.POSIXct(as.numeric(date_epoch), origin="1970-01-01", tz="GMT")) %>% 
    mutate(month = lubridate::month(date),
           year = lubridate::year(date)) %>% 
    dplyr::select(-date, -date_epoch, -clim_var) %>%
    # > correct value (convert L in mL)
    mutate(tp = tp*10^3)
  
  # > merge
  data_year_i <- t2m_year_i %>% 
    left_join(tp_year_i)
  
  # > add to previous years
  data_all_years <- rbind(data_all_years, data_year_i)
  
  cat("=")
  
}

# > compute the mean by grid-cell:
data_mean_2000_2023 <- data_all_years %>% 
  split(.$gridcode) %>% 
  map_dfr(., ~{
    
    data_i <- .x 
    
    if(unique(data_i$Planting_time) < unique(data_i$Harvesting_time))
    {
      
      data_i <- data_i %>% 
        mutate(growing_period = year) %>% 
        group_by(gridcode, growing_period) %>%
        mutate(to_keep = if_else(month >= Planting_time & month <= Harvesting_time, 1, 0)) %>% 
        filter(to_keep == 1) 
      #%>% 
      #  group_by(gridcode, set, x, y, growing_period) %>%  
      #  summarise(mean_t2m_year = mean(t2m),
      #            sum_tp_year = sum(tp)) 
      
    }
    
    if(unique(data_i$Planting_time) > unique(data_i$Harvesting_time))
    {
      
      data_i <- data_i %>% 
        mutate(growing_period = if_else(month <= Harvesting_time, year-1, year)) %>% 
        group_by(gridcode, growing_period) %>%
        mutate(to_keep = if_else(month >= Planting_time | month <= Harvesting_time, 1, 0)) %>% 
        filter(to_keep == 1, growing_period > 2000, growing_period < 2023) 
      #%>% 
      #  group_by(gridcode, set, x, y, growing_period) %>%  
      #  summarise(mean_t2m_year = mean(t2m),
      #            sum_tp_year = sum(tp)) 
        
      
    }
    
    data_i
    
  }, .id="gridcode")
  
  



data_mean_2000_2023_year <- data_mean_2000_2023 %>% 
  group_by(set, x, y, growing_period) %>%  
  summarise(mean_t2m_year = mean(t2m),
            sum_tp_year = sum(tp)) 
#%>% 
#  group_by(set, x, y) %>% 
#  summarise(mean_t2m = mean(mean_t2m_year),
#            mean_sum_tp = mean(sum_tp_year))
#


# > Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
} 





p1 <- ggplot() +
  geom_sf(data=eu27 %>% 
            filter(continent=="Europe"), fill="#F69C73FF") + 
  geom_sf(data=world %>% 
            filter(continent !="Antarctica"), color="grey50", fill="transparent")  + 
  geom_point(data = country_coords_expe_ma %>% 
               mutate(region =factor(region, levels = c("Europe", "Africa", "Asia", "Australia", "North America", "South America"))), 
             aes(x=x, y=y, fill = region, pch=region),size=2)+
  scale_fill_viridis_d(name="Continents", direction=-1,
                       guide = guide_legend(nrow=1), 
                       na.value = "red") +
  scale_shape_manual(values = c(21:25,21), 
                     name = "Continents", na.value = 3) + 
  theme_map() +
  theme(legend.position = "bottom"); p1

# > based on 2000-2023 conditions
tile_data_mean_2000_2023_year <- data_mean_2000_2023_year[which(data_mean_2000_2023_year$set=="eu"),]
density_eu <- get_density(x = tile_data_mean_2000_2023_year$mean_t2m_year, 
                          y = tile_data_mean_2000_2023_year$sum_tp_year, 
                          n=100)
tile_data_mean_2000_2023_year$density <- density_eu

p2 <- ggplot() +
  geom_point(data = tile_data_mean_2000_2023_year, 
             aes(x=mean_t2m_year, y=sum_tp_year, color=density), pch = 15) +
  geom_point(data = data_mean_2000_2023_year[which(data_mean_2000_2023_year$set != "eu"),] %>% 
               left_join(., country_coords_expe_ma %>% 
                           mutate(x=if_else(x < 0, x + 360, x)) %>% 
                           mutate(set="expe_ma"), 
                         by=c("x","y","set")) %>% 
               mutate(region = factor(region, 
                                      levels = c("Europe", "Africa", "Asia", "Australia", 
                                                 "North America", "South America"))) %>% 
               group_by(set, x, y, region) %>% 
               summarise(mean_t2m = mean(mean_t2m_year),
                         mean_sum_tp = mean(sum_tp_year)),  
             aes(x=mean_t2m, y=mean_sum_tp, 
                 fill=region, shape = region),
             size=2) +
  scale_color_viridis_c(direction=-1, name = "Density\n(for grid-cells with intercropping projections)", 
                        option = "rocket", end = 0.9,
                        guide = guide_colorbar(barwidth=10, barheight=0.5)) +
  scale_fill_viridis_d(name="Continents", direction=-1,
                       guide = guide_legend(nrow=3), na.value = "red") +
  scale_shape_manual(values = c(21:25,21), name = "Continents", na.value = 3) + 
  theme_cowplot() +
  theme(legend.position = "bottom",
        legend.title.position = "top",
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        panel.grid.major = element_line(color="grey92"),
        panel.border = element_rect(color="black", linewidth = 0.75),
        axis.line = element_blank()) +
  scale_y_log10() +
  labs(x = "Average temperature\n(°C)", y = "Average total precipitations\n(log. 10 scale, mL)") ; p2

plot_grid(p1+theme(legend.position="none"),
          p2,
          labels=c("a.", "b."), 
          vjust = c(1.5, 0.5),
          rel_heights = c(0.4,0.6),
          ncol=1)

ggsave("E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/05_NATURE_COMS/01_R2/FIGURES/MA_vs_data.png",
       width = 7, height = 9, dpi=300, bg="white")


