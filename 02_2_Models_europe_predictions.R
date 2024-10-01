# -------------------------------------------------------------------------
# 
#       02-2. Predict yield in Europe based on climate and irrigated portion
#       Author: M. Chen, Inrae, 2024 
# 
# -------------------------------------------------------------------------

# ----------------------------------
# > Packages
library(tidyverse) ; library(CCMHr)
# > machine learning
library(parallel) ; library(doParallel); library(foreach)
library(caret) ; library(ranger) ; 
library(terra) ; library(raster) ; library("rnaturalearth") ; library("rnaturalearthdata") ; library(sf) ; library(sp) ; library(rworldmap) ; library(rmapshaper) ; library(tidygeocoder) ; 
library(cowplot) ; library(metR)
# > others
library(hydroGOF)

# -----------------------------------
# Home-made functions
source("E:/POSTDOC INRAE/ANALYSES/00_TOOLS/00_Functions.R")

# Coordinates of pixels in EU
load("E:/POSTDOC INRAE/ANALYSES/B_OPTIMISATION/00_Data/00_dat_coords_EU27.rda")
dat_coords_eu27$gridcode = paste0(dat_coords_eu27$x, "_", dat_coords_eu27$y)

load("E:/POSTDOC INRAE/ANALYSES/B_OPTIMISATION/00_Data/00_dat_coords_EU42.rda")
dat_coords_EU <- dat_coords_eu42
dat_coords_eu42$gridcode = paste0(dat_coords_eu42$x, "_", dat_coords_eu42$y)

# ----------------------------------
# DATA

# Climate in Europe + irrigated portion
tab_soybean_EU <- loadRDa("E:/POSTDOC INRAE/ANALYSES/B_OPTIMISATION/00_Data/02_tab_eu_soybean.rda") %>% 
  # rename VPD PC 
  rename("PC1_month_vapor_pressure_deficit"="PC1_month_vpd_1",
         "PC2_month_vapor_pressure_deficit"="PC2_month_vpd_1",
         "PC3_month_vapor_pressure_deficit"="PC3_month_vpd_1",
         "PC4_month_vapor_pressure_deficit"="PC4_month_vpd_1",
         "PC5_month_vapor_pressure_deficit"="PC5_month_vpd_1",
         "PC6_month_vapor_pressure_deficit"="PC6_month_vpd_1",
         "PC7_month_vapor_pressure_deficit"="PC7_month_vpd_1") %>% 
  dplyr::select(-country_name) %>% 
  left_join(dat_coords_EU %>% 
              dplyr::select(x, y, country_name)) %>%
  # distinguish between those in EU27 and this in EU42
  mutate(id_eu27 = if_else(gridcode %in% unique(dat_coords_eu27$gridcode), 1, 0), 
         id_ext  = if_else(gridcode %in% unique(dat_coords_eu42$gridcode), 1, 0), 
         id_eu_ext = if_else(id_eu27 != id_ext & id_eu27 == 0, 1, 0))

# check 
tab_soybean_EU %>% 
  distinct(gridcode, id_eu27, id_ext, id_eu_ext) %>% 
  gather(key=id, value=value, -gridcode) %>%
  group_by(id, value) %>%
  count() %>% 
  filter(value==1)
# 2699 pixels in EU27
# 1493 pixels outside of EU27
# 4192 pixels in total

tab_maize_EU   <- loadRDa("E:/POSTDOC INRAE/ANALYSES/B_OPTIMISATION/00_Data/02_tab_eu_maize.rda") %>% 
  rename("PC1_month_vapor_pressure_deficit"="PC1_month_vpd_1",
         "PC2_month_vapor_pressure_deficit"="PC2_month_vpd_1",
         "PC3_month_vapor_pressure_deficit"="PC3_month_vpd_1",
         "PC4_month_vapor_pressure_deficit"="PC4_month_vpd_1",
         "PC5_month_vapor_pressure_deficit"="PC5_month_vpd_1",
         "PC6_month_vapor_pressure_deficit"="PC6_month_vpd_1",
         "PC7_month_vapor_pressure_deficit"="PC7_month_vpd_1",
         "PC8_month_vapor_pressure_deficit"="PC8_month_vpd_1") %>% 
  dplyr::select(-country_name) %>% 
  left_join(dat_coords_EU %>% 
              dplyr::select(x, y, country_name)) %>%
  # distinguish between those in EU27 and this in EU42
  mutate(id_eu27 = if_else(gridcode %in% unique(dat_coords_eu27$gridcode), 1, 0), 
         id_ext  = if_else(gridcode %in% unique(dat_coords_eu42$gridcode), 1, 0), 
         id_eu_ext = if_else(id_eu27 != id_ext & id_eu27 == 0, 1, 0))

# check 
tab_maize_EU %>% 
  distinct(gridcode, id_eu27, id_ext, id_eu_ext) %>% 
  gather(key=id, value=value, -gridcode) %>%
  group_by(id, value) %>%
  count() %>% 
  filter(value==1)
# 2699 pixels in EU27
# 1493 pixels outside of EU27
# 4192 pixels in total

# Surface of cropland per pixel
# > initial file
cropland_area <- rast("E:/POSTDOC INRAE/DATA/02_YIELDS/SASAM/Global-cropland-percentage-map.tif") ; cropland_area
#class       : SpatRaster 
#dimensions  : 30521, 84399, 1  (nrow, ncol, nlyr)
#resolution  : 0.004166667, 0.004166667  (x, y)
#extent      : -173.1125, 178.55, -55.97917, 71.19167  (xmin, xmax, ymin, ymax)
#coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#source      : Global-cropland-percentage-map.tif 
#name        : Global-cropland-percentage-map 

# > reduce resolution 
agg_cropland_area <- terra::aggregate(cropland_area, fact=120, fun="mean", cores=7) ; agg_cropland_area
#class       : SpatRaster 
#dimensions  : 255, 704, 1  (nrow, ncol, nlyr)
#resolution  : 0.5, 0.5  (x, y)
#extent      : -173.1125, 178.8875, -56.30833, 71.19167  (xmin, xmax, ymin, ymax)
#coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#source(s)   : memory
#name        : Global-cropland-percentage-map 
#min value   :                      0.0000000 
#max value   :                      0.9887098 

# Rough coordinates of European extent
ext_eur <- c(-14,51,34,71)

# Crop the temperature layer (bio1) to roughly European extent
agg_cropland_area_eu <- terra::crop(agg_cropland_area, ext_eur)
surf_agg_area_eu <- cellSize(agg_cropland_area_eu, unit="ha")
surf_agg_cropland_area_eu <- agg_cropland_area_eu*surf_agg_area_eu

# > realign cropland  data on yield data
# > load 1 initial yield file to resample era5 data 
yield_ref <- rast("E:/POSTDOC INRAE/DATA/02_YIELDS/GDHY_v1.3/gdhy_v1.2_v1.3_20190128/maize/yield_1981.nc4")

resample_agg_cropland_area <- resample(surf_agg_cropland_area_eu, yield_ref, method="med") ; resample_agg_cropland_area
#class       : SpatRaster 
#dimensions  : 360, 720, 1  (nrow, ncol, nlyr)
#resolution  : 0.5, 0.5  (x, y)
#extent      : 0, 360, -90, 90  (xmin, xmax, ymin, ymax)
#coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#source(s)   : memory
#name        : Global-cropland-percentage-map 
#min value   :                      0.0000000 
#max value   :                      0.9764748 

# > table with proportion of cropland per pixel
tab_cropland_area <- as.data.frame(resample_agg_cropland_area, xy=T) %>% 
  rename("cropland_area_ha"="Global-cropland-percentage-map") %>% 
  # > round x and y to 2 digits to be consistent with yield data 
  mutate(x=round(x,2),
         y=round(y,2)) %>% 
  mutate(x=if_else(x>180, x-360, x))

summary(tab_cropland_area$cropland_area_ha)

# Min.    1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.00000 0.00000 0.03737 0.00000 0.97647

# 0.00      0.00     49.27  26294.52  35461.88 212096.23 

plot(density(tab_cropland_area$cropland_area_ha))

dat_coords_EU %>% 
  left_join(tab_cropland_area) %>%
  ggplot(.) +
  geom_tile(aes(x=x, y=y, fill=cropland_area_ha)) +
  geom_sf(data=europe, fill="transparent") +
  scale_fill_gradientn(colors = c("red", viridis::viridis(99, direction=-1))) +
  theme_map() +
  lims(x=c(-14,51), y=c(34,71))

# > Check surfaces
# per country
dat_coords_eu27 %>% 
  left_join(tab_cropland_area) %>% 
  group_by(country_name) %>% 
  summarise(sum_cropland_area=sum(cropland_area_ha)) %>% 
  arrange(desc(sum_cropland_area))

#  country_name   sum_cropland_area
#1 France                 16885607.
#2 Spain                  13156454.
#3 Poland                 11522993.
#4 Germany                10786577.
#5 Romania                 8347285.
#6 Italy                   6843613.

# Total in EU27
dat_coords_eu27 %>% 
  left_join(tab_cropland_area) %>% 
  pull(cropland_area_ha) %>% sum(.) # 105629027 ha = 1056290 km2

# -> 1/4 cropland in EU = 105629027/4 = 26407257 ha (26.4 Mha)

# ----------------------------------
# Scenarios to predict

# Models list (top 3 of models based on previous analyses)
list_models <- list(
  pca.m.2      = list(name = "pca.m.2"),
  pca.m.3      = list(name = "pca.m.3"),
  avg.s        = list(name = "avg.s"),
  avg.m        = list(name = "avg.m")
)

# ----------------------------------
# PREDICTIONS FOR EUROPE

preds_eu <- list(soybean = list(crop = "soybean", tab_eu = tab_soybean_EU, list_models = list_models),
                 maize   = list(crop = "maize",   tab_eu = tab_maize_EU,   list_models = list_models)) %>% 
  map_dfr(., ~{
    
    # Crop 
    crop <- .x$crop
    
    # Data for predictions
    tab_eu <- .x$tab_eu 
    
    # Set irrigation to 0% 
    tab_eu$irrigated_portion_perc_init <-  tab_eu$irrigated_portion_perc
    tab_eu$irrigated_portion_perc <- 0
    
    # Used the model for each crop and 
    # predict yields for europe
    preds_crop_eu <- list_models %>% 
      map_dfr(., ~{
        
        # > load the model
        mod_i <- loadRDa(paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_mods/", crop, "_", .x$name, "_train.rda"))
        
        # > predict for EU 
        preds_eu <- predict(mod_i, 
                            data = tab_eu, 
                            type = "response", 
                            seed = 101,
                            num.trees = mod_i$num.trees)
        
        # > predictions in the dataset 
        tab_eu$Ya_pred <- as.numeric(as.character(preds_eu$predictions))
        
        # > output
        tab_eu %>% 
          ungroup() %>% 
          dplyr::select(x, y, year, Ya_pred, country_name, id_eu27, id_ext, id_eu_ext) %>% 
          # > Merge with cropland area per pixel 
          left_join(., tab_cropland_area, by = c("x", "y")) 
        
        
      }, .id = "model")
    
    # out
    preds_crop_eu
    
  }, .id = "crop") 


# > Distribution of predicted yields in Europe 
preds_eu %>% 
  group_by(crop, model) %>% 
  summarise(mean = mean(Ya_pred), 
            sd   = sd(Ya_pred),
            min  = min(Ya_pred),
            max  = max(Ya_pred)) 

#   crop    model    mean    sd    min     max
# 1 maize   avg.m    5.83  2.94    0       14.7 
# 2 maize   avg.s    5.27  3.21    0       15.3 
# 3 maize   pca.m.2  5.28  2.36    0.0795  10.4 
# 4 maize   pca.m.3  5.47  2.38    0.109   10.8 
# 5 soybean avg.m    1.84  1.18    0        5.93
# 6 soybean avg.s    1.54  1.19    0        6.09
# 7 soybean pca.m.2  1.76  0.963   0.00501  4.84
# 8 soybean pca.m.3  1.84  1.05    0.0117   5.27

# > Number of pixels with yields < 1 t ha
preds_eu %>% 
  group_by(crop, model, x, y) %>% 
  summarise(mean = mean(Ya_pred)) %>%
  mutate(id_low=if_else(mean<1, "lower than 1", "higher than 1")) %>%
  group_by(crop, model, id_low) %>% 
  count() %>% 
  spread(key=id_low, value=n)

#  crop    model   `higher than 1` `lower than 1`
#1 maize   avg.m              3668            524
#2 maize   avg.s              3576            616
#3 maize   pca.m.2            3723            469   <- 
#4 maize   pca.m.3            3771            421
#5 soybean avg.m              3108           1084
#6 soybean avg.s              2640           1552
#7 soybean pca.m.2            3240            952   <- 
#8 soybean pca.m.3            3254            938

# > Temporal variation of yield per year 
preds_eu %>% 
  ggplot(., aes(x=year, Ya_pred)) +
  geom_boxplot(outlier.shape = NA) + theme_bw() +
  facet_grid(model~crop, scales = "free")

# > Geographical distribution
# > gradient blue - yellow - red
pal <- wesanderson::wes_palette("Zissou1", 6, type="continuous")

# > countries maps
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
europe <- rnaturalearth::ne_countries(country = c('Finland', 'Sweden', 'Estonia', 'Latvia', 'Denmark', 'Lithuania', 
                                                  'Ireland', 'Germany', 'Poland', 'Netherlands', 'Belgium', 'France', 
                                                  'Czech Republic', 'Luxembourg', 'Slovakia', 'Austria', 'Hungary', 
                                                  'Romania', 'Italy', 'Slovenia', 'Croatia', 'Bulgaria', 'Spain', 'Portugal', 
                                                  'Greece', 'Northern Cyprus', 'Cyprus', 
                                                  'Norway', 'United Kingdom', 'Belarus', 'Ukraine', 'Moldova', 'Switzerland', 
                                                  'Republic of Serbia', 'Bosnia and Herzegovina', 'Montenegro', 'Kosovo', 'Georgia', 
                                                  'Albania', 'Macedonia', 'Turkey', 'Azerbaijan', 'Armenia'), scale = "medium", returnclass = "sf")
eu27 <- rnaturalearth::ne_countries(country = c('Finland', 'Sweden', 'Estonia', 'Latvia', 'Denmark', 'Lithuania', 
                                                'Ireland', 'Germany', 'Poland', 'Netherlands', 'Belgium', 'France', 
                                                'Czech Republic', 'Luxembourg', 'Slovakia', 'Austria', 'Hungary', 
                                                'Romania', 'Italy', 'Slovenia', 'Croatia', 'Bulgaria', 'Spain', 'Portugal', 
                                                'Greece', 'Northern Cyprus', 'Cyprus'), scale = "medium", returnclass = "sf")
not_eu27 <- rnaturalearth::ne_countries(country = c('Norway', 'United Kingdom', 'Belarus', 'Ukraine', 'Moldova', 'Switzerland', 
                                                    'Republic of Serbia', 'Bosnia and Herzegovina', 'Montenegro', 'Kosovo', 'Georgia', 
                                                    'Albania', 'Macedonia', 'Turkey', 'Azerbaijan', 'Armenia'), scale = "medium", returnclass = "sf")


# > maps

# Breaks for geom_contour_fill, geom_contour, geom_text_contour 
breaks_plot <- c(0, seq(0.5, 3.5, by=0.5), seq(4, 10, by=1), 12, 14, 16)
breaks_labels <- seq(1,10, by=1)

plot <- preds_eu %>% 
  group_by(crop, model, x, y) %>% 
  summarise(mean_Ya_pred = mean(Ya_pred)) %>%
  split(.$crop) %>% 
  purrr::map(., ~{
    .x %>% 
      split(.$model) %>% 
      purrr::map(., ~{ 
        
        .x %>%
          # > plot
          ggplot(.) + 
          geom_contour_fill(aes(x=x, y=y, z=mean_Ya_pred), 
                            breaks = breaks_plot,
                            na.fill = T, 
                            global.breaks = F,
                            clip = europe) +
          geom_sf(data=world, fill="transparent") +
          #geom_contour(aes(x=x, y=y, z=mean_Ya_pred), 
          #             color = "white", 
          #             linewidth = 0.1,
          #             breaks = breaks_plot) +
          facet_grid(.~crop) + 
          theme_map() + 
          lims(x = c(-11,51), y=c(33,71)) + 
          theme(legend.position = "bottom",, 
                legend.title = element_text(size=12), 
                legend.text = element_text(size=10),
                strip.text = element_text(face = "bold", hjust = 0.1, size=15)) +
          scale_fill_gradientn(colours = c("transparent", c("red", viridis::viridis(direction=-1, n=100))),
                               breaks = breaks_plot, 
                               labels = breaks_plot,
                               guide = guide_colorbar(barwidth = 20, barheight = 0.5, title.position = "top", title = "Mean yield (t/ha)"))
        
        
        })
    
  })
  

# -----------------
# Plots 

# > Plot per model
# avg.m
p_avg.m <- plot_grid(plot$maize$avg.m + ggtitle("Predictions based on monthly averages of climate data"), 
          plot$soybean$avg.m, 
          ncol = 2, axis = "lrtb", align = "hv") ; p_avg.m

ggsave(p_avg.m, filename = "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_preds/04_EU/plot_pred_EU_avg.m.png", 
       width=9, height=5, bg = "white")

# avg.s
p_avg.s <- plot_grid(plot$maize$avg.s + ggtitle("Predictions based on seasonal averages of climate data"), 
          plot$soybean$avg.s, 
          ncol = 2, axis = "lrtb", align = "hv") ; p_avg.s

ggsave(p_avg.s, filename = "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_preds/04_EU/plot_pred_EU_avg.s.png", 
       width=9, height=5, bg = "white")

# pca.m.2
p_pca.m.2 <- plot_grid(plot$maize$pca.m.2 + ggtitle("Predictions based on two first principal components\nfrom monthly averages of climate data"), 
          plot$soybean$pca.m.2, 
          ncol = 2, axis = "lrtb", align = "hv") ; p_pca.m.2

ggsave(p_pca.m.2, filename = "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_preds/04_EU/plot_pred_EU_pca.m.2.png", 
       width=9, height=5, bg = "white")

# pca.m.3
p_pca.m.3 <- plot_grid(plot$maize$pca.m.3 + ggtitle("Predictions based on three first principal components\nfrom monthly averages of climate data"), 
          plot$soybean$pca.m.3, 
          ncol = 2, axis = "lrtb", align = "hv") ; p_pca.m.3

ggsave(p_pca.m.3, filename = "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_preds/04_EU/plot_pred_EU_pca.m.3.png", 
       width=9, height=5, bg = "white")

# -----------------
# > Plot per crop 
# soybean 
p_soybean <- plot_grid(plot$soybean$avg.m   + ggtitle("Soybean", subtitle = "Model with monthly averages of climate data"), 
                       plot$soybean$avg.s   + ggtitle("",        subtitle = "Model with annual averages of climate data"), 
                       plot$soybean$pca.m.2 + ggtitle("",        subtitle = "Model with 2 scores from month-based PCA"),
                       plot$soybean$pca.m.3 + ggtitle("",        subtitle = "Model with 3 scores from month-based PCA"), 
                       nrow = 1, axis = "lrtb", align = "hv") ; p_soybean

ggsave(p_soybean, filename = "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_preds/04_EU/plot_pred_EU_soybean.png", 
       width=16, height=8, bg = "white")

# maize 
p_maize <- plot_grid(plot$maize$avg.m   + ggtitle("Maize", subtitle = "Model with monthly averages of climate data"), 
                     plot$maize$avg.s   + ggtitle("",      subtitle = "Model with annual averages of climate data"), 
                     plot$maize$pca.m.2 + ggtitle("",      subtitle = "Model with 2 scores from month-based PCA"),
                     plot$maize$pca.m.3 + ggtitle("",      subtitle = "Model with 3 scores from month-based PCA"), 
                     nrow = 1, axis = "lrtb", align = "hv") ; p_maize

ggsave(p_maize, filename = "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_preds/04_EU/plot_pred_EU_maize.png", 
       width=16, height=8, bg = "white")

# -----------------
# > Estimation of the quantity produced
prod_preds_eu <- preds_eu %>% 
  # > Potential yield per pixel depending on the frequency 
  #   of crop cultivation 
  mutate(s33 = 0.33*cropland_area_ha*Ya_pred,     # 33% of cropland ~ 1 year on 3
         s25 = 0.25*cropland_area_ha*Ya_pred,     # 25% of cropland ~ 1 year on 4
         s20 = 0.20*cropland_area_ha*Ya_pred,     # 20% of cropland ~ 1 year on 5
         s16 = 0.16*cropland_area_ha*Ya_pred) %>% # 16% of cropland ~ 1 year on 6
  # > Format of data 
  gather(key="crop_frequency", value="prod_Ya_pred", s33, s25, s20, s16) %>% 
  # > Labels
  mutate(crop_frequency_lab = recode(crop_frequency, 
                                     "s33"="1 year in 3",
                                     "s25"="1 year in 4", 
                                     "s20"="1 year in 5",
                                     "s16"="1 year in 6"))

# > Production (Mt) of soybean and maize in each scenario
prod_preds_eu %>% 
  group_by(crop, model, crop_frequency_lab) %>% 
  summarise(total_prod_Ya_pred = sum(prod_Ya_pred)/10^6) %>% 
  spread(key=crop_frequency_lab, value=total_prod_Ya_pred)

#  crop    model   `1 year in 3` `1 year in 4` `1 year in 5` `1 year in 6`
#1 maize   avg.m          11179.         8469.         6775.         5420.
#2 maize   avg.s           9918.         7514.         6011.         4809.
#3 maize   pca.m.2        10161.         7698.         6158.         4927.
#4 maize   pca.m.3        10393.         7874.         6299.         5039.
#5 soybean avg.m           3896.         2952.         2361.         1889.
#6 soybean avg.s           3267.         2475.         1980.         1584.
#7 soybean pca.m.2         3639.         2757.         2205.         1764.
#8 soybean pca.m.3         3735.         2830.         2264.         1811.

# -----------------
# Save 
# > Split by crop, model, and crop frequency and save predictions for optimization
preds_eu_wider <- prod_preds_eu %>% 
  split(.$crop) %>% 
  map_dfr(.,~{
    
    .x %>% 
      split(.$model) %>% 
      map_dfr(.,~{
        
        .x %>% 
          split(.$crop_frequency) %>% 
          map_dfr(., ~{
            
            dat_i <- .x 
            
            mod_i  <- unique(dat_i$model)
            crop_i <- unique(dat_i$crop)
            freq_i <- unique(dat_i$crop_frequency)
            
            tab_preds_eu <- dat_i %>% 
              dplyr::select(x, y, year, area, cropland_area_perc, cropland_area_ha, prod_Ya_pred) %>% 
              ## > compute mean and SD of production
              group_by(x, y) %>%
              mutate(mean_Ya_pred = mean(prod_Ya_pred),
                     sd_Ya_pred = sd(prod_Ya_pred)) %>%
              pivot_wider(names_from = year, values_from = prod_Ya_pred)
            
            # > save 
            write.csv2(tab_preds_eu, 
                       paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_preds/04_EU/", crop_i, "_pred_2000_2023_", mod_i, "_", freq_i,".csv"))
            # > out
            tab_preds_eu
            
          }, .id = "crop_frequency") 
        
      }, .id = "model")
      
    }, .id="crop")


# > Save the raw predictions 
Ya_pred_eu <- preds_eu

save(Ya_pred_eu, file = "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_preds/04_EU/Ya_pred_eu_2000_2023.rda")


# -----------------
# Prediction for Alberto

# > Estimation of the quantity produced
prod_preds_eu_Alberto <- preds_eu %>% 
  # > Constraint production surface to min(0.2*cropland; 0.01*area)
  mutate(area_1 = 0.01*area,
         cropland_area_20 = 0.2*cropland_area_ha) %>% 
  group_by(x, y, country_name, id_eu27, id_ext, id_eu_ext, year, model, crop) %>% 
  mutate(soybean_area = min(c(area_1, cropland_area_20))) %>% 
  # > Soybean production 
  mutate(prod_Ya_pred = soybean_area*Ya_pred) 

# > Production (Mt) of soybean and maize in each scenario
prod_preds_eu_Alberto %>% 
  group_by(crop, model) %>% 
  summarise(total_prod_Ya_pred = sum(prod_Ya_pred)/10^6) %>% 
  spread(key=crop, value=total_prod_Ya_pred)

#  1% du pixel 
#   model   maize soybean
# 1 avg.m   1105.    365.
# 2 avg.s    992.    305.
# 3 pca.m.2  994.    339.
# 4 pca.m.3 1027.    356.

# > Estimation of the quantity produced when 20% of the cropland is allocated
prod_preds_eu_Alberto_20 <- preds_eu %>% 
  # > Constraint production surface to 20% of the cropland is allocated
  mutate(soybean_area = 0.2*cropland_area_ha) %>% 
  # > Soybean production 
  mutate(prod_Ya_pred = soybean_area*Ya_pred) 

# > Production (Mt) of soybean and maize in each scenario
prod_preds_eu_Alberto_20 %>% 
  group_by(crop, model) %>% 
  summarise(total_prod_Ya_pred = sum(prod_Ya_pred)/10^6) %>% 
  spread(key=crop, value=total_prod_Ya_pred)

# 20% du pixel
#  model   maize soybean
# 1 avg.m   6775.   2361.
# 2 avg.s   6011.   1980.
# 3 pca.m.2 6158.   2205.
# 4 pca.m.3 6299.   2264.

# Save 
# > Split by crop, model, and crop frequency and save predictions for optimization
# > 1 file for continental EU
preds_eu_wider <- prod_preds_eu_Alberto %>% 
  split(.$crop) %>% 
  map_dfr(.,~{
    
    .x %>% 
      split(.$model) %>% 
      map_dfr(.,~{
        
        dat_i <- .x 
            
            mod_i  <- unique(dat_i$model)
            crop_i <- unique(dat_i$crop)
            
            tab_preds_eu <- dat_i %>% 
              dplyr::select(x, y, year, soybean_area, prod_Ya_pred) %>% 
              ## > compute mean and SD of production
              group_by(x, y) %>%
              mutate(mean_Ya_pred = mean(prod_Ya_pred),
                     sd_Ya_pred = sd(prod_Ya_pred)) %>%
              pivot_wider(names_from = year, values_from = prod_Ya_pred)
            
            # > save 
            write.csv2(tab_preds_eu, 
                       paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_preds/04_EU/", crop_i, "_pred_2000_2023_", mod_i, "_new_1perc.csv"))
            
            # > out
            tab_preds_eu
        
      }, .id = "model")
    
  }, .id="crop")

# > 1 for EU27
preds_eu27_wider <- prod_preds_eu_Alberto %>% 
  split(.$crop) %>% 
  map_dfr(.,~{
    
    .x %>% 
      split(.$model) %>% 
      map_dfr(.,~{
        
        dat_i <- .x 
        
        mod_i  <- unique(dat_i$model)
        crop_i <- unique(dat_i$crop)
        
        tab_preds_eu <- dat_i %>% 
          filter(id_eu27==1) %>%
          dplyr::select(x, y, year, soybean_area, prod_Ya_pred) %>% 
          ## > compute mean and SD of production
          group_by(x, y) %>%
          mutate(mean_Ya_pred = mean(prod_Ya_pred),
                 sd_Ya_pred = sd(prod_Ya_pred)) %>%
          pivot_wider(names_from = year, values_from = prod_Ya_pred)
        
        # > save 
        write.csv2(tab_preds_eu, 
                   paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_preds/04_EU/", crop_i, "_pred_2000_2023_", mod_i, "_new_1perc_eu27.csv"))
        
        # > out
        tab_preds_eu
        
      }, .id = "model")
    
  }, .id="crop")

# save
# 20%
preds_eu_wider20 <- prod_preds_eu_Alberto_20 %>% 
  split(.$crop) %>% 
  map_dfr(.,~{
    
    .x %>% 
      split(.$model) %>% 
      map_dfr(.,~{
        
        dat_i <- .x 
        
        mod_i  <- unique(dat_i$model)
        crop_i <- unique(dat_i$crop)
        
        tab_preds_eu <- dat_i %>% 
          dplyr::select(x, y, year, soybean_area, prod_Ya_pred) %>% 
          ## > compute mean and SD of production
          group_by(x, y) %>%
          mutate(mean_Ya_pred = mean(prod_Ya_pred),
                 sd_Ya_pred = sd(prod_Ya_pred)) %>%
          pivot_wider(names_from = year, values_from = prod_Ya_pred)
        
        # > save 
        write.csv2(tab_preds_eu, 
                   paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_preds/04_EU/", crop_i, "_pred_2000_2023_", mod_i, "_new_20perc.csv"))
        
        # > out
        tab_preds_eu
        
      }, .id = "model")
    
  }, .id="crop")

# > 1 for EU27
preds_eu27_wider20 <- prod_preds_eu_Alberto_20 %>% 
  split(.$crop) %>% 
  map_dfr(.,~{
    
    .x %>% 
      split(.$model) %>% 
      map_dfr(.,~{
        
        dat_i <- .x 
        
        mod_i  <- unique(dat_i$model)
        crop_i <- unique(dat_i$crop)
        
        tab_preds_eu <- dat_i %>% 
          filter(id_eu27==1) %>%
          dplyr::select(x, y, year, soybean_area, prod_Ya_pred) %>% 
          ## > compute mean and SD of production
          group_by(x, y) %>%
          mutate(mean_Ya_pred = mean(prod_Ya_pred),
                 sd_Ya_pred = sd(prod_Ya_pred)) %>%
          pivot_wider(names_from = year, values_from = prod_Ya_pred)
        
        # > save 
        write.csv2(tab_preds_eu, 
                   paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_preds/04_EU/", crop_i, "_pred_2000_2023_", mod_i, "_new_20perc_eu27.csv"))
        
        # > out
        tab_preds_eu
        
      }, .id = "model")
    
  }, .id="crop")


