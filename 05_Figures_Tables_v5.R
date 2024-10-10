# -------------------------------------------------------------------------
#
#       04 - Graphics and tables
#       Author: M. Chen, Inrae, 2024
#         
# -------------------------------------------------------------------------

# ----------------------------------------
# Packages & tools
library(tidyverse)
library(stringr)
library(lubridate)
library(terra) ; library(rnaturalearth)
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

# ----------------------------------------
# Functions 

# > Function for allocation
#source("E:/POSTDOC INRAE/ANALYSES/B_OPTIMISATION/00_0_Function_for_allocation_2.R")

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

# > Home-made function to format correctly the results from allocation
format_alloc <- function(results_allocations)
{
  
  formatted_result <- results_allocations %>% 
    # format
    separate(col = "scenario", into = c("pLER_s", "pLER_m", "freq_crop", "target_soybean", "max_surf"), sep = "_", remove = F) %>% 
    mutate(pLER_s    = as.numeric(as.character(pLER_s)),
           pLER_m    = as.numeric(as.character(pLER_m)),
           freq_crop = as.numeric(as.character(freq_crop)),
           target_soybean = as.numeric(as.character(target_soybean))/10^6,
           max_surf       = as.numeric(as.character(max_surf))/10^6) %>% 
    ungroup() %>%
    # Format and labels 
    # > pLERS
    mutate(pLER_lab = paste0(pLER_s, " - ", pLER_m)) %>% 
    # > Crop frequencies
    mutate(freq_crop_lab = recode(freq_crop, 
                                  "0.14"="1 year in 7", 
                                  "0.16"="1 year in 6", 
                                  "0.20"="1 year in 5", 
                                  "0.25"="1 year in 4", 
                                  "0.33"="1 year in 3", 
                                  "0.5"="1 year in 2"),
           freq_crop_lab = factor(freq_crop_lab, levels = c("1 year in 7", "1 year in 6", "1 year in 5", 
                                                            "1 year in 4", "1 year in 3", "1 year in 2"))) %>% 
    # > Production target
    mutate(target_soybean_lab = recode(target_soybean, 
                                       "9.075"="25%", 
                                       "18.150"="50%", 
                                       "27.225"="75%", 
                                       "36.300"="100%"),
           target_soybean_lab = factor(target_soybean_lab, c("25%", "50%", "75%", "100%"))) %>% 
    # > Maximum surface allocated
    mutate(max_surf_lab = "25% EU cropland area")
  
  return(formatted_result)
  
}

# ----------------------------------------
# Data

# PC
# > Predicted yields - Europe
load("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_preds/04_EU/Ya_pred_eu_2000_2023.rda")
# > EU coordinates
load("E:/POSTDOC INRAE/ANALYSES/B_OPTIMISATION/00_Data/00_dat_coords_EU.rda")
# > MAIN ANALYSIS : Results of simulations
#   of soybean and maize allocations in Europe
#   in EU27, on pixels with mean yield > 1 t/ha
#   (previously computed using the 03_Allocation_europe_full_simulations.R script)
allocations <- CCMHr::loadRDa("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/08_allocations/allocations_with_max_surf/allocations_soybean_maize_eu_restricted_min_rdt_check.rda")
#load("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/08_allocations/allocations_with_max_surf/allocations_soybean_maize_eu_restricted.rda")

# CIRAD
# > Predicted yields - Europe
load("D:/Mes Donnees/POSTDOC INRAE/ANALYSES/02_Maize_Soybean_intercropping_Europe/00_DATA/Ya_pred_eu_2000_2023.rda")
# > EU coordinates
load("D:/Mes Donnees/POSTDOC INRAE/ANALYSES/02_Maize_Soybean_intercropping_Europe/00_DATA/00_dat_coords_EU.rda")
# > MAIN ANALYSIS : Results of simulations
load("D:/Mes Donnees/POSTDOC INRAE/ANALYSES/02_Maize_Soybean_intercropping_Europe/00_DATA/allocations_soybean_maize_eu_restricted_min_rdt_check.rda")
allocations <- sensi_allocations ; rm(sensi_allocations)

# ----------------------------------------
# ------------- MAIN ANALYSIS ------------ 
# ----------------------------------------

# > Shapping allocations results
# 1. Allocation by pixel 
res1 <- format_alloc(results_allocations =  allocations %>% map_dfr(., ~{.x$data_res}, .id="scenario"))

# 2. Surface required to cover different shares of 
#    EU's consumption of soybean 
res2 <- format_alloc(results_allocations =  allocations %>% map_dfr(., ~{.x$res1}, .id="scenario"))

# 3. Surface required to produce the same production (maize+soybean)
#    as intercropping 
res3 <- format_alloc(results_allocations =  allocations %>% map_dfr(., ~{.x$res3}, .id="scenario"))

# 4. Performance of intercropping vs sole crops 
res4 <- format_alloc(results_allocations =  allocations %>% map_dfr(., ~{.x$res2}, .id="scenario"))

# ----------------------------------------
# ----------------------------------------
# Main figures 

# --------------------------
# Figure 1. Yield projections in Europe 

# Breaks for geom_contour_fill, geom_contour, geom_text_contour 
breaks_plot <- c(0, seq(0.5, 3.5, by=0.5), seq(4, 10, by=1))
breaks_labels <- seq(0,10, by=1)

# Map of average yield predictions in Europe
mean_Ya_pred_eu <- Ya_pred_eu %>% 
  filter(model=="pca.m.2", id_eu27==1) %>%
  # > label according to crops
  mutate(crop = if_else(crop=="maize", "b.", "a.")) %>% 
  # > long format 
  gather(key=year, value=Ya_pred, starts_with("X2")) %>% 
  # > for each pixel, recompute yields (in t/ha)
  group_by(crop, x, y) %>% 
  mutate(Ya_pred_t_ha = Ya_pred/cropland_area_ha) %>%
  summarise(mean_Ya_pred = mean(Ya_pred, na.rm=T)) %>%
  mutate(mean_Ya_pred = ifelse(is.na(mean_Ya_pred)==T, 0, mean_Ya_pred))

# > plot
p1 <- ggplot() + 
  geom_sf(data=eu27, fill="grey94") +
  geom_contour_fill(data=mean_Ya_pred_eu,
                    aes(x=x, y=y, z=mean_Ya_pred), 
                    breaks = breaks_plot,
                    na.fill = T, 
                    global.breaks = F,
                    clip = eu27) +
  geom_sf(data=eu27_ext, fill="transparent") +
  geom_contour(aes(x=x, y=y, z=mean_Ya_pred), 
               color = "white", 
               linewidth = 0.1,
               breaks = breaks_plot) +
  geom_text_contour(aes(x=x, y=y, z=mean_Ya_pred), 
                    color="black",
                    size=3, 
                    breaks = breaks_labels) +
  facet_grid(.~crop) + 
  theme_map() + 
  lims(x = c(-11,35), y=c(33,71)) + 
  #lims(x = c(-10,35), y=c(33,70)) + 
  theme(legend.position = "bottom",, 
        legend.title = element_text(size=12), 
        legend.text = element_text(size=10),
        strip.text = element_text(face = "bold", hjust = 0.1, size=15)) +
  scale_fill_gradientn(colours = c("transparent", 
                                   #viridis::plasma(n=100)[50:100], viridis::viridis(direction=-1, n=100)[1:75]),
                                   viridis::viridis(direction=-1, n=100)),
                       breaks = breaks_labels, 
                       labels = breaks_labels,
                       guide = guide_colorbar(barwidth = 20, barheight = 0.5, title.position = "top", title = expression("Mean yield t " ~ ha^-1))) ; p1

ggsave(p1, 
       filename = "E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/01_pred_maps.png", 
       height=18, width=22, bg="white", units = "cm", dpi = 300)

ggsave(p1, 
       filename = "E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/01_pred_maps.pdf", 
       height=18, width=22, bg="white", units = "cm", dpi = 300)

# Yield distribution 
mean_Ya_pred_eu %>% 
  mutate(crop=recode(crop,"a."="soybean","b."="maize")) %>% 
  spread(key=crop, value=mean_Ya_pred) %>%
  ungroup() %>% dplyr::select(-x,-y) %>% 
  summary(.)
#     maize          soybean     
#Min.   :0.223   Min.   :0.113  
#1st Qu.:4.305   1st Qu.:1.032  
#Median :6.220   Median :2.080  
#Mean   :5.261   Mean   :1.803  
#3rd Qu.:6.986   3rd Qu.:2.570  
#Max.   :9.002   Max.   :4.164

quantiles_soybean <- quantile(mean_Ya_pred_eu[which(mean_Ya_pred_eu$crop=="a."),]$mean_Ya_pred, probs=c(0.25,0.5,0.75)) ; quantiles_soybean
#      25%      50%      75% 
# 1.032215 2.080199 2.569741 
quantiles_maize <- quantile(mean_Ya_pred_eu[which(mean_Ya_pred_eu$crop=="b."),]$mean_Ya_pred, probs=c(0.25,0.5,0.75)) ; quantiles_maize
#      25%      50%      75% 
# 4.304674 6.220117 6.986286

# Area with given level of productivity for each crop
mean_Ya_pred_eu_cat <- mean_Ya_pred_eu %>% 
  mutate(crop=recode(crop,"a."="soybean","b."="maize")) %>% 
  spread(crop, mean_Ya_pred) %>% 
  # > soybean yield categories
  mutate(mean_Ya_soybean_cat_1   = if_else(soybean >= quantiles_soybean[1], 1, 0),
         mean_Ya_soybean_cat_2   = if_else(soybean >= quantiles_soybean[2], 1, 0),
         mean_Ya_soybean_cat_3   = if_else(soybean >= quantiles_soybean[3], 1, 0)) %>% 
  # > maize yield categories
  mutate(mean_Ya_maize_cat_1   = if_else(maize >= quantiles_maize[1], 1, 0),
         mean_Ya_maize_cat_2   = if_else(maize >= quantiles_maize[2], 1, 0),
         mean_Ya_maize_cat_3   = if_else(maize >= quantiles_maize[3], 1, 0)) %>%
  left_join(Ya_pred_eu %>% 
              dplyr::select(x,y,cropland_area_ha) %>% 
              distinct(), by=c("x","y"))

# total per cat
# > soybean
mean_Ya_pred_eu_cat %>%
  dplyr::select(-starts_with("mean_Ya_maize")) %>%
  gather(key=mean_Ya_soybean_cat, value=soybean_value, starts_with("mean_Ya_soybean")) %>%
  # > total area (106 Mha)
  group_by(mean_Ya_soybean_cat) %>% 
  mutate(total_cropland = sum(cropland_area_ha)) %>%
  ungroup() %>% 
  filter(soybean_value==1) %>%
  group_by(mean_Ya_soybean_cat, total_cropland) %>% 
  reframe(n_pixel=n(),
          cropland_area_sum = sum(cropland_area_ha)/10^6) %>% 
  mutate(cropland_area_prop = 100*(cropland_area_sum/(total_cropland/10^6))) 
#   mean_Ya_soybean_cat   total_cropland n_pixel cropland_area_sum cropland_area_prop
# 1 mean_Ya_soybean_cat_1     105629027.    2024             103.                97.9
# 2 mean_Ya_soybean_cat_2     105629027.    1350              81.8               77.5
# 3 mean_Ya_soybean_cat_3     105629027.     675              36.2               34.3

# > maize
mean_Ya_pred_eu_cat %>%
  dplyr::select(-starts_with("mean_Ya_soybean")) %>%
  gather(key=mean_Ya_maize_cat, value=maize_value, starts_with("mean_Ya_maize")) %>% 
  # > total area (106 Mha)
  group_by(mean_Ya_maize_cat) %>% 
  mutate(total_cropland = sum(cropland_area_ha)) %>%
  ungroup() %>% 
  filter(maize_value==1) %>%
  group_by(mean_Ya_maize_cat, total_cropland) %>% 
  reframe(n_pixel=n(),
          cropland_area_sum = sum(cropland_area_ha)/10^6) %>% 
  mutate(cropland_area_prop = 100*(cropland_area_sum/(total_cropland/10^6))) 
#   mean_Ya_maize_cat   total_cropland n_pixel cropland_area_sum cropland_area_prop
# 1 mean_Ya_maize_cat_1     105629027.    2024              99.1               93.8
# 2 mean_Ya_maize_cat_2     105629027.    1350              61.7               58.4
# 3 mean_Ya_maize_cat_3     105629027.     675              26.7               25.3

# cross tables
mean_Ya_pred_eu_cat %>% 
  gather(key=mean_Ya_soybean_cat, value=soybean_value, starts_with("mean_Ya_soybean"))%>%
  gather(key=mean_Ya_maize_cat, value=maize_value, starts_with("mean_Ya_maize")) %>%
  filter(soybean_value==1) %>% 
  group_by(mean_Ya_soybean_cat, mean_Ya_maize_cat) %>% 
  summarise(cropland_area_sum = sum(cropland_area_ha*(maize_value))/10^6) %>%
  spread(key=mean_Ya_maize_cat, value=cropland_area_sum)

#   mean_Ya_soybean_cat   mean_Ya_maize_cat_1 mean_Ya_maize_cat_2 mean_Ya_maize_cat_3
# 1 mean_Ya_soybean_cat_1                98.9                61.7                26.7
# 2 mean_Ya_soybean_cat_2                80.8                56.1                25.6
# 3 mean_Ya_soybean_cat_3                36.2                24.1                17.2

# --------------------------
# Figure 2. Maximum levels of soybean self-sufficiency reached with intercropping 
# according to varying values of pLER of soybean and crop return frequency

# > is it possible to cover EU27 soybean 
# consumption using intercropping?

p2 <- res2 %>% 
  # > Keep results for the plot
  filter(target_soybean == 36.3,        # soybean target production = 100% consumption
         strategy == "Intercropping", # intercropping
         crop == "Soybean",           # soybean
         pLER_s %in% c(0.3, 0.4, 0.5, 0.56, 0.6, 0.7),
         pLER_m == 0.79               # pLER maize is set as reference value
  ) %>%
  # > Self-sufficiency coverage for each scenario 
  mutate(perc_eu_supply=(production/(36.3*10^6))*100) %>%
  # > Reference values as reference for the dotted lines
  mutate(pLER_lab = case_when(
    pLER_lab == "References values for pLERs"~0, 
    TRUE~1)) %>% 
  # > Plot
  ggplot(., aes(x = freq_crop_lab, 
                y = perc_eu_supply,
                group=pLER_s)) +
  geom_hline(yintercept = c(25, 50, 75, 100), 
             linetype=2, 
             color = "grey90") +
  geom_line(aes(color=as.factor(pLER_s)),
            linewidth=1.2) +
  geom_line(aes(alpha=as.factor(pLER_lab)), 
            linewidth=1, color="white", linetype=3) +
  geom_point(aes(color=as.factor(pLER_s), shape=as.factor(pLER_s)), 
             size=3) +
  # current level of self-sufficiency (cake+grain)
  geom_hline(yintercept = 16, color='black', lty=2) +
  annotate("text", x = 0.5, y = 11.5, hjust=0, 
           label = "Current level of soybean (cake + grains)\nself-sufficiency in the EU: 7.4%", 
           color="black", fontface = 'italic') +
  scale_color_manual(values = c(viridis::viridis(6)[4:6], 
                                viridis::inferno(direction = -1, 6)[2:4])) +
  scale_alpha_manual(values = c(0,1), guide=guide_none()) +
  scale_y_continuous(breaks = c(25,50,75,100)) +
  guides(color = guide_legend(title = "Partial land equivalent ratio for soybean", reverse = T, nrow=1),
         shape = guide_legend(title = "Partial land equivalent ratio for soybean", reverse = T, nrow=1)) +
  theme_cowplot() +
  theme(strip.text.y = element_text(angle=0, size=11),
        strip.text.x = element_text(size=11),
        axis.text = element_text(size=11),
        axis.title = element_text(size=11),
        #legend.position = c(-0.0,.8),
        legend.position = "bottom",
        legend.title = element_text(size=11),
        legend.text = element_text(size=10),
        legend.background = element_rect(fill="white")#,plot.margin = unit(c(0,3,0,0), "cm")
        ) +
  #lims(y=c(0,101)) +
  coord_cartesian(clip = "off") +
  labs(x = "\nReturn frequency", y = "Soybean self-sufficiency level (%)\n") ; p2

# Save
ggsave(p2, 
       filename = "E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/02_selfsufficiency_perc_v2.png", 
       width=22, height=18, bg="white", units = "cm", dpi=300)

ggsave(p2, 
       filename = "E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/02_selfsufficiency_perc_v2.pdf", 
       width=22, height=18, bg="white", units = "cm", dpi=300)

# Associated table
tab_p2 <- res2 %>% 
  # > Keep results for the plot
  filter(target_soybean_lab == "100%",        # soybean target production = 100% consumption
         strategy == "Intercropping", # intercropping
         crop == "Soybean",           # soybean
         pLER_s %in% c(0.3, 0.4, 0.5, 0.56, 0.6, 0.7),
         pLER_m == 0.79               # pLER maize is set as reference value
  ) %>%
  # > Self-sufficiency coverage for each scenario 
  mutate(perc_eu_supply=(production/(target_soybean*10^6))*100,
         perc_eu_supply=if_else(perc_eu_supply>100, 100, perc_eu_supply), 
         surface = surface/10^6,
         production=production/10^6) %>%
  dplyr::select(freq_crop_lab, pLER_s, surface, production, perc_eu_supply) %>%
  gather(key=var, value=val, -freq_crop_lab, -pLER_s) %>%
  mutate(val=round(val,1)) %>%
  spread(key=freq_crop_lab, value=val) %>%
  arrange(var, pLER_s)

save(tab_p2, file = "E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/02_selfsufficiency_perc_v2.rda")

# Try with another type of presentation
# PC
load("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/08_allocations/allocations_with_max_surf/sensi/sensi_4_res1.rda")
# CIRAD
#load("D:/Mes Donnees/POSTDOC INRAE/ANALYSES/02_Maize_Soybean_intercropping_Europe/00_DATA/sensi/sensi_4.rda")

sensi_4_res1%>% 
  # > Keep results for the plot
  filter(target_soybean == 36.3,        # soybean target production = 100% consumption
         strategy == "Intercropping", # intercropping
         pLER_s %in% c(0.3, 0.4, 0.5, 0.56, 0.6, 0.7),
         pLER_m %in% c(0.5, 0.6, 0.7, 0.79, 0.8, 0.9)
  ) %>%
  mutate(pLER=if_else(crop=="Soybean", pLER_s, pLER_m)) %>% 
  mutate(keep=case_when(crop=="Soybean" & pLER_m==0.79~1,
                        crop=="Maize"&pLER_s==0.56~1,
                        TRUE~0)) %>%
  filter(keep==1) %>%
  # > Self-sufficiency coverage for each scenario 
  mutate(target_soybean = if_else(crop=="Soybean", target_soybean, 85.1), 
         perc_eu_supply=(production/(target_soybean*10^6))*100) %>% 
  mutate(crop=recode(crop,"Soybean"="a. Soybean", "Maize"="b. Maize")) %>% 
  mutate(freq_crop_lab=recode(freq_crop_lab, 
                              "1 year in 7"="one-in-seven",
                              "1 year in 6"="one-in-six",
                              "1 year in 5"="one-in-five",
                              "1 year in 4"="one-in-four",
                              "1 year in 2"="one-in-three",
                              "1 year in 3"="one-in-two")) %>%
  # > Reference values as reference for the dotted lines
  mutate(pLER_lab = case_when(
    pLER_lab == "0.56 - 0.79"~1, 
    TRUE~0)) %>%
  # > Plot
  ggplot(., aes(x = freq_crop_lab, 
                y = perc_eu_supply)) +
  # current level of self-sufficiency 
  geom_hline(aes(yintercept = target_soybean), color='black', lty=2, linewidth=1) +
  # refs
  geom_hline(yintercept = c(25,50,75,100,125,150), 
             linetype=2, 
             color = "grey90") +
  geom_path(aes(color=as.factor(pLER), group=interaction(crop,pLER)), 
            linewidth=1) +
  geom_point(aes(color=as.factor(pLER), shape=as.factor(pLER_lab)),
             size=2) +
  scale_color_manual(values = c(viridis::viridis(6), 
                                viridis::inferno(direction = -1, 7)[2:7])) +
  scale_shape_manual(values=c(3,15)) +
  scale_y_continuous(breaks = c(25,50,75,100,125,150)) +
  guides(color = guide_legend(title = "Partial land equivalent ratio (pLER):", nrow=1),
         shape = guide_none()) +
  theme_cowplot() +
  theme(strip.text.y = element_text(angle=0, size=11),
        strip.text.x = element_text(size=11, face = "bold", hjust = 0),
        strip.background = element_blank(),
        axis.text.x = element_text(size=10, angle=45, 
                                   vjust = 0.85, hjust = 0.9),
        axis.text.y = element_text(size=10),
        axis.title = element_text(size=11),
        panel.border = element_rect(color="black", linewidth = 0.1),
        axis.line = element_line(color="black", linewidth = 0.1),
        #legend.position = c(-0.0,.8),
        legend.position = "bottom",
        legend.title = element_text(size=11),
        legend.title.position = "top",
        legend.text = element_text(size=10),
        legend.background = element_rect(fill="white")#,plot.margin = unit(c(0,3,0,0), "cm")
  ) +
  #lims(y=c(0,101)) +
  facet_grid(.~crop) +
  coord_cartesian(clip = "off") +
  labs(x = "\nReturn frequency", y = "Self-sufficiency level (%)\n")

ggsave(filename = "E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/02_selfsufficiency_perc_v2_maize.png", 
       width=22, height=18, bg="white", units = "cm", dpi=300)


# -------------------------
# Figure 3 - Area requirements to reach the same production
# of maize and soybean in intercropping and sole cropping

# > Max surface = 1/4 of EU27 cropland ~ 25 Mha
num_eu_cropland25<-25

# Maps - maize and soybean allocation in EU when grown in 
# intercropping or as sole crops
alloc_maps <- res1 %>% 
  # > Keep the data for the figure (reference values)
  filter(pLER_s == 0.56,
         pLER_m == 0.79,          
         freq_crop == 0.25) %>% 
  split(.$target_soybean) %>% 
  map(., ~{
    
    tab_i <- .x %>% 
      filter(pixel_reach_equal_prod==1) %>% 
      # > pixel where maize as sole crop could be grown
      mutate(pixel_maize_solecropping = if_else(pixel_solecropping == 0 & pixel_landsaving == 0, 1, 0)) %>%
      # > keep usefull columns 
      dplyr::select(x, y, target_soybean_lab, pixel_intercropping, "pixel_soybean_solecropping"="pixel_solecropping", pixel_maize_solecropping, pixel_landsaving) %>% 
      gather(key = pixel_type, value = pixel_color, starts_with("pixel_")) %>%
      mutate(Crop_Design_lab = case_when(
        pixel_type == "pixel_intercropping"         ~ "Int", 
        pixel_type == "pixel_soybean_solecropping"  ~ "Soybean as SC", 
        pixel_type == "pixel_maize_solecropping"    ~ "Maize as SC", 
        pixel_type == "pixel_landsaving"            ~ "Additional area required as maize in SC\nto reach the same co-production as Int")) %>% 
      mutate(Crop_Design_lab = factor(Crop_Design_lab, levels = c("Int", 
                                                                  "Soybean as SC", 
                                                                  "Maize as SC",
                                                                  "Additional area required as maize in SC\nto reach the same co-production as Int"))) %>%
      filter(Crop_Design_lab %in% c("Int", 
                                    "Soybean as SC", 
                                    "Maize as SC",
                                    "Additional area required as maize in SC\nto reach the same co-production as Int")) %>% 
      mutate(Crop_Design_lab2 = if_else(Crop_Design_lab == "Int", "Intercropping (Int)", "Sole crops (SC)")) %>% 
      mutate(target_soybean_lab = paste0("Soybean\nself-sufficiency\nlevel: ", target_soybean_lab))
    
    # Maps surface
    ggplot() +
      #geom_text(data = tab_i, x = -5, y = 68, aes(label = surfaces_lab), size= 3) +
      geom_sf(data=eu27, fill="grey94", color="transparent") +
      geom_tile(data = tab_i, 
                aes(x=x, 
                    y=y, 
                    fill=Crop_Design_lab, 
                    alpha=as.factor(pixel_color))) +
      scale_fill_manual(values = c("#3CBC75FF", "#2D718EFF", "#FDE725FF", "purple"), 
                        guide=guide_legend("", nrow = 1, position = "bottom")) +
      scale_alpha_manual(values = c(0,1), guide=guide_none()) +
      geom_sf(data = not_eu27, color = NA, fill = "white",size = 0.2) +
      geom_sf(data = ocean, color = NA, fill = "white",size = 0.2) +
      geom_sf(data = europe, fill="transparent") +
      theme_map() + 
      theme(strip.text.y.left = element_text(angle=0, size=25),
            strip.background = element_rect(fill="lightgrey", color="transparent"),
            strip.text.x = element_text(size=25),
            legend.text = element_text(size=25)) +
      facet_grid(target_soybean_lab ~ Crop_Design_lab2, switch = "y") +
      lims(x = c(-11,35), y=c(33,70)) 
    
  }) ; alloc_maps[[2]]


# Maps - maize and soybean allocation in EU when grown in 
# intercropping or as sole crops 
# surface displayed per pixel (rather than presence/absence)
alloc_maps_surf_pixel <- res1 %>%
  # > Keep the data for the figure (reference values)
  filter(pLER_s == 0.56,
         pLER_m == 0.79,          
         freq_crop == 0.25) %>% 
  split(.$target_soybean_lab) %>% 
  map(., ~{
    
    tab_i <- .x %>% 
      filter(pixel_reach_equal_prod==1) %>% 
      # > pixel where maize as sole crop could be grown
      mutate(pixel_maize_solecropping = if_else(pixel_solecropping == 0 & pixel_landsaving == 0, 1, 0)) %>%
      # > keep useful columns 
      dplyr::select(x, y, area, target_soybean_lab, pixel_intercropping, "pixel_soybean_solecropping"="pixel_solecropping", pixel_maize_solecropping, pixel_landsaving) %>% 
      gather(key = pixel_type, value = pixel_color, starts_with("pixel_")) %>%
      mutate(Crop_Design_lab = case_when(
        pixel_type == "pixel_intercropping"         ~ "Maize-soybean intercropping", 
        pixel_type == "pixel_soybean_solecropping"  ~ "Soybean as sole crop", 
        pixel_type == "pixel_maize_solecropping"    ~ "Maize as sole crop", 
        pixel_type == "pixel_landsaving"            ~ "Land saved by intercropping")) %>% 
      mutate(Crop_Design_lab = factor(Crop_Design_lab, levels = c("Maize-soybean intercropping",
                                                                  "Soybean as sole crop", 
                                                                  "Maize as sole crop", 
                                                                  "Land saved by intercropping"))) %>%
      mutate(Crop_Design_lab2 = if_else(Crop_Design_lab == "Maize-soybean intercropping", "Intercropping", "Sole crops")) %>% 
      mutate(target_soybean_lab = paste0("Soybean self-sufficiency level: ", target_soybean_lab)) %>% 
      filter(pixel_color==1)
    
    # Maps surface
    p <- plot_grid(
      # > Intercropping
      ggplot() +
        geom_sf(data=eu27, fill="grey35", color="transparent") +
        geom_tile(data = tab_i %>% filter(pixel_type == "pixel_intercropping"), 
                  aes(x=x, y=y, fill = pixel_color*area/10^3)) +
        scale_fill_gradientn(colours = c("white", "#0D0887FF"), guide=guide_colorbar("Area (1000 ha)", nrow = 1, position = "bottom", barwidth=15, barheight=0.7), limits = c(0,180)) +
        geom_sf(data = not_eu27, color = NA, fill = "white",size = 0.2) + geom_sf(data = ocean, color = NA, fill = "white",size = 0.2) + geom_sf(data = europe, fill="transparent", color="white") +
        theme_map() + 
        theme(strip.background = element_rect(fill="lightgrey", color="transparent"),
              strip.text.x = element_text(size=25),
              legend.text = element_text(size=25),
              legend.title = element_text(size=25),
              legend.title.position = "top",
              title = element_text(size=25)) +
        facet_grid(. ~ Crop_Design_lab, switch = "y") +
        lims(x = c(-11,35), y=c(33,70)) +
        #ggtitle(paste0(unique(tab_i$target_soybean_lab))) +
        ggtitle("a. Geographical allocation"),
      # > Sole soybean
      ggplot() +
        geom_sf(data=eu27, fill="grey35", color="transparent") +
        geom_tile(data = tab_i %>% filter(pixel_type == "pixel_soybean_solecropping"), 
                  aes(x=x, y=y, fill = pixel_color*area/10^3)) +
        scale_fill_gradientn(colours = c("white", "#8405A7FF"), guide=guide_colorbar("Area (1000 ha)", nrow = 1, position = "bottom", barwidth=15, barheight=0.7), limits = c(0,180)) +
        geom_sf(data = not_eu27, color = NA, fill = "white",size = 0.2) + geom_sf(data = ocean, color = NA, fill = "white",size = 0.2) + geom_sf(data = europe, fill="transparent", color="white") +
        theme_map() + 
        theme(strip.background = element_rect(fill="lightgrey", color="transparent"),
              strip.text.x = element_text(size=25),
              legend.text = element_text(size=25),
              legend.title = element_text(size=25),
              legend.title.position = "top") +
        facet_grid(. ~ Crop_Design_lab, switch = "y") +
        lims(x = c(-11,35), y=c(33,70)), 
      # > Sole maize
      ggplot() +
        geom_sf(data=eu27, fill="grey35", color="transparent") +
        geom_tile(data = tab_i %>% filter(pixel_type == "pixel_maize_solecropping"), 
                  aes(x=x, y=y, fill = pixel_color*area/10^3)) +
        scale_fill_gradientn(colours = c("white", "#DF6263FF"), guide=guide_colorbar("Area (1000 ha)", nrow = 1, position = "bottom", barwidth=15, barheight=0.7), limits = c(0,180)) +
        geom_sf(data = not_eu27, color = NA, fill = "white",size = 0.2) + geom_sf(data = ocean, color = NA, fill = "white",size = 0.2) + geom_sf(data = europe, fill="transparent", color="white") +
        theme_map() + 
        theme(strip.background = element_rect(fill="lightgrey", color="transparent"),
              strip.text.x = element_text(size=25),
              legend.text = element_text(size=25),
              legend.title = element_text(size=25),
              legend.title.position = "top") +
        facet_grid(. ~ Crop_Design_lab, switch = "y") +
        lims(x = c(-11,35), y=c(33,70)),
      # > Land saving
      ggplot() +
        geom_sf(data=eu27, fill="grey35", color="transparent") +
        geom_tile(data = tab_i %>% filter(pixel_type == "pixel_landsaving"), 
                  aes(x=x, y=y, fill = pixel_color*area/10^3)) +
        scale_fill_gradientn(colours = c("white", "#F0F921FF"), guide=guide_colorbar("Area (1000 ha)", nrow = 1, position = "bottom", barwidth=15, barheight=0.7), limits = c(0,180)) +
        geom_sf(data = not_eu27, color = NA, fill = "white",size = 0.2) + geom_sf(data = ocean, color = NA, fill = "white",size = 0.2) + geom_sf(data = europe, fill="transparent", color="white") +
        theme_map() + 
        theme(strip.background = element_rect(fill="lightgrey", color="transparent"),
              strip.text.x = element_text(size=25),
              legend.text = element_text(size=25),
              legend.title = element_text(size=25), legend.title.position = "top") +
        facet_grid(. ~ Crop_Design_lab, switch = "y") +
        lims(x = c(-11,35), y=c(33,70)),
      nrow=1, align="hv"
      
    ) 
    
  }) ; alloc_maps_surf_pixel[["50%"]]


# Total surface 
alloc_surf <- res3 %>%
  filter(pLER_s == 0.56,
         pLER_m == 0.79,          
         #target_soybean %in% c(22.50, 33.75), 
         freq_crop == 0.25) %>% 
  split(.$target_soybean_lab) %>%
  map(., ~{
    
    .x %>% 
      # > keep only 1 line for intercropping (same surface for both crop)
      mutate(to_keep = if_else(strategy=="intercrop" & crop =="crop 2", 0, 1)) %>% 
      filter(to_keep == 1) %>% 
      # > re-label crop and strategy
      mutate(crop = if_else(strategy=="intercrop", "Maize-soybean intercropping", crop),
             crop = recode(crop, "landsaving"="Land saved by intercropping", "crop 1"="Soybean as sole crop", "crop 2" = "Maize as sole crop")) %>% 
      mutate(strategy = recode(strategy, "intercrop"="Int", "sole crop"="SC")) %>% 
      ggplot(data=.) +
      # > total surface covered by each strategy
      geom_col(aes(x=strategy , y = total_surface/10^6, fill=crop),
               col="transparent", width=0.75) +
      # > total surface available for production (25% of European cropland)
      geom_hline(yintercept=num_eu_cropland25, 
                 color="black", linetype = 2, lwd=2) +
      annotate("text", x = 2.7, y = num_eu_cropland25, 
               label = "25% of\ncroplands\nin the EU", hjust=0, size=9,
               color="black", fontface = 'italic') + 
      coord_cartesian(clip = "off", xlim = c(1, 2)) + 
      theme_cowplot() + 
      theme(strip.background = element_blank(),
            strip.text = element_blank(),
            legend.position = "bottom",
            legend.text = element_text(size=25),
            panel.border = element_rect(color="black"),
            title = element_text(size=25),
            axis.text = element_text(size=25),
            axis.title = element_text(size=25),
            plot.margin = unit(c(0,4,0,0), "cm")) +
      facet_wrap(.~target_soybean_lab, ncol=1) +
      scale_fill_manual(values = c("#0D0887FF", "#8405A7FF", "#DF6263FF","#F0F921FF"), 
                        limits=c("Maize-soybean intercropping", "Soybean as sole crop", "Maize as sole crop", "Land saved by intercropping"), 
                        guide=guide_legend("", nrow = 1, position = "bottom")) +
      ggtitle("b. Area") +
      labs(x = "", y = "Total (Mha)") +
      lims(y=c(0,27))
    
  }) ; alloc_surf$`50%`

# Total co-production of maize and soybean in each strategy
alloc_prod <- res3 %>%
  filter(pLER_s == 0.56,
         pLER_m == 0.79,          
         #target_soybean %in% c(22.50, 33.75), 
         freq_crop == 0.25) %>% 
  split(.$target_soybean_lab) %>%
  map(., ~{
    
    # Total production
    .x %>% 
      # > re-label crop and strategy 
      mutate(crop2 = crop,
             crop2 = factor(crop2, levels = c("intercrop", "landsaving", "crop 2", "crop 1"))) %>% 
      mutate(crop = if_else(strategy=="intercrop", "Int.", crop),
             crop = factor(crop, levels = c("Int.", "landsaving", "crop 2", "crop 1"))) %>% 
      mutate(strategy = recode(strategy, "intercrop"="Int", "sole crop"="SC")) %>% 
      # > compute total production for intercropping
      group_by(scenario, target_soybean_lab, target_soybean, crop, crop2, strategy) %>% 
      summarize(sum=sum(total_production)) %>% 
      ggplot(data=.) +
      # > production by crop and by strategy
      geom_col(aes(x=strategy , y = sum/10^6, fill=crop2, col=crop),
               width=0.75, linewidth=1) + 
      # > production target for soybean 
      geom_hline(aes(yintercept=target_soybean), color="black", linetype = 2, lwd=2) + 
      annotate("text", x = 2.7, y = as.numeric(as.character(unique(.x$target_soybean))), hjust=0, size=9,
               label = paste0(unique(.x$target_soybean_lab), "\nsoybean\nself-sufficiency"), 
               color="black", fontface = 'italic') + 
      coord_cartesian(clip = "off", xlim = c(1, 2)) + 
      theme_cowplot() + 
      theme(strip.background = element_blank(),
            strip.text = element_blank(),
            legend.position = "none",
            panel.border = element_rect(color="black"),
            title = element_text(size=25),
            axis.text = element_text(size=25),
            axis.title = element_text(size=25),
            plot.margin = unit(c(0,4,0,0), "cm")) +
      facet_wrap(.~target_soybean_lab, scales = "free", ncol=1) +
      scale_color_manual(values = c("#0D0887FF", "#F0F921FF", "#DF6263FF","#8405A7FF"), 
                        guide=guide_legend("", ncol = 2, position = "bottom")) +
      scale_fill_manual(values = rev(c("#8405A7FF", "#DF6263FF","#F0F921FF")), 
                        guide=guide_legend("", ncol = 2, position = "bottom")) +
      ggtitle("c. Production") +
      labs(x = "", y = "Total (Mt)") +
      lims(y=c(0,170))
    
  }) ; alloc_prod$`50%`


# Merge plots
#leg_p4 <- get_legend(plot = p4_maps$'22.5')

leg_p3 = cowplot::get_plot_component(alloc_surf[[2]], 'guide-box-bottom', return_all = TRUE)

p3 <- plot_grid(plot_grid(alloc_maps_surf_pixel[[2]], 
                          ggplot() + theme_void(),
                          plot_grid(ggplot() + theme_void(),
                                    alloc_surf[[2]] + theme(legend.position = "none"),   
                                    ggplot() + theme_void(),
                                    alloc_prod[[2]] + theme(legend.position = "none"), 
                                    ggplot() + theme_void(),
                                    nrow=1, rel_widths = c(0.17, 0.25,0.15,0.25,0.18)),  
                          ncol=1, align="hv", rel_heights = c(0.425,0.05,0.425),
                          axis = "btrl"), 
                #cowplot::ggdraw(leg_p3),
                plot_grid(ggplot() + theme_void(), cowplot::ggdraw(leg_p3), ggplot() + theme_void(), nrow=1, rel_widths = c(0.1,0.75,0.15)),
                ncol=1, rel_heights = c(0.94,0.06), align="hv", axis = "btrl"
                )
# Save
ggsave(p3,
       filename = "E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/03_landsaving.png", 
       width=22, height=18, bg="white", dpi=300)

ggsave(p3,
       filename = "E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/03_landsaving.pdf", 
       width=18, height=22, bg="white", dpi=300)

# Automatic saving figures for the 3 self-sufficiency levels
path_fig <- "D:/Mes Donnees/POSTDOC INRAE/PAPERS/02_PNAS/02_Figures/"

for(i in unique(res1$target_soybean_lab))
{
  
  # legend 
  leg_p3 = cowplot::get_plot_component(alloc_surf[[paste0(i)]], 'guide-box-bottom', return_all = TRUE)
  # plot
  p3 <- plot_grid(plot_grid(alloc_maps_surf_pixel[[paste0(i)]], 
                            ggplot() + theme_void(),
                            plot_grid(ggplot() + theme_void(),
                                      alloc_surf[[paste0(i)]] + theme(legend.position = "none"),   
                                      ggplot() + theme_void(),
                                      alloc_prod[[paste0(i)]] + theme(legend.position = "none"), 
                                      ggplot() + theme_void(),
                                      nrow=1, rel_widths = c(0.17, 0.25,0.15,0.25,0.18)),  
                            ncol=1, align="hv", rel_heights = c(0.425,0.05,0.425),
                            axis = "btrl"), 
                  #cowplot::ggdraw(leg_p3),
                  plot_grid(ggplot() + theme_void(), cowplot::ggdraw(leg_p3), ggplot() + theme_void(), nrow=1, rel_widths = c(0.1,0.75,0.15)),
                  ncol=1, rel_heights = c(0.94,0.06), align="hv", axis = "btrl"
  )
  # save
  ggsave(p3,
         filename = paste0(path_fig, "03_landsaving_", gsub("%", "", i), ".png"), 
         width=22, height=18, bg="white", dpi=300)
  
}




p4 <- plot_grid(
  plot_grid(p4_surf_barplots[[1]] + annotate("text", x = 2.2, y = num_eu_cropland25-4, 
                                             label = "25% of\ncroplands\nin the EU", hjust=1,
                                             color="darkorange", fontface = 'italic', size=7),  
            p4_prod_barplots[[1]]  + 
              annotate("text", x = 2.7, y = 11.25+2, hjust=0, size=7,
                       label = "Soybean\nself-\nsufficiency\ntarget", 
                       color="red", fontface = 'italic') + 
              coord_cartesian(clip = "off", xlim = c(1, 2)),  
            nrow=1, align="hv", axis = "btrl"),
  plot_grid(p4_surf_barplots[[2]],  p4_prod_barplots[[2]],  nrow=1, align="hv", axis = "btrl"),
  plot_grid(p4_surf_barplots[[3]],  p4_prod_barplots[[3]],  nrow=1, align="hv", axis = "btrl"), 
  plot_grid(p4_surf_barplots[[4]],  p4_prod_barplots[[4]],  nrow=1, align="hv", axis = "btrl"), 
  plot_grid(cowplot::ggdraw(leg_p4), ggplot() + theme_void(), nrow=1),
  ncol=1, align="hv", axis = "btrl", #labels = c("a.", "b.", "c.", ""), 
  rel_heights = c(0.18, 0.18, 0.18, 0.18, 0.06)
) ; p4

# Save
ggsave(p4,
       filename = "E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/04_landsaving_v2_V2_2.png", 
       width=18, height=22, bg="white", dpi=300)

# Associated table 
tab_p4 <- res3 %>%
  filter(pLER_s == 0.56,
         pLER_m == 0.79, 
         freq_crop == 0.25) %>% 
  dplyr::select(crop, strategy, target_soybean_lab, total_production, total_surface) %>% 
  mutate(total_production = total_production/10^6,
         total_surface = total_surface/10^6) %>%
  unite("lab", crop:strategy) %>% 
  gather(key=metric, value=value, total_production, total_surface) %>%
  spread(key=lab, value=value) %>% 
  arrange(metric, target_soybean_lab) %>%
  dplyr::select(metric,target_soybean_lab, `crop 1_intercrop`, `crop 1_sole crop`, `crop 2_intercrop`, `crop 2_sole crop`, `landsaving_sole crop`) %>%
  mutate(TOTAL_Intercropping = `crop 1_intercrop`+ `crop 2_intercrop`, 
         TOTAL_Sole = `crop 1_sole crop` + `crop 2_sole crop` + `landsaving_sole crop`) %>%
  mutate(diff = TOTAL_Intercropping - TOTAL_Sole)

save(tab_p4, file="E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/04_landsaving_v2.rda")

# ------------------------------------------------------#

# ------------------------------------------------------#
# -------------   {Sensitivity analyses}   -------------#
# ------------------------------------------------------#

# ------------------------------------------------------#

# -------------------------
# PC
path_alloc <- "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/08_allocations/allocations_with_max_surf/"

# CIRAD
path_alloc <- "D:/Mes Donnees/POSTDOC INRAE/ANALYSES/02_Maize_Soybean_intercropping_Europe/00_DATA/"

# -------------------------
# > sensitivity analysis with different yield threshold 

# > Soybean allocations based on yields projeted from various models 
load(paste0(path_alloc, "sensi/sensi_1.rda"))
sensi_1 <- sensi_allocations ; rm(sensi_allocations)

load(paste0(path_alloc, "sensi/sensi_2_1.rda"))
sensi_2_1 <- sensi_allocations ; rm(sensi_allocations)

load(paste0(path_alloc, "sensi/sensi_2_2.rda"))
sensi_2_2 <- sensi_allocations ; rm(sensi_allocations)

# > Allocations and area requierements
alloc_sensi_1 <- list("Areas with soybean mean yield equal or higher than 1.0 t.ha-1 (EU)"          = allocations,
                      "Areas with soybean mean yield equal or higher than 2.8 t.ha-1 (EU)"          = sensi_1,
                      "Areas with soybean mean yield equal or higher than 1.0 t.ha-1 (EU extended)" = sensi_2_1,
                      "Areas with soybean mean yield equal or higher than 2.8 t.ha-1 (EU extended)" = sensi_2_2)

# > Shapping allocations results
# Allocations
sensi_1_data_res <- format_alloc(results_allocations =  alloc_sensi_1 %>% 
                                   map_dfr(., ~{ 
                                     
                                     map_dfr(.x, ~{ .x$data_res }, .id="scenario")
                                     
                                   }, .id="sensi_set")) %>% 
  # > Set of sensitivity analyses
  mutate(sensi_set = recode(sensi_set, 
                            "Areas with soybean mean yield equal or higher than 1.0 t.ha-1 (EU)"          = "Mean yield\nequal or higher\nthan 1.0 t.ha-1\n(EU)"         ,
                            "Areas with soybean mean yield equal or higher than 2.8 t.ha-1 (EU)"          = "Mean yield\nequal or higher\nthan 2.8 t.ha-1\n(EU)"         ,
                            "Areas with soybean mean yield equal or higher than 1.0 t.ha-1 (EU extended)" = "Mean yield\nequal or higher\nthan 1.0 t.ha-1\n(EU extended)",
                            "Areas with soybean mean yield equal or higher than 2.8 t.ha-1 (EU extended)" = "Mean yield\nequal or higher\nthan 2.8 t.ha-1\n(EU extended)" )) %>%
  mutate(sensi_set = factor(sensi_set, levels=c("Mean yield\nequal or higher\nthan 1.0 t.ha-1\n(EU)"         ,
                                                "Mean yield\nequal or higher\nthan 2.8 t.ha-1\n(EU)"         ,
                                                "Mean yield\nequal or higher\nthan 1.0 t.ha-1\n(EU extended)",
                                                "Mean yield\nequal or higher\nthan 2.8 t.ha-1\n(EU extended)" )))

# Surface required to cover different shares of 
#    EU's consumption of soybean 
sensi_2_res1 <- format_alloc(results_allocations =  alloc_sensi_1 %>% 
                               map_dfr(., ~{ 
                                 
                                 map_dfr(.x, ~{ .x$res1 }, .id="scenario")
                                 
                               }, .id="sensi_set")) %>%
  # > Set of sensitivity analyses
  mutate(sensi_set_lab = recode(sensi_set, 
                            "Areas with soybean mean yield equal or higher than 1.0 t.ha-1 (EU)"          = "Mean yield\nequal or higher\nthan 1.0 t.ha-1\n(EU)"         ,
                            "Areas with soybean mean yield equal or higher than 2.8 t.ha-1 (EU)"          = "Mean yield\nequal or higher\nthan 2.8 t.ha-1\n(EU)"         ,
                            "Areas with soybean mean yield equal or higher than 1.0 t.ha-1 (EU extended)" = "Mean yield\nequal or higher\nthan 1.0 t.ha-1\n(EU extended)",
                            "Areas with soybean mean yield equal or higher than 2.8 t.ha-1 (EU extended)" = "Mean yield\nequal or higher\nthan 2.8 t.ha-1\n(EU extended)" )) %>%
  mutate(sensi_set_lab = factor(sensi_set_lab, levels=c("Mean yield\nequal or higher\nthan 1.0 t.ha-1\n(EU)"         ,
                                                "Mean yield\nequal or higher\nthan 2.8 t.ha-1\n(EU)"         ,
                                                "Mean yield\nequal or higher\nthan 1.0 t.ha-1\n(EU extended)",
                                                "Mean yield\nequal or higher\nthan 2.8 t.ha-1\n(EU extended)" )))

# Surface required to produce the same production (maize+soybean)
#    as intercropping 
sensi_2_res3 <- format_alloc(results_allocations =  alloc_sensi_1 %>% 
                               map_dfr(., ~{ 
                                 
                                 map_dfr(.x, ~{ .x$res3 }, .id="scenario")
                                 
                               }, .id="sensi_set")) %>%
  mutate(sensi_set = recode(sensi_set, 
                            "Areas with soybean mean yield equal or higher than 1.0 t.ha-1 (EU)"          = "Mean yield\nequal or higher\nthan 1.0 t.ha-1\n(EU)"         ,
                            "Areas with soybean mean yield equal or higher than 2.8 t.ha-1 (EU)"          = "Mean yield\nequal or higher\nthan 2.8 t.ha-1\n(EU)"         ,
                            "Areas with soybean mean yield equal or higher than 1.0 t.ha-1 (EU extended)" = "Mean yield\nequal or higher\nthan 1.0 t.ha-1\n(EU extended)",
                            "Areas with soybean mean yield equal or higher than 2.8 t.ha-1 (EU extended)" = "Mean yield\nequal or higher\nthan 2.8 t.ha-1\n(EU extended)" )) %>%
  mutate(sensi_set = factor(sensi_set, levels=c("Mean yield\nequal or higher\nthan 1.0 t.ha-1\n(EU)"         ,
                                                "Mean yield\nequal or higher\nthan 2.8 t.ha-1\n(EU)"         ,
                                                "Mean yield\nequal or higher\nthan 1.0 t.ha-1\n(EU extended)",
                                                "Mean yield\nequal or higher\nthan 2.8 t.ha-1\n(EU extended)" )))

# > Compare with main analysis
sensi_2_res1 %>% 
  filter(pLER_s==0.56, pLER_m==0.79, freq_crop == 0.25, target_soybean_lab == "75%") %>% 
  dplyr::select(sensi_set, crop, strategy, surface) %>% 
  spread(key = strategy, value = surface)

sensi_2_res3 %>% 
  filter(pLER_s==0.56, pLER_m==0.79, freq_crop == 0.25, target_soybean_lab == "75%") %>% 
  dplyr::select(sensi_set, crop, strategy, total_production) %>% 
  spread(key = strategy, value = total_production)

sensi_2_res3 %>% 
  filter(pLER_s==0.56, pLER_m==0.79, freq_crop == 0.25, target_soybean_lab == "75%") %>% 
  dplyr::select(sensi_set, crop, strategy, total_surface) %>% 
  spread(key = strategy, value = total_surface)

# > % of self-sufficiency reached in each scenario
#   depending on pLER and crop frequency
sensi_2_res1 %>% 
  # > Keep results for the plot
  filter(target_soybean_lab == "100%",        # soybean target production = 100% consumption
         strategy == "Intercropping", # intercropping
         crop == "Soybean",           # soybean
         pLER_s %in% c(0.3, 0.4, 0.5, 0.56, 0.6, 0.7),
         pLER_m == 0.79               # pLER maize is set as reference value
  ) %>%
  # > Self-sufficiency coverage for each scenario 
  mutate(perc_eu_supply=(production/(36.3*10^6))*100) %>%
  # > Reference values as reference for the dotted lines
  mutate(pLER_lab = case_when(
    pLER_lab == "References values for pLERs"~0, 
    TRUE~1)) %>% 
  mutate(pLER_s=paste0("Soybean pLER = ", pLER_s)) %>%
  # > Plot
  ggplot(., aes(x = freq_crop_lab, 
                y = perc_eu_supply,
                group=sensi_set)) +
  # current level of self-sufficiency (cake+grain)
  geom_hline(yintercept = 16, color='red', lty=2) +
  geom_hline(yintercept = c(25, 50, 75, 100), 
             linetype=2, 
             color = "grey90") +
  geom_line(aes(color=as.factor(sensi_set)), #linetype = 2,
            linewidth=0.75) +
  geom_point(aes(color=as.factor(sensi_set), shape=as.factor(sensi_set)), 
             size=2) +
  scale_color_manual(values = pal[c(1,3,4,5)])+
  scale_linetype_manual(values=c(2,1), guide=guide_none()) +
  scale_y_continuous(breaks=c(25,50,75,100))+
  guides(color = guide_legend(title = "", ncol=1),
         shape=guide_legend(title = "", ncol=1)) +
  theme_cowplot() +
  theme(strip.text.y = element_text(angle=0, size=11),
        strip.text.x = element_text(size=11),
        axis.text.y = element_text(size=11), 
        axis.text.x = element_text(size=11, 
                                   angle=45, 
                                   vjust = 0.85, hjust = 0.9),
        axis.title = element_text(size=11),
        legend.position = "bottom",
        legend.title = element_text(size=11),
        legend.text = element_text(size=10)) +
  #lims(y=c(0,101)) +
  labs(x = "\nReturn frequency", y = "Soybean self-sufficiency level (%)\n") +
  facet_wrap(.~pLER_s, nrow=1) 

# other presentation
sensi_2_res1 %>% 
  # > Keep results for the plot
  filter(target_soybean_lab == "100%",        # soybean target production = 100% consumption
         strategy == "Intercropping", # intercropping
         crop == "Soybean") %>%
  # > Self-sufficiency coverage for each scenario 
  mutate(perc_eu_supply=(production/(36.3*10^6))*100) %>%
  # > Plot
  ggplot(., aes(x = freq_crop_lab, 
                y = perc_eu_supply,
                group=pLER_s)) +
  # current level of self-sufficiency (cake+grain)
  geom_hline(yintercept = 16, color='red', lty=2) +
  geom_hline(yintercept = c(25, 50, 75, 100), 
             linetype=2, 
             color = "grey90") +
  geom_line(aes(color=as.factor(pLER_s)), #linetype = 2,
            linewidth=0.75) +
  geom_point(aes(color=as.factor(pLER_s), shape=as.factor(pLER_s)), 
             size=2) +
  scale_color_manual(values = c(viridis::viridis(6)[4:6], 
                                viridis::inferno(direction = -1, 6)[2:4])) +
  scale_linetype_manual(values=c(2,1), guide=guide_none()) +
  scale_y_continuous(breaks=c(25,50,75,100))+
  guides(color = guide_legend(title = "", ncol=1),
         shape=guide_legend(title = "", ncol=1)) +
  theme_cowplot() +
  theme(strip.text.y = element_text(angle=0, size=11),
        strip.text.x = element_text(size=11),
        axis.text.y = element_text(size=11), 
        axis.text.x = element_text(size=11, 
                                   angle=45, 
                                   vjust = 0.85, hjust = 0.9),
        axis.title = element_text(size=11),
        #legend.position = "bottom",
        legend.title = element_text(size=11),
        legend.text = element_text(size=10)) +
  #lims(y=c(0,101)) +
  labs(x = "\nReturn frequency", y = "Soybean self-sufficiency level (%)\n") +
  facet_wrap(.~sensi_set_lab, nrow=1) 

sensi_2_res1 %>% 
  # > Keep results for the plot
  filter(target_soybean_lab == "100%",        # soybean target production = 100% consumption
         strategy == "Intercropping", # intercropping
         crop == "Soybean") %>%
  # > Self-sufficiency coverage for each scenario 
  mutate(perc_eu_supply=(production/(36.3*10^6))*100) %>% 
  filter(pLER_s==0.56, freq_crop==0.25) %>% 
  dplyr::select(sensi_set, production, surface, perc_eu_supply)

#                                                                    sensi_set production  surface perc_eu_supply
#1          Areas with soybean mean yield equal or higher than 1.0 t.ha-1 (EU)   33131294 25001535       91.27078
#2          Areas with soybean mean yield equal or higher than 2.8 t.ha-1 (EU)    5436812  3147557       14.97744
#3 Areas with soybean mean yield equal or higher than 1.0 t.ha-1 (EU extended)   32839809 25004410       90.46779
#4 Areas with soybean mean yield equal or higher than 2.8 t.ha-1 (EU extended)    6473288  3755520       17.83275






# -------------------------
# > sensitivity analyses with different yield predictive models 

# > Projected yields 
# Breaks for geom_contour_fill, geom_contour, geom_text_contour 
breaks_plot <- c(0, seq(0.5, 3.5, by=0.5), seq(4, 10, by=1))
breaks_labels <- seq(0,10, by=1)

# Map of average yield predictions in continental Europe
Ya_pred_eu %>% 
  filter(id_eu27==1) %>%
  # > label according to crops
  mutate(crop = if_else(crop=="maize", "Maize", "Soybean")) %>% 
  mutate(crop = factor(crop, levels=c("Soybean", "Maize"))) %>% 
  mutate(model = recode(model,
                        "avg.m"  ="Model incorportaing\nmonthly averages",  
                        "avg.s"  ="Model incorporating\nseasonal averages",  
                        "pca.m.2"="Model incorportaing\ntwo first principal components\nof each climate variable\n(main analysis)",
                        "pca.m.3"="Model incorportaing\nthree first principal components\nof each climate variable" ) ) %>%
  mutate(model = factor(model, levels=c("Model incorportaing\ntwo first principal components\nof each climate variable\n(main analysis)",
                                        "Model incorportaing\nmonthly averages",  
                                        "Model incorporating\nseasonal averages",  
                                        "Model incorportaing\nthree first principal components\nof each climate variable"))) %>%
  # > long format 
  gather(key=year, value=Ya_pred, starts_with("X2")) %>% 
  # > for each pixel, recompute yields (in t/ha)
  group_by(model, crop, x, y) %>% 
  mutate(Ya_pred_t_ha = Ya_pred/cropland_area_ha) %>%
  summarise(mean_Ya_pred = mean(Ya_pred, na.rm=T)) %>%
  mutate(mean_Ya_pred = ifelse(is.na(mean_Ya_pred)==T, 0, mean_Ya_pred)) %>% 
  # > plot
  ggplot(.) + 
  geom_sf(data=eu27, fill="grey94") +
  geom_contour_fill(aes(x=x, y=y, z=mean_Ya_pred), 
                    breaks = breaks_plot,
                    na.fill = T, 
                    global.breaks = F,
                    clip = eu27) +
  geom_sf(data=eu27_ext, fill="transparent") +
  geom_contour(aes(x=x, y=y, z=mean_Ya_pred), 
               color = "white", 
               linewidth = 0.05,
               breaks = breaks_plot) +
  #geom_text_contour(aes(x=x, y=y, z=mean_Ya_pred), 
  #                  color="black",
  #                  size=1.5, 
  #                  breaks = breaks_labels) +
  facet_grid(crop~model, switch = "y") + 
  theme_map() + 
  lims(x = c(-11,51), y=c(33,71)) +
  theme(legend.position = "bottom",
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        strip.text.x = element_text(hjust = 0.1, face = "bold"),
        strip.text.y.left = element_text(hjust = 0.1, angle = 0, face = "bold")) +
  scale_fill_gradientn(colours = c("transparent", viridis::plasma(n=100)[50:100], viridis::viridis(direction=-1, n=100)[1:75]),
                       breaks = breaks_labels, 
                       labels = breaks_labels,
                       guide = guide_colorbar(barwidth = 20, barheight = 0.5, title.position = "left", title = expression("Mean yield t " ~ ha^-1)))

ggsave(filename = "E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/Sensitivity_analyses/sensi_1_pred_maps_eu_models.png", 
       height=8, width=14, bg="white", dpi = 300)

# > Soybean allocations based on yields projeted from various models 
sensi_1_avg_m   <- loadRDa(paste0(path_alloc, "sensi/sensi_1_avg_m.rda"))
sensi_1_avg_s   <- loadRDa(paste0(path_alloc, "sensi/sensi_1_avg_s.rda"))
sensi_1_pca_m_3 <- loadRDa(paste0(path_alloc, "sensi/sensi_1_pca_m_3.rda"))

# > Allocations and area requierements
alloc_sensi_1 <- list("Main analysis: model incoporating\nthe two first principal components\nfor each climate variable (N=2047)"            = allocations,
                      "Sensitivity analysis 1: model incoporating\nmonthly averages (N=1984)"                                                = sensi_1_avg_m,
                      "Sensitivity analysis 2: model incoporating\nseasonal averages (N=1778)"                                               = sensi_1_avg_s,
                      "Sensitivity analysis 3: model incoporating\nthe three first principal components\nfor each climate variable (N=2042)" = sensi_1_pca_m_3)


# > Shapping allocations results
# Allocations
sensi_1_data_res <- format_alloc(results_allocations =  alloc_sensi_1 %>% 
                               map_dfr(., ~{ 
                                 
                                 map_dfr(.x, ~{ .x$data_res }, .id="scenario")
                                 
                               }, .id="sensi_set")) %>% 
  # > Set of sensitivity analyses
  mutate(sensi_set = recode(sensi_set, 
                            "Main analysis: model incoporating\nthe two first principal components\nfor each climate variable (N=2047)"           ="Model incoporating\nthe two first\nprincipal components\nfor each climate\nvariable\n(main analysis)",
                            "Sensitivity analysis 1: model incoporating\nmonthly averages (N=1984)"                                               ="Model incoporating\nmonthly averages",
                            "Sensitivity analysis 2: model incoporating\nseasonal averages (N=1778)"                                              ="Model incoporating\nseasonal averages",
                            "Sensitivity analysis 3: model incoporating\nthe three first principal components\nfor each climate variable (N=2042)"="Model incoporating\nthe three first\nprincipal components\nfor each climate\nvariable")) %>%
  mutate(sensi_set = factor(sensi_set, levels=c("Model incoporating\nthe two first\nprincipal components\nfor each climate\nvariable\n(main analysis)",
                                                "Model incoporating\nmonthly averages",
                                                "Model incoporating\nseasonal averages",
                                                "Model incoporating\nthe three first\nprincipal components\nfor each climate\nvariable")))

# Surface required to cover different shares of 
#    EU's consumption of soybean 
sensi_1_res1 <- format_alloc(results_allocations =  alloc_sensi_1 %>% 
                               map_dfr(., ~{ 
                                 
                                 map_dfr(.x, ~{ .x$res1 }, .id="scenario")
                                 
                               }, .id="sensi_set")) %>%
  # > Set of sensitivity analyses
  mutate(sensi_set = recode(sensi_set, 
                            "Main analysis: model incoporating\nthe two first principal components\nfor each climate variable (N=2047)"           ="Main\nanalysis\n",
                            "Sensitivity analysis 1: model incoporating\nmonthly averages (N=1984)"                                               ="Sensitivity\nanalysis\n1",
                            "Sensitivity analysis 2: model incoporating\nseasonal averages (N=1778)"                                              ="Sensitivity\nanalysis\n2",
                            "Sensitivity analysis 3: model incoporating\nthe three first principal components\nfor each climate variable (N=2042)"="Sensitivity\nanalysis\n3"))

# Surface required to produce the same production (maize+soybean)
#    as intercropping 
sensi_1_res3 <- format_alloc(results_allocations =  alloc_sensi_1 %>% 
                       map_dfr(., ~{ 
                         
                         map_dfr(.x, ~{ .x$res3 }, .id="scenario")
                         
                       }, .id="sensi_set")) %>%
  # > Set of sensitivity analyses
  mutate(sensi_set = recode(sensi_set, 
                            "Main analysis: model incoporating\nthe two first principal components\nfor each climate variable (N=2047)"           ="Main\nanalysis\n",
                            "Sensitivity analysis 1: model incoporating\nmonthly averages (N=1984)"                                               ="Sensitivity\nanalysis\n1",
                            "Sensitivity analysis 2: model incoporating\nseasonal averages (N=1778)"                                              ="Sensitivity\nanalysis\n2",
                            "Sensitivity analysis 3: model incoporating\nthe three first principal components\nfor each climate variable (N=2042)"="Sensitivity\nanalysis\n3"))


# > Compare with main analysis
sensi_1_res1 %>% 
  filter(pLER_s==0.56, pLER_m==0.79, freq_crop == 0.25, target_soybean_lab == "75%") %>% 
  dplyr::select(sensi_set, crop, strategy, surface) %>% 
  spread(key = strategy, value = surface)

sensi_1_res3 %>% 
  filter(pLER_s==0.56, pLER_m==0.79, freq_crop == 0.25, target_soybean == 33.75) %>% 
  dplyr::select(sensi_set, crop, strategy, total_production) %>% 
  spread(key = strategy, value = total_production)

sensi2_res3 %>% 
  filter(pLER_s==0.56, pLER_m==0.79, freq_crop == 0.25, target_soybean == 33.75) %>% 
  dplyr::select(sensi_set, crop, strategy, total_surface) %>% 
  spread(key = strategy, value = total_surface)

# > % of self-sufficiency reached 
sensi_1_res1 %>% 
  # > Keep results for the plot
  filter(target_soybean_lab == "100%",        # soybean target production = 100% consumption
         strategy == "Intercropping", # intercropping
         crop == "Soybean",           # soybean
         pLER_s %in% c(0.3, 0.4, 0.5, 0.56, 0.6, 0.7),
         pLER_m == 0.79               # pLER maize is set as reference value
  ) %>%
  # > Self-sufficiency coverage for each scenario 
  mutate(perc_eu_supply=(production/(58.4*10^6))*100) %>%
  # > Reference values as reference for the dotted lines
  mutate(pLER_lab = case_when(
    pLER_lab == "References values for pLERs"~0, 
    TRUE~1)) %>% 
  mutate(pLER_s=paste0("Soybean pLER = ", pLER_s)) %>%
  # > Set of sensitivity analyses
  mutate(sensi_set = recode(sensi_set, 
                            "Main\nanalysis\n"="Main analysis: Model incoporating the two first principal components for each climate variable",
                            "Sensitivity\nanalysis\n1"="Model incoporating monthly averages",
                            "Sensitivity\nanalysis\n2"="Model incoporating seasonal averages",
                            "Sensitivity\nanalysis\n3"="Model incoporating the three first principal components for each climate variable")) %>%
  mutate(sensi_set = factor(sensi_set, levels=c("Main analysis: Model incoporating the two first principal components for each climate variable",
                                                "Model incoporating monthly averages",
                                                "Model incoporating seasonal averages",
                                                "Model incoporating the three first principal components for each climate variable"))) %>%  
  # > Plot
  ggplot(., aes(x = freq_crop_lab, 
                y = perc_eu_supply,
                group=sensi_set)) +
  # current level of self-sufficiency (cake+grain)
  geom_hline(yintercept = 42.3, color='red', lty=2) +
  geom_hline(yintercept = c(25, 50, 75, 100), 
             linetype=2, 
             color = "grey90") +
  geom_line(aes(color=as.factor(sensi_set)), #linetype = 2,
            linewidth=0.75) +
  geom_point(aes(color=as.factor(sensi_set), shape=as.factor(sensi_set)), 
             size=2) +
  scale_color_manual(values = c(viridis::viridis(4)[2:3], viridis::inferno(5)[3:4])) +
  scale_linetype_manual(values=c(2,1), guide=guide_none()) +
  scale_y_continuous(breaks=c(25,50,75,100))+
  guides(color = guide_legend(title = "", ncol=1),
         shape=guide_legend(title = "", ncol=1)) +
  theme_cowplot() +
  theme(strip.text.y = element_text(angle=0, size=11),
        strip.text.x = element_text(size=11),
        axis.text.y = element_text(size=11), 
        axis.text.x = element_text(size=11, 
                                   angle=45, 
                                   vjust = 0.85, hjust = 0.9),
        axis.title = element_text(size=11),
        legend.position = "bottom",
        legend.title = element_text(size=11),
        legend.text = element_text(size=10)) +
  #lims(y=c(0,101)) +
  labs(x = "\nReturn frequency", y = "Soybean self-sufficiency level (%)\n") +
  facet_wrap(.~pLER_s, nrow=1) 

# Save
ggsave(filename = "E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/Sensitivity_analyses/sensi_1_selfsufficiency_perc.png", 
       width=14, height=6, bg="white", dpi=300)

# > Maps of allocation 
sensi_1_data_res %>% 
  # > Keep the data for the figure (reference values)
  filter(pLER_s == 0.56,
         pLER_m == 0.79, 
         freq_crop == 0.25) %>% 
  split(.$target_soybean_lab) %>% 
  map(., ~{
    
    tab_i <- .x %>% 
      # > keep usefull columns 
      dplyr::select(x, y, target_soybean_lab, pixel_intercropping, pixel_solecropping, pixel_intercropping_free, sensi_set) %>% 
      gather(key = pixel_type, value = pixel_color, starts_with("pixel_")) %>%
      mutate(Crop_Design_lab = case_when(
        pixel_type == "pixel_intercropping"      ~ "Intercropping", 
        pixel_type == "pixel_solecropping"       ~ "Soybean as sole crop", 
        pixel_type == "pixel_intercropping_free" ~ "Maize as sole crop")) %>% 
      mutate(Crop_Design_lab = factor(Crop_Design_lab, levels = c("Intercropping", 
                                                                  "Soybean as sole crop", 
                                                                  "Maize as sole crop"))) %>%
      filter(Crop_Design_lab %in% c("Intercropping", 
                                    "Soybean as sole crop", 
                                    "Maize as sole crop")) %>% 
      mutate(Crop_Design_lab2 = if_else(Crop_Design_lab == "Intercropping", "Intercropping", "Sole crops"))
    
    # Maps surface
    plot_i <- ggplot() +
      geom_sf(data=eu27, fill="grey94", color="transparent") +
      geom_tile(data = tab_i, 
                aes(x=x, 
                    y=y, 
                    fill=Crop_Design_lab, 
                    alpha=as.factor(pixel_color))) +
      scale_fill_manual(values = c("#3CBC75FF", "#2D718EFF", "#FDE725FF"), 
                        guide=guide_legend("", nrow = 1)) +
      scale_alpha_manual(values = c(0,1), guide=guide_none()) +
      geom_sf(data = ocean, color = NA, fill = "white",size = 0.2) +
      geom_sf(data = not_eu27, fill = "transparent",size = 0.2) +
      geom_sf(data = europe, fill="transparent") +
      theme_map() + 
      theme(strip.text.y.left = element_text(angle=0),
            strip.background = element_rect(fill="lightgrey", color="transparent"),
            legend.position = "bottom"
      ) +
      facet_grid(sensi_set ~ Crop_Design_lab2, switch = "y") +
      lims(x = c(-11,51), y=c(33,71)) +
      ggtitle(paste0("EU's soybean self-sufficiency level targeted: ", unique(.x$target_soybean_lab), " (", unique(.x$target_soybean), " Mt)"))
    
    # > Save
    ggsave(plot_i,
           filename = paste0("E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/Sensitivity_analyses/sensi_1_map", gsub("%", "", unique(.x$target_soybean_lab)), ".png"), 
           width=11, height=12, bg="white", dpi=300)
    
  })

# > Barplots - area requirements and co-production of maize and soybean in each strategy
sensi_1_res1 %>%
  filter(pLER_s == 0.56,
         pLER_m == 0.79, 
         freq_crop == 0.25)  %>% 
  split(.$target_soybean_lab) %>%
  map(., ~{
    
    # > Area requirements
    surf_barplots_i <- .x %>% 
      # > keep only 1 line for intercropping (same surface occupied by both crop)
      mutate(to_keep = if_else(strategy == "Sole crop" | strategy == "Intercropping" & crop == "Soybean", 1, 0)) %>% 
      filter(to_keep == 1) %>% 
      # > re-label crop and strategy
      mutate(crop = if_else(strategy == "Intercropping", "Intercropping", paste0(crop, " as sole crop"))) %>% 
      mutate(crop=factor(crop, levels =c("Intercropping", "Maize as sole crop", "Soybean as sole crop"))) %>%
      #mutate(strategy = recode(strategy, "Sole crop"="SC", "Intercropping"="Int")) %>%
      ggplot(data=.) +
      # > total surface intercropping and by crop in sole crops
      geom_col(aes(x=sensi_set , y = surface/10^6, fill=crop),
               col="transparent", width=0.75) +
      # > total surface available (when considering crop frequency = 0.25)
      geom_hline(aes(yintercept=max_surf), 
                 color="darkorange", linetype = 2) + 
      theme_cowplot() + 
      theme(legend.position = "bottom",
            panel.border = element_rect(color="black"),
            axis.title = element_text(size=11),
            axis.text = element_text(size=8.5)
      ) +
      scale_fill_manual(values = c("#3CBC75FF", "#FDE725FF", "#2D718EFF"), name="") +
      facet_wrap(strategy~., scales = "free", ncol=2) +
      ggtitle("") +
      labs(x = "", y = "Total area (Mha)") +
      lims(y=c(0,50))
    
    # > Associated co-productions
    prod_barplot_i <- .x %>%
      # > re-label crop and strategy 
      mutate(crop2 = crop,
             crop2 = factor(crop2, levels = c("Maize", "Soybean"))) %>%
      mutate(crop = if_else(strategy=="Intercropping", "Int.", crop),
             crop = factor(crop, levels = c("Int.", "Maize", "Soybean"))) %>% 
      #mutate(strategy = recode(strategy, "Intercropping"="Int", "Sole crop"="SC")) %>%  
      # > compute sum production especially for intercropping
      group_by(sensi_set, scenario, target_soybean_lab, target_soybean, crop, crop2, strategy) %>% 
      summarize(sum=sum(production)) %>%
      ggplot(data=.) +
      # > total production for each strategy
      geom_col(aes(x=sensi_set , y = sum/10^6, fill=crop2, col=crop), width=0.75, linewidth=1) +
      # > soybean production target
      geom_hline(aes(yintercept=target_soybean), color="red", linetype = 2) + 
      theme_cowplot() + 
      theme(legend.position = "none",
            panel.border = element_rect(color="black"),
            axis.title = element_text(size=11),
            axis.text = element_text(size=8.5)
      ) +
      facet_wrap(.~strategy, scales = "free", ncol=2) +
      #scale_fill_manual(values = c("#3CBC75FF","gray","gray20")) +
      scale_fill_manual(values = c("#FDE725FF","#2D718EFF")) +
      scale_color_manual(values = c("#3CBC75FF","#FDE725FF","#2D718EFF")) +
      ggtitle("") +
      labs(x = "", y = "Total co-production (Mt)") +
      lims(y=c(0,300))
    
    # > Save
    p <- plot_grid(
      plot_grid(surf_barplots_i + 
                  theme(legend.position = "none") +
                  ggtitle(paste0("EU's soybean self-sufficiency level targeted: ", unique(.x$target_soybean_lab)," (", unique(.x$target_soybean), " Mt)"),
                          "a. Area"),  
                prod_barplot_i + 
                  ggtitle(label = "", subtitle = "b. Production"),  nrow=2, align="hv", axis = "btrl"),
      cowplot::ggdraw(cowplot::get_plot_component(surf_barplots_i, 'guide-box-bottom', return_all = TRUE)), 
      ncol=1, align="hv", axis = "btrl", rel_heights = c(0.9, 0.1)
    )
    
    ggsave(p, 
           filename = paste0("E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/Sensitivity_analyses/sensi_1_barplots_", gsub("%", "", unique(.x$target_soybean_lab)), ".png"), 
           width=7, height=7, bg="white", dpi=300)
    
    
    
  })

# -------------------------
# > sensitivity analysis with different yield threshold

# Sets of pixels with productivity higher than a certain threshold
# > pixels with yields >1 t/ha
#sensi_2_pixels <- Ya_pred_eu %>% 
#  filter(model=="pca.m.2", crop=="soybean") %>% 
#  group_by(x, y, id_eu27) %>% 
#  summarise(mean_Ya_pred = mean(Ya_pred)) %>%
#  mutate(cat_pixel_main  = if_else(mean_Ya_pred >= 1 & id_eu27 == 1, 1, 0),
#         cat_pixel_sensi = if_else(mean_Ya_pred >= 2 & id_eu27 == 1, 1, 0))
#
#plot_grid(
#  ggplot() +
#    geom_sf(data = eu27, fill="grey94") +
#    geom_tile(data=sensi_2_pixels, 
#              aes(x=x,y=y, fill=as.factor(cat_pixel_main))) +
#    geom_sf(data = ocean, color = NA, fill = "white",size = 0.2) +
#    geom_sf(data = russia, fill="white", color="white") +
#    geom_sf(data = not_eu27, fill="white") +
#    geom_sf(data = eu27, fill="transparent") +
#    theme_map() + #theme(plot.title = element_text(size=10)) +
#    lims(x = c(-11,51), y=c(33,71)) +
#    scale_fill_manual(values = c("transparent", viridis::viridis(4)[2]), guide=guide_none()) + 
#    ggtitle("a. Main analysis:\ngrid-cells with soybean mean yield > 1 t ha-1 (N=2047)"),
#  ggplot() +
#    geom_sf(data=eu27, fill="grey94") +
#    geom_tile(data=sensi_2_pixels, 
#              aes(x=x,y=y, fill=as.factor(cat_pixel_sensi))) +
#    geom_sf(data = ocean, color = NA, fill = "white",size = 0.2) +
#    geom_sf(data = russia, fill="white", color="white") +
#    geom_sf(data = not_eu27, fill="white") +
#    geom_sf(data = eu27_ext, fill="transparent") +
#    theme_map() + #theme(plot.title = element_text(size=10)) +
#    lims(x = c(-11,51), y=c(33,71)) +
#    scale_fill_manual(values = c("transparent", viridis::inferno(5)[4]), guide=guide_none()) + 
#    ggtitle("b. Sensitivity analysis\ngrid-cells with soybean mean yield > 2 t ha-1 (N=1387)"),
#  nrow=1
#  
#)
#
#ggsave(filename = "E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/Sensitivity_analyses/sensi_2_pixels.png",
#       width = 12, height = 7, dpi=300, bg="white")

# > Soybean allocations based on yields projeted from various models 
load(paste0(path_alloc, "sensi/sensi_1.rda"))
sensi_2 <- sensi_allocations ; rm(sensi_allocations)

# > Allocations and area requierements
alloc_sensi_2 <- list("Areas with soybean mean yield equal or higher than 1.0 t ha-1" = allocations,
                      "Areas with soybean mean yield equal or higher than 2.8 t ha-1" = sensi_2)


# > Shapping allocations results
# Allocations
sensi_2_data_res <- format_alloc(results_allocations =  alloc_sensi_2 %>% 
                                   map_dfr(., ~{ 
                                     
                                     map_dfr(.x, ~{ .x$data_res }, .id="scenario")
                                     
                                   }, .id="sensi_set")) %>% 
  # > Set of sensitivity analyses
  mutate(sensi_set = recode(sensi_set, 
                            "Areas with soybean mean yield equal or higher than 1.0 t ha-1"="Areas with\nsoybean mean yield\n> 1.0 t ha-1",
                            "Areas with soybean mean yield equal or higher than 2.8 t ha-1"="Areas with\nsoybean mean yield\n> 2.8 t ha-1")) %>%
  mutate(sensi_set = factor(sensi_set, levels=c("Areas with\nsoybean mean yield\n> 1.0 t ha-1",
                                                "Areas with\nsoybean mean yield\n> 2.8 t ha-1")))

# Surface required to cover different shares of 
#    EU's consumption of soybean 
sensi_2_res1 <- format_alloc(results_allocations =  alloc_sensi_2 %>% 
                               map_dfr(., ~{ 
                                 
                                 map_dfr(.x, ~{ .x$res1 }, .id="scenario")
                                 
                               }, .id="sensi_set")) %>%
  # > Set of sensitivity analyses
  mutate(sensi_set = recode(sensi_set, 
                            "Areas with soybean mean yield equal or higher than 1.0 t ha-1" = "Areas with\nsoybean mean yield\n> 1.0 t ha-1",
                            "Areas with soybean mean yield equal or higher than 2.8 t ha-1" = "Areas with\nsoybean mean yield\n> 2.8 t ha-1")) %>%
  mutate(sensi_set = factor(sensi_set, levels=c("Areas with\nsoybean mean yield\n> 1.0 t ha-1",
                                                "Areas with\nsoybean mean yield\n> 2.8 t ha-1")))

# Surface required to produce the same production (maize+soybean)
#    as intercropping 
sensi_2_res3 <- format_alloc(results_allocations =  alloc_sensi_2 %>% 
                               map_dfr(., ~{ 
                                 
                                 map_dfr(.x, ~{ .x$res3 }, .id="scenario")
                                 
                               }, .id="sensi_set")) %>%
  mutate(sensi_set = recode(sensi_set, 
                            "Areas with soybean mean yield equal or higher than 1.0 t ha-1"="Areas with\nsoybean mean yield\n> 1.0 t ha-1",
                            "Areas with soybean mean yield equal or higher than 2.8 t ha-1"="Areas with\nsoybean mean yield\n> 2.8 t ha-1")) %>%
  mutate(sensi_set = factor(sensi_set, levels=c("Areas with\nsoybean mean yield\n> 1.0 t ha-1",
                                                "Areas with\nsoybean mean yield\n> 2.8 t ha-1")))

# > Compare with main analysis
sensi_2_res1 %>% 
  filter(pLER_s==0.56, pLER_m==0.79, freq_crop == 0.25, target_soybean_lab == "75%") %>% 
  dplyr::select(sensi_set, crop, strategy, surface) %>% 
  spread(key = strategy, value = surface)

sensi_2_res3 %>% 
  filter(pLER_s==0.56, pLER_m==0.79, freq_crop == 0.25, target_soybean_lab == "75%") %>% 
  dplyr::select(sensi_set, crop, strategy, total_production) %>% 
  spread(key = strategy, value = total_production)

sensi_2_res3 %>% 
  filter(pLER_s==0.56, pLER_m==0.79, freq_crop == 0.25, target_soybean_lab == "75%") %>% 
  dplyr::select(sensi_set, crop, strategy, total_surface) %>% 
  spread(key = strategy, value = total_surface)

# > % of self-sufficiency reached 
sensi_2_res1 %>% 
  # > Keep results for the plot
  filter(target_soybean_lab == "100%",        # soybean target production = 100% consumption
         strategy == "Intercropping", # intercropping
         crop == "Soybean",           # soybean
         pLER_s %in% c(0.3, 0.4, 0.5, 0.56, 0.6, 0.7),
         pLER_m == 0.79               # pLER maize is set as reference value
  ) %>%
  # > Self-sufficiency coverage for each scenario 
  mutate(perc_eu_supply=(production/(58.4*10^6))*100) %>%
  # > Reference values as reference for the dotted lines
  mutate(pLER_lab = case_when(
    pLER_lab == "References values for pLERs"~0, 
    TRUE~1)) %>% 
  mutate(pLER_s=paste0("Soybean pLER = ", pLER_s)) %>%
  # > Labels
  # > Set of sensitivity analyses
  mutate(sensi_set = recode(sensi_set, 
                            "Grid-cells with\nsoybean mean yield\n> 1 t ha-1\n(main analysis)" = "Main analysis: grid-cells with soybean mean yield > 1 t ha-1"       ,
                            "Grid-cells with\nsoybean mean yield\n> 2 t ha-1"                  = "Sensitivity analysis: grid-cells with soybean mean yield > 2 t ha-1")) %>%
  # > Plot
  ggplot(., aes(x = freq_crop_lab, 
                y = perc_eu_supply,
                group=sensi_set)) +
  # current level of self-sufficiency (cake+grain)
  geom_hline(yintercept = 16, color='red', lty=2) +
  geom_hline(yintercept = c(25, 50, 75, 100), 
             linetype=2, 
             color = "grey90") +
  geom_line(aes(color=as.factor(sensi_set)), #linetype = 2,
            linewidth=0.75) +
  geom_point(aes(color=as.factor(sensi_set), shape=as.factor(sensi_set)), 
             size=2) +
  scale_color_manual(values = c(viridis::viridis(4)[2], viridis::inferno(5)[4])) +
  scale_linetype_manual(values=c(2,1), guide=guide_none()) +
  scale_y_continuous(breaks=c(25,50,75,100))+
  guides(color = guide_legend(title = "", ncol=1),
         shape=guide_legend(title = "", ncol=1)) +
  theme_cowplot() +
  theme(strip.text.y = element_text(angle=0, size=11),
        strip.text.x = element_text(size=11),
        axis.text.y = element_text(size=11), 
        axis.text.x = element_text(size=11, 
                                   angle=45, 
                                   vjust = 0.85, hjust = 0.9),
        axis.title = element_text(size=11),
        legend.position = "bottom",
        legend.title = element_text(size=11),
        legend.text = element_text(size=10)) +
  #lims(y=c(0,101)) +
  labs(x = "\nReturn frequency", y = "Soybean self-sufficiency level (%)\n") +
  facet_wrap(.~pLER_s, nrow=1) 

# Save
ggsave(filename = "E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/Sensitivity_analyses/sensi_2_selfsufficiency_perc.png", 
       width=14, height=6, bg="white", dpi=300)


# Check with sensitivity analysis with intercropped soybean 
# restricted to areas where soybean yield > mean in EU (2.8 t.ha-1)
res2_sensi_p2 %>% 
  # > Keep results for the plot
  filter(target_soybean_lab == "100%",        # soybean target production = 100% consumption
         strategy == "Intercropping", # intercropping
         crop == "Soybean",           # soybean
         pLER_s %in% c(0.3, 0.4, 0.5, 0.56, 0.6, 0.7),
         pLER_m == 0.79               # pLER maize is set as reference value
  ) %>%
  # > Self-sufficiency coverage for each scenario 
  mutate(perc_eu_supply=(production/(target_soybean*10^6))*100,
         perc_eu_supply=if_else(perc_eu_supply>100, 100, perc_eu_supply), 
         surface = surface/10^6,
         production=production/10^6) %>%
  dplyr::select(freq_crop_lab, pLER_s, surface, production, perc_eu_supply) %>%
  gather(key=var, value=val, -freq_crop_lab, -pLER_s) %>%
  mutate(val=round(val,1)) %>%
  spread(key=freq_crop_lab, value=val) %>%
  arrange(var, pLER_s)

#    pLER_s            var 1 year in 4
# 1   0.56 perc_eu_supply        15.0
# 2   0.56     production         5.4
# 3   0.56        surface         3.1

res2_sensi_p2 %>% 
  # > Keep results for the plot
  filter(target_soybean == 36.3,        # soybean target production = 100% consumption
         strategy == "Intercropping", # intercropping
         crop == "Soybean",           # soybean
         pLER_s %in% c(0.3, 0.4, 0.5, 0.56, 0.6, 0.7),
         pLER_m == 0.79               # pLER maize is set as reference value
  ) %>%
  # > Self-sufficiency coverage for each scenario 
  mutate(perc_eu_supply=(production/(36.3*10^6))*100) %>%
  # > Reference values as reference for the dotted lines
  mutate(pLER_lab = case_when(
    pLER_lab == "References values for pLERs"~0, 
    TRUE~1)) %>% 
  # > Plot
  ggplot(., aes(x = freq_crop_lab, 
                y = perc_eu_supply,
                group=pLER_s)) +
  geom_hline(yintercept = c(25, 50, 75, 100), 
             linetype=1, 
             color = "grey90") +
  geom_line(aes(color=as.factor(pLER_s)),
            linewidth=1.2) +
  geom_line(aes(alpha=as.factor(pLER_lab)), 
            linewidth=1, color="white", linetype=3) +
  geom_point(aes(color=as.factor(pLER_s), shape=as.factor(pLER_s)), 
             size=3) +
  # current level of self-sufficiency (cake+grain)
  geom_hline(yintercept = 16, color='black', lty=2) +
  annotate("text", x = 6.75, y = 16, hjust=0, 
           label = "Current\nself-sufficiency\nrate in\nthe EU: 16%", 
           color="black", fontface = 'italic') +
  scale_color_manual(values = c(viridis::viridis(6)[4:6], 
                                viridis::inferno(direction = -1, 6)[2:4])) +
  scale_alpha_manual(values = c(0,1), guide=guide_none()) +
  scale_y_continuous(breaks = c(25,50,75,100)) +
  guides(color = guide_legend(title = "Partial land equivalent ratio for soybean", reverse = T, nrow=2),
         shape = guide_legend(title = "Partial land equivalent ratio for soybean", reverse = T, nrow=2)) +
  theme_cowplot() +
  theme(strip.text.y = element_text(angle=0, size=11),
        strip.text.x = element_text(size=11),
        axis.text = element_text(size=11),
        axis.title = element_text(size=11),
        panel.border = element_rect(color="black"),
        legend.position = "bottom",
        legend.title = element_text(size=11),
        legend.text = element_text(size=10),
        legend.background = element_rect(fill="white"),
        plot.margin = unit(c(1,4,1,1), "cm")
  ) +
  #lims(y=c(0,101)) +
  coord_cartesian(clip = "off", xlim=c(1,6)) +
  labs(x = "\nReturn frequency", y = "Soybean self-sufficiency level (%)\n")

# Save
ggsave(filename = "E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/02_selfsufficiency_perc_v2_mean_yield.png", 
       width=22, height=18, bg="white", units = "cm", dpi=300)



# > Maps of allocation 
sensi_2_data_res %>% 
  # > Keep the data for the figure (reference values)
  filter(pLER_s == 0.56,
         pLER_m == 0.79, 
         freq_crop == 0.25) %>% 
  split(.$target_soybean_lab) %>% 
  map(., ~{
    
    tab_i <- .x %>% 
      # > keep usefull columns 
      dplyr::select(x, y, target_soybean_lab, pixel_intercropping, pixel_solecropping, pixel_intercropping_free, sensi_set) %>% 
      gather(key = pixel_type, value = pixel_color, starts_with("pixel_")) %>%
      mutate(Crop_Design_lab = case_when(
        pixel_type == "pixel_intercropping"      ~ "Intercropping", 
        pixel_type == "pixel_solecropping"       ~ "Soybean as sole crop", 
        pixel_type == "pixel_intercropping_free" ~ "Maize as sole crop")) %>% 
      mutate(Crop_Design_lab = factor(Crop_Design_lab, levels = c("Intercropping", 
                                                                  "Soybean as sole crop", 
                                                                  "Maize as sole crop"))) %>%
      filter(Crop_Design_lab %in% c("Intercropping", 
                                    "Soybean as sole crop", 
                                    "Maize as sole crop")) %>% 
      mutate(Crop_Design_lab2 = if_else(Crop_Design_lab == "Intercropping", "Intercropping", "Sole crops"))
    
    # Maps surface
    plot_i <- ggplot() +
      geom_sf(data=eu27, fill="grey94", color="transparent") +
      geom_tile(data = tab_i, 
                aes(x=x, 
                    y=y, 
                    fill=Crop_Design_lab, 
                    alpha=as.factor(pixel_color))) +
      scale_fill_manual(values = c("#3CBC75FF", "#2D718EFF", "#FDE725FF"), 
                        guide=guide_legend("", nrow = 1)) +
      scale_alpha_manual(values = c(0,1), guide=guide_none()) +
      geom_sf(data = ocean, color = NA, fill = "white",size = 0.2) +
      geom_sf(data = not_eu27, fill = "transparent",size = 0.2) +
      geom_sf(data = europe, fill="transparent") +
      theme_map() + 
      theme(strip.text.y.left = element_text(angle=0),
            strip.background = element_rect(fill="lightgrey", color="transparent"),
            legend.position = "bottom"
      ) +
      facet_grid(sensi_set ~ Crop_Design_lab2, switch = "y") +
      lims(x = c(-11,51), y=c(33,71)) +
      ggtitle(paste0("EU's soybean self-sufficiency level targeted: ", unique(.x$target_soybean_lab), " (", unique(.x$target_soybean), " Mt)"),
              "a. Crop allocation")
    
    # > Save
    #ggsave(plot_i,
    #       filename = paste0("E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/Sensitivity_analyses/sensi_2_map", gsub("%", "", unique(.x$target_soybean_lab)), ".png"), 
    #       width=9, height=8, bg="white", dpi=300)
    
    plot_i
    
  })

# > Barplots - area requirements and co-production of maize and soybean in each strategy
p_sensi_2_1 <- sensi_2_res1 %>%
  filter(pLER_s == 0.56,
         pLER_m == 0.79, 
         freq_crop == 0.25)  %>% 
  split(.$target_soybean_lab) %>%
  map(., ~{
    
    # > Area requirements
    surf_barplots_i <- .x %>% 
      # > keep only 1 line for intercropping (same surface occupied by both crop)
      mutate(to_keep = if_else(strategy == "Sole crop" | strategy == "Intercropping" & crop == "Soybean", 1, 0)) %>% 
      filter(to_keep == 1) %>% 
      # > re-label crop and strategy
      mutate(crop = if_else(strategy == "Intercropping", "Intercropping", paste0(crop, " as sole crop"))) %>% 
      mutate(crop=factor(crop, levels =c("Intercropping", "Maize as sole crop", "Soybean as sole crop"))) %>%
      #mutate(strategy = recode(strategy, "Sole crop"="SC", "Intercropping"="Int")) %>%
      ggplot(data=.) +
      # > total surface intercropping and by crop in sole crops
      geom_col(aes(x=sensi_set , y = surface/10^6, fill=crop),
               col="transparent", width=0.75) +
      # > total surface available (when considering crop frequency = 0.25)
      geom_hline(aes(yintercept=max_surf), 
                 color="darkorange", linetype = 2) + 
      theme_cowplot() + 
      theme(legend.position = "bottom",
            panel.border = element_rect(color="black"),
            axis.title = element_text(size=11),
            axis.text = element_text(size=8.5)
      ) +
      scale_fill_manual(values = c("#3CBC75FF", "#FDE725FF", "#2D718EFF"), name="") +
      facet_wrap(strategy~., scales = "free", ncol=2) +
      ggtitle("") +
      labs(x = "", y = "Total area (Mha)") +
      lims(y=c(0,50))
    
    # > Associated co-productions
    prod_barplot_i <- .x %>%
      # > re-label crop and strategy 
      mutate(crop2 = crop,
             crop2 = factor(crop2, levels = c("Maize", "Soybean"))) %>%
      mutate(crop = if_else(strategy=="Intercropping", "Int.", crop),
             crop = factor(crop, levels = c("Int.", "Maize", "Soybean"))) %>% 
      #mutate(strategy = recode(strategy, "Intercropping"="Int", "Sole crop"="SC")) %>%  
      # > compute sum production especially for intercropping
      group_by(sensi_set, scenario, target_soybean_lab, target_soybean, crop, crop2, strategy) %>% 
      summarize(sum=sum(production)) %>%
      ggplot(data=.) +
      # > total production for each strategy
      geom_col(aes(x=sensi_set , y = sum/10^6, fill=crop2, col=crop), width=0.75, linewidth=1) +
      # > soybean production target
      geom_hline(aes(yintercept=target_soybean), color="red", linetype = 2) + 
      theme_cowplot() + 
      theme(legend.position = "none",
            panel.border = element_rect(color="black"),
            axis.title = element_text(size=11),
            axis.text = element_text(size=8.5)
      ) +
      facet_wrap(.~strategy, scales = "free", ncol=2) +
      #scale_fill_manual(values = c("#3CBC75FF","gray","gray20")) +
      scale_fill_manual(values = c("#FDE725FF","#2D718EFF")) +
      scale_color_manual(values = c("#3CBC75FF","#FDE725FF","#2D718EFF")) +
      ggtitle("") +
      labs(x = "", y = "Total co-production (Mt)") +
      lims(y=c(0,300))
    
    # > Save
    #p <- plot_grid(
    #  plot_grid(surf_barplots_i + 
    #              theme(legend.position = "none") +
    #              ggtitle(paste0("EU's soybean self-sufficiency level targeted: ", unique(.x$target_soybean_lab)," (", unique(.x$target_soybean), " Mt)"),
    #                      "a. Area"),  
    #            prod_barplot_i + 
    #              ggtitle(label = "", subtitle = "b. Production"),  nrow=2, align="hv", axis = "btrl"),
    #  cowplot::ggdraw(cowplot::get_plot_component(surf_barplots_i, 'guide-box-bottom', return_all = TRUE)), ncol=1, align="hv", axis = "btrl", rel_heights = c(0.9, 0.1)
    #)
    
    p <- plot_grid(surf_barplots_i + 
                     theme(legend.position = "none") +
                     ggtitle(label=paste0(unique(.x$target_soybean_lab)), subtitle = "b. Area"),  
                   prod_barplot_i + 
                     ggtitle(label = "", subtitle = "c. Production"),  
                   nrow=2, align="hv", axis = "btrl")
    
    #ggsave(p, 
    #       filename = paste0("E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/Sensitivity_analyses/sensi_2_barplots_", gsub("%", "", unique(.x$target_soybean_lab)), ".png"), 
    #       width=6, height=8, bg="white", dpi=300)
    p
  })

p_sensi_2_1[[3]]

# -------------------------
# > sensitivity analysis including countries outside of EU

# Sets of pixels in the EU members states or outside of the EU
sensi_3_pixels <- Ya_pred_eu %>% 
  filter(model=="pca.m.2", crop=="soybean") %>% 
  group_by(x, y, id_eu27) %>% 
  summarise(mean_Ya_pred = mean(Ya_pred)) %>%
  mutate(cat_pixel_main     = if_else(mean_Ya_pred >= 1 & id_eu27 == 1, 1, 0),
         cat_pixel_sensi = if_else(mean_Ya_pred >= 1, 1, 0))

plot_grid(
  ggplot() +
    geom_sf(data = eu27, fill="grey94") +
    geom_tile(data=sensi_3_pixels, 
              aes(x=x,y=y, fill=as.factor(cat_pixel_main))) +
    geom_sf(data = ocean, color = NA, fill = "white",size = 0.2) +
    geom_sf(data = russia, fill="white", color="white") +
    geom_sf(data = not_eu27, fill="white") +
    geom_sf(data = eu27, fill="transparent") +
    theme_map() + #theme(plot.title = element_text(size=10)) +
    lims(x = c(-11,51), y=c(33,71)) +
    scale_fill_manual(values = c("transparent", viridis::viridis(4)[2]), guide=guide_none()) + 
    ggtitle("a. Main analysis:\nEU member states (N=2047)"),
  ggplot() +
    geom_sf(data=eu27_ext, fill="grey94") +
    geom_tile(data=sensi_3_pixels, 
              aes(x=x,y=y, fill=as.factor(cat_pixel_sensi))) +
    geom_sf(data = ocean, color = NA, fill = "white",size = 0.2) +
    geom_sf(data = russia, fill="white", color="white") +
    geom_sf(data = eu27_ext, fill="transparent") +
    theme_map() + #theme(plot.title = element_text(size=10)) +
    lims(x = c(-11,51), y=c(33,71)) +
    scale_fill_manual(values = c("transparent", viridis::inferno(5)[4]), guide=guide_none()) + 
    ggtitle("b. Sensitivity analysis\nContinental Europe (N=3240)"),
  nrow=1
  
)
    
ggsave(filename = "E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/Sensitivity_analyses/sensi_3_pixels.png",
       width = 12, height = 7, dpi=300, bg="white")

# > Projected yields 
# Breaks for geom_contour_fill, geom_contour, geom_text_contour 
breaks_plot <- c(0, seq(0.5, 3.5, by=0.5), seq(4, 10, by=1))
breaks_labels <- seq(0,10, by=1)

# Map of average yield predictions in continental Europe
Ya_pred_eu %>% 
  filter(model=="pca.m.2") %>%
  # > label according to crops
  mutate(crop = if_else(crop=="maize", "b. Maize", "a. Soybean")) %>% 
  # > long format 
  gather(key=year, value=Ya_pred, starts_with("X2")) %>% 
  # > for each pixel, recompute yields (in t/ha)
  group_by(crop, x, y) %>% 
  mutate(Ya_pred_t_ha = Ya_pred/cropland_area_ha) %>%
  summarise(mean_Ya_pred = mean(Ya_pred, na.rm=T)) %>%
  mutate(mean_Ya_pred = ifelse(is.na(mean_Ya_pred)==T, 0, mean_Ya_pred)) %>% 
  # > plot
  ggplot(.) + 
  geom_sf(data=eu27_ext, fill="grey94") +
  geom_contour_fill(aes(x=x, y=y, z=mean_Ya_pred), 
                    breaks = breaks_plot,
                    na.fill = T, 
                    global.breaks = F,
                    clip = eu27_ext) +
  geom_sf(data=eu27_ext, fill="transparent") +
  geom_contour(aes(x=x, y=y, z=mean_Ya_pred), 
               color = "white", 
               linewidth = 0.01,
               breaks = breaks_plot) +
  geom_text_contour(aes(x=x, y=y, z=mean_Ya_pred), 
                    color="black",
                    size=3, 
                    breaks = breaks_labels) +
  facet_grid(.~crop) + 
  theme_map() + 
  lims(x = c(-11,51), y=c(33,71)) +
  theme(legend.position = "bottom",, 
        legend.title = element_text(size=12), 
        legend.text = element_text(size=10),
        strip.text = element_text(face = "bold", hjust = 0.1, size=15)) +
  scale_fill_gradientn(colours = c("transparent", viridis::plasma(n=100)[50:100], viridis::viridis(direction=-1, n=100)[1:75]),
                       breaks = breaks_labels, 
                       labels = breaks_labels,
                       guide = guide_colorbar(barwidth = 20, barheight = 0.5, title.position = "top", title = expression("Mean yield t " ~ ha^-1)))

ggsave(filename = "E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/Sensitivity_analyses/sensi_3_pred_maps_eu_ext.png", 
       height=7, width=11, bg="white", dpi = 300)



# > Load allocations
sensi_3 <- loadRDa(paste0(path_alloc, "sensi/sensi_3_eu42.rda"))

# > Allocations and area requierements
alloc_sensi_3 <- list("Main analysis: EU member states (N=2047)"          = allocations,
                      "Sensitivity analysis: Continental Europe (N=3240)" = sensi_3)


# > Shapping allocations results
# Allocations
sensi_3_data_res <- format_alloc(results_allocations =  alloc_sensi_3 %>% 
                                   map_dfr(., ~{ 
                                     
                                     map_dfr(.x, ~{ .x$data_res }, .id="scenario")
                                     
                                   }, .id="sensi_set")) %>% 
  # > Set of sensitivity analyses
  mutate(sensi_set = recode(sensi_set, 
                            "Main analysis: EU member states (N=2047)"         ="EU member states\n(main analysis)",
                            "Sensitivity analysis: Continental Europe (N=3240)"="Continental Europe")) %>%
  mutate(sensi_set = factor(sensi_set, levels=c("EU member states\n(main analysis)",
                                                "Continental Europe")))

# Surface required to cover different shares of 
#    EU's consumption of soybean 
sensi_3_res1 <- format_alloc(results_allocations =  alloc_sensi_3 %>% 
                               map_dfr(., ~{ 
                                 
                                 map_dfr(.x, ~{ .x$res1 }, .id="scenario")
                                 
                               }, .id="sensi_set")) %>%
  # > Set of sensitivity analyses
  mutate(sensi_set = recode(sensi_set, 
                            "Main analysis: EU member states (N=2047)"         ="Main analysis:\nEU member states"         ,
                            "Sensitivity analysis: Continental Europe (N=3240)"="Sensitivity analysis:\nContinental Europe")) %>%
  mutate(sensi_set = factor(sensi_set, levels=c("Main analysis:\nEU member states",
                                                "Sensitivity analysis:\nContinental Europe")))

# Surface required to produce the same production (maize+soybean)
#    as intercropping 
sensi_3_res3 <- format_alloc(results_allocations =  alloc_sensi_3 %>% 
                               map_dfr(., ~{ 
                                 
                                 map_dfr(.x, ~{ .x$res3 }, .id="scenario")
                                 
                               }, .id="sensi_set")) %>%
  # > Set of sensitivity analyses
  mutate(sensi_set = recode(sensi_set, 
                            "Main analysis: EU member states (N=2047)"         ="Main analysis:\nEU member states"         ,
                            "Sensitivity analysis: Continental Europe (N=3240)"="Sensitivity analysis:\nContinental Europe")) %>%
  mutate(sensi_set = factor(sensi_set, levels=c("Main analysis:\nEU member states",
                                                "Sensitivity analysis:\nContinental Europe")))

# > Compare with main analysis
sensi_3_res1 %>% 
  filter(pLER_s==0.56, pLER_m==0.79, freq_crop == 0.25, target_soybean_lab == "75%") %>% 
  dplyr::select(sensi_set, crop, strategy, surface) %>% 
  spread(key = strategy, value = surface)

sensi_3_res1 %>% 
  filter(pLER_s==0.56, pLER_m==0.79, freq_crop == 0.25, sensi_set=="Sensitivity analysis: Continental Europe (N=3240)") %>% 
  dplyr::select(sensi_set, crop, strategy, surface) %>% 
  spread(key = strategy, value = surface)

sensi_3_res3 %>% 
  filter(pLER_s==0.56, pLER_m==0.79, freq_crop == 0.25, target_soybean == 33.75) %>% 
  dplyr::select(sensi_set, crop, strategy, total_production) %>% 
  spread(key = strategy, value = total_production)

sensi_3_res3 %>% 
  filter(pLER_s==0.56, pLER_m==0.79, freq_crop == 0.25, target_soybean == 33.75) %>% 
  dplyr::select(sensi_set, crop, strategy, total_surface) %>% 
  spread(key = strategy, value = total_surface)

# > % of self-sufficiency reached 
sensi_3_res1 %>% 
  # > Keep results for the plot
  filter(target_soybean_lab == "100%",        # soybean target production = 100% consumption
         strategy == "Intercropping", # intercropping
         crop == "Soybean",           # soybean
         pLER_s %in% c(0.3, 0.4, 0.5, 0.56, 0.6, 0.7),
         pLER_m == 0.79               # pLER maize is set as reference value
  ) %>%
  # > Self-sufficiency coverage for each scenario 
  mutate(perc_eu_supply=(production/(58.4*10^6))*100) %>%
  # > Reference values as reference for the dotted lines
  mutate(pLER_lab = case_when(
    pLER_lab == "References values for pLERs"~0, 
    TRUE~1)) %>% 
  mutate(pLER_s=paste0("Soybean pLER = ", pLER_s)) %>%
  # > Labels
  # > Set of sensitivity analyses
  mutate(sensi_set = recode(sensi_set, 
                            "Main analysis:\nEU member states"="Main analysis: EU member states (N=2047)",
                            "Sensitivity analysis:\nContinental Europe"="Sensitivity analysis: Continental Europe (N=3240)")) %>%
  # > Plot
  ggplot(., aes(x = freq_crop_lab, 
                y = perc_eu_supply,
                group=sensi_set)) +
  # current level of self-sufficiency (cake+grain)
  geom_hline(yintercept = 42.3, color='red', lty=2) +
  geom_hline(yintercept = c(25, 50, 75, 100), 
             linetype=2, 
             color = "grey90") +
  geom_line(aes(color=as.factor(sensi_set)), #linetype = 2,
            linewidth=0.75) +
  geom_point(aes(color=as.factor(sensi_set), shape=as.factor(sensi_set)), 
             size=2) +
  scale_color_manual(values = c(viridis::viridis(4)[2], viridis::inferno(5)[4])) +
  scale_linetype_manual(values=c(2,1), guide=guide_none()) +
  scale_y_continuous(breaks=c(25,50,75,100))+
  guides(color = guide_legend(title = "", ncol=1),
         shape=guide_legend(title = "", ncol=1)) +
  theme_cowplot() +
  theme(strip.text.y = element_text(angle=0, size=11),
        strip.text.x = element_text(size=11),
        axis.text.y = element_text(size=11), 
        axis.text.x = element_text(size=11, 
                                   angle=45, 
                                   vjust = 0.85, hjust = 0.9),
        axis.title = element_text(size=11),
        legend.position = "bottom",
        legend.title = element_text(size=11),
        legend.text = element_text(size=10)) +
  #lims(y=c(0,101)) +
  labs(x = "\nReturn frequency", y = "Soybean self-sufficiency level (%)\n") +
  facet_wrap(.~pLER_s, nrow=1) 

# Save
ggsave(filename = "E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/Sensitivity_analyses/sensi_3_selfsufficiency_perc.png", 
       width=14, height=6, bg="white", dpi=300)

# > Maps of allocation 
sensi_3_data_res %>% 
  # > Keep the data for the figure (reference values)
  filter(pLER_s == 0.56,
         pLER_m == 0.79, 
         freq_crop == 0.25) %>% 
  split(.$target_soybean_lab) %>% 
  map(., ~{
    
    tab_i <- .x %>% 
      # > keep usefull columns 
      dplyr::select(x, y, target_soybean_lab, pixel_intercropping, pixel_solecropping, pixel_intercropping_free, sensi_set) %>% 
      gather(key = pixel_type, value = pixel_color, starts_with("pixel_")) %>%
      mutate(Crop_Design_lab = case_when(
        pixel_type == "pixel_intercropping"      ~ "Intercropping", 
        pixel_type == "pixel_solecropping"       ~ "Soybean as sole crop", 
        pixel_type == "pixel_intercropping_free" ~ "Maize as sole crop")) %>% 
      mutate(Crop_Design_lab = factor(Crop_Design_lab, levels = c("Intercropping", 
                                                                  "Soybean as sole crop", 
                                                                  "Maize as sole crop"))) %>%
      filter(Crop_Design_lab %in% c("Intercropping", 
                                    "Soybean as sole crop", 
                                    "Maize as sole crop")) %>% 
      mutate(Crop_Design_lab2 = if_else(Crop_Design_lab == "Intercropping", "Intercropping", "Sole crops"))
    
    # Maps surface
    plot_i <- ggplot() +
      geom_sf(data=eu27_ext, fill="grey94", color="transparent") +
      geom_tile(data = tab_i, 
                aes(x=x, 
                    y=y, 
                    fill=Crop_Design_lab, 
                    alpha=as.factor(pixel_color))) +
      scale_fill_manual(values = c("#3CBC75FF", "#2D718EFF", "#FDE725FF"), 
                        guide=guide_legend("", nrow = 1)) +
      scale_alpha_manual(values = c(0,1), guide=guide_none()) +
      geom_sf(data = ocean, color = NA, fill = "white",size = 0.2) +
      geom_sf(data = not_eu27, fill = "transparent",size = 0.2) +
      geom_sf(data = europe, fill="transparent") +
      theme_map() + 
      theme(strip.text.y.left = element_text(angle=0),
            strip.background = element_rect(fill="lightgrey", color="transparent"),
            legend.position = "bottom"
      ) +
      facet_grid(sensi_set ~ Crop_Design_lab2, switch = "y") +
      lims(x = c(-11,51), y=c(33,71)) +
      ggtitle(paste0("EU's soybean self-sufficiency level targeted: ", unique(.x$target_soybean_lab), " (", unique(.x$target_soybean), " Mt)"),
              "a. Crop allocation")
    
    # > Save
    ggsave(plot_i,
           filename = paste0("E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/Sensitivity_analyses/sensi_3_map", gsub("%", "", unique(.x$target_soybean_lab)), ".png"), 
           width=9, height=8, bg="white", dpi=300)
    
  })

# > Barplots - area requirements and co-production of maize and soybean in each strategy
sensi_3_res1 %>%
  filter(pLER_s == 0.56,
         pLER_m == 0.79, 
         freq_crop == 0.25)  %>% 
  split(.$target_soybean_lab) %>%
  map(., ~{
    
    # > Area requirements
    surf_barplots_i <- .x %>% 
      # > keep only 1 line for intercropping (same surface occupied by both crop)
      mutate(to_keep = if_else(strategy == "Sole crop" | strategy == "Intercropping" & crop == "Soybean", 1, 0)) %>% 
      filter(to_keep == 1) %>% 
      # > re-label crop and strategy
      mutate(crop = if_else(strategy == "Intercropping", "Intercropping", paste0(crop, " as sole crop"))) %>% 
      mutate(crop=factor(crop, levels =c("Intercropping", "Maize as sole crop", "Soybean as sole crop"))) %>%
      #mutate(strategy = recode(strategy, "Sole crop"="SC", "Intercropping"="Int")) %>%
      ggplot(data=.) +
      # > total surface intercropping and by crop in sole crops
      geom_col(aes(x=sensi_set , y = surface/10^6, fill=crop),
               col="transparent", width=0.75) +
      # > total surface available (when considering crop frequency = 0.25)
      geom_hline(aes(yintercept=max_surf), 
                 color="darkorange", linetype = 2) + 
      theme_cowplot() + 
      theme(legend.position = "bottom",
            panel.border = element_rect(color="black"),
            axis.title = element_text(size=11),
            axis.text = element_text(size=8.5)
      ) +
      scale_fill_manual(values = c("#3CBC75FF", "#FDE725FF", "#2D718EFF"), name="") +
      facet_wrap(strategy~., scales = "free", ncol=2) +
      ggtitle("") +
      labs(x = "", y = "Total area (Mha)") +
      lims(y=c(0,50))
    
    # > Associated co-productions
    prod_barplot_i <- .x %>%
      # > re-label crop and strategy 
      mutate(crop2 = crop,
             crop2 = factor(crop2, levels = c("Maize", "Soybean"))) %>%
      mutate(crop = if_else(strategy=="Intercropping", "Int.", crop),
             crop = factor(crop, levels = c("Int.", "Maize", "Soybean"))) %>% 
      #mutate(strategy = recode(strategy, "Intercropping"="Int", "Sole crop"="SC")) %>%  
      # > compute sum production especially for intercropping
      group_by(sensi_set, scenario, target_soybean_lab, target_soybean, crop, crop2, strategy) %>% 
      summarize(sum=sum(production)) %>%
      ggplot(data=.) +
      # > total production for each strategy
      geom_col(aes(x=sensi_set , y = sum/10^6, fill=crop2, col=crop), width=0.75, linewidth=1) +
      # > soybean production target
      geom_hline(aes(yintercept=target_soybean), color="red", linetype = 2) + 
      theme_cowplot() + 
      theme(legend.position = "none",
            panel.border = element_rect(color="black"),
            axis.title = element_text(size=11),
            axis.text = element_text(size=8.5)
      ) +
      facet_wrap(.~strategy, scales = "free", ncol=2) +
      #scale_fill_manual(values = c("#3CBC75FF","gray","gray20")) +
      scale_fill_manual(values = c("#FDE725FF","#2D718EFF")) +
      scale_color_manual(values = c("#3CBC75FF","#FDE725FF","#2D718EFF")) +
      ggtitle("") +
      labs(x = "", y = "Total co-production (Mt)") +
      lims(y=c(0,300))
    
    # > Save
    #p <- plot_grid(
    #  plot_grid(surf_barplots_i + 
    #              theme(legend.position = "none") +
    #              ggtitle(paste0("EU's soybean self-sufficiency level targeted: ", unique(.x$target_soybean_lab)," (", unique(.x$target_soybean), " Mt)"),
    #                      "a. Area"),  
    #            prod_barplot_i + 
    #              ggtitle(label = "", subtitle = "b. Production"),  nrow=2, align="hv", axis = "btrl"),
    #  cowplot::ggdraw(cowplot::get_plot_component(surf_barplots_i, 'guide-box-bottom', return_all = TRUE)), ncol=1, align="hv", axis = "btrl", rel_heights = c(0.9, 0.1)
    #)
    
    p <- plot_grid(surf_barplots_i + 
                     theme(legend.position = "none") +
                     ggtitle(label="", subtitle = "b. Area"),  
                   prod_barplot_i + 
                     ggtitle(label = "", subtitle = "c. Production"),  
                   nrow=2, align="hv", axis = "btrl")
    
    ggsave(p, 
           filename = paste0("E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/Sensitivity_analyses/sensi_3_barplots_", gsub("%", "", unique(.x$target_soybean_lab)), ".png"), 
           width=6, height=8, bg="white", dpi=300)
    
  })



# -------------------------
# > comparison of intercropping and pure stands stategies 
# according pLER values and crop return frequencies

# > data loading or extracting+load
#load_data <- T  # to extract and shape data from raw allocations
load_data <- F # to only load the data 

# > extract and shape data from raw allocations
if(load_data==T)
{
  
  # Load allocation in all scenarios (N=40344)
  sensi_4_list <- list()
  
  for(file_i in list.files("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/08_allocations/allocations_with_max_surf/sensi_pLERs", full.names = T))
  {
    
    # load rda
    sensi_i <- loadRDa(paste0(file_i))
    
    # store data in a list
    sensi_4_list[[paste0(file_i)]] <- sensi_i %>% 
      map_dfr(., ~{
        
        .x$res2
        
      }, .id="scenario") 
    
  }
  
  # put all data in the same table
  sensi_4_res2 <- format_alloc(sensi_4_list %>% map_dfr(., ~{ .x }, .id="file_i")) 
  
  # save
  save(sensi_4_res2, file = "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/08_allocations/allocations_with_max_surf/sensi/sensi_4.rda")
  
}

# > just load
if(load_data==F)
{
  
  sensi_4_res2 <- loadRDa("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/08_allocations/allocations_with_max_surf/sensi/sensi_4.rda")
  
}

# Teta values 
summary(sensi_4_res2$Teta_1) 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  1.000   1.021   1.032   1.036   1.053   1.113 

summary(sensi_4_res2$Teta_2)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's  
#  0.960   1.009   1.025   1.029   1.052   1.084    5043

# > 

plot_grid(
  sensi_4_res2 %>% 
    filter(pLER_m == 0.79) %>% 
    ggplot(., aes(x=freq_crop_lab, y=Teta_1)) +
    # > Teta variability by crop frequency and soybean target
    geom_boxplot(fill="grey84", outlier.shape = NA) +
    geom_jitter(aes(color=pLER_s), 
                width=0.25, size=0.5) +
    # > Teta=1 means that yields are homogenous
    geom_hline(yintercept = 1, color="black", lty=2)+
    # > Teta of the reference pLERs for soybean 
    geom_point(data=sensi_4_res2 %>% 
                 filter(pLER_m == 0.79, pLER_s==0.56),
               aes(x=freq_crop_lab, y=Teta_1),
               color="black", fill='yellow', size=2, shape=21) +
    scale_color_gradientn(
      colours = RColorBrewer::brewer.pal(9,"BuPu"),
      name = "Partial land equivalent ratio for soybean           ",
      guide = guide_colorbar(barwidth = 10, barheight = 0.5, title.position = "left")
    ) +
    facet_grid(.~target_soybean_lab) +
    theme_cowplot() +
    theme(legend.position = "bottom",
          legend.text = element_text(size=9),
          legend.title = element_text(size=9),
          axis.text = element_text(size=9),
          axis.title = element_text(size=10),
          strip.text = element_text(size=10),
          plot.title = element_text(size=11),
          plot.subtitle = element_text(size=11)) +
    lims(y=c(0.94,1.12))+
    labs(x="",y="Coefficient of homogeneity of yields in intercropping vs sole crops areas") +
    ggtitle("a. Soybean", "Variability according to partial land equivalent ratio for soybean\n(partial land equivalent ratio for maize fixed at 0.79)") +
    coord_flip(),
  
  sensi_4_res2 %>% 
      filter(pLER_s == 0.56) %>% 
      ggplot(., aes(x=freq_crop_lab, y=Teta_1)) +
      # > Teta variability by crop frequency and soybean target
      geom_boxplot(fill="grey84",outlier.shape = NA) +
      geom_jitter(aes(color=pLER_m), 
                  width=0.25, size=0.5) +
      # > Teta=1 means that yields are homogenous
      geom_hline(yintercept = 1, color="black", lty=2)+
      # > Teta of the reference pLERs for soybean 
      geom_point(data=sensi_4_res2 %>% 
                   filter(pLER_m == 0.79, pLER_s==0.56),
                 aes(x=freq_crop_lab, y=Teta_1),
                 color="black", fill='yellow', size=2, shape=21) +
      scale_color_gradientn(
        colours = RColorBrewer::brewer.pal(9,"BuGn"),
        name = "Partial land equivalent ratio for maize           \n",
        guide = guide_colorbar(barwidth = 10, barheight = 0.5, title.position = "left")
      ) +
      facet_grid(.~target_soybean_lab) +
      theme_cowplot() +
      theme(legend.position = "bottom",
            legend.text = element_text(size=9),
            legend.title = element_text(size=9),
            axis.text = element_text(size=9),
            axis.title = element_text(size=10),
            strip.text = element_text(size=10),
            plot.title = element_text(size=11),
            plot.subtitle = element_text(size=11)) +
      lims(y=c(0.94,1.12))+
      labs(x="",y="Coefficient of homogeneity of yields in intercropping vs sole crops areas") +
      ggtitle("", "Variability according partial land equivalent ratio for maize\n(partial land equivalent ratio for soybean fixed at 0.56)") +
      coord_flip(),
  
  sensi_4_res2 %>% 
    filter(pLER_m == 0.79) %>% 
    mutate(Teta_2 = if_else(is.na(Teta_2), min(sensi_4_res2$Teta_2, na.rm = T)-0.01, Teta_2)) %>% 
    ggplot(., aes(x=freq_crop_lab, y=Teta_2)) +
    # > Teta variability by crop frequency and soybean target
    geom_boxplot(fill="grey84",outlier.shape = NA) +
    geom_jitter(aes(color=pLER_s), 
                width=0.25, size=0.5) +
    # > Teta=1 means that yields are homogenous
    geom_hline(yintercept = 1, color="black", lty=2)+
    # > Teta of the reference pLERs for soybean 
    geom_point(data=sensi_4_res2 %>% 
                 filter(pLER_m == 0.79, pLER_s==0.56) %>%
                 mutate(Teta_2 = if_else(is.na(Teta_2), min(sensi_4_res2$Teta_2, na.rm = T)-0.01, Teta_2)),
               aes(x=freq_crop_lab, y=Teta_2),
               color="black", fill='yellow', size=2, shape=21) +
    scale_color_gradientn(
      colours = RColorBrewer::brewer.pal(9,"BuPu"),
      name = "Partial land equivalent ratio for soybean           ",
      guide = guide_colorbar(barwidth = 10, barheight = 0.5, title.position = "left")
    ) +
    facet_grid(.~target_soybean_lab) +
    theme_cowplot() +
    theme(legend.position = "bottom",
          legend.text = element_text(size=9),
          legend.title = element_text(size=9),
          axis.text = element_text(size=9),
          axis.title = element_text(size=10),
          strip.text = element_text(size=10),
          plot.title = element_text(size=11),
          plot.subtitle = element_text(size=11)) +
    lims(y=c(0.94,1.12))+
    labs(x="",y="Coefficient of homogeneity of yields in intercropping vs sole crops areas") +
    ggtitle("b. Maize", "Variability according partial land equivalent ratio for soybean\n(partial land equivalent ratio for maize fixed at 0.79)") +
    coord_flip(),
  
  sensi_4_res2 %>% 
    filter(pLER_s == 0.56) %>% 
    mutate(Teta_2 = if_else(is.na(Teta_2), min(sensi_4_res2$Teta_2, na.rm = T)-0.01, Teta_2)) %>% 
    ggplot(., aes(x=freq_crop_lab, y=Teta_2)) +
    # > Teta variability by crop frequency and soybean target
    geom_boxplot(fill="grey84",outlier.shape = NA) +
    geom_jitter(aes(color=pLER_m), 
                width=0.25, size=0.5) +
    # > Teta=1 means that yields are homogenous
    geom_hline(yintercept = 1, color="black", lty=2)+
    # > Teta of the reference pLERs for soybean 
    geom_point(data=sensi_4_res2 %>% 
                 filter(pLER_m == 0.79, pLER_s==0.56) %>%
                 mutate(Teta_2 = if_else(is.na(Teta_2), min(sensi_4_res2$Teta_2, na.rm = T)-0.01, Teta_2)),
               aes(x=freq_crop_lab, y=Teta_2),
               color="black", fill='yellow', size=2, shape=21) +
    scale_color_gradientn(
      colours = RColorBrewer::brewer.pal(9,"BuGn"),
      name = "Partial land equivalent ratio for maize           ",
      guide = guide_colorbar(barwidth = 10, barheight = 0.5, title.position = "left")
    ) +
    facet_grid(.~target_soybean_lab) +
    theme_cowplot() +
    theme(legend.position = "bottom",
          legend.text = element_text(size=9),
          legend.title = element_text(size=9),
          axis.text = element_text(size=9),
          axis.title = element_text(size=10),
          strip.text = element_text(size=10),
          plot.title = element_text(size=11),
          plot.subtitle = element_text(size=11)) +
    lims(y=c(0.94,1.12))+
    labs(x="",y="Coefficient of homogeneity of yields in intercropping vs sole crops areas") +
    ggtitle("", "Variability according to partial land equivalent ratio for maize\n(partial land equivalent ratio for soybean fixed at 0.56)") +
    coord_flip(),
  
  ncol=2
  
)

ggsave(filename = "E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/Sensitivity_analyses/sensi_4_teta_boxplots.png", 
       dpi=400, width=16, height=8.5, bg="white")


plot_grid(
  sensi_4_res2 %>%
    filter(#freq_crop == '1 year in 4', 
      pLER_m == 0.79) %>% 
    mutate(id_shape = if_else(pLER_s==0.56, 1, 0)) %>% 
    ggplot(., aes(x=pLER_s, y = Teta_1, 
                  color=target_soybean_lab)) + 
    # Teta=1 means that yields are homogenous
    geom_hline(yintercept=1, color="black", lty=2) +
    geom_point(aes(alpha=id_shape)) +
    geom_path() +
    theme_cowplot() + 
    facet_grid(.~freq_crop) +
    theme(strip.text.y = element_text(angle=0, size=11),
          strip.text.x = element_text(size=11),
          axis.text = element_text(size=11),
          axis.title = element_text(size=11),
          legend.position = "none") +
    scale_alpha_continuous(range = c(0,1)) +
    scale_color_manual(values = viridis::inferno(5)[1:4]) +
    labs(y = "Ratio between mean yields\nin sole crop and intercrop areas", 
         x = "Partial land equivalent ratio for soybean", 
         color = "EU's soybean self-sufficency target") +
    lims(y = c(0.94,1.15)) +
    ggtitle("a. Soybean"),
  sensi_4_res2 %>%
    filter(pLER_m == 0.79) %>% 
    # some NAs values 
    mutate(id_shape1 = if_else(is.na(Teta_2), 1, 0),
           id_shape2 = if_else(pLER_s==0.56, 1, 0), 
           Teta_2 = if_else(is.na(Teta_2), 0.95, Teta_2)) %>% 
    ggplot(., aes(x=pLER_s, y = Teta_2, 
                  color=target_soybean_lab)) + 
    # Teta=1 means that yields are homogenous
    geom_hline(yintercept=1, color="black", lty=2) +
    geom_point(aes(shape=as.factor(id_shape1),alpha=id_shape2)) +
    geom_path() +
    theme_cowplot() + 
    theme(strip.text.y = element_text(angle=0, size=11),
          strip.text.x = element_text(size=11),
          axis.text = element_text(size=11),
          axis.title = element_text(size=11),
          legend.position = "bottom",
          legend.title = element_text(size=11),
          legend.text = element_text(size=10)) +
    labs(y = "Ratio between mean yields\nin sole crop and intercrop areas", 
         x = "Partial land equivalent ratio for soybean", 
         color = "EU's soybean self-sufficency levels targeted") +
    facet_grid(.~freq_crop) +
    scale_shape_manual(values=c(16, 4)) +
    scale_color_manual(values = viridis::inferno(5)[1:4]) +
    scale_alpha_continuous(range = c(0,1)) +
    guides(shape=guide_none(), alpha=guide_none()) +
    lims(y = c(0.94,1.15)) +
    ggtitle("b. Maize"),
  ncol=1, rel_heights = c(0.45,0.55)
)

ggsave(filename = "E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/Sensitivity_analyses/sensi_4_teta_lines.png", 
       dpi=400, width=12, height=7, bg="white")

# > Heatmap
summary(sensi_4_res2$Ratio_2)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.759   1.045   1.163   1.163   1.279   1.633    5043 

limits.min.plots <- 0.7
limits.max.plots <- 1.7
midpoint <- 1

sensi_4_htp <- sensi_4_res2 %>%
  filter(freq_crop_lab == '1 year in 4', target_soybean_lab=="50%") %>% 
  mutate(freq_crop_lab = factor(freq_crop_lab, levels = c('1 year in 2', '1 year in 3', '1 year in 4', '1 year in 5', '1 year in 6', '1 year in 7'))) %>% 
  ggplot(.) + 
  geom_raster(aes(y=pLER_s, x = pLER_m, fill=Ratio_2)) + 
  geom_point(y = 0.56, x = 0.79, pch = 15, size=2, color="white") + 
  geom_linerange(y=0.56, xmin=0, xmax=0.79, linetype=2, color="white") + 
  geom_linerange(ymin=0, ymax=0.56, x=0.79, linetype=2, color="white") + 
  theme_cowplot() + 
  theme(strip.text.y = element_text(angle=0, size=11),
        strip.text.x = element_text(size=11),
        axis.text = element_text(size=11),
        axis.title = element_text(size=11),
        legend.position = "bottom",
        legend.title = element_text(size=11),
        legend.text = element_text(size=10)) +
  #facet_grid(freq_crop_lab~target_soybean_lab, scales = "free") + 
  labs(y = "Partial land equivalent ratio for\nsoybean", x = "Partial land equivalent for\nmaize", 
       fill = "Maize production ratio between intercrop\nand sole crops strategies") + 
  scale_fill_gradientn(
    colours = rev(c(viridis::mako(30), "white", viridis::rocket(30, direction = -1))),
    values = c(0, (midpoint - limits.min.plots)/(limits.max.plots - limits.min.plots), 1), 
    limits = c(limits.min.plots, limits.max.plots),
    breaks = c(seq(limits.min.plots, limits.max.plots, by = 0.1)),
    na.value = "white", 
    guide = guide_colorbar(barwidth = 15, barheight = 0.5, position = "bottom")) ; sensi_4_htp

#ggsave(plot = sensi_4_htp,
#       filename = "E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/Sensitivity_analyses/sensi_4_heat_map_Ratio2.png", 
#       dpi=400, width=9, height=11, bg="white")

ggsave(plot = sensi_4_htp,
       filename = "E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/Sensitivity_analyses/sensi_4_heat_map_Ratio2_small.png", 
       dpi=400, width=8, height=8, bg="white")

# ----------------------------------------
# ----------------------------------------
# Supplementary figures 

load("E:/POSTDOC INRAE/ANALYSES/B_OPTIMISATION/00_Data/00_tab_soybean.rda")
load("E:/POSTDOC INRAE/ANALYSES/B_OPTIMISATION/00_Data/00_tab_maize.rda")

# > Maps 
# Analyses on the "global" datasets 
tab_soybean %>%
  distinct(x,y,country_name) %>% 
  mutate(desert=if_else(country_name=="Desert", 1, 0)) %>%
  group_by(desert) %>%
  count()
#  desert     n
#1      0  2626
#2      1   660

tab_maize %>%
  distinct(x,y,country_name) %>% 
  mutate(desert=if_else(country_name=="Desert", 1, 0)) %>%
  group_by(desert) %>%
  count()
#  desert     n
#1      0  1660
#2      1   404


# > Soybean
g_s <- ggplot() + 
  geom_sf(data = world) +
  geom_point(data = tab_soybean %>% 
               distinct(x, y, country_name) %>% 
               mutate(desert=if_else(country_name=="Desert", 
                                     "Selected sites in geographical areas unsuitable for soybean production (yield = 0 t/ha) (N=660)", 
                                     "Selected sites in major soybean producing areas (N=2626)")) %>% 
               mutate(desert=factor(desert, levels = c("Selected sites in major soybean producing areas (N=2626)", "Selected sites in geographical areas unsuitable for soybean production (yield = 0 t/ha) (N=660)"))), 
             aes(x=x, y=y, color=as.factor(desert)), size=0.1) + 
  theme_map() + 
  theme(legend.position = "bottom", 
        legend.text = element_text(size=9)) +
  scale_color_manual(values = c("darkgreen", "black")) +
  guides(color=guide_legend(title = "", override.aes = list(size=2), nrow = 2)) +  
  ggtitle("a. Soybean")

# > Maize
g_m <- ggplot() + 
  geom_sf(data = world) +
  geom_point(data = tab_maize %>% 
               distinct(x, y, country_name) %>% 
               mutate(desert=if_else(country_name=="Desert", 
                                     "Selected sites in geographical areas unsuitable for maize production (yield = 0 t/ha) (N=404)", 
                                     "Selected sites in maize producing areas in North Hemisphere (N=1660)")) %>% 
               mutate(desert=factor(desert, levels = c("Selected sites in maize producing areas in North Hemisphere (N=1660)",
                                                       "Selected sites in geographical areas unsuitable for maize production (yield = 0 t/ha) (N=404)"))), 
             aes(x=x, y=y, color=as.factor(desert)), size=0.1) + 
  theme_map() + 
  theme(legend.position = "bottom", 
        legend.text = element_text(size=9)) +
  scale_color_manual(values = c("darkorange", "black")) +
  guides(color=guide_legend(title = "", override.aes = list(size=2), nrow = 2)) +  
  ggtitle("b. Maize")

sp1 <- plot_grid(g_s, g_m, nrow=2, axis="tbrl", align="hv")

ggsave(sp1, 
       filename = "E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/Supporting information/01_maps.png",
       width=7, height=9, bg = "white")

# ----------------------------------------
# Distribution of climate predictors 

tab_soybean %>% 
  group_by(x,y) %>% 
  summarise(mean_PC1_month_max_temp = mean(PC1_month_max_temp)) %>% 
  ggplot(.) + 
  geom_sf(data=world) + 
  geom_point(aes(x=x,y=y,color=mean_PC1_month_max_temp), size=0.1)

# ---------------------------------------
# Distribution of observed yields - world

paste0(round(mean(tab_soybean$Ya), 2), " (", round(sd(tab_soybean$Ya), 2), ")\n[", round(min(tab_soybean$Ya), 2), " - ", round(max(tab_soybean$Ya), 2), "]")
# "2.04 (1.38) [0 - 6.7]"

paste0(round(mean(tab_maize$Ya), 2), " (", round(sd(tab_maize$Ya), 2), ")\n[", round(min(tab_maize$Ya), 2), " - ", round(max(tab_maize$Ya), 2), "]")
# "7.58 (5.08) [0 - 19.97]"

# Global and by country
sp2 <- plot_grid(
  plot_grid(
    ggplot() + 
      geom_density(data=tab_soybean %>% filter(country_name != "Desert"), aes(x=Ya)) + 
      geom_point(data=tab_soybean %>% filter(country_name != "Desert") %>% group_by(country_name) %>% summarise(Ya_med = median(Ya)), aes(x=Ya_med, y=0, color=country_name)) + 
      theme_cowplot() + 
      theme(legend.position = "none",
            axis.title.x = element_blank(),
            plot.margin = unit(c(0,0,0,0), "cm")) + 
      ggtitle("a. Soybean") + 
      labs(x = "Yield (tons/ha)", y = "Density"),
    
    ggplot() + 
      geom_boxplot(data=tab_soybean %>% 
                     filter(country_name != "Desert"), 
                   aes(x=reorder(country_name, Ya, median), y=Ya), outlier.shape = NA) + 
      geom_point(data=tab_soybean %>% filter(country_name != "Desert") %>% 
                   group_by(country_name) %>% 
                   summarise(Ya_med = median(Ya)), 
                 aes(x=reorder(country_name, Ya_med, max), y=Ya_med, color=country_name)) + 
      theme_cowplot() + 
      theme(legend.position = "none",
            plot.margin = unit(c(0,0,0,0), "cm")) + 
      labs(y = "Yield (tons/ha)", x="") +
      coord_flip(),
    ggplot() +
      theme_void(),
    ncol=1, align ="hv", axis="tbrl", rel_heights = c(0.33,0.33,0.33) 
  ),
  
  plot_grid(
    ggplot() + 
      geom_density(data=tab_maize %>% filter(country_name != "Desert"), aes(x=Ya)) + 
      geom_point(data=tab_maize %>% filter(country_name != "Desert") %>% group_by(country_name) %>% summarise(Ya_med = median(Ya)), aes(x=Ya_med, y=0, color=country_name)) + 
      theme_cowplot() + 
      theme(legend.position = "none",
            axis.title.x = element_blank(),
            plot.margin = unit(c(0,0,0,0), "cm")) + 
      ggtitle("b. Maize") + 
      labs(x = "Yield (tons/ha)", y = "Density"),
    
    ggplot() + 
      geom_boxplot(data=tab_maize %>% 
                     filter(country_name != "Desert"), 
                   aes(x=reorder(country_name, Ya, median), y=Ya), outlier.shape = NA) + 
      geom_point(data=tab_maize %>% filter(country_name != "Desert") %>% 
                   group_by(country_name) %>% 
                   summarise(Ya_med = median(Ya)), 
                 aes(x=reorder(country_name, Ya_med, max), y=Ya_med, color=country_name)) + 
      theme_cowplot() + 
      theme(legend.position = "none",
            plot.margin = unit(c(0,0,0,0), "cm")) + 
      labs(y = "Yield (tons/ha)", x="") +
      coord_flip(),
    ncol=1, align ="hv", axis="tbrl", rel_heights = c(0.33,0.66) 
  ),
  ncol=2, align ="hv", axis="tbrl" 
)

ggsave(plot = sp2, 
       paste0("E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/Supporting information/02_density_yields.png"),
       #paste0("E:/POSTDOC INRAE/ANALYSES/B_OPTIMISATION/00_graph_evaluation_predictive_models/",crop_i, "_", mod_i, ".png"),
       width=10, height=8, bg="white")


# Mean yield by year and country
# > soybean
ggplot(tab_soybean %>% 
         filter(country_name != "Desert") %>% 
         group_by(country_name, year) %>% 
         summarise(median_Ya = median(Ya)), 
       aes(x=year, y=median_Ya, color=country_name, group=country_name)) + 
  geom_point() + geom_line() +
  theme_bw() + 
  theme(legend.position = "bottom")
# > maize  
ggplot(tab_maize %>% 
         filter(country_name != "Desert") %>% 
         group_by(country_name, year) %>% 
         summarise(median_Ya = median(Ya)), 
       aes(x=year, y=median_Ya, color=country_name, group=country_name)) + 
  geom_point() + geom_line() +
  theme_bw() + 
  theme(legend.position = "bottom")


tab_soybean %>% 
  filter(country_name != "Desert") %>% 
  group_by(country_name) %>% 
  mutate(median_Ya = median(Ya)) %>% 
  mutate(id = if_else(median_Ya > 2.5, "high yields", "low yields")) %>% 
  ggplot(., aes(x=Ya, color=id)) + geom_density()

tab_maize %>% 
  filter(Ya < 20) %>% 
  filter(country_name != "Desert") %>% 
  group_by(country_name) %>% 
  mutate(median_Ya = median(Ya)) %>% 
  mutate(id = if_else(median_Ya > 7.5, "high yields", "low yields")) %>% 
  ggplot(., aes(x=Ya, color=id)) + geom_density()

# ---------------------------------------
# Models performances

load("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_preds/05_GLOBAL/tab_perf_models.rda")

# > Compute mean performance 
tab_perf_labelled_i_mean <- tab_perf_labelled %>% 
  group_by(Outcome, Crop, Model, pred_perf, pred_perf_lab) %>% 
  # MEAN OF BOTH CROSS-VALIDATION
  summarise(mean_pred_perf = mean(pred_perf_value)) %>% 
  ungroup() %>%
  mutate(Type_cv_lab = "c. Mean of cross-validation results")

# > Performances
# NSE
rbind(tab_perf_labelled %>% ungroup() %>% dplyr::select(Outcome, Crop, Model, pred_perf, pred_perf_lab, pred_perf_value, Type_cv_lab),
      tab_perf_labelled_i_mean %>% dplyr::select(Outcome, Crop, Model, pred_perf, pred_perf_lab, 'pred_perf_value'='mean_pred_perf', Type_cv_lab)) %>% 
  mutate(Type_cv_lab = recode(Type_cv_lab, 
                              "Cross-validation on years"="a. Cross-validation on years", 
                              "Cross-validation on sites"="b. Cross-validation on sites")) %>% 
  filter(pred_perf_lab=="NSE\n(higher is better)") %>% 
  spread(key="Type_cv_lab", value="pred_perf_value") %>%
  arrange(desc(Crop)) %>% 
  dplyr::select(-Outcome, -pred_perf, -pred_perf_lab)

#  Crop    Model     `a. Cross-validation on years` `b. Cross-validation on sites` `c. Mean of cross-validation results`
#  soybean pca.m.3                           0.910                          0.939                                 0.925
#  soybean pca.m.2                           0.919                          0.932                                 0.925
#  soybean avg.m                             0.900                          0.942                                 0.921
#  soybean avg.s                             0.913                          0.900                                 0.906
#  maize   pca.m.3                           0.902                          0.911                                 0.907
#  maize   pca.m.2                           0.911                          0.905                                 0.908
#  maize   avg.m                             0.866                          0.917                                 0.892
#  maize   avg.s                             0.920                          0.872                                 0.896

rbind(tab_perf_labelled %>% ungroup() %>% dplyr::select(Outcome, Crop, Model, pred_perf, pred_perf_lab, pred_perf_value, Type_cv_lab),
      tab_perf_labelled_i_mean %>% dplyr::select(Outcome, Crop, Model, pred_perf, pred_perf_lab, 'pred_perf_value'='mean_pred_perf', Type_cv_lab)) %>% 
  mutate(Type_cv_lab = recode(Type_cv_lab, 
                              "Cross-validation on years"="a. Cross-validation on years", 
                              "Cross-validation on sites"="b. Cross-validation on sites")) %>% 
  filter(pred_perf_lab=="RMSEP\n(lower is better)") %>% 
  spread(key="Type_cv_lab", value="pred_perf_value") %>%
  arrange(desc(Crop)) %>% 
  dplyr::select(-Outcome, -pred_perf, -pred_perf_lab)

#  Crop    Model     `a. Cross-validation on years` `b. Cross-validation on sites` `c. Mean of cross-validation results`
#  soybean pca.m.3                          0.413                          0.342                                 0.377
#  soybean pca.m.2                          0.394                          0.361                                 0.377
#  soybean avg.m                            0.437                          0.332                                 0.384
#  soybean avg.s                            0.407                          0.436                                 0.422
#  maize   pca.m.3                          1.59                           1.52                                  1.55 
#  maize   pca.m.2                          1.52                           1.57                                  1.54 
#  maize   avg.m                            1.86                           1.46                                  1.66 
#  maize   avg.s                            1.44                           1.82                                  1.63

# ---------------------------------------
# Models evaluation 

load("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_preds/05_GLOBAL/tab_preds.rda")

# > Observed vs predicted - world
tab_preds_density <- tab_preds %>% 
  # > compute residuals
  mutate(Residuals = Predicted - Observed) %>% 
  # > compute RMSEP and NSE
  group_by(Type_cv, Crop, Model) %>% 
  mutate(
    R2    = caret::R2(obs = Observed,      pred = Predicted),
    RMSEP = caret::RMSE(obs = Observed,      pred = Predicted),
    NSE   = hydroGOF::NSE(obs = Observed,    sim = Predicted)
  ) %>% 
  mutate(Type_cv = recode(Type_cv, 
                          "01_YEARS"="Cross-validation on years",
                          "02_GEO"="Cross-validation on sites"))

plots <- tab_preds_density %>% 
  filter(Model == "pca.m.2") %>% 
  split(.$Crop) %>% 
  map(., ~{ 
    
    .x %>% 
      split(.$Model) %>% 
      map(., ~{ 
        
        # Compute density of points per combination of pred and obs values for each cross-validation strategy
        tab_preds_density_i <- .x %>% 
          split(.$Type_cv) %>% 
          map_dfr(., ~{
            
            density <- get_density(x = .x$Observed, y = .x$Predicted, n=100)
            
            .x$density <- density
            
            .x
            
          }, .id = "Type_cv") %>% 
          dplyr::select(-RMSEP, -NSE, -R2) %>% 
          left_join(tab_preds_density %>% 
                      distinct(Type_cv, Crop, Model, Outcome, RMSEP, NSE, R2), 
                    by=c("Crop", "Model","Outcome", "Type_cv")) %>% 
          mutate(label=paste0("R²=" , round(R2,2), "\nRMSEP=", round(RMSEP,2), " t ha -1"))
        
        mod_i <- unique(tab_preds_density_i$Model)
        crop_i <- unique(tab_preds_density_i$Crop)
        
        
        # > Pred vs Obs
        p.pred.obs <- tab_preds_density_i %>%
          ggplot(data = .) + 
          geom_point(aes(x = Observed, y = Predicted, color = density), 
                     size = 0.5) + 
          geom_abline(color = "red") + 
          geom_text(x=0.5, y=max(tab_preds_density_i$Predicted)-0.75, 
                   aes(label=label), hjust=0) + 
          theme_cowplot() + 
          theme(legend.position = "bottom",
                strip.text = element_text(size=10), 
                #strip.background = element_blank(),
                title = element_text(size=10),
                axis.text = element_text(size=10),
                axis.title.x = element_text(size=10),
                legend.text = element_text(size=10)
          ) + 
          facet_wrap(.~Type_cv, ncol=1) +
          ggtitle("") +
          scale_color_viridis_c(guide = guide_colorbar(barwidth = 10, barheight = 0.5)) +
          labs(x = "Observed yields (t ha -1)", y = "Predicted yields (t ha -1)", color = "Density per combination of\n observed/predicted yield values") 
        
        # > Residuals
        p.res <- tab_preds_density_i %>%
          # mean residuals 
          group_by(Crop, Model, Type_cv) %>% 
          mutate(mean_Residuals = mean(Residuals)) %>% 
          ggplot(.) + 
          geom_histogram(aes(x=Residuals), 
                         bins = 100, fill = "darkgrey", color = "transparent") + 
          geom_vline(aes(xintercept = mean_Residuals), color = "black", 
                     linetype = 2, linewidth = 0.5) +
          theme_cowplot() + 
          theme(legend.position = "none",
                strip.text = element_text(size=10), 
                #strip.background = element_blank(),
                title = element_text(size=10),
                axis.text = element_text(size=10),
                axis.title.x = element_text(size=10)) + 
          facet_wrap(.~Type_cv, ncol=1) +
          ggtitle("") +
          labs(x = "Predicted - observed yields (t ha -1)", y = "Counts") 
        
        list("p.pred.obs"=p.pred.obs,
             "p.res"=p.res)
        
        })
    
    })


sp4_1 <- plot_grid(plots$soybean$pca.m.2$p.pred.obs+ggtitle("a. Soybean"), 
          plots$maize$pca.m.2$p.pred.obs+ggtitle("b. Maize"), 
          nrow=1)

ggsave(plot = sp4_1, 
       "E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/Supporting information/04_model_perf_pca.m.2_1.png",
       width=10, height=10, bg="white")

sp4_2 <- plot_grid(plots$maize$pca.m.2$p.res+ggtitle("a. Soybean"),
                   plots$soybean$pca.m.2$p.res+ggtitle("b. Maize"), 
                   nrow=1)

ggsave(plot = sp4_2, 
       "E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/Supporting information/04_model_perf_pca.m.2_2.png",
       width=10, height=9.5, bg="white")

# --------------------------
# Variables importance

tab_imp <- list(#"mod_maize_avg.m"=mod_maize_avg.m, 
                #"mod_soybean_avg.m"=mod_soybean_avg.m, 
                #"mod_maize_avg.s"=mod_maize_avg.s, 
                #"mod_soybean_avg.s"=mod_soybean_avg.s, 
                #"mod_maize_pca.m.3"=mod_maize_pca.m.3,
                #"mod_soybean_pca.m.3"=mod_soybean_pca.m.3,
                "mod_maize_pca.m.2"=mod_maize_pca.m.2, 
                "mod_soybean_pca.m.2"=mod_soybean_pca.m.2) %>%  
  map_dfr(., ~{
    
    data.frame(Predictor = names(.x$variable.importance), 
               Importance = .x$variable.importance)
    
  }, .id="Model") %>%
  separate(Model, c("mod", "crop", "model"), sep="_", remove=T) %>%
  dplyr::select(-mod) %>% 
  group_by(model, crop) %>%
  arrange(desc(Importance)) %>%
  mutate(Rank = 1:n()) %>%
  mutate(Predictor_lab = recode(Predictor, 
                                "irrigated_portion_perc"="Crop irrigated portion",
                                'PC1_month_prec'='Score 1 - total precipitation', 
                                'PC1_month_min_temp'='Score 1 - min. temperature', 
                                'PC1_month_max_temp'='Score 1 - max. temperature', 
                                'PC1_month_vapor_pressure_deficit'='Score 1 - vapor pressure deficit', 
                                'PC1_month_rad'='Score 1 - radiation', 
                                'PC1_month_et0'='Score 1 - ref. evapotranspiration', 
                                'irrigated_portion'='Irrigation', 
                                'PC2_month_et0'='Score 2 - ref. evapotranspiration', 
                                'PC2_month_rad'='Score 2 - radiation', 
                                'PC2_month_min_temp'='Score 2 - min. temperature', 
                                'PC2_month_vapor_pressure_deficit'='Score 2 - vapor pressure deficit', 
                                'PC2_month_max_temp'='Score 2 - max. temperature', 
                                'PC2_month_prec'='Score 2 - total precipitation',
                                'PC3_month_et0'='Score 3 - ref. evapotranspiration', 
                                'PC3_month_rad'='Score 3 - radiation', 
                                'PC3_month_min_temp'='Score 3 - min. temperature', 
                                'PC3_month_vapor_pressure_deficit'='Score 3 - vapor pressure deficit', 
                                'PC3_month_max_temp'='Score 3 - max. temperature', 
                                'PC3_month_prec'='Score 3 - total precipitation')) %>% 
  mutate(crop = recode(crop, "soybean"="Soybean", "maize"="Maize"), 
         crop=factor(crop, levels=c("Soybean", "Maize")))
  
# Plot
plot_imp <- tab_imp %>% 
  split(.$crop) %>% 
  map(., ~{
    
    .x %>% 
      split(.$model) %>% 
      map(., ~{
        
        # title plot
        mod_i <- unique(.x$model)
        crop_i <- unique(.x$crop)
        
        # plot
        .x %>%
          ggplot(.) + 
          geom_col(aes(x = Importance, y = reorder(Predictor_lab, Rank, min, decreasing = T)), width=0.75) +
          theme_cowplot() + 
          theme(axis.title.y = element_blank(),
                axis.text = element_text(size=10),
                panel.grid.major = element_line(color="lightgrey", linewidth=0.5),
                panel.grid.minor= element_line(color="lightgrey", linewidth=0.25)) +
          ggtitle("")
        
        
      })
    
  })

# Merge 
sp5 <- plot_grid(
  plot_imp$Soybean$pca.m.2,
  plot_imp$Maize$pca.m.2,
  ncol=2, labels = c("a. Soybean", "b. Maize")
)

# Save
ggsave(sp5, 
       filename = paste0("E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/Supporting information/05_variable_imp.png"),
       width = 10, height=5)

# --------------------------------------------------------------------------
# Partial dependency + correlation with monthly averages

# ----------- To do only once --------------
# > name variables and abbrevations
vars_names <- data.frame(clim.var = c("max_2m_temperature", "min_2m_temperature", 
                                      "et0", "surface_net_solar_radiation", 
                                      "total_precipitation", "vapor_pressure_deficit")) %>% 
  mutate(clim.var_abb = recode(clim.var, 
                               "min_2m_temperature"         ="min_temp",
                               "max_2m_temperature"         ="max_temp",
                               "et0"                        ="et0",
                               "surface_net_solar_radiation"="rad",
                               "total_precipitation"        ="prec",
                               "vapor_pressure_deficit"     ="vpd"))  %>% 
  mutate(clim.var_lab = recode(clim.var, 
                               "min_2m_temperature"         ="Minimum temperature (°C)",
                               "max_2m_temperature"         ="Maximum temperature (°C)",
                               "et0"                        ="Evapotransp. ref (mm/day)",
                               "surface_net_solar_radiation"="Net solar radiations (MJ/m²)",
                               "total_precipitation"        ="Total precipitations (mm)",
                               "vapor_pressure_deficit"     ="Vapor pressure deficit")) %>% 
  mutate(clim.var_lab2 = recode(clim.var, 
                                "min_2m_temperature"         ="minimum temperature\n(°C)",
                                "max_2m_temperature"         ="maximum temperature\n(°C)",
                                "et0"                        ="evapotranspiration of reference\n(mm/day)",
                                "surface_net_solar_radiation"="net solar radiations\n(MJ/m²)",
                                "total_precipitation"        ="total precipitations\n(mm)",
                                "vapor_pressure_deficit"     ="vapor pressure deficit\n"))

# > Computation pdp for each predictor of each model
library(pdp) ; library(ranger)
library(parallel) ; library(doParallel); library(foreach)
data_path <- "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/"

# > Soybean
variable.pdp <- list()

# > setting for parallelization
n.cores <- parallel::detectCores() - 3 ; n.cores # 7 on the computer, 11 on the laptop

# >>> create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

# register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

for(var in names(mod_soybean_pca.m.2$variable.importance))
{ 
  # > compute pdp for each variable 
  # (takes long time!)
  variable.pdp[[paste0(var)]] <- mod_soybean_pca.m.2 %>% 
    pdp::partial(pred.var = var, train=tab_soybean) %>% 
    dplyr::rename("value"=1) %>%
    mutate(variable = var)
}

# > transform list into dataframe
variable.pdp.tab <- plyr::ldply(variable.pdp, data.frame, .id = "variable") %>% 
  mutate(clim.var_abb = if_else(variable != "irrigated_portion", substr(variable, 11, nchar(variable)), variable)) %>%
  mutate(clim.var_abb = if_else(clim.var_abb == "vapor_pressure_deficit", "vpd", clim.var_abb)) %>% 
  left_join(., vars_names, by = "clim.var_abb")

# >>> stop cluster//
stopCluster(my.cluster)

save(variable.pdp.tab, file = paste0(data_path, "06_pdp_imp/best_03_WORLD_soybean.rda"))

rm(variable.pdp, variable.pdp.tab)

# > Maize
variable.pdp <- list()

# >>> create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

# register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

for(var in names(mod_maize_pca.m.2$variable.importance))
{ 
  # > compute pdp for each variable 
  # (takes long time!)
  variable.pdp[[paste0(var)]] <- mod_maize_pca.m.2 %>% 
    pdp::partial(pred.var = var, train=tab_maize) %>% 
    dplyr::rename("value"=1) %>%
    mutate(variable = var)
}

# > transform list into dataframe
variable.pdp.tab <- plyr::ldply(variable.pdp, data.frame, .id = "variable") %>% 
  mutate(clim.var_abb = if_else(variable != "irrigated_portion", substr(variable, 11, nchar(variable)), variable)) %>%
  mutate(clim.var_abb = if_else(clim.var_abb == "vapor_pressure_deficit", "vpd", clim.var_abb)) %>% 
  left_join(., vars_names, by = "clim.var_abb")

# >>> stop cluster//
stopCluster(my.cluster)

save(variable.pdp.tab, file = paste0(data_path, "06_pdp_imp/best_03_WORLD_maize.rda"))
rm(variable.pdp, variable.pdp.tab)

# ----------- Start here if PDP already computed -----------
# Correlation between monthly average and PCA scores

# > load PCA scores 
pca_soybean <- loadRDa("C:/Users/benni/Documents/Post doc/ERA5_daily/soybean/pca_soybean.rda")
pca_maize   <- loadRDa("C:/Users/benni/Documents/Post doc/ERA5_daily/maize/pca_maize.rda")

# > load PDP
variable.pdp.tab_soybean <- loadRDa("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/06_pdp_imp/best_03_WORLD_soybean.rda")
variable.pdp.tab_maize   <- loadRDa("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/06_pdp_imp/best_03_WORLD_maize.rda")

# Eigenvalues of the scores associated 
# with each principal component
variable.eig.tab_soybean <- pca_soybean$list_pca_per_variable %>% 
  map_dfr(., ~{
    
    # > Eigenvalues
    data.frame(CP = paste0("Dim.", 1:length(.x$pca.eig[,2])), 
               variance.percent = .x$pca.eig[,2])
    
  }, .id='var') %>%
  mutate(CP_lab = paste0("Score ", substr(CP, nchar(CP), nchar(CP)))) %>% 
  mutate(var = if_else(var == "vapor_pressure_deficit", "vpd", var)) %>% 
  left_join(., vars_names, by = c("var"="clim.var_abb")) 
variable.eig.tab_maize <- pca_maize$list_pca_per_variable %>% 
  map_dfr(., ~{
    
    # > Eigenvalues
    data.frame(CP = paste0("Dim.", 1:length(.x$pca.eig[,2])), 
               variance.percent = .x$pca.eig[,2])
    
  }, .id='var') %>%
  mutate(CP_lab = paste0("Score ", substr(CP, nchar(CP), nchar(CP)))) %>% 
  mutate(var = if_else(var == "vapor_pressure_deficit", "vpd", var)) %>% 
  left_join(., vars_names, by = c("var"="clim.var_abb")) 

# Contribution and correlation of averages climate data
# with the scores associated with principal components
variable.contrib.tab_soybean <- pca_soybean$list_pca_per_variable %>% 
  map_dfr(., ~{
    
    # > Contribution + correlation of each variable to the first dimensions
    .x$pca.var %>% 
      mutate(metric = rownames(.)) %>%
      gather(key = "pca_feature_full", value = "value", -metric) %>%
      separate(col = "pca_feature_full", into = c("pca_feature", "B", "Dim"), remove = T) %>% 
      unite("CP", B:Dim, sep = ".") %>%
      spread(key = "pca_feature", value = "value") 
    
  }, .id='var') %>% 
  mutate(var = if_else(var == "vapor_pressure_deficit", "vpd", var)) %>% 
  mutate(month_lab = paste0("Month ", substr(metric, nchar(metric), nchar(metric)))) %>% 
  mutate(CP_lab = paste0("Score ", substr(CP, nchar(CP), nchar(CP)))) %>% 
  left_join(., vars_names, by = c("var"="clim.var_abb"))

variable.contrib.tab_maize <- pca_maize$list_pca_per_variable %>% 
  map_dfr(., ~{
    
    # > Contribution + correlation of each variable to the first dimensions
    .x$pca.var %>% 
      mutate(metric = rownames(.)) %>%
      gather(key = "pca_feature_full", value = "value", -metric) %>%
      separate(col = "pca_feature_full", into = c("pca_feature", "B", "Dim"), remove = T) %>% 
      unite("CP", B:Dim, sep = ".") %>%
      spread(key = "pca_feature", value = "value") 
    
  }, .id='var') %>% 
  mutate(var = if_else(var == "vapor_pressure_deficit", "vpd", var)) %>% 
  mutate(month_lab = paste0("Month ", substr(metric, nchar(metric), nchar(metric)))) %>% 
  mutate(CP_lab = paste0("Score ", substr(CP, nchar(CP), nchar(CP)))) %>% 
  left_join(., vars_names, by = c("var"="clim.var_abb"))

# Variable importance
variable.imp.tab_soybean <- tab_imp %>%
  filter(model == "pca.m.2") %>% 
  filter(crop == "Soybean") %>%
  mutate(Score_lab = case_when(Predictor %in% c("PC1_month_prec", "PC1_month_min_temp", "PC1_month_max_temp", "PC1_month_vapor_pressure_deficit", "PC1_month_rad", "PC1_month_et0") ~ "Score 1",
                               Predictor %in% c("PC2_month_prec", "PC2_month_min_temp", "PC2_month_max_temp", "PC2_month_vapor_pressure_deficit", "PC2_month_rad", "PC2_month_et0") ~ "Score 2", 
                               Predictor == "irrigated_portion" ~ "Crop irrigated fraction")) 
variable.imp.tab_maize <- tab_imp %>%
  filter(model == "pca.m.2") %>% 
  filter(crop == "Maize") %>%
  mutate(Score_lab = case_when(Predictor %in% c("PC1_month_prec", "PC1_month_min_temp", "PC1_month_max_temp", "PC1_month_vapor_pressure_deficit", "PC1_month_rad", "PC1_month_et0") ~ "Score 1",
                               Predictor %in% c("PC2_month_prec", "PC2_month_min_temp", "PC2_month_max_temp", "PC2_month_vapor_pressure_deficit", "PC2_month_rad", "PC2_month_et0") ~ "Score 2", 
                               Predictor == "irrigated_portion" ~ "Crop irrigated fraction")) 

variable.imp.tab_soybean %>%
  dplyr::select(model, crop, Predictor_lab, Score_lab, Importance, Rank)
variable.imp.tab_maize %>%
  dplyr::select(model, crop, Predictor_lab, Score_lab, Importance, Rank)

# > Link pdp and importance 
variable.pdp.tab_soybean2 <- variable.pdp.tab_soybean %>% 
  # importance values 
  left_join(variable.imp.tab_soybean, by = c("variable" = "Predictor")) %>% 
  # label for irrigation 
  mutate(clim.var = if_else(Predictor_lab == "Crop irrigated portion", "Crop irrigated portion", clim.var),
         clim.var_lab = if_else(Predictor_lab == "Crop irrigated portion", "Crop irrigated portion", clim.var_lab),
         clim.var_lab2 = if_else(Predictor_lab == "Crop irrigated portion", "Crop irrigated portion", clim.var_lab2)) 
variable.pdp.tab_maize2 <- variable.pdp.tab_maize %>% 
  # importance values 
  left_join(variable.imp.tab_maize, by = c("variable" = "Predictor")) %>% 
  # label for irrigation 
  mutate(clim.var = if_else(Predictor_lab == "Crop irrigated portion", "Crop irrigated portion", clim.var),
         clim.var_lab = if_else(Predictor_lab == "Crop irrigated portion", "Crop irrigated portion", clim.var_lab),
         clim.var_lab2 = if_else(Predictor_lab == "Crop irrigated portion", "Crop irrigated portion", clim.var_lab2)) 

# Graphics 
# Soybean
explain_plots_soybean <- variable.pdp.tab_soybean2 %>% 
  mutate(clim.var_lab = recode(clim.var_lab, 
                               "Maximum temperature (°C)"="Maximum temperature",
                               "Minimum temperature (°C)" = "Minimum temperature",
                               "Evapotransp. ref (mm/day)" = "Evapotransp. ref",
                               "Net solar radiations (MJ/m²)"="Net solar radiation",
                               "Total precipitations (mm)" = "Total precipitation")) %>% 
  split(.$clim.var_abb) %>% 
  map(., ~{
    
    .x %>% 
      split(.$Score_lab) %>% 
      map(., ~{
        
        var_i <- unique(.x$clim.var_abb)
        lab_var_i <- unique(.x$clim.var_lab)
        score_i <- unique(.x$Score_lab)
        pred_i <- unique(.x$Predictor_lab)
        imp_i <- round(unique(.x$Importance), 1)
        rank_i <- unique(.x$Rank)
        
        # Contribution and correlation of averages climate data
        # with the scores associated with principal components
        p1 <- variable.contrib.tab_soybean %>%
          filter(var == var_i, 
                 CP_lab == score_i) %>%
          # Add explained variance into the graph
          left_join(variable.eig.tab_soybean, by = c("CP", "CP_lab", "var", "clim.var",  "clim.var_lab", "clim.var_lab2")) %>%
          mutate(CP_lab2 = paste0(CP_lab, " (", round(variance.percent, 1), "%)"),
                 contrib_lab = round(contrib, 1),
                 cor_lab = round(cor, 2)) %>%
          # Lab with contribution 
          mutate(nudge_x_lab = if_else(cor < 0, -1, 1),
                 x_lab = cor + 0.1*nudge_x_lab) %>%
          # reverse y lab order
          mutate(month_lab = factor(month_lab, levels = rev(c('Month 1', 'Month 2', 'Month 3', 'Month 4', 'Month 5', 'Month 6', 'Month 7')))) %>% 
          # recode
          mutate(clim.var_lab = recode(clim.var_lab, 
                                       "Maximum temperature (°C)"="Maximum\ntemperature",
                                       "Minimum temperature (°C)" = "Minimum\ntemperature",
                                       "Evapotransp. ref (mm/day)" = "Evapotransp.\nref",
                                       "Net solar radiations (MJ/m²)"="Net solar\nradiation",
                                       "Total precipitations (mm)" = "Total\nprecipitation",
                                       "Vapor pressure deficit"="Vapor pressure\ndeficit"), 
                 clim.var_lab = paste0(clim.var_lab, "\n\nImportance score:\n", imp_i)) %>% 
          # Plot
          ggplot(., aes(x = month_lab, y = cor, 
                        label = cor_lab)) + 
          geom_hline(yintercept = 0, color="darkgrey", linetype = 2) + 
          geom_col(fill = "lightgrey", color="black", width = 0.75) + 
          geom_text(aes(y = x_lab), color = "black", size=2.5) + 
          theme_cowplot() +
          theme(panel.grid = element_blank(),
                legend.position = "bottom",
                strip.placement = "outside", 
                strip.text.y.left = element_text(angle = 0, size=10),
                axis.title.y = element_blank(),
                axis.title.x = element_text(size=10),
                axis.text = element_text(size=10)) +
          coord_flip() + 
          facet_grid(clim.var_lab~., switch = "y") +
          labs(y = paste0("Correlation between monthly averages and ", score_i)) + 
          lims(y = c(-1.1, 1.1))
        
        # plot PDP
        p2 <- .x %>% 
          ggplot(.) +
          geom_line(aes(x = value, y = yhat)) + 
          theme_cowplot() + 
          theme(axis.title = element_text(size=10),
                axis.text = element_text(size=10), 
                panel.grid.major = element_line(color="lightgrey", linewidth=0.5),
                panel.grid.minor= element_line(color="lightgrey", linewidth=0.25)) +
          labs(x = paste0(score_i), y = "Predicted yield\n(t/ha)") + 
          lims(y = c(1, 2.5))
        
        list("contrib" = p1, 
             "pdp" = p2)
        
        
      })
    
  })

# order in which the plots are displayed
variable.imp.tab_soybean %>% ungroup() %>%
  dplyr::select(Predictor_lab, Rank)

#1 Score 1 - average precipitation       1
#2 Score 1 - vapor pressure deficit      2
#3 Score 1 - min. temperature            3
#4 Score 1 - max. temperature            4
#5 Score 1 - radiation                   5
#6 Score 1 - ref. evapotranspiration     6
#7 Crop irrigated portion                7
#8 Score 2 - ref. evapotranspiration     8
#9 Score 2 - radiation                   9
#10 Score 2 - min. temperature           10
#11 Score 2 - vapor pressure deficit     11
#12 Score 2 - max. temperature           12
#13 Score 2 - average precipitation      13

# > scores 1
sp6 <- plot_grid(
  #plot_grid(
  # contrib
  plot_grid(explain_plots_soybean$prec$`Score 1`$contrib + 
              theme(legend.position = "none") +
              ggtitle("a."), 
            explain_plots_soybean$vpd$`Score 1`$contrib + theme(legend.position = "none"), 
            explain_plots_soybean$min_temp$`Score 1`$contrib + theme(legend.position = "none"), 
            explain_plots_soybean$max_temp$`Score 1`$contrib + theme(legend.position = "none"), 
            explain_plots_soybean$rad$`Score 1`$contrib + theme(legend.position = "none"), 
            explain_plots_soybean$et0$`Score 1`$contrib + theme(legend.position = "none"), 
            ncol=1, rel_heights = c(0.17, 0.15, 0.15, 0.15, 0.15, 0.15)),
  # pdp
  plot_grid(explain_plots_soybean$prec$`Score 1`$pdp +
              ggtitle("b."), 
            explain_plots_soybean$vpd$`Score 1`$pdp, 
            explain_plots_soybean$min_temp$`Score 1`$pdp, 
            explain_plots_soybean$max_temp$`Score 1`$pdp, 
            explain_plots_soybean$rad$`Score 1`$pdp, 
            explain_plots_soybean$et0$`Score 1`$pdp,
            ncol=1, rel_heights = c(0.17, 0.15, 0.15, 0.15, 0.15, 0.15)), 
  # settings
  align="hv", axis="trlb", 
  rel_widths = c(0.55,0.45)
  # legend
  #, plot_grid(get_legend(explain_plots$prec$`Score 1`$contrib), ggplot() + theme_void(), nrow=1),
  #ncol=1, rel_heights = c(0.95, 0.05)
)

ggsave(sp6, 
       filename=paste0("E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/Supporting information/06_pdp_soybean_1.png"), 
       width = 11.5, height=13, dpi=300, bg="white")

# ---------------------------------------------------------
# > score 2
sp7 <- plot_grid(
  #plot_grid(
  # contrib
  plot_grid(explain_plots_soybean$et0$`Score 2`$contrib + 
              theme(legend.position = "none") +
              ggtitle("a."), 
            explain_plots_soybean$rad$`Score 2`$contrib + theme(legend.position = "none"), 
            explain_plots_soybean$min_temp$`Score 2`$contrib + theme(legend.position = "none"), 
            explain_plots_soybean$vpd$`Score 2`$contrib + theme(legend.position = "none"), 
            explain_plots_soybean$max_temp$`Score 2`$contrib + theme(legend.position = "none"), 
            explain_plots_soybean$prec$`Score 2`$contrib + theme(legend.position = "none"), 
            ncol=1, rel_heights = c(0.17, 0.15, 0.15, 0.15, 0.15, 0.15)),
  # pdp
  plot_grid(explain_plots_soybean$et0$`Score 2`$pdp +
              ggtitle("b."), 
            explain_plots_soybean$rad$`Score 2`$pdp, 
            explain_plots_soybean$min_temp$`Score 2`$pdp, 
            explain_plots_soybean$vpd$`Score 2`$pdp, 
            explain_plots_soybean$max_temp$`Score 2`$pdp, 
            explain_plots_soybean$prec$`Score 2`$pdp,
            ncol=1, rel_heights = c(0.17, 0.15, 0.15, 0.15, 0.15, 0.15)), 
  # settings
  align="hv", axis="trlb", 
  rel_widths = c(0.55,0.45)
  # legend
  #, plot_grid(get_legend(explain_plots$prec$`Score 1`$contrib), ggplot() + theme_void(), nrow=1),
  #ncol=1, rel_heights = c(0.95, 0.05)
)

ggsave(sp7, 
       filename=paste0("E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/Supporting information/07_pdp_soybean_2.png"), 
       width = 11.5, height=13, dpi=300, bg="white")

# ---------------------------------------------------------
# Maize
explain_plots_maize <- variable.pdp.tab_maize2 %>% 
  mutate(clim.var_lab = recode(clim.var_lab, 
                               "Maximum temperature (°C)"="Maximum temperature",
                               "Minimum temperature (°C)" = "Minimum temperature",
                               "Evapotransp. ref (mm/day)" = "Evapotransp. ref",
                               "Net solar radiations (MJ/m²)"="Net solar radiation",
                               "Total precipitations (mm)" = "Total precipitation")) %>% 
  split(.$clim.var_abb) %>% 
  map(., ~{
    
    .x %>% 
      split(.$Score_lab) %>% 
      map(., ~{
        
        var_i <- unique(.x$clim.var_abb)
        lab_var_i <- unique(.x$clim.var_lab)
        score_i <- unique(.x$Score_lab)
        pred_i <- unique(.x$Predictor_lab)
        imp_i <- round(unique(.x$Importance), 1)
        rank_i <- unique(.x$Rank)
        
        # Contribution and correlation of averages climate data
        # with the scores associated with principal components
        p1 <- variable.contrib.tab_maize %>%
          filter(var == var_i, 
                 CP_lab == score_i) %>%
          # Add explained variance into the graph
          left_join(variable.eig.tab_maize, by = c("CP", "CP_lab", "var", "clim.var",  "clim.var_lab", "clim.var_lab2")) %>%
          mutate(CP_lab2 = paste0(CP_lab, " (", round(variance.percent, 1), "%)"),
                 cor_lab = round(cor, 2)) %>%
          # Lab with contribution 
          mutate(nudge_x_lab = if_else(cor < 0, -1, 1),
                 x_lab = cor + 0.1*nudge_x_lab) %>%
          # reverse y lab order
          mutate(month_lab = factor(month_lab, levels = rev(c('Month 1', 'Month 2', 'Month 3', 'Month 4', 'Month 5', 'Month 6', 'Month 7', 'Month 8')))) %>% 
          # recode
          mutate(clim.var_lab = recode(clim.var_lab, 
                                       "Maximum temperature (°C)"="Maximum\ntemperature",
                                       "Minimum temperature (°C)" = "Minimum\ntemperature",
                                       "Evapotransp. ref (mm/day)" = "Evapotransp.\nref",
                                       "Net solar radiations (MJ/m²)"="Net solar\nradiation",
                                       "Total precipitations (mm)" = "Total\nprecipitation",
                                       "Vapor pressure deficit"="Vapor pressure\ndeficit"), 
                 clim.var_lab = paste0(clim.var_lab, "\n\nImportance score:\n", imp_i)) %>% 
          # Plot
          ggplot(., aes(x = month_lab, y = cor, 
                        label = cor_lab)) + 
          geom_hline(yintercept = 0, color="darkgrey", linetype = 2) + 
          geom_col(fill = "lightgrey", color="black", width = 0.75) + 
          geom_text(aes(y = x_lab), color = "black", size=2.5) + 
          theme_cowplot() +
          theme(panel.grid = element_blank(),
                legend.position = "bottom",
                strip.placement = "outside", 
                strip.text.y.left = element_text(angle = 0, size=10),
                axis.title.y = element_blank(),
                axis.title.x = element_text(size=10),
                axis.text = element_text(size=10)) +
          coord_flip() + 
          facet_grid(clim.var_lab~., switch = "y") +
          labs(y = paste0("Correlation between monthly averages and ", score_i)) + 
          lims(y = c(-1.1, 1.1))
        
        # plot PDP
        p2 <- .x %>% 
          ggplot(.) +
          geom_line(aes(x = value, y = yhat)) + 
          theme_cowplot() + 
          theme(axis.title = element_text(size=10),
                axis.text = element_text(size=10), 
                panel.grid.major = element_line(color="lightgrey", linewidth=0.5),
                panel.grid.minor= element_line(color="lightgrey", linewidth=0.25)) +
          labs(x = paste0(score_i), y = "Predicted yield\n(t/ha)") + 
          lims(y = c(4.5, 8.5))
        
        list("contrib" = p1, 
             "pdp" = p2)
        
        
      })
    
  })


# Soybean
variable.imp.tab_maize %>% ungroup() %>%
  dplyr::select(Predictor_lab, Rank)

# 1 Score 1 - average precipitation       1
# 2 Crop irrigated portion                2
# 3 Score 1 - min. temperature            3
# 4 Score 1 - ref. evapotranspiration     4
# 5 Score 1 - max. temperature            5
# 6 Score 1 - vapor pressure deficit      6
# 7 Score 1 - radiation                   7
# 8 Score 2 - ref. evapotranspiration     8
# 9 Score 2 - min. temperature            9
#10 Score 2 - radiation                  10
#11 Score 2 - average precipitation      11
#12 Score 2 - vapor pressure deficit     12
#13 Score 2 - max. temperature           13

# > scores 1
sp8 <- plot_grid(
  #plot_grid(
  # contrib
  plot_grid(explain_plots_maize$prec$`Score 1`$contrib + 
              theme(legend.position = "none") +
              ggtitle("a."), 
            explain_plots_maize$min_temp$`Score 1`$contrib + theme(legend.position = "none"), 
            explain_plots_maize$et0$`Score 1`$contrib + theme(legend.position = "none"), 
            explain_plots_maize$max_temp$`Score 1`$contrib + theme(legend.position = "none"), 
            explain_plots_maize$vpd$`Score 1`$contrib + theme(legend.position = "none"), 
            explain_plots_maize$rad$`Score 1`$contrib + theme(legend.position = "none"), 
            ncol=1, rel_heights = c(0.17, 0.15, 0.15, 0.15, 0.15, 0.15)),
  # pdp
  plot_grid(explain_plots_maize$prec$`Score 1`$pdp +
              ggtitle("b."), 
            explain_plots_maize$min_temp$`Score 1`$pdp, 
            explain_plots_maize$et0$`Score 1`$pdp,
            explain_plots_maize$max_temp$`Score 1`$pdp, 
            explain_plots_maize$vpd$`Score 1`$pdp, 
            explain_plots_maize$rad$`Score 1`$pdp, 
            ncol=1, rel_heights = c(0.17, 0.15, 0.15, 0.15, 0.15, 0.15)), 
  # settings
  align="hv", axis="trlb", 
  rel_widths = c(0.55,0.45)
  # legend
  #, plot_grid(get_legend(explain_plots$prec$`Score 1`$contrib), ggplot() + theme_void(), nrow=1),
  #ncol=1, rel_heights = c(0.95, 0.05)
)

ggsave(sp8, 
       filename=paste0("E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/Supporting information/08_pdp_maize_1.png"), 
       width = 11.5, height=13, dpi=300, bg="white")

# > score 2
# 8 Score 2 - ref. evapotranspiration     8
# 9 Score 2 - min. temperature            9
#10 Score 2 - radiation                  10
#11 Score 2 - average precipitation      11
#12 Score 2 - vapor pressure deficit     12
#13 Score 2 - max. temperature           13

# ---------------------------------------------------------
sp9 <- plot_grid(
  #plot_grid(
  # contrib
  plot_grid(explain_plots_maize$et0$`Score 2`$contrib + 
              theme(legend.position = "none") +
              ggtitle("a."), 
            explain_plots_maize$min_temp$`Score 2`$contrib + theme(legend.position = "none"), 
            explain_plots_maize$rad$`Score 2`$contrib + theme(legend.position = "none"), 
            explain_plots_maize$prec$`Score 2`$contrib + theme(legend.position = "none"), 
            explain_plots_maize$vpd$`Score 2`$contrib + theme(legend.position = "none"), 
            explain_plots_maize$max_temp$`Score 2`$contrib + theme(legend.position = "none"), 
            ncol=1, rel_heights = c(0.17, 0.15, 0.15, 0.15, 0.15, 0.15)),
  # pdp
  plot_grid(explain_plots_maize$et0$`Score 2`$pdp +
              ggtitle("b."), 
            explain_plots_maize$min_temp$`Score 2`$pdp, 
            explain_plots_maize$rad$`Score 2`$pdp, 
            explain_plots_maize$prec$`Score 2`$pdp,
            explain_plots_maize$vpd$`Score 2`$pdp, 
            explain_plots_maize$max_temp$`Score 2`$pdp, 
            ncol=1, rel_heights = c(0.17, 0.15, 0.15, 0.15, 0.15, 0.15)), 
  # settings
  align="hv", axis="trlb", 
  rel_widths = c(0.55,0.45)
  # legend
  #, plot_grid(get_legend(explain_plots$prec$`Score 1`$contrib), ggplot() + theme_void(), nrow=1),
  #ncol=1, rel_heights = c(0.95, 0.05)
)

ggsave(sp9, 
       filename=paste0("E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/Supporting information/09_pdp_maize_2.png"), 
       width = 11.5, height=13, dpi=300, bg="white")

# ---------------------------------------------------------
# Spattial patterns of soybean and maize yields
# Breaks for geom_contour_fill, geom_contour, geom_text_contour 
breaks_plot_ratio <- c(0, seq(1, 10, by=0.5), 12, 14)
breaks_labels_ratio <- seq(0,14, by=2)

# Map of average yield predictions in Europe
Ya_pred_eu %>% 
  filter(id_eu27==1) %>%
  # > long format 
  gather(key=year, value=Ya_pred, starts_with("X2")) %>% 
  # > for each pixel, recompute yields (in t/ha)
  group_by(crop, x, y) %>% 
  mutate(Ya_pred_t_ha = Ya_pred/cropland_area_ha) %>%
  summarise(mean_Ya_pred = mean(Ya_pred, na.rm=T)) %>%
  mutate(mean_Ya_pred = ifelse(is.na(mean_Ya_pred)==T, 0, mean_Ya_pred)) %>% 
  spread(crop, mean_Ya_pred) %>% 
  # > compute the ratio between maize yield and soybean yields
  ungroup() %>% 
  mutate(ratio=maize/soybean) %>% 
  # > plot
  ggplot(.) + 
  geom_sf(data=eu27, fill="grey94") +
  geom_contour_fill(aes(x=x, y=y, z=ratio), 
                    breaks = breaks_plot_ratio,
                    na.fill = T, 
                    global.breaks = F,
                    clip = eu27) +
  geom_sf(data=eu27_ext, fill="transparent") +
  theme_map() + 
  lims(x = c(-11,35), y=c(33,71)) + 
  #lims(x = c(-10,35), y=c(33,70)) + 
  theme(legend.position = "bottom",, 
        legend.title = element_text(size=12), 
        legend.text = element_text(size=10),
        strip.text = element_text(face = "bold", hjust = 0.1, size=15)) +
  scale_fill_gradientn(colours = c("transparent", viridis::inferno(direction=-1, n=100)),
                       breaks = breaks_labels_ratio, 
                       labels = breaks_labels_ratio,
                       guide = guide_colorbar(barwidth = 20, barheight = 0.5, title.position = "top", title = "Ratio between maize and soybean yields"))


ggsave(filename = "E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/Supporting information/10_pred_maps_ratio.png", 
       height=10, width=10, bg="white", dpi = 300)






# Allocation depending on other models 

# Run function for allocation
source("E:/POSTDOC INRAE/ANALYSES/B_OPTIMISATION/00_0_Function_for_allocation.R")

# Data with predictions from each model
pred_eu <- Ya_pred_eu %>% 
  spread(key="crop", value="Ya_pred") %>% 
  rename("M_s_i"="soybean", "M_m_i"="maize") %>%
  left_join(dat_coords_EU, by=c("x", "y"))

# Partial Land Equivalent Ratio (pLER) of soybean and maize
# based on the meta-analysis of Xu et al., 2020
num_pLER_s <- 0.56
num_pLER_m <- 0.79

# European Union 27 soybean supply
# Based on FAOSTATS - 2017-2022
num_eu_supply <- 45*10^6

# Simulation plan 
simul_1 <- expand.grid(pLER_s_i = num_pLER_s, 
                       pLER_m_i = num_pLER_m,          
                       freqs_i = 0.25,
                       target_i = c(num_eu_supply*0.5, num_eu_supply*0.75, num_eu_supply)) %>% 
  unite(col = "scenario", c(pLER_s_i, pLER_m_i, freqs_i, target_i), sep = "_", remove = F) %>% 
  mutate(pLER_s_i = as.numeric(as.character(pLER_s_i)),
         pLER_m_i = as.numeric(as.character(pLER_m_i)))

dim(simul_1) # 3 scenarios

# Lists to store outputs
list_alloc_mod_i <- list()
list_maps_alloc_l <- list()

for(mod_i in unique(pred_eu$model))
{
  
  # Prediction data 
  # > first only prediction from the model using the 2 first PC 
  #   from PCA applied on monthly data
  data_simul_1 <- pred_eu %>% 
    filter(model == mod_i)
  
  # Allocation based on each scenario 
  alloc_1 <- simul_1 %>% 
    split(.$scenario) %>% 
    map_dfr(., ~{
      
      # > Retrieve values from the selected scenario
      pLER_s_i_x <- unique(.x$pLER_s_i)
      pLER_m_i_x <- unique(.x$pLER_m_i)
      freqs_i_x  <- unique(.x$freqs_i)
      target_i_x <- unique(.x$target_i)
      
      # > Allocation
      alloc_i  <- alloc_2(data = data_simul_1,
                          #target_supply=num_eu_supply*0.5,  
                          target_supply=target_i_x,  
                          pLER_s = pLER_s_i_x, 
                          pLER_m = pLER_m_i_x, 
                          freq_s = freqs_i_x, freq_m = freqs_i_x, freq_i = freqs_i_x,
                          eu27ext = 0)
      
      # > Out
      alloc_i
      
    }, .id = "scenario")
  list_alloc_mod_i[[paste0(mod_i)]] <- alloc_1
  
  # Table for plot results
  alloc_1_for_plot <- alloc_1 %>% 
    separate(col = "scenario", into = c("pLER_s", "pLER_m", "freq_crop", "target_soybean"), sep = "_", remove = F) %>% 
    # Format labels 
    # > pLER soybean and maize 
    mutate(pLER_lab = paste0(pLER_s, " - ", pLER_m))  %>% 
    # > Crop frequency
    mutate(freq_crop = recode(freq_crop, "0.14"="1 year in 7", "0.16"="1 year in 6", "0.2"="1 year in 5", "0.25"="1 year in 4", "0.33"="1 year in 3", "0.5"="1 year in 2"),
           freq_crop = factor(freq_crop, levels = c("1 year in 7", "1 year in 6", "1 year in 5", "1 year in 4", "1 year in 3", "1 year in 2"))) %>% 
    # > Production target
    mutate(target_soybean_lab = recode(target_soybean, "22500000"="50%", "33750000"="75%", "4.5e+07"="100%"),
           target_soybean_lab = factor(target_soybean_lab, c("50%", "75%", "100%")))
  
  # Maps of soybean allocation in EU 
  maps_alloc_l <- alloc_1_for_plot %>% 
    split(.$scenario) %>%
    map_dfr(., ~{
      
      # Compute surface
      total_A_i_km2 <- round(sum(.x[which(.x$target_supply_intercrop==1),"A_i_i"])*0.01)
      total_A_s_km2 <- round(sum(.x[which(.x$target_supply_monocrop==1),"A_s_i"])*0.01)
      total_A_m_km2 <- round(sum(.x[which(.x$id_saved==1),"A_m_i"])*0.01)
      
      # Set surfaces in a table
      surfaces <- data.frame(
        Crop_Design_lab = c("Intercropping", "Sole soybean", "Sole maize"),
        surfaces        = c(total_A_i_km2,   total_A_s_km2,  total_A_m_km2),
        freq_crop       = rep(unique(.x$freq_crop), 3))
      
      # Join with allocation
      .x %>% 
        gather(key="Crop_Design", value = "Allocated", 
               target_supply_intercrop, target_supply_monocrop, id_saved) %>%
        mutate(Color_dot = case_when(
          Crop_Design == "target_supply_intercrop" & Allocated == 1 ~ "Intercropping",
          Crop_Design == "target_supply_monocrop" & Allocated == 1 ~ "Sole soybean",
          Crop_Design == "id_saved" & Allocated == 1 ~ "Sole maize",
          TRUE ~ "Not allocated")) %>% 
        mutate(Color_dot = factor(Color_dot, levels = c("Intercropping", "Sole soybean", "Sole maize", "Not allocated"))) %>% 
        mutate(Crop_Design_lab = recode(Crop_Design, 
                                        "target_supply_intercrop"="Intercropping",
                                        "target_supply_monocrop" ="Sole soybean",
                                        "id_saved"               ="Sole maize"), 
               Crop_Design_lab = factor(Crop_Design_lab, levels=c("Intercropping", 
                                                                  "Sole soybean", 
                                                                  "Sole maize"))) %>%
        left_join(surfaces, by=c("freq_crop", "Crop_Design_lab")) %>%
        mutate(surfaces_lab = round(as.numeric(as.character(surfaces))/1000,1)) 
      
      
    }, .id = "scenario") %>% 
    split(.$pLER_lab) %>%
    map(., ~{  
      
      .x %>% 
        split(.$target_soybean) %>% 
        map(., ~{
          
          tab_i <- .x %>% 
            mutate(Crop_Design_lab = recode(Crop_Design_lab, 
                                            "Intercropping"="Soybean-maize\nintercopping", 
                                            "Sole soybean" ="Soybean\nin pure stands", 
                                            "Sole maize"="Not allocated\nin pure stands")) %>% 
            mutate(Crop_Design_lab = factor(Crop_Design_lab, levels = c("Soybean-maize\nintercopping", "Soybean\nin pure stands", "Not allocated\nin pure stands")))
          
          # Maps surface
          ggplot() +
            geom_sf(data=world, fill="grey94") +
            geom_point(data = tab_i, 
                       aes(x=x,y=y,color=Color_dot), size=0.1) +
            #geom_text(data = tab_i, x = -5, y = 68, aes(label = surfaces_lab), size= 3) +
            theme_map() + 
            theme(strip.text.y.left = element_text(angle=0, size=10),
                  strip.background = element_rect(fill="lightgrey", color="transparent"),
                  strip.text.x = element_text(size=10)) +
            facet_grid(target_soybean_lab ~ Crop_Design_lab, switch = "y") +
            scale_color_manual(values = c("#0C7BDC", "darkgreen", "darkorange", "transparent"), 
                               guide = guide_none()) +
            scale_alpha_binned(range = c(0,1)) + 
            lims(x = c(-15,35), y=c(33,70)) 
          
          
        })
      
    })
  
  # 3 plots
  sp_map_i <- plot_grid(
    maps_alloc_l$`0.56 - 0.79`$"22500000",
    maps_alloc_l$`0.56 - 0.79`$"33750000"+theme(strip.background.x = element_rect(fill="white"), strip.text.x = element_text(color="transparent")),
    maps_alloc_l$`0.56 - 0.79`$`4.5e+07`+theme(strip.background.x = element_rect(fill="white"), strip.text.x = element_text(color="transparent")),
    ncol=1, labels = c("a.", "b.", "c.")
  )
  
  list_maps_alloc_l[[paste0(mod_i)]] <- maps_alloc_l
  
  # Save
  ggsave(sp_map_i,
         filename = paste0("E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/Supporting information/12_area_req_", mod_i, ".png"), 
         width=9, height=11, bg="white")
  
}


# Barplots total surfaces and productions
# Compute production & area for both crops 
# in intercropping & monocropping 
comp_1 <- list_alloc_mod_i %>% 
  map_dfr(., ~{
    
    .x %>% 
      split(.$scenario) %>% 
      map_dfr(., ~{
        
        
        # > Dataset with intercropping area
        alloc_i <- .x %>% filter(target_supply_intercrop == 1)
        
        # > Dataset monocropping 
        alloc_m <- .x %>% filter(target_supply_monocrop == 1)
        
        # > Dataset maize in monocropping
        alloc_m2 <- .x %>% filter(target_supply_intercrop == 1, target_supply_monocrop == 0)
        
        # - Production target of soybean (should be at least 22500000 t)
        # in intercropping
        Pi_1 <- sum(alloc_i$mean_prod_Y_s_i)
        # in monocropping 
        Pm_1 <- sum(alloc_m$mean_prod_M_s_i)
        
        # - Production maize
        # in intercropping
        Pi_2 <- sum(alloc_i$mean_prod_Y_m_i)
        # in monocropping 
        Pm_2 <- sum(alloc_m2$mean_prod_M_m_i)
        
        # - Surface of soybean needed to reach target production
        # in intercropping
        Ai_1 <- sum(alloc_i$A_i_i)
        # in monocropping
        Am_1 <- sum(alloc_m$A_s_i)
        
        # - Surface of maize 
        # in intercropping 
        Ai_2 <- sum(alloc_i$A_i_i)
        # in monocropping
        Am_2 <- sum(alloc_m2$A_m_i)
        
        # - Average monocropped soybean yield (t/ha) 
        #   ratio between total production / total surface required
        # in intercropping
        Y_1 <- sum(alloc_i$mean_prod_M_s_i)/ sum(alloc_i$A_i_i)
        # in monocropping
        teta_1_Y_1 <- sum(alloc_m$mean_prod_M_s_i)/ sum(alloc_m$A_s_i)
        
        # - Average monocropped maize yield (t/ha) 
        #   in the area allocated to soybean
        # in intercropping
        Y_2 <- sum(alloc_i$mean_prod_M_m_i)/sum(alloc_i$A_i_i)
        # in monocropping
        teta_2_Y_2 <- sum(alloc_m2$mean_prod_M_m_i)/ sum(alloc_m2$A_m_i)
        
        # Factors specifing the difference in average yield of crop 1 between area A_monocropping and A_intercropping
        # teta_1 and teta_2 != 1: yields vary in space 
        #teta_1 <- teta_1_Y_1/Y_1
        #teta_2 <- teta_2_Y_2/Y_2
        
        # > out
        data.frame(productions = c(Pi_1, Pm_1, Pi_2, Pm_2),
                   areas       = c(Ai_1, Am_1, Ai_2, Am_2),
                   yields      = c(Y_1, teta_1_Y_1, Y_2, teta_2_Y_2), 
                   crop        = c("Soybean", "Soybean", "Maize", "Maize"),
                   crop_system = c("Intercropping", "Monocropping", "Intercropping", "Monocropping"))
        
        
      }, .id = "scenario")  %>% 
      separate(col = "scenario", into = c("pLER_s", "pLER_m", "freq_crop", "target_soybean"), sep = "_", remove = F) %>% 
      # Format labels 
      # > pLER soybean and maize 
      mutate(pLER_lab = paste0(pLER_s, " - ", pLER_m))  %>% 
      # > Crop frequency
      mutate(freq_crop = recode(freq_crop, "0.14"="1 year in 7", "0.16"="1 year in 6", "0.2"="1 year in 5", "0.25"="1 year in 4", "0.33"="1 year in 3", "0.5"="1 year in 2"),
             freq_crop = factor(freq_crop, levels = c("1 year in 7", "1 year in 6", "1 year in 5", "1 year in 4", "1 year in 3", "1 year in 2"))) %>% 
      # > Production target
      mutate(target_soybean_lab = recode(target_soybean, "22500000"="50%", "33750000"="75%", "4.5e+07"="100%"),
             target_soybean_lab = factor(target_soybean_lab, c("50%", "75%", "100%")))
    
    
    
  }, .id="model")


# Graphics
sp13 <- plot_grid(
  comp_1 %>% 
    mutate(model = factor(model, levels = c("pca.m.2", "pca.m.3", "avg.m", "avg.s"))) %>% 
    mutate(crop = factor(crop, levels = rev(c("Soybean", "Maize")))) %>% 
    mutate(crop_system = recode(crop_system, "Monocropping"="Pure stands")) %>% 
    mutate(target_soybean=as.numeric(as.character(target_soybean))) %>% 
    ggplot() +
    geom_col(aes(x=model, y = productions/10^6, fill=crop),
             col="black", width=0.5) +
    geom_hline(aes(yintercept=target_soybean/10^6), 
               lty=2, color="darkred") +
    theme_cowplot() + 
    theme(panel.grid.major.y = element_line(linewidth = 0.25, color="lightgrey"),
          legend.position = "bottom",
          strip.placement = "outside", 
          strip.text.y.left = element_text(angle = 0, size=10),
          axis.title = element_text(size=10),
          axis.text = element_text(size=10)) +
    facet_grid(target_soybean_lab~crop_system) +
    scale_fill_manual(values = rev(c("darkgreen", "darkorange")), 
                      guide = guide_legend(reverse = TRUE)) +
    ggtitle("") +
    labs(x = "Soybean self-sufficiency", y = "Total co-production (Mt)", fill="Crop") + 
    lims(y=c(0,210)),
  
  ggplot() + theme_void(),
  
  comp_1 %>% 
    mutate(model = factor(model, levels = c("pca.m.2", "pca.m.3", "avg.m", "avg.s"))) %>% 
    mutate(check=paste0(crop, " - ", crop_system)) %>%
    filter(check!="Maize - Intercropping") %>% 
    mutate(crop = if_else(crop_system == "Intercropping", "Intercropping soybean-maize", crop)) %>%
    mutate(crop = factor(crop, levels = rev(c("Intercropping soybean-maize", "Soybean", "Maize")))) %>% 
    mutate(crop_system = recode(crop_system, "Monocropping"="Pure stands")) %>% 
    ggplot() +
    geom_col(aes(x=model, y = areas/100000, fill=crop),
             col="black", width=0.5) +
    theme_cowplot() + 
    theme(panel.grid.major.y = element_line(linewidth = 0.25, color="lightgrey"),
          legend.position = "bottom",
          strip.placement = "outside", 
          strip.text.y.left = element_text(angle = 0, size=10),
          axis.title = element_text(size=10),
          axis.text = element_text(size=10)) +
    facet_grid(target_soybean_lab~crop_system) +
    scale_fill_manual(values = rev(c("#0C7BDC", "darkgreen", "darkorange")),
                      guide = guide_legend(reverse = TRUE)) +
    labs(x = "Soybean self-sufficiency", y = expression(paste("Area requirements (", 10^{5}, " ha)")), fill="Crop"),
  
  ncol=1, labels = c("a.", "", "b."), rel_heights = c(0.45, 0.1, 0.45)
)

# Save
ggsave(sp13,
       filename = "E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/Supporting information/13_production_surface.png", 
       width=8.5, height=10, bg="white")


# Comparison of intercropping and pure stands stategies 
# in sensitivity analyses 

# load allocation according to 
# all scenarios
load("E:/POSTDOC INRAE/ANALYSES/B_OPTIMISATION/00_Data/sensi.optim.eu.rda")

comp_2 <- alloc_2 %>% 
  #filter(scenario == "0.56_0.79_0.2") %>% 
  split(.$scenario) %>% 
  map_dfr(., ~{
    
    
    # > Dataset with intercropping area
    alloc_i <- .x %>% filter(target_supply_intercrop == 1)
    
    # > Dataset monocropping 
    alloc_m <- .x %>% filter(target_supply_monocrop == 1)
    
    # > Dataset maize in monocropping
    alloc_m2 <- .x %>% filter(target_supply_intercrop == 1, target_supply_monocrop == 0)
    
    # - Production target of soybean (should be at least 22500000 t)
    # in intercropping
    Pi_1 <- sum(alloc_i$mean_prod_Y_s_i)
    # in monocropping 
    Pm_1 <- sum(alloc_m$mean_prod_M_s_i)
    
    # - Production maize
    # in intercropping
    Pi_2 <- sum(alloc_i$mean_prod_Y_m_i)
    # in monocropping 
    Pm_2 <- sum(alloc_m2$mean_prod_M_m_i)
    
    # - Surface of soybean needed to reach target production
    # in intercropping
    Ai_1 <- sum(alloc_i$A_i_i)
    # in monocropping
    Am_1 <- sum(alloc_m$A_s_i)
    
    # - Surface of maize 
    # in intercropping 
    Ai_2 <- sum(alloc_i$A_i_i)
    # in monocropping
    Am_2 <- sum(alloc_m2$A_m_i)
    
    # - Average monocropped soybean yield (t/ha) 
    #   ratio between total production / total surface required
    # in intercropping
    Y_1 <- sum(alloc_i$mean_prod_M_s_i) / sum(alloc_i$A_i_i)
    # in monocropping
    teta_1_Y_1 <- sum(alloc_m$mean_prod_M_s_i)/ sum(alloc_m$A_s_i)
    
    # - Average monocropped maize yield (t/ha) 
    #   in the area allocated to soybean
    # in intercropping
    Y_2 <- sum(alloc_i$mean_prod_M_m_i)/sum(alloc_i$A_i_i)
    # in monocropping
    teta_2_Y_2 <- sum(alloc_m2$mean_prod_M_m_i)/ sum(alloc_m2$A_m_i)
    
    # Factors specifing the difference in average yield of crop 1 between area A_monocropping and A_intercropping
    # teta_1 and teta_2 != 1: yields vary in space 
    #teta_1 <- teta_1_Y_1/Y_1
    #teta_2 <- teta_2_Y_2/Y_2
    
    # > out
    data.frame(productions = c(Pi_1, Pm_1, Pi_2, Pm_2),
               areas       = c(Ai_1, Am_1, Ai_2, Am_2),
               yields      = c(Y_1, teta_1_Y_1, Y_2, teta_2_Y_2), 
               crop        = c("Soybean", "Soybean", "Maize", "Maize"),
               crop_system = c("Intercropping", "Monocropping", "Intercropping", "Monocropping"))
    
    
  }, .id = "scenario") %>% 
  # some scenarios lead to 0 maize production from the sole crop strategy
  mutate(yields = ifelse(is.na(yields)==T, 0, yields)) %>% 
  # format
  separate(col = "scenario", into = c("pLER_s", "pLER_m", "freq_crop", "target_soybean"), sep = "_", remove = F) %>% 
  mutate(pLER_s    = as.numeric(as.character(pLER_s)),
         pLER_m    = as.numeric(as.character(pLER_m)),
         freq_crop = as.numeric(as.character(freq_crop)),
         target_soybean = as.numeric(as.character(target_soybean))/10^6) %>% 
  ungroup() %>%
  # Format and labels 
  # > pLERS
  mutate(pLER_lab = paste0(pLER_s, " - ", pLER_m),
         pLER_lab = if_else(pLER_lab == "0.56 - 0.79", "References values for pLERs\n(soybean: 0.56, maize: 0.79)", pLER_lab)) %>% 
  # > Crop frequencies
  mutate(freq_crop = recode(freq_crop, "0.14"="1 year in 7", "0.16"="1 year in 6", "0.20"="1 year in 5", "0.25"="1 year in 4", "0.33"="1 year in 3", "0.5"="1 year in 2"),
         freq_crop = factor(freq_crop, levels = c("1 year in 7", "1 year in 6", "1 year in 5", "1 year in 4", "1 year in 3", "1 year in 2"))) %>% 
  # > Production target
  mutate(target_soybean_lab = recode(target_soybean, "22.50"="50%", "33.75"="75%", "45.00"="100%"),
         target_soybean_lab = factor(target_soybean_lab, c("50%", "75%", "100%")))

# Heat map of teta values according to pLERs, crop frequency, soybean target production
summary(R_2$teta_2)
limits.min.plots <- 0.7
limits.max.plots <- 1.3

# teta 1
sp10 <- R_2 %>% 
  filter(freq_crop != "1 year in 2") %>% 
  mutate(freq_crop = factor(freq_crop, levels = c('1 year in 2', '1 year in 3', '1 year in 4', '1 year in 5', '1 year in 6', '1 year in 7'))) %>% 
  ggplot(.) + 
  geom_raster(aes(y=pLER_s, x = pLER_m, fill=teta_1)) + 
  geom_point(y = 0.56, x = 0.79, pch = 15, size=2, color="white") + 
  #geom_text(x = 0.5, y = 0.81, label = "Reference values of pLER", check_overlap = T, color="white") + 
  geom_linerange(y=0.56, xmin=0, xmax=0.79, linetype=2, color="white") + 
  geom_linerange(ymin=0, ymax=0.56, x=0.79, linetype=2, color="white") + 
  theme_cowplot() + 
  theme(strip.text.y = element_text(angle=0, size=11),
        strip.text.x = element_text(size=11),
        axis.text = element_text(size=11),
        axis.title = element_text(size=11),
        legend.position = "bottom",
        legend.title = element_text(size=11),
        legend.text = element_text(size=10)) +
  facet_grid(freq_crop~target_soybean_lab, scales = "free") + 
  labs(y = "Partial land equivalent ratio for soybean", 
       x = "Partial land equivalent for maize", 
       fill = "Coefficient of spatial variability of yields\nin intercropping vs pure stands (teta)\nfor soybean") + 
  scale_fill_gradientn(
    colours = c(viridis::mako(10), "white", viridis::rocket(10, direction = -1)),
    values = c(0, (0 - limits.min.plots)/(limits.max.plots - limits.min.plots), 1), 
    limits = c(limits.min.plots, limits.max.plots),
    breaks = c(seq(0.5, 1.5, by = 0.1)),
    na.value = "white",
    guide = guide_colorbar(barwidth = 10, barheight = 0.5))

ggsave(sp10, 
       filename = "E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/Supporting information/10_heat_map_teta1.png", 
       dpi=400, width=9, height=11, bg="white")

# teta 2
sp11 <- R_2 %>% 
  filter(freq_crop != "1 year in 2") %>% 
  mutate(freq_crop = factor(freq_crop, levels = c('1 year in 2', '1 year in 3', '1 year in 4', '1 year in 5', '1 year in 6', '1 year in 7'))) %>% 
  ggplot(.) + 
  geom_raster(aes(y=pLER_s, x = pLER_m, fill=teta_2)) + 
  geom_point(y = 0.56, x = 0.79, pch = 15, size=2, color="white") + 
  geom_linerange(y=0.56, xmin=0, xmax=0.79, linetype=2, color="white") + 
  geom_linerange(ymin=0, ymax=0.56, x=0.79, linetype=2, color="white") + 
  theme_cowplot() + 
  theme(strip.text.y = element_text(angle=0, size=11),
        strip.text.x = element_text(size=11),
        axis.text = element_text(size=11),
        axis.title = element_text(size=11),
        legend.position = "bottom",
        legend.title = element_text(size=11),
        legend.text = element_text(size=10)) +
  facet_grid(freq_crop~target_soybean_lab, scales = "free") + 
  labs(y = "Partial land equivalent ratio for soybean", 
       x = "Partial land equivalent for maize", 
       fill = "Coefficient of spatial variability of yields\nin intercropping vs pure stands (teta)\nfor maize") + 
  scale_fill_gradientn(
    colours = c(viridis::mako(10), "white", viridis::rocket(10, direction = -1)),
    values = c(0, (0 - limits.min.plots)/(limits.max.plots - limits.min.plots), 1), 
    limits = c(limits.min.plots, limits.max.plots),
    breaks = c(seq(0.5, 1.5, by = 0.1)),
    na.value = "white",
    guide = guide_colorbar(barwidth = 10, barheight = 0.5))

ggsave(sp11, 
       filename = "E:/POSTDOC INRAE/PAPERS/03_OPTIMIZATION/01_PNAS/FIGURES/Supporting information/11_heat_map_teta2.png", 
       dpi=400, width=9, height=11, bg="white")

















#ggsave(filename = "E:/POSTDOC INRAE/DATA/02_YIELDS/SASAM/plot_area_cropland_SASAM.png", 
#       width=5, height=5, bg = "white")

# > Geographical distribution
world <- ne_countries(scale = "medium", returnclass = "sf")
plot <- Ya_pred_eu %>% 
  group_by(crop, model, x, y) %>% 
  summarise(mean_Ya_pred = mean(Ya_pred)) %>%
  split(.$crop) %>% 
  purrr::map(., ~{
    .x %>% 
      split(.$model) %>% 
      purrr::map(., ~{ 
        
        .x %>% 
          ggplot(.) + 
          geom_sf(data=world, fill="transparent") +
          geom_point(aes(x=x, y=y, 
                         color = mean_Ya_pred),
                     size=0.4, shape=15) + 
          #facet_grid(.~crop_frequency_lab) + 
          theme_map() + 
          lims(x = c(-15,55), y=c(30,76)) + 
          theme(legend.position = "bottom", 
                legend.text = element_text(size=7)) +
          scale_color_viridis_c(direction = -1, option = "A",
                                name = "Yield (tons/hectare)", 
                                guide = guide_colorbar(barwidth = 15, barheight = 0.5, title.position = "top")) + 
          ggtitle("", paste0(unique(.x$crop)))
        
        
      })
    
  })

t <- read_csv2("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_preds/04_EU/soybean_pred_2000_2023_avg.m_s20.csv")
ggplot(data=t) + 
  geom_sf(data=world, fill="transparent") +
  geom_point(aes(x=x, y=y, 
                 color = mean_Ya_pred),
             size=0.4, shape=15) + 
  #facet_grid(.~crop_frequency_lab) + 
  theme_map() + 
  lims(x = c(-15,55), y=c(30,76)) + 
  theme(legend.position = "bottom", 
        legend.text = element_text(size=7)) +
  scale_color_viridis_c(direction = -1, option = "A",
                        name = "   ", 
                        guide = guide_colorbar(barwidth = 15, barheight = 0.5, title.position = "top"))


# > plot for the paper 
plot_grid(plot_cropland, 
          plot$soybean$pca.m.2 + ggtitle("b. Predicted soybean yields (t/ha)", ""),
          ncol=2)

ggsave(filename = "E:/POSTDOC INRAE/PAPERS/02_GECCO/Figures/plot_area_pred_Ya.png", 
       width=9, height=5, bg = "white")


ggsave(plot$soybean$pca.m.2 + ggtitle("", ""), # + ggtitle("Mean of predicted soybean yields\nin Europe between 2000 and 2023", ""), 
       filename = "E:/POSTDOC INRAE/PAPERS/02_GECCO/Figures/plot_pred_Ya.png", 
       width=7, height=6.5, bg = "white")
# > Mean over 2000-2010
plot_2 <- Ya_pred_eu %>% 
  split(.$crop) %>% 
  purrr::map(., ~{
    .x %>% 
      split(.$model) %>% 
      purrr::map(., ~{ 
        
        .x %>% 
          mutate(year_2 = if_else(year %in% 2000:2011, "2000-2011", "2012-2023")) %>%
          group_by(crop, model, x, y, year_2) %>% 
          summarise(mean_Ya_pred = mean(Ya_pred)) %>%
          ggplot(.) + 
          geom_sf(data=world, fill="transparent") +
          geom_point(aes(x=x, y=y, 
                         color = mean_Ya_pred),
                     size=0.4, shape=15) + 
          stat_contour(aes(x=x, y=y, 
                           z = mean_Ya_pred), color="white", 
                       breaks=c(2.5, 5,7.5,10,12.5,15,17.5,20)) +
          facet_grid(.~year_2) + 
          theme_map() + 
          lims(x = c(-15,55), y=c(30,76)) + 
          theme(legend.position = "bottom", 
                legend.text = element_text(size=7)) +
          scale_color_viridis_c(direction = -1, option = "A",
                                name = "   ", 
                                guide = guide_colorbar(barwidth = 15, barheight = 0.5, title.position = "top")) + 
          ggtitle(paste0(unique(.x$crop)))
        
        
      })
    
  })

# > plot for the paper 
plot_grid(plot_2$maize$pca.m.2,
          plot_2$soybean$pca.m.2,
          ncol=2)

ggsave(filename = "E:/POSTDOC INRAE/PAPERS/02_GECCO/Figures/plot_area_pred_Ya_decades.png", 
       width=10, height=5, bg = "white")



