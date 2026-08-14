# -------------------------------------------------------------------------
#
#       05_1 - Script to produce main figures and tables
#       Author: M. Chen, CIRAD, 2024-2026
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

# ----------------------------------------
# Countries maps (data from Natural Earth public dataset)
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
# > Home-made function to correctly format the results from allocations scenarios 
#   (the function is provided in the "00_0_Function_for_allocation.R" script)
format_alloc <- function(results_allocations)
{
  formatted_result <- results_allocations %>% 
    # format
    separate(col = "scenario", into = c("pLER_s", "pLER_m", "freq_crop", "target_soybean", "max_surf"), sep = "_", remove = F) %>% 
    mutate(pLER_s         = as.numeric(as.character(pLER_s)),
           pLER_m         = as.numeric(as.character(pLER_m)),
           freq_crop      = as.numeric(as.character(freq_crop)),
           target_soybean = as.numeric(as.character(target_soybean))/10^6,
           max_surf       = as.numeric(as.character(max_surf))/10^6) %>% 
    ungroup() %>%
    # Format and labels 
    # > pLERS
    mutate(pLER_lab = paste0(pLER_s, " - ", pLER_m)) %>% 
    # > Crop frequencies
    mutate(freq_crop_lab = recode(freq_crop, "0.14" = "1 year in 7", "0.16" = "1 year in 6", "0.20" = "1 year in 5", "0.25" = "1 year in 4", "0.33" = "1 year in 3", "0.5"  = "1 year in 2"),
           freq_crop_lab = factor(freq_crop_lab, levels = c("1 year in 7", "1 year in 6", "1 year in 5", "1 year in 4", "1 year in 3", "1 year in 2"))) %>% 
    # > Production target
    mutate(target_soybean_lab = recode(target_soybean, "9.075"  = "25%", "18.150" = "50%", "27.225" = "75%", "36.300" = "100%"),
           target_soybean_lab = factor(target_soybean_lab, c("25%", "50%", "75%", "100%"))) %>% 
    # > Maximum surface allocated
    mutate(max_surf_lab = "25% EU cropland area")
  
  return(formatted_result)
  
}

# ----------------------------------------
# Data
# // update based on your architecture \\
path_project <- "..."

# > Predicted yields - Europe
load(paste0(path_project, "/00_DATA/Ya_pred_eu_2000_2023.rda"))

# > Coordinates of sites in the EU
load(paste0(path_project, "/00_DATA/00_dat_coords_EU.rda"))

# > Results of maize-soybean intercropping and sole cropping allocations (from the "03_1_Allocation_europe_main_analysis.R")
load(paste0(path_project, "00_DATA/allocations_soybean_maize_eu_restricted_min_rdt_check.rda"))
# change the name of the object
allocations <- sensi_allocations ; rm(sensi_allocations)

# > Shapping allocations results
# 1. Allocation by pixel 
res1 <- format_alloc(results_allocations =  allocations %>% map_dfr(., ~{.x$data_res}, .id="scenario"))

# 2. Surface required to cover different shares of 
#    EU's consumption of soybean 
res2 <- format_alloc(results_allocations =  allocations %>% map_dfr(., ~{.x$res1}, .id="scenario"))

# 3. Surface required to produce the same production (maize+soybean)
#    as intercropping 
res3 <- format_alloc(results_allocations =  allocations %>% map_dfr(., ~{.x$res3}, .id="scenario"))

# > Results of sensitivity analysis with varying pLERs and crop frequency values
load(paste0(path_project, "00_DATA/sensi/sensi_4_res1.rda"))

# ----------------------------------------
# ---------------------------------------- 
# Figure 1. 
# Legend: Soybean (a) and maize (b) rainfed yield projections (2000-2023 averages), and maize/soybean ratio (c) in the European Union. 
# Subtitle: Projections were obtained from random forest models based on climate inputs (derived from the ERA5-land dataset) covering the crop growing seasons (April-November for soybean and April-December for maize) and irrigated fraction (estimated by SPAM2010, and set to zero for the yield projections). Climate inputs included monthly minimum and maximum temperatures, total precipitation, solar radiation, reference evapotranspiration, and vapor pressure deficit. For each climate variable, the first two scores derived from a principal component analysis were used as climate predictors in the model. For panels (a) and (b), density graphics on the top show the variability in average projected yields, with 25th, 50th (median), and 75th quantiles indicated as vertical yellow, red, and black lines, respectively. Below, the maps represent the spatial distribution of projected average yields, with darker shades indicating regions with higher mean yields, whereas lighter shades representing lower yields (e.g., areas with average crop yield ~ 0 t.ha-1 are displayed in light grey). Regions corresponding to the 25th, 50th (median), and 75th yield quantiles are further delineated by yellow, red, and black contours, respectively. For panel (c), darker shades indicate high ratio between maize and soybean average projected yields, for example in the areas close to the Alps characterized by low soybean yields (i.e., <1 t.ha-1) but moderate-to-high maize yields (i.e., >6 t.ha-1). Base map based on Natural Earth data, created using the R package rnaturalearth.

# Define colors palette 
density_pal <- c("#fcfdbf", "#b73779", "#000004")

# SOYBEAN 
# 1. Distribution of projected yields in the EU
quantiles_soybean <- quantile(mean_Ya_pred_eu[which(mean_Ya_pred_eu$crop=="a."),]$mean_Ya_pred, probs=c(0.25,0.5,0.75)) ; quantiles_soybean
#      25%      50%      75% 
# 1.032215 2.080199 2.569741 
quantiles_crop <- quantiles_soybean
quantiles_crop_lab <- paste(c("25th\nquantile:\n", "Median:\n", "75th quantile:\n"), format(round(quantiles_crop, digits = 1), nsmall = 1))

pa1 <- mean_Ya_pred_eu %>% 
  filter(crop=="a.") %>% 
  ggplot(.) +
  geom_density(aes(x=mean_Ya_pred)) +
  geom_linerange(ymin = 0, ymax = 0.45, x = quantiles_crop[1], color = density_pal[1], linewidth = 0.75) + 
  geom_text(x = quantiles_crop[1], y = 0.5, label = quantiles_crop_lab[1], check_overlap = T, hjust = 0.9, vjust = 0, size=2.5) +
  geom_linerange(ymin = 0, ymax = 0.55, x = quantiles_crop[2], color = density_pal[2], linewidth = 0.75) +
  geom_text(x = quantiles_crop[2], y = 0.575, label = quantiles_crop_lab[2], check_overlap = T, hjust = 0.5, vjust = 0, size=2.5) +
  geom_linerange(ymin = 0, ymax = 0.55, x = quantiles_crop[3], color = density_pal[3], linewidth = 0.75) + 
  geom_text(x = quantiles_crop[3], y = 0.575, label = quantiles_crop_lab[3], check_overlap = T, hjust = 0, vjust = 0, size=2.5) +
  #annotate(geom = "text", x = quantiles_soybean, y = 0.575, label = quantiles_soybean_lab, vjust=0) +
  coord_cartesian(ylim = c(0, 0.7)) +
  scale_y_continuous(breaks = seq(0,0.5,by=0.1)) + 
  theme_cowplot(font_size = 10) +
  theme(plot.margin = margin(1, 0, 0, 0)) +
  labs(x = expression("Mean projected yield t." ~ ha^-1), y = "Density"); pa1

# 2. Spatial variability of yields 
breaks_plot <- c(0, seq(0.5, 3.5, by=0.5), seq(4, 10, by=1))
breaks_labels <- quantile(mean_Ya_pred_eu[which(mean_Ya_pred_eu$crop=="a."),]$mean_Ya_pred, probs=c(0.25,0.5,0.75))

pa2 <- mean_Ya_pred_eu[which(mean_Ya_pred_eu$crop == "a."),] %>% 
  ggplot(data=.) + 
  geom_sf(data=eu27, fill="grey94") +
  geom_contour_fill(aes(x=x, y=y, z=mean_Ya_pred), 
                    breaks = breaks_plot,
                    na.fill = T, 
                    global.breaks = F,
                    clip = eu27) +
  geom_sf(data=eu27_ext, fill="transparent") +
  geom_contour(
    data = data_for_plot,
    aes(x, y, z = mean_Ya_pred),
    breaks = quantiles_crop[1],
    colour = density_pal[1]
  ) +
  geom_contour(
    data = data_for_plot,
    aes(x, y, z = mean_Ya_pred),
    breaks = quantiles_crop[2],
    colour = density_pal[2]
  ) +
  geom_contour(
    data = data_for_plot,
    aes(x, y, z = mean_Ya_pred),
    breaks = quantiles_crop[3],
    colour = density_pal[3]
  ) +
  theme_cowplot(font_size = 10) + 
  lims(x = c(-11,35), y=c(33,71)) + 
  theme(legend.position = "bottom",, 
        legend.title = element_text(size=10), 
        legend.text = element_text(size=10),
        axis.title = element_blank(), 
        axis.text = element_text(size=6),
        plot.margin = margin(0, 0, 0, 0)) +
  scale_fill_gradientn(colours = c("transparent", 
                                   rev(grey.colors(100))),
                       breaks = seq(0, 4.5, by = 1), 
                       labels = seq(0, 4.5, by = 1),
                       guide = guide_colorbar(barwidth = 6, barheight = 0.5, 
                                              title.position = "top", 
                                              title = expression("Mean projected yield t." ~ ha^-1))) ; pa2

# MAIZE
# 1. Distributions
quantiles_maize <- quantile(mean_Ya_pred_eu[which(mean_Ya_pred_eu$crop=="b."),]$mean_Ya_pred, probs=c(0.25,0.5,0.75)) ; quantiles_maize
#      25%      50%      75% 
# 4.304674 6.220117 6.986286
quantiles_crop <- quantiles_maize
quantiles_crop_lab <- paste(c("25th\nquantile:\n", "Median:\n", "75th\nquantile:\n"), format(round(quantiles_crop, digits = 1), nsmall = 1))

pb1 <- mean_Ya_pred_eu %>% 
  filter(crop=="b.") %>% 
  ggplot(.) +
  geom_density(aes(x=mean_Ya_pred)) +
  geom_linerange(ymin = 0, ymax = 0.35, x = quantiles_crop[1], color = density_pal[1], linewidth = 0.75) + 
  geom_text(x = quantiles_crop[1], y = 0.375, label = quantiles_crop_lab[1], check_overlap = T, hjust = 0.9, vjust = 0, size=2.5) +
  geom_linerange(ymin = 0, ymax = 0.35, x = quantiles_crop[2], color = density_pal[2], linewidth = 0.75) +
  geom_text(x = quantiles_crop[2], y = 0.375, label = quantiles_crop_lab[2], check_overlap = T, hjust = 0.9, vjust = 0, size=2.5) +
  geom_linerange(ymin = 0, ymax = 0.45, x = quantiles_crop[3], color = density_pal[3], linewidth = 0.75) + 
  geom_text(x = quantiles_crop[3], y = 0.475, label = quantiles_crop_lab[3], check_overlap = T, hjust = 0, vjust = 0, size=2.5) +
  #annotate(geom = "text", x = quantiles_soybean, y = 0.575, label = quantiles_soybean_lab, vjust=0) +
  coord_cartesian(ylim = c(0, 0.7)) +
  scale_y_continuous(breaks = seq(0,0.5,by=0.1)) + 
  theme_cowplot(font_size = 10) +
  theme(plot.margin = margin(1, 0, 0, 0)) +
  labs(x = expression("Mean projected yield t." ~ ha^-1), y = "Density"); pb1

# 2. Spatial variability
pb2 <- mean_Ya_pred_eu[which(mean_Ya_pred_eu$crop == "b."),] %>% 
  ggplot(data=.) + 
  geom_sf(data=eu27, fill="grey94") +
  geom_contour_fill(aes(x=x, y=y, z=mean_Ya_pred), 
                    breaks = breaks_plot,
                    na.fill = T, 
                    global.breaks = F,
                    clip = eu27) +
  geom_sf(data=eu27_ext, fill="transparent") +
  geom_contour(
    data = data_for_plot,
    aes(x, y, z = mean_Ya_pred),
    breaks = quantiles_crop[1],
    colour = density_pal[1]
  ) +
  geom_contour(
    data = data_for_plot,
    aes(x, y, z = mean_Ya_pred),
    breaks = quantiles_crop[2],
    colour = density_pal[2]
  ) +
  geom_contour(
    data = data_for_plot,
    aes(x, y, z = mean_Ya_pred),
    breaks = quantiles_crop[3],
    colour = density_pal[3]
  ) +
  theme_cowplot(font_size = 10) + 
  lims(x = c(-11,35), y=c(33,71)) + 
  theme(legend.position = "bottom",, 
        legend.title = element_text(size=10), 
        legend.text = element_text(size=10),
        axis.title = element_blank(), 
        axis.text = element_text(size=6),
        plot.margin = margin(0, 0, 0, 0)) +
  scale_fill_gradientn(colours = c("transparent", 
                                   rev(grey.colors(100))),
                       breaks = seq(0, 8.5, by = 1), 
                       labels = seq(0, 8.5, by = 1),
                       guide = guide_colorbar(barwidth = 10, barheight = 0.5, 
                                              title.position = "top", 
                                              title = expression("Mean projected yield t." ~ ha^-1))) ; pb2


# PROJECTED YIELD MAIZE / PROJECTED YIELD SOYBEAN RATIOS 
# Spatial variability of yields ratios 
breaks_plot_ratio <- c(0, seq(1, 8, by=0.5), 12)
breaks_labels_ratio <- c(seq(0, 8, by=2), 12)

p1c <- Ya_pred_eu %>% 
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
  mutate(ratio = case_when(
    soybean < 0.5 & maize < 0.5 ~ 0,
    #soybean < 1 & maize > 5 ~ 0, 
    TRUE ~ maize/soybean)) %>% 
  #mutate(ratio=maize/soybean) %>% 
  mutate(crop="c.") %>%
  # > plot
  ggplot(.) + 
  geom_sf(data=eu27, fill="grey94") +
  geom_contour_fill(aes(x=x, y=y, z=ratio), 
                    breaks = breaks_plot_ratio,
                    na.fill = T, 
                    global.breaks = F,
                    clip = eu27) +
  geom_sf(data=eu27_ext, fill="transparent") +
  theme_cowplot(font_size = 10) + 
  lims(x = c(-11,35), y=c(33,71)) + 
  theme(legend.position = "bottom",, 
        legend.title = element_text(size=10), 
        legend.text = element_text(size=10),
        axis.title = element_blank(), 
        axis.text = element_text(size=6),
        plot.margin = margin(0, 0, 0, 0)) +
  scale_fill_gradientn(colours = c("transparent", viridis::inferno(direction=-1, n=100), "darkgrey"),
                       breaks = breaks_labels_ratio, 
                       labels = breaks_labels_ratio,
                       guide = guide_colorbar(barwidth = 10, 
                                              barheight = 0.5, 
                                              title.position = "top", 
                                              title = "Maize yield /Soybean yield ratio")) ; p1c


# Wrap the 3 panels into 1 
plot_grid(
  plot_grid(pa1, pb1, ggplot() + theme_void() + theme(plot.margin = margin(1, 0, 0, 0)), 
            nrow = 1, labels = c("a.", "b.", "c.")), 
  plot_grid(pa2, pb2, p1c, nrow=1, align="hv"), 
  nrow=2, rel_heights = c(0.3, 0.7))

# Save plot 
ggsave(filename = paste0(path_project, "FIGURES/Figure1.pdf"), 
       height=18, width=21, bg="white", units = "cm", dpi = 300)


# ----------------------------------------
# ---------------------------------------- 
# Table 1. 
# Legend: Soybean and maize coproduction and area requirement for 25, 50, 75, and 100% soybean self-sufficiency in the European Union (EU) achieved by intercropping or sole cropping. 

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

# Save table
save(tab_p4, = paste0(path_project, "FIGURES/Table_1.rda"))

# ----------------------------------------
# ---------------------------------------- 
# Figure 3. 
# Legend: Levels of soybean (a) and maize (b) self-sufficiency in the European Union (EU) achieved from intercropping for different assumptions of partial land equivalent ratio (pLER) and crop return frequency. 
# Subtitle: Crop frequencies of one year in seven, six, five, four, three, or two correspond to allocating maize-soybean intercropping on 14%, 16%, 20%, 25%, 33%, or 50% of cropland area in each grid-cell. Intercropping was first allocated to soybean highest-yielding grid-cells, and in all cases total intercropping area was constrained to not exceed 25% of croplands in the EU (i.e., 25 Mha). Squares represent the results obtained with average productive performances of maize-soybean intercropping systems (i.e., pLER of 0.56 for soybean and 0.79 for maize) estimated in a previous meta-analysis (Xu et al., 2020).

# Format of the results
tab_p <- sensi_4_res1 %>% 
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
  mutate(target_soybean = if_else(crop == "Soybean", target_soybean, 85.1), 
         perc_eu_supply = (production / (target_soybean*10^6)) * 100) %>% 
  # > Change the labels for the plot 
  mutate(crop = recode(crop, "Soybean" = "a.", "Maize" = "b.")) %>% 
  mutate(freq_crop_lab=recode(freq_crop_lab, "1 year in 7"="one-in-seven","1 year in 6"="one-in-six","1 year in 5"="one-in-five","1 year in 4"="one-in-four","1 year in 3"="one-in-three","1 year in 2"="one-in-two")) %>%
  # > Reference values as reference for the dotted lines
  mutate(pLER_lab = case_when(pLER_lab == "0.56 - 0.79" ~ 1, TRUE ~ 0))

# > Plot
plot_grid(
  # > SOYBEAN
  ggplot(tab_p[which(tab_p$crop=="a."),], aes(x = freq_crop_lab, y = perc_eu_supply)) +
    # refs
    geom_hline(yintercept = c(25,50,75,100,125,150), 
               linetype=2, lwd=0.5, color = "grey90") +
    geom_path(aes(color=as.factor(pLER),
                  lty = as.factor(pLER_lab),
                  group=interaction(crop,pLER)), 
              linewidth=0.75) +
    geom_point(aes(color=as.factor(pLER), shape=as.factor(pLER)),
               size=2.5) +
    # current level of self-sufficiency 
    geom_hline(yintercept = 16, color='black', lty=2, linewidth=0.5) +
    annotate("text", x = 6, y = 11.5, hjust=1, 
             label = "Current level in the EU: 16%", 
             color="black", fontface = 'italic') +
    scale_color_manual(values = viridis::magma(13, direction = -1)[2:7]) +
    #scale_color_grey() +
    scale_linetype(guide = guide_none()) +
    scale_y_continuous(breaks = c(25,50,75,100,125,150), limits = c(10,166)) +
    guides(color = guide_legend(title = "\nPartial land equivalent ratio for soybean:", ncol=3),
           shape = guide_none()) +
    facet_wrap(.~crop)+
    theme_cowplot(font_size = 10) +
    theme(strip.text.y.left = element_text(angle=0, size=13, face = "bold"),
          strip.text.x = element_text(size=10, face = "bold", hjust = 0),
          strip.background = element_blank(),
          axis.text = element_text(size=9),
          axis.title = element_text(size=10),
          axis.line = element_line(color="black", linewidth = 0.1),
          panel.border = element_rect(color="black"),
          plot.caption = element_text(hjust = 0),
          legend.position = "bottom",
          legend.title = element_text(size=10),
          legend.title.position = "top",
          legend.text = element_text(size=9),
          legend.background = element_rect(fill="white")) +
    labs(x = "\nYears with intercopping", y = "Self-sufficiency level (%)\n", 
         caption = "Note: partial land equivalent ratio\nfor maize fixed at 0.79"),
  
  # > MAIZE
  ggplot(tab_p[which(tab_p$crop!="a."),], aes(x = freq_crop_lab, y = perc_eu_supply)) +
    # refs
    geom_hline(yintercept = c(25,50,75,100,125,150), 
               linetype=2, lwd=0.5,
               color = "grey90") +
    geom_path(aes(color=as.factor(pLER),
                  lty = as.factor(pLER_lab),
                  group=interaction(crop,pLER)), 
              linewidth=0.75) +
    geom_point(aes(color=as.factor(pLER), shape=as.factor(pLER)),
               size=2) +
    # current level of self-sufficiency 
    geom_hline(yintercept = 81, color='black', lty=2, linewidth=0.5) +
    annotate("text", x = 6, y = 69.5, hjust=1, 
             label = "Current level\n in the EU: 81%", 
             color="black", fontface = 'italic') +
    scale_color_manual(values = viridis::magma(13, direction = -1)[8:13]) +
    scale_linetype(guide = guide_none()) +
    scale_y_continuous(breaks = c(25,50,75,100,125,150), 
                       limits = c(10,166)) +
    guides(color = guide_legend(title = "\nPartial land equivalent ratio for maize:", ncol=3),
           shape = guide_none()) +
    facet_wrap(.~crop)+
    theme_cowplot(font_size = 10) +
    theme(strip.text.y.left = element_text(angle=0, size=13, face = "bold"),
          strip.text.x = element_text(size=10, face = "bold", hjust = 0),
          strip.background = element_blank(),
          axis.text.x = element_text(size=9),
          axis.text.y = element_blank(),
          axis.title.x = element_text(size=10),
          axis.title.y = element_blank(),
          axis.line = element_line(color="black", linewidth = 0.1),
          panel.border = element_rect(color="black"),
          plot.caption = element_text(hjust = 0),
          legend.position = "bottom",
          legend.title = element_text(size=10),
          legend.title.position = "top",
          legend.text = element_text(size=9),
          legend.background = element_rect(fill="white")
    ) +
    labs(x = "\nYears with intercopping", y = "Self-sufficiency level (%)\n", 
         caption = "Note: partial land equivalent ratio\nfor soybean fixed at 0.56"),
  nrow=1, rel_widths = c(0.55, 0.45)
)


ggsave(filename = paste0(path_project, "FIGURES/Figure3.pdf"),
       width=18, height=18, bg="white", units = "cm", dpi=300)

