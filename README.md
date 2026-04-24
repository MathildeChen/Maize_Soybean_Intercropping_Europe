# High potential contribution of intercropping to soybean and maize self-sufficiency in Europe

## Authors: 

- Mathilde Chen (https://orcid.org/0000-0002-5982-2143) $1,2,3$ 

- Nicolas Guilpart (https://orcid.org/0000-0003-3804-0211), ${4}$  

- David Makowski (https://orcid.org/0000-0001-6385-3703) ${1}$


Affiliations:

$1$  Université Paris-Saclay, INRAE, AgroParisTech, UMR MIA PS, 91120 Palaiseau, France

$2$ CIRAD, UMR PHIM, F-34398 Montpellier, France

$3$ PHIM, Univ Montpellier, CIRAD, INRAE, Institut Agro, IRD, Montpellier, France

$4$ Université Paris-Saclay, AgroParisTech, INRAE, UMR Agronomie, 91120 Palaiseau, France


## Summary

This repository contains scripts supporting a paper aiming to examine the potential for maize-soybean intercropping to improve the European Union’s self-sufficiency in soybean and maize without cropland expansion and by minimizing crop substitution.  

Briefly, we developped several models forecasting soybean and maize yields, respectively, from climate predictors and irrigation practices (see Chen et al., 2024 for more details on the full procedure). The models are trained at the global scale and compared based on their predictive accuracy. For each crop, the best model is choosen and used to project crop yields in the European Union. 

Then, yields projections in the EU are used to simulate several scenarios of soybean and maize allocation, either in sole cropping or in intercropping. For the latter, yields projections are weighted using estimates of partial land equivalent ratios (pLERs) estimated from a global meta-analysis (Xu et al., 2020). Basically these pLERs are the ratios of intercropped and monocropped yield for each species. 

Finally, these scenarios are used to assess the efficiency of intercropping over sole cropping to reach various levels of soybean self-sufficiency in the EU, while not reducing maize self-sufficiency, and limiting the area used. Several sets of sensitivity analyses are performed. 


## Packages required: 
- for data management: *tidyverse*, *stringr*, *lubridate*, *CCMHr*
- for graphics: *cowplot* (plot arrangement), *wesanderson* and *RColorBrewer* (color palettes)
- for spatial data management and maps production: *terra*, *raster*, *rnaturalearth*, *rnaturalearthdata*, *sf*, *sp*, *rworldmap*, *rmapshaper*, *tidygeocoder*, *metR*
- for parallelization of the analyses: *parallel*, *doParallel*, *foreach*
- for random forest models fitting: *ranger*

All analyses were done using R version 4.3.1. 

## Analyses step: 

### Step 0: Prior the analysis: 
- **00_0_Function_for_allocation_3.R**: Home-made functions used to perform the allocation of soybean and maize, as shown in the Figure 2 and Supplementary Material 2 of the paper. 
- **00_0_Functions.R**: home-made functions to 1) estimate the climate predictors based on PCA of the monthly climate variables; 2) to perform the evaluation of the predictive models 
- **00_1_Load_coordinates_EU.R**: selection of the sites where to project soybean and maize in the EU
  

### Step 1: Train the models at the global scale

The scripts used to perform this step are: 
- **01_0_Load_climate_data.R**: derive monthly averages of climate variables (minimum and maximum temperatures, total precipitations, vapor pressure deficit, evapotranspiration, and solar radiations) from the ERA5-land dataset, for the sites with enought soybean production (>1% of the surface dedicated to soybean) and those with 0 soybean (located in regions with unfavorable conditions for crop production) used in the train dataset (at the global scale); 
- **01_1a_Dataset_preparation.R**: merging climate, irrigation, and yield datasets together for the training dataset;
- **01_1b_Dataset_preparation_irrigation_SPAM2020.R**: add gridded irrigation fraction (dataset: SPAM)
- **01_1c_Dataset_preparation_Nrate.R**: add gridded N fertilizer rates (dataset: NPKGRIDS)
- **01_3_Models_world_train.R**: train the models; 
- **01_4_Models_world_evaluation.R**: evaulate the models; Models' prediction accuracy were evaluated using two metrics were used (root mean square error and an equivalent of R², the model efficiency) through two cross-validation procedure unsuring good model transferability in time and in space. 

At the end of this step, the best model for each crop (i.e., soybean and maize) is used to project crop productivity the EU. 

_Note: The data (not provided in this repository) are from the global dataset of historical yields for major crops (GDHY, Iisumi and Sakai, 2020), ERA5-land database (Hersbach et al., 2023), and SPAM dataset on irrigation practices (Yu et al., 2020) which are all freely availble (see references below). The data were on 3447 sites worldwide from 1981 to 2016._ 

References and access to the dataset: 
- Historical yields: Iizumi, T., Sakai, T. The global dataset of historical yields for major crops 1981–2016. Sci Data 7, 97 (2020). https://doi.org/10.1038/s41597-020-0433-7 ; the dataset is accessible here: https://doi.pangaea.de/10.1594/PANGAEA.909132
- Historical climate data: Hersbach, H., Bell, B., Berrisford, P., Biavati, G., Horányi, A., Muñoz Sabater, J., Nicolas, J., Peubey, C., Radu, R., Rozum, I., Schepers, D., Simmons, A., Soci, C., Dee, D., Thépaut, J-N. (2023): ERA5 hourly data on single levels from 1940 to present. Copernicus Climate Change Service (C3S) Climate Data Store (CDS), DOI: 10.24381/cds.adbb2d47 (Accessed in February 2023) ; the data is accessible here: https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=download  
- Irrigation fraction: Yu, Q., You, L., Wood-Sichra, U., Ru, Y., Joglekar, A. K. B., Fritz, S., Xiong, W., Lu, M., Wu, W., and Yang, P.: A cultivated planet in 2010 – Part 2: The global gridded agricultural-production maps, Earth Syst. Sci. Data, 12, 3545–3572, https://doi.org/10.5194/essd-12-3545-2020, 2020. SPAM2010 can be downloaded via an open-data repository (DOI: https://doi.org/10.7910/DVN/PRFF8V; IFPRI, 2019).
- Crop calendars provided by the Agricultural Market Information System: www.amis-outlook.org/amis-about/calendars/soybeancal/en/

## Step 2: Use the best models to projet soybean and maize yield in the EU 

- **02_0_Load_climate_data_EU.R**: derive monthly averages of climate variables (minimum and maximum temperatures, total precipitations, vapor pressure deficit, evapotranspiration, and solar radiations) from the ERA5-land dataset, for the sites in the EU (selected in the 00_1_Load_coordinates_EU.R);
- **02_1_Dataset_preparation_EU.R**: merging climate, irrigation, and yield datasets together for the EU;
- **02_2_Models_europe_predictions.R**: projection of soybean and maize yields in the EU; 

The output of this step is the yield projections in the EU, which are further used for soybean and maize allocation in the EU. 

## Step 3: Scenarios of soybean and maize allocation in the EU

- **03_1_Allocation_europe_main_analyses.R**: this script performs the allocation of soybean and maize in sole or intercropping assuming given pLERs values, a production target of soybean (in t), a maximum surface allocated (25 Mha), the frequence of crop in the rotation (range between 0 and 1, 0 = never and 1 = every year). Here, the scenarios considered included various pLERs for soybean (from 0.3 to 0.7 by 0.1 increment) and for maize (from 0.5 to 0.9 by 0.1 increment), and of cropping frequencies (from 1/7 to 1/2).
- **03_2_Allocation_europe_full_simulations.R**: this script performs sensitivity analyses with finer grid of pLERs (with 0.01 increment); 

## Step 4: Figures production 
These scripts listed below produce the figures and supplementary figures of the paper: 
- **04_Figures.R**: main Figures and supplementary figures

## Step 5: new analyses required during the revision process: 

Following the revision of the article by 3 reviewers, 2 additional scripts have been added: 
- **05_1_Allocation_europe_sensi_analyses_N_use.R**: performs the sensitivity analysis in less favorable conditions for intercropping (i.e., higher N inputs and a null temporal niche differentiation);
- **05_2_Representativity_Xu_et_al_vs_EU.R**: this scripts compares the climate conditions in the sites used to project soybean and maize in the EU and the sites included in the meta-analysis of Xu et al. (2020). 


All these scripts will be shared publicly once the paper has been accepted for publication in a scientific review. 
