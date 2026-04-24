# -------------------------------------------------------------------------
# 
#       01-2. Train predictive models on full dataset
#       Author: M. Chen, Inrae, 2024  
# 
# -------------------------------------------------------------------------

# ----------------------------------
# > Packages
library(tidyverse) ; library(CCMHr)
# > machine learning
library(parallel) ; library(doParallel); library(foreach)
library(caret) ; library(ranger) ; 
# > others
library(hydroGOF)

# -------------------------------------------------------------------------

# ----------------------------------
# Data 
# Soybean
# > add specific path to the project
load(".../data/00_tab_soybean.rda")

# Maize
# > add specific path to the project
load(".../data/00_tab_maize.rda")

# ----------------------------------
# Data to use for models fitting
# > select the crop/dataset 

crop <- "soybean"
#crop <- "maize"

dat_pred <- tab_soybean
#dat_pred <- tab_maize

tab_sites <- dat_pred %>% 
  distinct(x, y, gridcode)

length(unique(tab_sites$gridcode)) 
# Soybean: 3286
# Maize: 2064
 
# ---------------------------------- 
# Train the model on full data set 

# > FIT MODELS ON FULL WORLD DATASET

# Models list 
# top 4 of models based on previous analyses in Guilpart et al. (2022, 10.1038/s43016-022-00481-3) and Chen et al. (2024, 10.1088/1748-9326/ad42b5)
list_models <- list(
  pca.m.2      = list(name = "pca.m.2", formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PC1_month"), starts_with("PC2_month"))), collapse = " + ")),
  pca.m.3      = list(name = "pca.m.3", formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PC1_month"), starts_with("PC2_month"), starts_with("PC3_month"))), collapse = " + ")),
  avg.s        = list(name = "avg.s",   formula = "year_min_2m_temperature + year_max_2m_temperature + year_surface_net_solar_radiation + year_et0 + year_total_precipitation + year_vapor_pressure_deficit"),
  avg.m        = list(name = "avg.m",   formula = paste0(names(dat_pred %>% dplyr::select(starts_with("monthly_"))), collapse = " + "))
)

# > setting for parallelization
n.cores <- parallel::detectCores() - 3 ; n.cores # 7 on the computer, 11 on the laptop

# >>> create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

# register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

# Models fitting and cross-validation on years (full data provided)
list_models %>% 
  map(., ~{
    
    # model formula
    model_formula <- paste0("Ya ~ irrigated_portion_perc + ", .x$formula)
    model_name    <- .x$name
    
    # > fit
    set.seed(101)
    mod  <- ranger(as.formula(model_formula),
                   data=dat_pred, 
                   num.tree=500,
                   importance="impurity") 
    
    # > indicate the correct path where to save the outputs
    save(mod, file = paste0(".../", crop, "_", model_name, "_train.rda"))

})

# >>> stop cluster//
stopCluster(my.cluster)








