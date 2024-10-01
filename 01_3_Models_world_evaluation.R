# -------------------------------------------------------------------------
# 
#       01-3. Evaluate predictive models on full dataset
#       Author: M. Chen, Inrae, 2024  
# 
# -------------------------------------------------------------------------

# > Packages
library(tidyverse) ; library(stringr) ; library(lubridate) ; library(CCMHr)
# > machine learning
library(parallel) ; library(doParallel); library(foreach)
library(caret) ; library(ranger) ; library(fastshap) ; library(boot) ; library(coxed)
# > evaluation metrics
library(hydroGOF)

# ---------------------------------- 
# Functions 
# > Load function to cross-validate the models (year-by-year cross-validation)
source("E:/POSTDOC INRAE/ANALYSES/B_OPTIMISATION/00_0_Functions.R")

# ----------------------------------
# Data 
# Soybean
load("E:/POSTDOC INRAE/ANALYSES/B_OPTIMISATION/00_Data/00_tab_soybean.rda")

# Maize
load("E:/POSTDOC INRAE/ANALYSES/B_OPTIMISATION/00_Data/00_tab_maize.rda")

# ----------------------------------
# CROSS-VALIDATION ON YEARS

# > setting for parallelization
n.cores <- parallel::detectCores() - 3 ; n.cores # 7 on the computer, 11 on the laptop

# >>> create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

# register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

system.time({  # estimate run time
  list(
    soybean = list(tab_world = tab_soybean, crop ="soybean"),
    maize   = list(tab_world = tab_maize,   crop ="maize")) %>% 
    map(., ~{
      
      # > data to use for models fitting, name of the analysis, and bootstrap procedure
      dat_pred <- .x$tab_world
      crop     <- .x$crop
      tab_sites <- dat_pred %>% 
        distinct(x, y, gridcode)
      
      # --------------------------------------
      # > FIT MODELS & ASSESS PREDICTIVE PERFORMANCES
      #   BASED ON CROSS-VALIDATION ON YEARS (N years = 35)
      
      # Models list
      list_models <- list(
        pca.m.2      = list(name = "pca.m.2", recompute_scores = TRUE),
        pca.m.3      = list(name = "pca.m.3", recompute_scores = TRUE),
        avg.s        = list(name = "avg.s",   recompute_scores = FALSE),
        avg.m        = list(name = "avg.m",   recompute_scores = FALSE)
      )
      
      # --------------------------------------
      
      # Models fitting and cross-validation on years (full data provided)
      list_fit_cv <- list_models %>% 
        map_dfr(., ~{
          
          # load model 
          mod_i     <- loadRDa(paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_mods/", crop, "_", .x$name, "_train.rda"))
          predictors_i <- paste(names(mod_i$variable.importance), collapse = " + ")
          
          # model formula
          model_formula <- paste0("Ya ~ ", predictors_i)
          model_name    <- .x$name
          
          # recompute scores for pca models 
          recompute_scores_i <- .x$recompute_scores
          
          # > cross-validation & save prediction for each model and site-year
          mod_cv <- function_cv_year(model            = model_formula, 
                                     outcome          = "Ya",
                                     data             = dat_pred, 
                                     model_name       = model_name, 
                                     recompute_scores = recompute_scores_i,
                                     res              = "perf",
                                     save             = TRUE, 
                                     path_save        = paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_preds/05_GLOBAL/", crop, "/01_YEARS"))
          
          
        }) 
      
    })
  
  # >>> stop cluster//
  stopCluster(my.cluster)
  
})

# ----------------------------------
# CROSS-VALIDATION ON SITES

# > setting for parallelization
n.cores <- parallel::detectCores() - 3 ; n.cores # 7 on the computer, 11 on the laptop

# >>> create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

# register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

system.time({  # estimate run time
  list(
    soybean = list(tab_world = tab_soybean, crop ="soybean"),
    maize   = list(tab_world = tab_maize,   crop ="maize")
  ) %>% 
    map(., ~{
      
      # > data to use for models fitting, name of the analysis, and bootstrap procedure
      dat_pred <- .x$tab_world
      crop     <- .x$crop
      tab_sites <- dat_pred %>% 
        distinct(x, y, gridcode)
      
      # --------------------------------------
      # > FIT MODELS & ASSESS PREDICTIVE PERFORMANCES
      #   BASED ON CROSS-VALIDATION ON YEARS (N years = 35)
      
      # Models list
      list_models <- list(
        pca.m.2      = list(name = "pca.m.2", recompute_scores = TRUE),
        pca.m.3      = list(name = "pca.m.3", recompute_scores = TRUE),
        avg.s        = list(name = "avg.s",   recompute_scores = FALSE),
        avg.m        = list(name = "avg.m",   recompute_scores = FALSE)
      )
      
      # --------------------------------------
      
      # Models fitting and cross-validation on years (full data provided)
      list_fit_cv <- list_models %>% 
        map_dfr(., ~{
          
          # load model 
          mod_i     <- loadRDa(paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_mods/", crop, "_", .x$name, "_train.rda"))
          predictors_i <- paste(names(mod_i$variable.importance), collapse = " + ")
          
          # model formula
          model_formula <- paste0("Ya ~ ", predictors_i)
          model_name    <- .x$name
          
          # recompute scores for pca models 
          recompute_scores_i <- .x$recompute_scores
          
          # > cross-validation & save prediction for each model and site-year
          mod_cv <- function_cv_geo(model            = model_formula, 
                                     outcome          = "Ya",
                                     data             = dat_pred, 
                                     model_name       = model_name, 
                                     recompute_scores = recompute_scores_i,
                                     res              = "perf",
                                     save             = TRUE, 
                                     path_save        = paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_preds/05_GLOBAL/", crop, "/02_GEO"))
          
          
        }) 
      
    })
  
  # >>> stop cluster//
  stopCluster(my.cluster)
  
})

# ----------------------------------
# MERGE PREDICTIONS 
data_path <- "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/"

# > Function to read rda into list
rda2list <- function(file) {
  e <- new.env()
  load(file, envir = e)
  as.list(e)
}

# -------------------------------------------------------------------------
# > Load predictions for each model, each country, and both outcomes 
tab_preds <- list(
  soybean = list(name="soybean"),
  maize   = list(name="maize")) %>% 
  map_dfr(., ~{
    
    # -----------------------
    # -----------------------
    # Predictions of yield 
    # > Cross-validated on years
    # > folder where preds are stored
    folder <- paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_preds/05_GLOBAL/", .x$name, "/01_YEARS")
    files <- list.files(folder, pattern = ".rda$")
    
    # > merge all preds for 1 country
    tab_preds <- Map(rda2list, file.path(folder, files)) %>% 
      plyr::ldply(., data.frame, .id = "path") %>% 
      mutate(outcome="01_Ya", 
             type_cv = "01_YEARS") %>% 
      # rename for consistency among cv procedures
      dplyr::rename("preds.fold_cv"="preds.year")
    # -----------------------
    # > Cross-validated on sites
    # > folder where preds are stored
    folder <- paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_preds/05_GLOBAL/", .x$name, "/02_GEO")
    files <- list.files(folder, pattern = ".rda$")
    
    # > merge all preds for 1 country
    tab_preds_geo <- Map(rda2list, file.path(folder, files)) %>% 
      plyr::ldply(., data.frame, .id = "path") %>% 
      mutate(outcome="01_Ya", 
             type_cv = "02_GEO")
    
    # -----------------------
    rbind(tab_preds, tab_preds_geo)
    
  }, .id="Country") %>% 
  dplyr::select(-path) 

# > Rename columns' name
head(tab_preds)
colnames(tab_preds) <- c("Crop", "Fold_cv", "Model", "Site_year", "Observed", "Predicted", "N_predictors", "Outcome", "Type_cv")
names(tab_preds)

# > Compute RMSEP, NSE, R2 and Bias
# for each model in each country

tab_perf <- tab_preds %>% 
  #mutate(Observed = if_else(is.na(Observed) == T, 0, Observed)) %>% 
  group_by(Outcome, Crop, Model, Type_cv) %>% 
  summarise(
    RMSEP = caret::RMSE(obs = Observed,      pred = Predicted),
    NSE   = hydroGOF::NSE(obs = Observed,    sim = Predicted),
    R2    = caret::R2(obs = Observed,        pred = Predicted),
    Bias  = Metrics::bias(actual = Observed, predicted = Predicted),
    N_pred = max(N_predictors)
  )

# > Label models and country for plots
tab_perf_labelled <- tab_perf %>% 
  # set long format 
  pivot_longer(cols=c(RMSEP, NSE, R2, Bias), names_to = "pred_perf", values_to = "pred_perf_value") %>% 
  # > change order in Models
  mutate(Model = factor(Model, levels = rev(c("avg.s", "avg.m", 
                                              "pca.m.2",   "pca.m.3")))) %>% 
  # > type of data (daily, month, annual)
  mutate(data_type=case_when(
    Model %in% c("avg.s") ~ "Seasonal climatic predictors", 
    TRUE ~ "Monthly climatic predictors")) %>% 
  mutate(data_type = factor(data_type, levels = c("Monthly climatic predictors", "Seasonal climatic predictors"))) %>% 
  # > dimension reduction techniques
  mutate(dimension_reduction=case_when(
    Model %in% c('pca.m.2',   'pca.m.3')  ~"PCA",
    TRUE ~ "Averages")) %>% 
  mutate(dimension_reduction = factor(dimension_reduction, levels = c("Averages", "PCA"))) %>% 
  # > crop label 
  mutate(Crop_lab = case_when(
    Crop == "soybean" ~ "Soybean\n(N=122121)", 
    Crop == "maize" ~ "Maize\n(N=75526)"
  )) %>% 
  mutate(Crop_lab = factor(Crop_lab, levels = c("Soybean\n(N=122121)", 
                                                "Maize\n(N=75526)"))) %>% 
  # > outcome
  mutate(Outcome_lab = recode(Outcome, "01_Ya" = "Yield")) %>% 
  # > type cv 
  mutate(Type_cv_lab = recode(Type_cv, "01_YEARS" = "Cross-validation on years", "02_GEO" = "Cross-validation on sites"),
         Type_cv_lab = factor(Type_cv_lab, levels = c("Cross-validation on years", "Cross-validation on sites"))) %>% 
  # > best model per country and indicator
  group_by(Type_cv_lab, Outcome_lab, Crop, pred_perf) %>% 
  mutate(best_pred_perf_value = if_else(pred_perf %in% c("RMSEP", "Bias"), min(abs(pred_perf_value)), max(pred_perf_value))) %>% 
  filter(pred_perf %in% c("RMSEP", "NSE")) %>%
  mutate(pred_perf_lab = case_when(
    pred_perf == "RMSEP" ~ "RMSEP\n(lower is better)", 
    pred_perf == "Bias"  ~ "Bias\n(lower is better)", 
    pred_perf == "NSE"   ~ "NSE\n(higher is better)",
    pred_perf == "R2"    ~ "R squared\n(higher is better)"),
    pred_perf_lab = factor(pred_perf_lab, levels = c("RMSEP\n(lower is better)", 
                                                     "NSE\n(higher is better)",
                                                     "R squared\n(higher is better)",
                                                     "Bias\n(lower is better)"))) %>% 
  # > negative values
  mutate(
    negative_pred_perf_value = if_else(pred_perf_value < 0, 1, 0),
    pred_perf_value_lab = if_else(pred_perf_value < 0, 0, pred_perf_value)) 

# > Save predictions and model performances
save(tab_preds,         file = paste0(data_path, "05_preds/05_GLOBAL/tab_preds.rda"))
save(tab_perf_labelled, file = paste0(data_path, "05_preds/05_GLOBAL/tab_perf_models.rda"))




