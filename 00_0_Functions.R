# -------------------------------------------------------------------------
# 
#       00-0. FUNCTIONS 
#       Author: M. Chen, Inrae, 2023  
# 
# -------------------------------------------------------------------------

# ----------------------------------
# > Packages
library(tidyverse) ; library(stringr) ; library(lubridate) ; library(CCMHr)
# > plots
library(cowplot) ; 
# > maps 
library(terra) ; library(raster) ; library("rnaturalearth") ; library("rnaturalearthdata") ; library(sf) ; library(sp) ; library(rworldmap) ; library(rmapshaper) ; library(tidygeocoder) ; 
# > exploratory analyses
library(arsenal) ; library(mgcv) ; library(multcomp) ; library(corrplot)
# > PCA
library(FactoMineR) ; library(factoextra) 
# > machine learning
library(parallel) ; library(doParallel); library(foreach)
library(caret) ; library(ranger) ; library(fastshap)
# > functional analysis & partial least square
library(fda) ; library(MFPCA) ; library(pls) ; library(fda.usc)
# > others
library(hydroGOF)


# ----------------------------------
# Functions to run PCA
# Argument: 
# - data_for_pca: 

# > Function to run PCA and store results 
do.pca <- function(data_for_pca){
  
  # > Check if data are standardized
  #   Mean of means should be = 0, mean of sd should be = 1
  mean_data <- mean(round(apply(data_for_pca, MARGIN = 2, FUN = mean)))
  sd_data <- mean(round(apply(data_for_pca, MARGIN = 2, FUN = sd)))
  testthat::expect_equal(mean_data, 0)
  testthat::expect_equal(sd_data, 1)
  
  # > Check if any NA
  #   There should be no na
  na_data <- sum(apply(is.na(data_for_pca),2,sum))
  testthat::expect_equal(na_data, 0)
  
  # > Compute PCA
  pca <- prcomp(x = data_for_pca, tol = 0.05)
  
  # > Get results from PCA
  # Eigen values
  pca.eig <- get_eigenvalue(pca)
  
  # Variables results
  res.var <- suppressWarnings(get_pca_var(pca))
  pca.var <- data.frame(
    coord = res.var$coord,        # Coordinates
    contrib = res.var$contrib,    # Axis contribution
    cos2 = res.var$cos2,          # Representation quality 
    cor = res.var$cor             # Correlation
  )
  
  # Individuals results
  res.ind <- get_pca_ind(pca)
  pca.ind <- data.frame(
    ind.coord = res.ind$coord,        # Coordinates
    ind.contrib = res.ind$contrib,    # Axis contribution
    ind.cos2 = res.ind$cos2           # Representation quality
  )
  
  # > Store results
  pca.res <- list(pca = pca,
                  pca.eig = pca.eig,
                  pca.var = pca.var,
                  pca.ind = pca.ind)
  
  return(pca.res)
  
}

# ----------------------------------
# PRINCIPAL COMPONENT ANALYSIS 
# - type_data   : "D" for daily averages, "M" for monthly averages
# - data        : table including climatic features for each site-year 
#                 if type_data = "M", a data.frame with 1 line per site-year and 1 column per climatic monthly averages 
#                 if type_data = "D", a list containing 1 data.frame per climatic feature, in which 1 line corresponds to 1 site-year*date (long-format table)
# - load_data   : should data need to be loaded? if load_data=TRUE, provide to data argument a data.frame containing the site-year ID for checking whether all site-year have data
# - vars_names  : table including the name, full name, and abbreviations of climatic features
# - cum_clim.var: do climatic data need to be accumulated over the growing season (typically for daily data)
# - scale       : do data need to be scaled (default=TRUE)

function_pca <- function(type_data, 
                         data,
                         load_data    = F,
                         vars_names   = vars_names,
                         cum_clim.var = T,
                         scale        = T){
  
  # -------------------------
  # ------------------------- 
  
  # > Object to store the PCA and outputs
  list_pca        <- list()
  list_scores_pca <- list()
  
  # -------------------------
  # -------------------------
  # > Monthly data 
  if(type_data == "M")
  {
    
    # > For each variable 
    for(var_i in unique(vars_names$clim.var))
    {
      
      # > variable abbreviation
      var_i_abb <- vars_names[which(vars_names$clim.var==var_i), "clim.var_abb"]
      message(paste0("pca on ", var_i_abb))
      
      # -------------------------
      # > DATA SELECTION & PREPARATION FOR PCA
      # > select data for the considered variable
      data_var_i <- data %>%  
        dplyr::select(site_year, starts_with(paste0("monthly_", var_i)))
      
      # > check if each site-year has data
      testthat::expect_equal(sort(unique(data_var_i$site_year)), sort(unique(data$site_year)))
      
      #message("ok data")
      
      # > SCALE DATA 
      if(scale == TRUE)
      {
        # > select data for the considered variable and scale it
        z_tab_pca_var_i <- data_var_i %>% 
          dplyr::select(-site_year) %>%
          scale(.) 
      }
      
      if(scale == FALSE)
      {
        # > do not scale data (not recommanded)
        z_tab_pca_var_i <- data_var_i %>% 
          dplyr::select(-site_year) 
        warning("Warning: data used for PCA is not scaled")
        
      }
      # > rownames are sites-year
      rownames(z_tab_pca_var_i) <- data_var_i$site_year
      
      #message("ok scale data")
      
      # -------------------------
      # > PCA
      # > do pca on data 
      pca_var_i <- do.pca(z_tab_pca_var_i)
      
      #message("ok pca")
      
      # > extract scores derived from PCA
      var_i_pca_scores           <- pca_var_i$pca$x
      colnames(var_i_pca_scores) <- paste0("PC", 1:ncol(var_i_pca_scores), "_month_", var_i_abb)
      rownames(var_i_pca_scores) <- data_var_i$site_year
      
      #message("ok scores")
      
      # > STORE PCA AND SCORES  
      list_pca[[paste0(var_i_abb)]]        <- pca_var_i
      list_scores_pca[[paste0(var_i_abb)]] <- var_i_pca_scores
      
    }
    
    # > Extract scores for all variables 
    tab_PCA_scores <- do.call(cbind, list_scores_pca) %>% as.data.frame(.)
    
    # > Return the PCA for each variable and the scores 
    res <- list("tab_PCA_scores"        = tab_PCA_scores,
                "list_pca_per_variable" = list_pca)
    
  }
  
  # -------------------------
  # -------------------------
  # > Daily data
  if(type_data == "D")
  {
    
    # > For each variable 
    for(var_i in unique(vars_names$clim.var)){
      
      # > variable abbreviation
      var_i_abb <- vars_names[which(vars_names$clim.var==var_i), "clim.var_abb"]
      message(paste0("pca on ", var_i_abb))
      
      # -------------------------
      # > LOAD DATA (if load_data=T)
      if(load_data == T)
      {
        
        # > Load climatic data
        era5daily_var_i <- loadRDa(paste0("C:/Users/benni/Documents/Post doc/ERA5_daily/soybean/era5daily_data_", var_i, ".rda"))
        
        # > Check if each site-year has data
        testthat::expect_equal(sort(unique(era5daily_var_i$site_year)), sort(unique(data$site_year)))
        
      }
      if(load_data == F)
      {
        era5daily_var_i <- data[[paste0(var_i)]]
      }
      
      # -------------------------
      # > DATA SELECTION & PREPARATION FOR PCA
      # > select data for the considered variable
      # > if accumulated data over growing season is used
      if(cum_clim.var == TRUE)
      {
        data_var_i <- era5daily_var_i %>% 
          dplyr::select(site_year, day_of_year, clim.var, cum_clim.value) %>% 
          pivot_wider(names_from = c("clim.var", "day_of_year"), 
                      values_from = "cum_clim.value", 
                      names_prefix = "day_", 
                      names_sep = ".")
      }
      
      # > if raw data is used
      if(cum_clim.var == FALSE)
      {
        data_var_i <- era5daily_var_i %>% 
          dplyr::select(site_year, day_of_year, clim.var, clim.value) %>% 
          pivot_wider(names_from = c("clim.var", "day_of_year"), 
                      values_from = "clim.value", 
                      names_prefix = "day_", 
                      names_sep = ".")
        warning("Warning: data used for PCA is not accumulated over growing season")
      }
      
      # > SCALE DATA 
      if(scale == TRUE)
      {
        # > select data for the considered variable and scale it
        z_tab_pca_var_i <- data_var_i %>% 
          ungroup() %>% 
          dplyr::select(starts_with(paste0("day_"))) %>% 
          scale(.) 
      }
      
      if(scale == FALSE)
      {
        # > not scaled data (not recommended)
        z_tab_pca_var_i <- data_var_i %>% 
          ungroup() %>% 
          dplyr::select(starts_with(paste0("day_")))
        warning("Warning: data used for PCA is not scaled")
        
      }
      
      # > rownames are sites-year
      rownames(z_tab_pca_var_i) <- data_var_i$site_year
      
      # -------------------------
      # > PCA
      # > do pca on data 
      pca_var_i <- do.pca(z_tab_pca_var_i)
      
      # > extract scores derived from PCA
      var_i_pca_scores           <- pca_var_i$pca$x
      colnames(var_i_pca_scores) <- paste0("PC", 1:ncol(var_i_pca_scores), "_day_", var_i_abb)
      rownames(var_i_pca_scores) <- data_var_i$site_year
      
      # > STORE PCA AND SCORES  
      list_pca[[paste0(var_i_abb)]]        <- pca_var_i
      list_scores_pca[[paste0(var_i_abb)]] <- var_i_pca_scores
      
    }
    
    # > Extract scores for all variables 
    tab_PCA_scores <- do.call(cbind, list_scores_pca) %>% as.data.frame(.)
    
    # > Return the PCA for each variable and the scores 
    res <- list("tab_PCA_scores"        = tab_PCA_scores,
                "list_pca_per_variable" = list_pca)
  }
  
  # -------------------------
  # -------------------------
  
  return(res)
  
}

# Function to compute new score 
newscores_pca <- function(type_data, 
                          vars_names, 
                          init_data,
                          new_data, 
                          data_day = NULL){
  
  # > Check if data has PLS scores 
  test_name_init_data <- init_data %>% 
    dplyr::select(starts_with("PC"))
  test_name_new_data <- new_data %>% 
    dplyr::select(starts_with("PC"))
  
  if(ncol(test_name_init_data) != 0){ stop("Error: the initial data already contains PC scores.") }
  if(ncol(test_name_new_data) != 0){ stop("Error: the new data already contains PC scores.") }
  
  # ----- Monthly data ------
  if(type_data=="M")
  {
    
    # > PCA on initial data 
    init_pca <- function_pca(type_data  = "M", 
                             vars_names = vars_names,
                             data       = init_data)
    
    # > List to store the future scores 
    list_scores_pca <- list()
    
    # > Apply PCA loads on new data 
    for(var_j in unique(vars_names$clim.var)) 
    {
      
      # > Variable j 
      clim.var_abb_j <- vars_names[which(vars_names$clim.var==var_j),]$clim.var_abb
      
      # > Loads from PCA from initial data for variable j 
      loads_pca_j <- init_pca$list_pca_per_variable[[paste0(clim.var_abb_j)]]$pca$rotation
      
      # > Mean and sd of the original data to standardize new data from these values
      mu_init_data <- init_data %>% ungroup(.) %>%
        dplyr::select(starts_with(paste0("monthly_", var_j))) %>% 
        colMeans(.)
      sd_init_data <- init_data %>% ungroup(.) %>%
        dplyr::select(starts_with(paste0("monthly_", var_j))) %>% 
        apply(., 2, sd)
      
      # > Scale new data based on mean and sd of the original data 
      new_data_j <- new_data %>% ungroup(.) %>%
        dplyr::select(site_year, starts_with(paste0("monthly_", var_j)))
      z_new_data_j <- scale(new_data_j[,-1], center = mu_init_data, scale=sd_init_data)
      
      # > PCA scores from newdata
      list_scores_pca_j <- list()
      for(k in 1:ncol(loads_pca_j))
      {
        
        list_scores_pca_j[[paste0("PC", k, "_month_")]] <- 
          data.frame(site_year = new_data_j$site_year,
                     score = sapply(1:ncol(z_new_data_j), function(x) z_new_data_j[,x] * loads_pca_j[x,k] ) %>% apply(., 1, sum))
        
      }
      
      # > Store
      list_scores_pca[[paste0(clim.var_abb_j)]] <- plyr::ldply(list_scores_pca_j, data.frame, .id = "PC_id") 
      #colnames(tab_scores_plsr_j) <- paste0("PLS", 1:ncol(tab_scores_plsr_j),"_month_",clim.var_abb_j)
      
    }
    
    # > Scores for all variables 
    tab_scores_pca <- plyr::ldply(list_scores_pca, data.frame, .id="var_name") %>% 
      unite(col = "PC_name", c("PC_id", "var_name"), sep = "", remove=T) %>% 
      pivot_wider(names_from = c("PC_name"), values_from = "score")
    
    # > Merge with initial data 
    init_data_pca <- init_data  %>% 
      cbind(., init_pca$tab_PCA_scores)
    
    # > Merge with new data
    new_data_pca <- new_data %>% 
      cbind(., tab_scores_pca)
    
    # > Output 
    res <- list(init_data = init_data_pca, 
                new_data  = new_data_pca)
    
  }
  # ----- Daily data ------
  if(type_data=="D")
  {
    
    # > Detect data for days 
    if(is.null(init_data_day)==T){ stop("Error: No initial daily data provided (init_data_day = NULL)")}
    if(is.null(new_data_day)==T){  stop("Error: No new daily data provided (new_data_day = NULL)")}
    
    # > Daily train and test data 
    data_day_init <- data_day %>% 
      map(., ~ {
        .x %>% filter(site_year %in% unique(init_data$site_year))  
      })
    
    data_day_new <- data_day %>% 
      map(., ~ {
        .x %>% filter(site_year %in% unique(new_data$site_year))  
      })
    
    message("refit PCA scores for daily data")
    
    # > Apply PCA loads on new data 
    init_pca <- function_pca(type_data  = "D",
                             vars_names = vars_names,
                             data       = data_day_init)
    
    # > List to store the future scores 
    list_scores_pca <- list()
    
    # > Apply PLSR loads on new data 
    for(var_j in unique(vars_names$clim.var)) 
    {
      
      # > Variable j 
      clim.var_abb_j <- vars_names[which(vars_names$clim.var==var_j),]$clim.var_abb
      
      # > Loads from PCA from initial data for variable j 
      loads_pca_j <- init_pca$list_pca_per_variable[[paste0(clim.var_abb_j)]]$pca$rotation
      
      # > Select daily data for the var_j 
      init_data_j <- data_day_init[[paste0(var_j)]]
      new_data_j  <- data_day_new[[paste0(var_j)]]
      
      # > Prepare data for scaling 
      tab_init_data_j <- init_data_j %>% 
        dplyr::select(site_year, Ya, day_of_year, clim.var, cum_clim.value) %>% 
        pivot_wider(names_from   = c("clim.var", "day_of_year"), 
                    values_from  = "cum_clim.value", 
                    names_prefix = "day_", 
                    names_sep    = ".") %>% 
        ungroup(.)
      
      # > Mean and sd of the original data to standardize new data 
      mu_init_data <- tab_init_data_j %>% ungroup(.) %>%
        dplyr::select(starts_with(paste0("day_", var_j))) %>% 
        colMeans(.)
      sd_init_data <- tab_init_data_j %>% ungroup(.) %>%
        dplyr::select(starts_with(paste0("day_", var_j))) %>% 
        apply(., 2, sd)
      
      # > Scale new data based on mean and sd of the original data 
      tab_new_data_j <- new_data_j %>% 
        dplyr::select(site_year, day_of_year, clim.var, cum_clim.value) %>% 
        pivot_wider(names_from   = c("clim.var", "day_of_year"), 
                    values_from  = "cum_clim.value", 
                    names_prefix = "day_", 
                    names_sep    = ".") %>% 
        ungroup(.)
      z_tab_new_data_j <- scale(tab_new_data_j[,-1], center = mu_init_data, scale=sd_init_data)
      
      # > PCA score from scaled newdata
      list_scores_pca_j <- list()
      for(k in 1:ncol(loads_pca_j))
      {
        
        list_scores_pca_j[[paste0("PC", k, "_day_")]] <- 
          data.frame(site_year = tab_new_data_j$site_year,
                     score = sapply(1:ncol(z_tab_new_data_j), function(x) z_tab_new_data_j[,x] * loads_pca_j[x,k] ) %>% apply(., 1, sum))
        
      }
      
      # > Store
      list_scores_pca[[paste0(clim.var_abb_j)]] <- plyr::ldply(list_scores_pca_j, data.frame, .id = "PC_id")
      
    }
    
    # > Scores for all variables
    tab_scores_pca <- plyr::ldply(list_scores_pca, data.frame, .id="var_name") %>% 
      unite(col = "PC_name", c("PC_id", "var_name"), sep = "", remove=T) %>% 
      pivot_wider(names_from = c("PC_name"), values_from = "score")
    
    # > Merge with initial data 
    init_data_pca <- init_data  %>% 
      cbind(., init_pca$tab_PCA_scores)
    
    # > Merge with new data
    new_data_pca <- new_data %>% 
      cbind(., tab_scores_pca)
    
    # > Output 
    res <- list(init_data = init_data_pca, 
                new_data  = new_data_pca)
    
  }
  
  # > Out
  return(res)
  
}

# ----------------------------------
# WRAPPER FOR COMPUTING NEW SCORES (for models based on PCA)

newscores_wraper <- function(model_init    = NULL, 
                             init_data     = NULL, 
                             pred_data     = NULL, 
                             init_data_day = NULL, 
                             pred_data_day = NULL,
                             outcome       = NULL, 
                             vars_names    = NULL){
  
  # Checks
  if(is.null(model_init)) { stop("No model provided, fill 'model_init' argument") }
  if(is.null(init_data) | is.null(pred_data)){ stop("No initial or predicted data provided, fill arguments 'init_data' or 'pred_data'") }
  if(is.null(vars_names)){
    vars_names <- data.frame(clim.var = c("max_2m_temperature", "min_2m_temperature", 
                                          "et0", "surface_net_solar_radiation", 
                                          "total_precipitation", "vapor_pressure_deficit")) %>% 
      mutate(clim.var_abb = recode(clim.var, 
                                   "min_2m_temperature"         ="min_temp",
                                   "max_2m_temperature"         ="max_temp",
                                   "et0"                        ="et0",
                                   "surface_net_solar_radiation"="rad",
                                   "total_precipitation"        ="prec",
                                   "vapor_pressure_deficit_1"   ="vpd"))  %>% 
      mutate(clim.var_lab = recode(clim.var, 
                                   "min_2m_temperature"         ="Minimum temperature",
                                   "max_2m_temperature"         ="Maximum temperature",
                                   "et0"                        ="Evapotranspiration ref",
                                   "surface_net_solar_radiation"="Solar radiation",
                                   "total_precipitation"        ="Precipitation",
                                   "vapor_pressure_deficit_1"   ="Vapor pressure deficit"))
  }
  
  # --------------------------------------
  # Compute scores from train dataset and predict for test data set 
  # --------------------------------------
  # PCA: 
  # > monthly data
  if(model_init %in% c("pca.m.1", "pca.m.2", "pca.m.3", "pca.m.all"))
  {
    
    message("refit PCA scores for monthly data")
    
    # > Compute new scores 
    new_scores_pca_M <- newscores_pca(type_data = "M", 
                                      vars_names = vars_names, 
                                      init_data = init_data %>% dplyr::select(-starts_with("PC")),
                                      new_data =  pred_data  %>% dplyr::select(-starts_with("PC")), 
                                      data_day = NULL)
    
    # > Data train and test with new scores
    new_init_data  <- new_scores_pca_M$init_data
    new_pred_data  <- new_scores_pca_M$new_data
    
  }
  # > daily data
  if(model_init %in% c("pca.d.1", "pca.d.2", "pca.d.3", "pca.d.all"))
  {
    
    message("refit PCA scores for daily data")
    
    # > Compute new scores 
    new_scores_pca_D <- newscores_pca(type_data  = "M", 
                                      vars_names = vars_names, 
                                      init_data  = init_data %>% dplyr::select(-starts_with("PC")),
                                      new_data   = pred_data %>% dplyr::select(-starts_with("PC")),
                                      data_day   = init_data_day)
    
    # > Data train and test with new scores
    new_init_data  <- new_scores_pca_D$init_data
    new_pred_data  <- new_scores_pca_D$new_data
    
  }
  
  # --------------------------------------
  # PLSR: compute scores from train dataset and predict for test data set 
  # > monthly data 
  if(model_init %in% c("pls.m.1", "pls.m.2", "pls.m.3", "pls.m.all"))
  {
    
    if(is.null(outcome)) { stop("No outcome provided, fill 'outcome' argument (either Ya or Ya_ano)") }
    
    message("refit PLSR scores for monthly data")
    
    # > Compute new scores 
    new_scores_plsr2_M <- newscores_plsr2(type_data     = "M", 
                                          vars_names    = vars_names, 
                                          init_data     = init_data %>% dplyr::select(-starts_with("PLS")),
                                          new_data      = pred_data %>% dplyr::select(-starts_with("PLS")), 
                                          init_data_day = NULL, 
                                          new_data_day  = NULL, 
                                          outcome       = outcome)
    
    # > Data train and test with new scores
    new_init_data  <- new_scores_plsr2_M$init_data
    new_pred_data  <- new_scores_plsr2_M$new_data
    
  }
  # > daily data 
  if(model_init %in% c("pls.d.1", "pls.d.2", "pls.d.3", "pls.d.all"))
  {
    
    if(is.null(outcome)) { stop("No outcome provided, fill 'outcome' argument (either Ya or Ya_ano)") }
    
    message("refit PLSR scores for daily data")
    
    # > Compute new scores
    new_scores_plsr2_D <- newscores_plsr2(type_data     = "D", 
                                          vars_names    = vars_names, 
                                          init_data     = init_data %>% dplyr::select(-starts_with("PLS")),
                                          new_data      = pred_data %>% dplyr::select(-starts_with("PLS")),
                                          init_data_day = init_data_day, 
                                          new_data_day  = pred_data_day, 
                                          outcome       = outcome)
    
    # > Data train and test with new scores
    new_init_data  <- new_scores_plsr2_D$init_data
    new_pred_data  <- new_scores_plsr2_D$new_data
    
  }
  
  # Output
  res <- list(new_init_data = new_init_data, 
              new_pred_data = new_pred_data)
  return(res)
  
  
}


# ----------------------------------
# CROSS VALIDATION ON YEARS
# res: "preds" or "perf"
# model_gpe: "rf or "lm",
# outcome: name of the predicted variable "Ya" or "Ya_ano"
function_cv_year <- function(model      = NULL, 
                             outcome    = NULL,
                             model_gpe  = "rf",
                             recompute_scores = FALSE,
                             data, 
                             data_day,
                             seed=101, 
                             model_name = "Give_me_a_name",
                             save = F,
                             path_save = NULL, 
                             res = "preds"){
  
  # Data to store predicted values
  data_pred <- list()
  
  # CROSS VALIDATION ON YEARS
  for(y in unique(data$year))
  {
    
    # Name of the model 
    model_init <- model_name
    
    # ---------------------
    # Define test and train datasets
    data_train_0 <- data[which(data$year != y),]
    data_test_0 <- data[which(data$year == y),]
    
    # --------------------------------------
    # If necessary, compute scores from train dataset and predict for test data set 
    if(recompute_scores == TRUE)
    {
      
      new_scores <- newscores_wraper(model_init = model_init, 
                                     init_data  = data_train_0, 
                                     pred_data  = data_test_0, 
                                     #data_day   = data_day, 
                                     outcome    = outcome)
      
      # Test and train datasets with new scores
      data_train <- new_scores$new_init_data
      data_test  <- new_scores$new_pred_data
    }
    if(recompute_scores == FALSE)
    {
      
      # Test and train datasets provided 
      data_train <- data_train_0
      data_test  <- data_test_0
    }
    
    
    # ---------------------
    # Fit model on the train dataset
    # random forest
    if(model_gpe == "rf")
    {
      # > fit on train
      set.seed(seed)
      mod_train <- ranger(as.formula(model),
                          data=data_train, 
                          num.tree=500,
                          importance="impurity")   
      
      # > predict on test
      pred.test <- predict(mod_train, 
                           data = data_test, 
                           type = "response", 
                           seed = seed, 
                           num.trees = 500)
      # > retrieve predictions & nb of predictors
      pred.test_vec <- pred.test$predictions
      N_predictors <- mod_train$num.independent.variables
      
    }
    # linear regression (benchmark)
    if(model_gpe == "lm")
    {
      # > fit on train
      mod_train <- lm(as.formula(model),
                      data=data_train)   
      
      # > predict on test
      pred.test <- predict(mod_train,
                           newdata = data_test, 
                           type = "response")
      # > retrieve predictions
      pred.test_vec <- as.numeric(as.character(pred.test))
      N_predictors <- length(mod_train$coefficients) - 1
    }
    
    
    # ---------------------
    # Store the predictions
    data_pred[[paste0(y)]] <- data.frame(Model     = model_name,
                                         site_year = data_test$site_year,
                                         Ya_obs    = data_test[paste0(outcome)],
                                         Ya_pred   = pred.test_vec,
                                         N_predictors = N_predictors) %>% 
      dplyr::rename("Ya_obs" = 3)
    
  }
  
  # > predictions in a table
  preds <- plyr::ldply(data_pred, data.frame, .id = "year")
  
  # > save result
  if(save==TRUE)
  {
    save(preds, file = paste0(path_save, "/", model_name, ".rda"))
  }
  
  # RESULTS RETURNED
  # > return predictions
  if(res == "preds")
  {
    out <- preds
  }
  # > return prediction performance indicators
  if(res == "perf")
  {
    # > performance indicators
    RMSEP <- caret::RMSE(obs = preds$Ya_obs,      pred = preds$Ya_pred)
    NSE   <- hydroGOF::NSE(obs = preds$Ya_obs,    sim=preds$Ya_pred)
    #R2    <- caret::R2(obs = preds$Ya_obs,        pred = preds$Ya_pred)
    #Bias  <- Metrics::bias(actual = preds$Ya_obs, predicted = preds$Ya_pred)
    
    out <- data.frame("RMSEP"=RMSEP, "NSE"=NSE, #"R2"=R2, "Bias"=Bias, 
                      "N_predictors"=N_predictors)
  }
  
  return(out)
  
}


# ----------------------------------
# CROSS VALIDATION ON SITES
# res: "preds" or "perf"
# model_gpe: "rf or "lm",
# outcome: name of the predicted variable "Ya" or "Ya_ano"

function_cv_geo <- function(model, 
                            outcome,
                            model_gpe  = "rf",
                            recompute_scores = FALSE,
                            k_nb_folds=10,
                            data, 
                            data_day,
                            seed=101, 
                            model_name = "Give_me_a_name",
                            save = F,
                            path_save = NULL, 
                            res = "preds"){
  
  # Data to store predicted values
  data_pred <- list()
  
  # CROSS VALIDATION ON LOCALISATION
  # Splits sites in k_nb_folds groups
  data_sites <- data %>% 
    distinct(gridcode, x, y)
  
  set.seed(seed); folds.cv <- caret::createFolds(y = data_sites$gridcode, 
                                                 k = k_nb_folds, 
                                                 list = F)
  data_sites$fold_for_cv <- folds.cv
  
  # Merge with initial data 
  data_for_cv <- data %>% 
    left_join(., data_sites %>% dplyr::select(-gridcode), 
              by = c("x", "y"))
  
  # Remove the sites from the dataset
  for(i in 1:length(unique(data_for_cv$fold_for_cv)))
  {
    
    # --------------------------------------
    # Name of the model 
    model_init <- model_name
    
    # ---------------------
    # Define test and train datasets
    data_train_0 <- data_for_cv[which(data_for_cv$fold_for_cv != i),]
    data_test_0  <- data_for_cv[which(data_for_cv$fold_for_cv == i),]
    
    # --------------------------------------
    # If necessary, compute scores from train dataset and predict for test data set 
    if(recompute_scores == TRUE)
    {
      
      new_scores <- newscores_wraper(model_init = model_init, 
                                     init_data  = data_train_0, 
                                     pred_data  = data_test_0, 
                                     #data_day   = data_day, 
                                     outcome    = outcome)
      
      # Test and train datasets with new scores
      data_train <- new_scores$new_init_data
      data_test  <- new_scores$new_pred_data
    }
    if(recompute_scores == FALSE)
    {
      
      # Test and train datasets provided 
      data_train <- data_train_0
      data_test  <- data_test_0
    }
    
    # ---------------------
    # Fit model on the train dataset
    # random forest
    if(model_gpe == "rf")
    {
      # > fit on train
      set.seed(seed)
      mod_train <- ranger(as.formula(model),
                          data=data_train, 
                          num.tree=500,
                          importance="impurity")   
      
      # > predict on test
      pred.test <- predict(mod_train, 
                           data = data_test, 
                           type = "response", 
                           seed = seed, 
                           num.trees = 500)
      # > retrieve predictions & nb of predictors
      pred.test_vec <- pred.test$predictions
      N_predictors <- mod_train$num.independent.variables
      
    }
    # linear regression (benchmark)
    if(model_gpe == "lm")
    {
      # > fit on train
      mod_train <- lm(as.formula(model),
                      data=data_train)   
      
      # > predict on test
      pred.test <- predict(mod_train,
                           newdata = data_test, 
                           type = "response")
      # > retrieve predictions
      pred.test_vec <- as.numeric(as.character(pred.test))
      N_predictors <- length(mod_train$coefficients) - 1
    }
    
    
    # ---------------------
    # Store the predictions
    data_pred[[paste0(i)]] <- data.frame(Model     = model_name,
                                         site_year = data_test$site_year,
                                         Ya_obs    = data_test[,paste0(outcome)],
                                         Ya_pred   = pred.test_vec,
                                         N_predictors = N_predictors) %>% 
      dplyr::rename("Ya_obs" = 3)
    
  }
  
  # > predictions in a table
  preds <- plyr::ldply(data_pred, data.frame, .id = "fold_cv")
  
  # > save result
  if(save==TRUE)
  {
    save(preds, file = paste0(path_save, "/", model_name, ".rda"))
  }
  
  # RESULTS RETURNED
  # > return predictions
  if(res == "preds")
  {
    out <- preds
  }
  # > return prediction performance indicators
  if(res == "perf")
  {
    # > performance indicators
    RMSEP <- caret::RMSE(obs = preds$Ya_obs,      pred = preds$Ya_pred)
    NSE   <- hydroGOF::NSE(obs = preds$Ya_obs,    sim=preds$Ya_pred)
    #R2    <- caret::R2(obs = preds$Ya_obs,        pred = preds$Ya_pred)
    #Bias  <- Metrics::bias(actual = preds$Ya_obs, predicted = preds$Ya_pred)
    
    out <- data.frame("RMSEP"=RMSEP, "NSE"=NSE, #"R2"=R2, "Bias"=Bias, 
                      "N_predictors"=N_predictors)
  }
  
  return(out)
  
}
