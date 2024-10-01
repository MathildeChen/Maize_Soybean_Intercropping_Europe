# -------------------------------------------------------------------------
# 
#       00-2. Computation of maize and soybean pLER 
#       based on data from Li et al., 2022
#       Author: M. Chen, Inrae, 2024 
# 
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Loading packages 
library(tidyverse)
library(nlme)
library(broom.mixed)

# -------------------------------------------------------------------------
# Data 
tab_maize_soybean <- read.csv2("E:/POSTDOC INRAE/BIBLIO/INTERCROPPING/Li et al. 2023 - productivity of intercropping/Dataset.csv") %>% 
  filter(Scomb == "Maize/soybean") %>% 
  # > create column for Yield of maize and yield of soybean 
  mutate(Y_m = if_else(Latin.name_species.A == "Zea mays", Y1, Y2),
         Y_s = if_else(Latin.name_species.A == "Glycine max", Y1, Y2),
         M_m = if_else(Latin.name_species.A == "Zea mays", M1, M2),
         M_s = if_else(Latin.name_species.A == "Glycine max", M1, M2)) %>%
  # > enter yields as numeric 
  mutate(Y_m = as.numeric(as.character(Y_m)),
         Y_s = as.numeric(as.character(Y_s)),
         M_m = as.numeric(as.character(M_m)),
         M_s = as.numeric(as.character(M_s))) %>% 
  # > compute pLER for each crop
  mutate(pLER_m = Y_m/M_m,
         pLER_s = Y_s/M_s) %>% 
  # > compute LER
  mutate(LER = pLER_m+pLER_s) %>%
  # > compute log(pLER)
  mutate(log_pLER_m = log(pLER_m),
         log_pLER_s = log(pLER_s),
         log_LER    = log(LER))

# -------------------------------------------------------------------------
# Model to estimate the mean pLER and 95% confidence interval of pLER 

# > LER
Fit_M1<-lme(LER~1,random=~1|Study/Expt,data=tab_maize_soybean,method="ML")#choose Study/Expt as random effect
Fit_M1<-lme(LER~1,random=~1|Study/Expt,data=tab_maize_soybean)#refit the model
summary(Fit_M1)
tidy(Fit_M1,effects="fixed", conf.int=TRUE)[, c(3,8,9)]
# estimate conf.low conf.high
#     1.24     1.13      1.36

# > pLER maize
Fit_M2<-lme(pLER_m~1,random=~1|Study/Expt,data=tab_maize_soybean,method="ML")#choose Study/Expt as random effect
Fit_M2<-lme(pLER_m~1,random=~1|Study/Expt,data=tab_maize_soybean)#refit the model
tidy(Fit_M2,effects="fixed", conf.int=TRUE)[, c(3,8,9)]

#   estimate conf.low conf.high
# 1    0.784    0.709     0.859

# > pLER maize
Fit_M3<-lme(pLER_s~1,random=~1|Study/Expt,data=tab_maize_soybean,method="ML")#choose Study/Expt as random effect
Fit_M3<-lme(pLER_s~1,random=~1|Study/Expt,data=tab_maize_soybean)#refit the model
tidy(Fit_M3,effects="fixed", conf.int=TRUE)[, c(3,8,9)]

#   estimate conf.low conf.high
# 1    0.458    0.370     0.547


