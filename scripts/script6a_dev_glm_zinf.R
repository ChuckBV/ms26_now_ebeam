# The GLM_zinf model is as used

library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(car)
library(emmeans)
library(multcomp)

source("./scripts/themes.R")

##### Load data set ####################################################### #

df_dev <- readRDS("./data/dat_expt2_pupae.Rds")

# Create dose categories as in script5
df_dev <- df_dev %>% 
  mutate(dose2 = case_when(
    gy < 5 ~ 0,              # 0 only
    gy < 100 ~ 5,            # 6 only
    gy < 200 ~ 175,          # 176, 177
    gy < 250 ~ 210,          # 208, 211
    gy < 315 ~ 310,          # 310
    .default = 420             # 405, 445
  ))

df_dev$dose2 <- factor(df_dev$dose2)

df_dev
# A tibble: 349 × 5
#   Run      gy sex    pupae_larvae dose2
#   <fct> <dbl> <chr>         <dbl> <fct>
# 1 1         0 Female            1 0    
# 2 1         0 Female           57 0    
# 3 1         0 Female            2 0    
# 4 1         0 Female           33 0    
# 5 1         0 Female           15 0 


############### Use glmmTMB to compare neg. binomial and zero inflation ### #

m_nb <- glmmTMB(
  pupae_larvae ~ dose2 * sex + dose2 * Run + sex * Run,
  family = nbinom2,
  data = df_dev
)
# dropping columns from rank-deficient conditional model: dose2210:Run2, dose2310:Run2, dose2420:Run2
# Warning messages:
# 1: In finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old) :
#   Model convergence problem; non-positive-definite Hessian matrix. See vignette('troubleshooting')
# 2: In finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old) :
#   Model convergence problem; singular convergence (7). See vignette('troubleshooting'), help('diagnose')

m_zi <- glmmTMB(
  pupae_larvae ~ dose2 * sex + dose2 * Run + sex * Run,
  ziformula = ~ 1,
  family = nbinom2,
  data = df_dev
)
# dropping columns from rank-deficient conditional model: dose2210:Run2

# AIC comparison
# AIC(m_nb, m_zi)
#      df      AIC
# m_nb 17       NA
# m_zi 18 729.1308

anova(m_nb, m_zi)
# Data: df_dev2
# Models:
# m_nb: pupae_larvae ~ dose2 * sex + dose2 * Run + sex * Run, zi=~0, disp=~1
# m_zi: pupae_larvae ~ dose2 * sex + dose2 * Run + sex * Run, zi=~1, disp=~1
#      Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)
# m_nb 13 720.20 764.32 -347.10   694.20                        
# m_zi 14 721.13 768.64 -346.57   693.13 1.073      1     0.3003
    # No significant improvement with zero-inflation

# Residuals diagnostics
simulateResiduals(m_nb) %>% plot()
simulateResiduals(m_zi) %>% plot()
    # surprisingly, nb passes resids test. But AIC and Hessian rule it out

summary(m_zi)
# Family: nbinom2  ( log )
# Formula:          pupae_larvae ~ dose2 * sex + dose2 * Run + sex * Run
# Zero inflation:                ~1
# Data: df_dev
# 
# AIC       BIC    logLik -2*log(L)  df.resid 
# 729.1     798.5    -346.6     693.1       331 
# 
# 
# Dispersion parameter for nbinom2 family (): 0.36 
# 
# Conditional model:
# Estimate   Std. Error z value        Pr(>|z|)    
# (Intercept)            3.0037       0.4577   6.563 0.0000000000528 ***
# dose25                -4.2781       1.0591  -4.040 0.0000535532002 ***
# dose2175              -7.2265       1.6462  -4.390 0.0000113465452 ***
# dose2210             -42.2758 1320692.6725   0.000        0.999974    
# dose2310             -34.5309   24083.8705  -0.001        0.998856    
# dose2420             -25.6204   16529.2079  -0.002        0.998763    
# sexMale                0.5017       0.6083   0.825        0.409514    
# Run2                   8.3960       1.5602   5.382 0.0000000738441 ***
# dose25:sexMale         3.7044       1.2178   3.042        0.002351 ** 
# dose2175:sexMale       3.9469       1.6360   2.412        0.015844 *  
# dose2210:sexMale      37.0440 1320692.6725   0.000        0.999978    
# dose2310:sexMale       1.2402  140958.7532   0.000        0.999993    
# dose2420:sexMale      -0.7906   25044.2613   0.000        0.999975    
# dose25:Run2           -3.3594       1.2242  -2.744        0.006065 ** 
# dose2175:Run2         -4.6265       1.4021  -3.300        0.000968 ***
# dose2210:Run2              NA           NA      NA              NA    
# dose2310:Run2              NA           NA      NA              NA    
# dose2420:Run2              NA           NA      NA              NA    
# sexMale:Run2          -5.5715       1.2324  -4.521 0.0000061598867 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Zero-inflation model:
#             Estimate Std. Error z value Pr(>|z|)
# (Intercept)  -1.0290     0.8015  -1.284    0.199

# Likelihood ratio omnibus test

# Test the following using the drop1() function
# m_no_dose_sex <- update(m_full, . ~ . - dose2:sex)
# m_no_dose_run <- update(m_full, . ~ . - dose2:Run)
# m_no_sex_run  <- update(m_full, . ~ . - sex:Run)
# 
# anova(m_no_dose_sex, m_full)   # tests dose2:sex
# anova(m_no_dose_run, m_full)   # tests dose2:Run
# anova(m_no_sex_run,  m_full)   # tests sex:Run

drop1(m_zi, test = "Chisq")
# Single term deletions
# 
# Model:
# pupae_larvae ~ dose2 * sex + dose2 * Run + sex * Run
# Single term deletions
# 
# Model:
# pupae_larvae ~ dose2 * sex + dose2 * Run + sex * Run
#           Df    AIC    LRT   Pr(>Chi)    
# <none>       729.13                      
# dose2:sex  5 747.01 27.880 0.00003842 ***
# dose2:Run  2 737.19 12.059   0.002407 ** 
# sex:Run    1 745.71 18.580 0.00001629 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  
car::Anova(m_zi)
# Analysis of Deviance Table (Type II Wald chisquare tests)
#   
# Response: pupae_larvae
# Chisq Df Pr(>Chisq)    
# dose2     81.6395  5  3.808e-16 ***
# sex        0.8041  1   0.369879    
# Run        5.8671  1   0.015426 *  
# dose2:sex  9.7103  5   0.083873 .  
# dose2:Run 11.9141  2   0.002588 ** 
# sex:Run   20.4379  1  6.160e-06 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

  