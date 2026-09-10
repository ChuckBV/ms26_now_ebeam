# Examines mating success with dose as categorical, consistent with
# oviposition and pupal development analyses.

# Preliminary examination of the 0 Gy (stay at home) and 5 Gy (travel with 
# no eBeam) controls suggests a significant difference (lines 80 to 120)

# After examining both the models and the plots, it was decided that modeling
# dose as continuous is best for this response variable.


# Load libraries
library(tidyverse)
library(broom)          # Used with tidy and glance to get lm stats by group
library(DHARMa)
library(glmmTMB)

source("./scripts/themes.R")

##### Open and describe Experiment 2 mating ################# #

df_mating <- readRDS("./data/dat_expt2_mating.Rds")

df_mating %>% 
  group_by(run,Dose) %>% 
  summarise(nObs = n())
# # A tibble: 11 × 3
# # Groups:   run [2]
#  run  Dose  nObs
#  <dbl> <dbl> <int>
#  1     1     0    60
#  2     1     6    60
#  3     1   177    60
#  4     1   211    60
#  5     1   405    60
#  6     1   445    60
#  7     2     6    60
#  8     2   176    60
#  9     2   208    60
# 10     2   309    60
# 11     2   310    60

####### Make dose categorical ############################################# #

df_mating <- df_mating %>% 
  mutate(dose2 = case_when(
    Dose < 5 ~ 0,              # 0 only
    Dose < 100 ~ 5,            # 6 only
    Dose < 200 ~ 175,          # 176, 177
    Dose < 250 ~ 210,          # 208, 211
    Dose < 315 ~ 310,          # 310
    .default = 420             # 405, 445
  ))

# Make dose2 a factor
df_mating$dose2 <- factor(df_mating$dose2)


####### Get proportion mated by categorical dose2 ######################### #

df_pmated <- df_mating %>% 
  group_by(run,dose2,Sex,Mated) %>%
  summarise(nObs = n()) %>% 
  pivot_wider(names_from = Mated, values_from = nObs) %>% 
  rename(mated = Yes) %>% 
  mutate(trials = No + mated, prop = mated/trials) %>% 
  select(run,dose2,Sex,mated,trials,prop)

df_pmated
# # A tibble: 18 × 6
# # Groups:   run, dose2, Sex [18]
#   run dose2 Sex    mated trials  prop
#   <dbl> <fct> <fct>  <int>  <int> <dbl>
# 1     1 0     Male      23     30 0.767
# 2     1 0     Female    15     30 0.5  
# 3     1 5     Male      16     30 0.533
# 4     1 5     Female     8     30 0.267
# 5     1 175   Male      14     30 0.467

# Make run into a factor
df_pmated$run <- factor(df_pmated$run, levels = c(1,2))

######### Compare control doses 0 and 5 ################################### #

df_ctrl <- df_pmated %>% 
  filter(dose2 %in% c("0","5")) %>% 
  droplevels()

df_ctrl %>% 
  group_by(run,dose2) %>% 
  summarise(nObs = n(),
            mn = mean(prop, na.rm = T),
            sem = FSA::se(prop)
            )
# # A tibble: 3 × 5
# # Groups:   run [2]
#   run   dose2  nObs    mn    sem
#   <fct> <fct> <int> <dbl>  <dbl>
# 1 1     0         2 0.633 0.133 
# 2 1     5         2 0.4   0.133 
# 3 2     5         2 0.533 0.0667
    # Looks similar

# Examine with binomial GLM
m1 <- glm(cbind(mated, trials - mated) ~ dose2 + Sex + run,
          family = binomial,
          data = df_ctrl)

summary(m1)

# Test dispersion
dispersion <- sum(residuals(m1, type = "pearson")^2) / m1$df.residual
dispersion # 3.46

# Try quasibinomial
m_quasi  <- glm(cbind(mated, trials - mated) ~ dose2 + Sex + run,
                family = quasibinomial,
                data = df_ctrl)

dispersion_quasi <- sum(residuals(m_quasi, type = "pearson")^2) /
  m_quasi$df.residual

dispersion_quasi

summary(m_quasi)

# Try beta binomial

m1_bb <- glmmTMB(
  cbind(mated, trials - mated) ~ dose2 + Sex + run,
  family = betabinomial(),
  data = df_ctrl
)

AIC(m1_bb)

dispersion_m1bb <- sum(residuals(m1_bb, type = "pearson")^2) / m1_bb$df.residual

sim <- simulateResiduals(m1_bb, n = 1000)
plot(sim)

summary(m1_bb)
# Family: betabinomial  ( logit )
# Formula:          cbind(mated, trials - mated) ~ dose2 + Sex + run
# Data: df_ctrl
# 
# AIC       BIC    logLik -2*log(L)  df.resid 
# 39.5      38.4     -14.7      29.5         1 
# 
# 
# Dispersion parameter for betabinomial family ():  187 
# 
# Conditional model:
#             Estimate Std. Error z value Pr(>|z|)  
# (Intercept)   0.8354     0.3385   2.468   0.0136 *
# dose25       -0.9683     0.4064  -2.383   0.0172 *
# SexFemale    -0.5589     0.3301  -1.693   0.0904 .
# run2          0.5492     0.4014   1.368   0.1712  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


######### Plotted ################################### #

ggplot(df_pmated, aes(x = dose2, y = prop, color = run)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ylim(c(0,1)) +
  xlab("Dose (Gy)") +
  ylab("Proportion mated") +
  facet_grid(Sex ~ .) +
  theme_powerpoint()  

######### Model for complete data set ################################### #

m2_binomial <- glm(cbind(mated, trials - mated) ~ dose2*Sex + dose2*run + Sex*run,
                   family = binomial,
                   data = df_pmated)

dispersion <- sum(residuals(m2_binomial, type = "pearson")^2) / m2_binomial$df.residual
dispersion # [1] 4.798938

m2_bb <- glmmTMB(
  cbind(mated, trials - mated) ~ dose2 + Sex + run,
  family = betabinomial(),
  data = df_pmated
)

sim <- simulateResiduals(m2_bb, n = 1000)
plot(sim) # Failed miserably

m2_re <- glmmTMB(
  cbind(mated, trials - mated) ~ dose2*Sex + (1 | run),
  family = binomial,
  data = df_pmated
)

simulateResiduals(m2_re) %>% plot()

summary(m2_re)
# Family: binomial  ( logit )
# Formula:          cbind(mated, trials - mated) ~ dose2 + Sex + (1 | run)
# Data: df_pmated
# 
# AIC       BIC    logLik -2*log(L)  df.resid 
# 108.8     115.9     -46.4      92.8        10 
# 
# Random effects:
#   
#   Conditional model:
#   Groups Name        Variance  Std.Dev. 
# run    (Intercept) 8.775e-11 9.368e-06
# Number of obs: 18, groups:  run, 2
# 
# Conditional model:
#             Estimate Std. Error z value Pr(>|z|)   
# (Intercept)  0.73132    0.28141   2.599  0.00935 **
# dose25      -0.68551    0.32572  -2.105  0.03532 * 
# dose2175    -0.82133    0.32643  -2.516  0.01187 * 
# dose2210    -0.07166    0.32839  -0.218  0.82725   
# dose2310    -0.48366    0.32553  -1.486  0.13734   
# dose2420    -0.31459    0.32616  -0.965  0.33479   
# SexFemale   -0.36086    0.15821  -2.281  0.02256 * 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

car::Anova(m2_re)
# Analysis of Deviance Table (Type II Wald chisquare tests)
# 
# Response: cbind(mated, trials - mated)
#             Chisq Df Pr(>Chisq)  
# dose2     12.5147  5    0.02838 *
# Sex        5.0180  1    0.02508 *
# dose2:Sex  4.5003  5    0.47985  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


