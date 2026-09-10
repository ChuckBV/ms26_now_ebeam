# Objective--determine if there is a difference between the stay-at-home 
# controls and the sham controls that traveled to the eBeam facility but
# were not irradiated.

# Conclusion: sex significant (P < 0.01), run marginally sig.(P < 0.1),
# and sex ns. However, one of the stay-at-home controls was twice as high
# as anything else in the data set and was therefore likely and influential
# observation for linear regression.

library(tidyverse)
library(DHARMa)

####### Re-create data frames used in other analyses ##################### #

df_surv <- readRDS("./data/dat_expt1_longevity.Rds")

df_medians <- df_surv %>% 
  group_by(run, sex, cage, dose) %>% 
  summarise(nObs = n(),
            med = median(survival_time),
            n_surv = sum(status == "Censored"))

df_medians <- df_medians %>% 
  mutate(pct_14dsurv = 100*n_surv/nObs)

df_ctrl <- df_medians %>% 
  filter(dose < 10) %>% 
  mutate(control_type = if_else(dose <= 0,"stay","travel"))

####### Crude mean & SE for simple comparison and multi-factorial ########## #

# Percent survival to 14 days, simple comparison
df_ctrl %>%
  group_by(control_type) %>% 
  summarise(nObs = n(),
            mn = mean(pct_14dsurv),
            sem = FSA::se(pct_14dsurv))
# # A tibble: 2 × 4
#   control_type  nObs    mn   sem
#   <chr>        <int> <dbl> <dbl>
# 1 stay            12  29.4  6.69
# 2 travel          12  26.2  4.14

# Percent survival to 14 days, multi-factorial
dose_avg <- df_ctrl %>%
  group_by(run,sex,control_type) %>% 
  summarise(nObs = n(),
            mn = mean(pct_14dsurv),
            sem = FSA::se(pct_14dsurv))
dose_avg
# # A tibble: 8 × 6
# # Groups:   run, sex [4]
#     run sex   control_type  nObs    mn   sem
#   <dbl> <chr> <chr>        <int> <dbl> <dbl>
# 1     1 f     stay             3 31.1   5.88
# 2     1 f     travel           3 31.9   5.19
# 3     1 m     stay             3  6.67  6.67
# 4     1 m     travel           3 19.5   4.67
# 5     2 f     stay             3 53.3  16.7 
# 6     2 f     travel           3 30    10   
# 7     2 m     stay             3 26.7   8.82
# 8     2 m     travel           3 23.3  13.3 

####### GLM           ##################################################### #

m1 <- glm(cbind(n_surv, nObs - n_surv) ~ dose + sex + run,
          family = binomial,
          data = df_ctrl)

# Dispersion test
dispersion <- sum(residuals(m1, type = "pearson")^2) / m1$df.residual
dispersion # [1] 1.407268

# Residuals plot
sim <- simulateResiduals(m1, n = 1000)
plot(sim) # Good enough

summary(m1)
# 
# Call:
#   glm(formula = cbind(n_surv, nObs - n_surv) ~ dose + sex + run, 
#       family = binomial, data = df_ctrl)
# 
# Coefficients:
# Estimate Std. Error z value Pr(>|z|)   
# (Intercept) -1.34149    0.58363  -2.299  0.02153 * 
# dose        -0.01277    0.06003  -0.213  0.83151   
# sexm        -0.89437    0.30356  -2.946  0.00322 **
# run          0.53725    0.32281   1.664  0.09606 . 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 43.792  on 23  degrees of freedom
# Residual deviance: 31.200  on 20  degrees of freedom
# AIC: 93.213
# 
# Number of Fisher Scoring iterations: 4
