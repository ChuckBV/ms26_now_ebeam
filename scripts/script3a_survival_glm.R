#     Uses proportion surviving to 14 days for GLM w quasibonomial
#             This was used for the manuscript

# Load libraries
library(tidyverse)
library(DHARMa)
library(car)
library(emmeans)

source("./scripts/themes.R")

##### Open and describe Experiment 1 flight longevity data ################# #

df_surv <- readRDS("./data/dat_expt1_longevity.Rds")
df_surv 
# # A tibble: 716 × 7
#     run sex    dose  cage individual_moth survival_time status   
#   <dbl> <chr> <dbl> <dbl>           <dbl>         <dbl> <chr>    
# 1     1 m         0     1               1             1 Completed
# 2     1 m         0     1               2             5 Completed
# 3     1 m         0     1               3             5 Completed
# 4     1 m         0     1               4             7 Completed
# 5     1 m         0     1               5             8 Completed

##### Document the censored rate ################# #


# Histogram for raw distribution
hist(df_surv$survival_time)
# 14 is by far the most common frequency

# Group by cage
df_surv_cage <- df_surv %>% 
  group_by(run,sex,dose,cage) %>% 
  summarise(med = median(survival_time))
df_surv_cage

# Histogram for distribution of aggregated data
hist(df_surv_cage$med)
# more equal, but 14 still the most frequent obs

######### Examine median survival time and plot median vs dose ####### #

## Examine median survival time by dose
df_medians <- df_surv %>% 
  group_by(run, sex, cage, dose) %>% 
  summarise(nObs = n(),
            med = median(survival_time),
            n_surv = sum(status == "Censored"))

df_medians <- df_medians %>% 
  mutate(pct_14dsurv = 100*n_surv/nObs)

dose_avg <- df_medians %>%
  group_by(dose) %>% 
  summarise(nObs = n(),
            mn = mean(pct_14dsurv),
            sem = FSA::se(pct_14dsurv))

dose_avg
# # A tibble: 12 × 4
#  dose  nObs    mn   sem
#  <dbl> <int> <dbl> <dbl>
#  1    -1     6  40   10.3 
#  2     0     6  18.9  6.76
#  3     3     6  26.7  7.60
#  4     6     6  25.7  4.17
#  5   133     6  20    3.65
#  6   171     6  28.3  9.10
#  7   175     6  26.7  5.59
#  8   206     6  45    8.06
#  9   302     6  25   10.6 
# 10   324     6  38.3 11.1 
# 11   376     6  28.6  7.93
# 12   521     6  15    2.24

hist(df_medians$med)

# Proportion censored by Run and sex
df_surv %>%
  group_by(run,sex,status) %>%
  summarise(nObs = n()) %>%
  pivot_wider(names_from = status, values_from = nObs) %>%
  mutate(pct_censored = 100*Censored/(Censored + Completed))
# # A tibble: 4 × 5
# # Groups:   run, sex [4]
#     run sex   Censored Completed pct_censored
#   <dbl> <chr>    <int>     <int>        <dbl>
# 1     1 f           54       123         30.5
# 2     1 m           32       147         17.9
# 3     2 f           64       116         35.6
# 4     2 m           52       128         28.9  
    # Numerically more 14-d survivors among females than males, and
    # in run 2 compared to run 1

# # Ask how many survived 14 days
# df_surv %>% 
#   group_by(dose,status) %>% 
#   summarise(nObs = n()) %>% 
#   pivot_wider(names_from = status, values_from = nObs) %>% 
#   mutate(pct_surv = 100*Censored/(Censored + Completed))
#     # Mortality after 14 days is high but variable

df_survmeans <- df_medians %>% 
  group_by(sex,dose) %>% 
  summarise(nObs = n(),
            mn = mean(pct_14dsurv),
            sem = FSA::se(pct_14dsurv))
  
  

#### Crude ggplot of 14-day survival time ########################## #

# ggplot(df_medians,aes(x = dose, y = pct_14dsurv, color = factor(run))) +
#   geom_point(aes(shape = sex)) +
#   geom_smooth(method = "lm", aes(group = 1, color = "black"))

ggplot(df_survmeans,aes(x = dose, y = mn)) +
  geom_point() +
  geom_errorbar(aes(ymin = mn - sem, ymax = mn + sem)) +
  geom_smooth(method = "lm") +
  facet_grid(sex ~ .)


#### Fit df_medians to a GLM ########################## #

m1 <- glm(cbind(n_surv, nObs - n_surv) ~ dose*sex + dose*run + run*sex,
          family = binomial,
          data = df_medians)

dispersion <- sum(residuals(m1, type = "pearson")^2) / m1$df.residual # 1.65

sim_res <- simulateResiduals(fittedModel = m1, n = 1000)
plot(sim_res)  # passes

summary(m1)
# 
# Call:
#   glm(formula = cbind(n_surv, nObs - n_surv) ~ dose * sex + dose * 
#         run + run * sex, family = binomial, data = df_medians)
# 
# Coefficients:
#              Estimate Std. Error z value Pr(>|z|)   
# (Intercept) -1.437961   0.487022  -2.953  0.00315 **
# dose         0.001949   0.001945   1.002  0.31626   
# sexm        -1.541306   0.592856  -2.600  0.00933 **
# run          0.661667   0.297174   2.227  0.02598 * 
# dose:sexm    0.002658   0.001059   2.511  0.01205 * 
# dose:run    -0.002230   0.001117  -1.997  0.04580 * 
# sexm:run     0.362533   0.344715   1.052  0.29294   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 140.90  on 71  degrees of freedom
# Residual deviance: 114.68  on 65  degrees of freedom
# AIC: 292.04
# 
# Number of Fisher Scoring iterations: 4

m1_quasi <- glm(cbind(n_surv, nObs - n_surv) ~ dose*sex + dose*run + sex*run,
                family = quasibinomial,
                data = df_medians)

dispersion <- sum(residuals(m1_quasi, type = "pearson")^2) / m1_quasi$df.residual 
dispersion# 1.65

par(mfrow = c(2,2))
plot(m1_quasi)

summary(m1_quasi)
# 
# Call:
#   glm(formula = cbind(n_surv, nObs - n_surv) ~ dose * sex + dose * 
#         run + sex * run, family = quasibinomial, data = df_medians)
# 
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)  
# (Intercept) -1.437961   0.625008  -2.301   0.0246 *
# dose         0.001949   0.002496   0.781   0.4377  
# sexm        -1.541306   0.760827  -2.026   0.0469 *
# run          0.661667   0.381371   1.735   0.0875 .
# dose:sexm    0.002658   0.001359   1.956   0.0547 .
# dose:run    -0.002230   0.001433  -1.556   0.1245  
# sexm:run     0.362533   0.442381   0.820   0.4155  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for quasibinomial family taken to be 1.646925)
# 
# Null deviance: 140.90  on 71  degrees of freedom
# Residual deviance: 114.68  on 65  degrees of freedom
# AIC: NA
# 
# Number of Fisher Scoring iterations: 4

car::Anova(m1_quasi)
# Analysis of Deviance Table (Type II tests)
# 
# Response: cbind(n_surv, nObs - n_surv)
# LR Chisq Df Pr(>Chisq)  
# dose       0.6497  1    0.42021  
# sex        5.0564  1    0.02454 *
# run        3.7611  1    0.05246 .
# dose:sex   3.8692  1    0.04918 *
# dose:run   2.4391  1    0.11834  
# sex:run    0.6742  1    0.41158  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

emmeans(m1_quasi, ~ sex)

# Longevity by sex
df_medians %>% 
  group_by(sex) %>% 
  summarise(nObs = n(),
            mn = mean(pct_14dsurv),
            sem = FSA::se(pct_14dsurv))
# # A tibble: 2 × 4
#   sex    nObs    mn   sem
#   <chr> <int> <dbl> <dbl>
# 1 f        36  33.0  3.19
# 2 m        36  23.3  3.13

#### Calculate LM by sex ########################## #

df_medians_male <- df_medians %>% 
  filter(sex == "m")

df_medians_female <- df_medians %>% 
  filter(sex == "f")

# Linear regression, females

m_lm_fems <- lm(pct_14dsurv ~ dose, data = df_medians_female)
summary(m_lm_fems)
# 
# Call:
#   lm(formula = pct_14dsurv ~ dose, data = df_medians_female)
# 
# Residuals:
#    Min     1Q Median     3Q    Max 
# -28.33 -14.26  -3.35  11.55  47.69 
# 
# Coefficients:
#             Estimate Std. Error t value       Pr(>|t|)    
# (Intercept) 39.24830    4.68065   8.385 0.000000000864 ***
# dose        -0.03370    0.01899  -1.775         0.0849 .  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 18.6 on 34 degrees of freedom
# Multiple R-squared:  0.08478,	Adjusted R-squared:  0.05786 
# F-statistic:  3.15 on 1 and 34 DF,  p-value: 0.0849

# Convert 1 1-sided--extract t-value
t <- summary(m_lm_fems)$coefficients["dose", "t value"] # [1] -1.77

# Extract 1-sided P value
p_one_sided <- pt(t, df = df.residual(m_lm_fems)) # [1] 0.042

# Look for influential observations
car::influencePlot(m_lm_fems)

leverage <- hatvalues(m_lm_fems)
cooks_d <- cooks.distance(m_lm_fems)
dffits_vals <- dffits(m_lm_fems)
studentized_residuals <- rstudent(m_lm_fems)

summary(leverage)
summary(cooks_d)
summary(dffits_vals)
summary(studentized_residuals)

### None of the diagnostics suggest undue influence

# Linear regression, males

m_lm_males <- lm(pct_14dsurv ~ dose, data = df_medians_male)
# summary(m_lm_males)
# 
# Call:
#   lm(formula = pct_14dsurv ~ dose, data = df_medians_male)
# 
# Residuals:
#     Min      1Q  Median      3Q     Max 
# -25.618 -10.490  -3.070   6.249  44.382 
# 
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 19.74533    4.72757   4.177 0.000195 ***
# dose         0.01944    0.01918   1.014 0.317842    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 18.79 on 34 degrees of freedom
# Multiple R-squared:  0.02934,	Adjusted R-squared:  0.000793 
# F-statistic: 1.028 on 1 and 34 DF,  p-value: 0.3178

# Convert 1 1-sided--extract t-value
t <- summary(m_lm_males)$coefficients["dose", "t value"] # [1] 1.01

# Extract 1-sided P value
p_one_sided <- pt(t, df = df.residual(m_lm_males)) # [1] 0.84

# Look for influential observations
car::influencePlot(m_lm_males)

leverage <- hatvalues(m_lm_males)
cooks_d <- cooks.distance(m_lm_males)
dffits_vals <- dffits(m_lm_males)
studentized_residuals <- rstudent(m_lm_males)

summary(leverage)
summary(cooks_d)
summary(dffits_vals)
summary(studentized_residuals)

### None of the diagnostics suggest undue influence

