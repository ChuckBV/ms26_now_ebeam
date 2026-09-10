# Provides GLM with binomial and plot with linear smoothing.
# The GLM fit stats are good. With all interactions and car::Anova() type III,
# Sex is nominally significant (0.1 > P > 0.05), nothing else is close

# Load libraries
library(tidyverse)
library(broom)          # Used with tidy and glance to get lm stats by group

source("./scripts/themes.R")

##### Open and describe Experiment 2 mating ################# #

df_mating <- readRDS("./data/dat_expt2_mating.Rds")
df_mating
# # A tibble: 660 × 4
#     run  Dose Sex   Mated
#   <dbl> <dbl> <fct> <fct>
# 1     1     6 Male  No   
# 2     1     6 Male  Yes  
# 3     1     6 Male  No   
# 4     1     6 Male  Yes  
# 5     1     6 Male  No   


# Create contingency table
mating_table <- table(df_mating$Sex, df_mating$Mated)
print(mating_table)
# 
#        Yes  No
# Male   189 141
# Female 160 170

prop.table(mating_table, margin = 1)
# 
#              Yes        No
# Male   0.5727273 0.4272727
# Female 0.4848485 0.5151515

# Perform Chi-squared test
chi_test <- chisq.test(mating_table)
print(chi_test)
# 
# Pearson's Chi-squared test with Yates' continuity correction
# 
# data:  mating_table
# X-squared = 4.7673, df = 1, p-value = 0.029


## So 57% of males were mated, but only 48% of females
## For this sample size, this is a significant departure from a 1:1 sex ratio

###### Recode to y, n, and prop, then plot prop by sex, dose, and run ##### #

# Get successes, trials, and proportion success for rep*Dose*sex

df_pmated <- df_mating %>% 
  group_by(run,Dose,Sex,Mated) %>%
  summarise(nObs = n()) %>% 
  pivot_wider(names_from = Mated, values_from = nObs) %>% 
  rename(mated = Yes) %>% 
  mutate(trials = No + mated, prop = mated/trials) %>% 
  select(run,Dose,Sex,mated,trials,prop)

df_pmated
# # A tibble: 22 × 6
# # Groups:   rep, Dose, Sex [22]
#     rep  Dose Sex    mated trials  prop
#   <dbl> <dbl> <chr>  <int>  <int> <dbl>
# 1     2     0 Female    15     30 0.5  
# 2     2     0 Male      23     30 0.767
# 3     2     6 Female     8     30 0.267
# 4     2     6 Male      16     30 0.533
# 5     2   177 Female    15     30 0.5  

# Make run into a factor
df_pmated$run <- factor(df_pmated$run, levels = c(1,2))


# Plot proportion mated by vs dose by sex, with run color-coded

# p1 <- 
ggplot(df_pmated, aes(x = Dose, y = prop, color = run)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ylim(c(0,1)) +
  xlab("Dose (Gy)") +
  ylab("Proportion mated") +
  facet_grid(Sex ~ .) +
  theme_powerpoint()  

ggsave("./results/fig3_mating_by_sex_run.jpg",
       plot = p1,
       device = "jpg",
       width = 5.83,
       height = 5.83,
       units = "in",
       dpi = 300)


########## GLM with binomial distribution for run*sex ##################### #

data$rep <- data$rep - 1  # Change rep numbers from 2,3 to 1,2

data$rep = factor(data$rep, levels = c(1,2))

m <- glm(Mated ~ Dose*Sex + Dose*run + Sex*run,
         data = df_mating,
         family = binomial
         )

summary(m) 
# 
# Call:
#   glm(formula = Mated ~ Dose * Sex + Dose * run + Sex * run, family = binomial, 
#       data = df_mating)
# 
# Coefficients:
#                  Estimate Std. Error z value Pr(>|z|)  
# (Intercept)    -0.6617891  0.4923216  -1.344   0.1789  
# Dose           -0.0006611  0.0017052  -0.388   0.6982  
# SexFemale       0.9243058  0.5355307   1.726   0.0844 .
# run             0.2518892  0.3319529   0.759   0.4480  
# Dose:SexFemale -0.0005528  0.0010613  -0.521   0.6024  
# Dose:run        0.0004575  0.0012069   0.379   0.7047  
# SexFemale:run  -0.3129431  0.3149679  -0.994   0.3204  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 912.77  on 659  degrees of freedom
# Residual deviance: 904.34  on 653  degrees of freedom
# AIC: 918.34
# 
# Number of Fisher Scoring iterations: 4

# Model diagnostics
plot(m, which = 1)

disp <- sum(residuals(m, type = "pearson")^2) / m$df.null
disp # 1.0015, very good fit 
    # (1.5 would be likely overdispersion, 2 would be clear overdispersion)
    # passes diagnostics
    
car::Anova(m, type = 3)
# Analysis of Deviance Table (Type III tests)
# 
# Response: Mated
# LR Chisq Df Pr(>Chisq)  
# Dose      0.15034  1    0.69821  
# Sex       2.99563  1    0.08349 .
# run       0.57568  1    0.44801  
# Dose:Sex  0.27141  1    0.60239  
# Dose:run  0.14368  1    0.70465  
# Sex:run   0.98764  1    0.32032  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#### Examine linear regression stats for run*sex ########################## #

df_results <- df_pmated %>%
  group_by(Sex, Run) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(prop ~ Dose, data = .x)),
    tidy   = map(model, tidy),      # coefficients
    glance = map(model, glance)     # R², AIC, etc.
  )


df_results %>%
  select(Sex, Run, tidy) %>% 
  unnest(tidy)
# # A tibble: 8 × 7
# # Groups:   Sex, Run [4]
#   Sex    Run   term         estimate std.error statistic p.value
#   <fct>  <fct> <chr>           <dbl>     <dbl>     <dbl>   <dbl>
# 1 Male   1     (Intercept)  0.638     0.126        5.07  0.00716
# 2 Male   1     Dose        -0.000128  0.000466    -0.274 0.798  
# 3 Female 1     (Intercept)  0.413     0.0734       5.63  0.00489
# 4 Female 1     Dose         0.000365  0.000272     1.34  0.250  
# 5 Male   2     (Intercept)  0.437     0.0865       5.05  0.0150 
# 6 Male   2     Dose         0.000446  0.000375     1.19  0.320  
# 7 Female 2     (Intercept)  0.568     0.202        2.81  0.0671 
# 8 Female 2     Dose        -0.000435  0.000875    -0.497 0.653 

x <- df_results %>%
  select(Sex, Run, glance) %>%
  unnest(glance)
x
# # A tibble: 4 × 14
# # Groups:   Sex, Run [4]
#   Sex    Run   r.squared adj.r.squared  sigma statistic p.value    df logLik   AIC    BIC deviance
#   <fct>  <fct>     <dbl>         <dbl>  <dbl>     <dbl>   <dbl> <dbl>  <dbl> <dbl>  <dbl>    <dbl>
# 1 Male   1        0.0184       -0.227  0.198     0.0751   0.798     1   2.43  1.14  0.517   0.156 
# 2 Female 1        0.311         0.139  0.115     1.80     0.250     1   5.67 -5.34 -5.96    0.0531
# 3 Male   2        0.320         0.0936 0.0936    1.41     0.320     1   6.03 -6.05 -7.22    0.0263
# 4 Female 2        0.0760       -0.232  0.218     0.247    0.653     1   1.79  2.41  1.24    0.143 
# # ℹ 2 more variables: df.residual <int>, nobs <int>
