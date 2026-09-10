# Uses a GLM to assesses the impact of run, sex, and dose in 
# a single Analysis of Deviance 

# Load libraries and scripts
library(tidyverse)
library(DHARMa)
library(car)
library(detectseparation)


###### Expt1 flight cylinder data set ################################# #

df_flight <- readRDS("./data/dat_expt1_flt_cylinder.Rds")

## Plot of percent escaped by sex and run is produced in script 7
## and displayed in ./doc/results_overview.pdf

df_flight
# # A tibble: 72 × 5
#   Run   sex       y  dose     n
#   <fct> <chr> <dbl> <dbl> <dbl>
# 1 1     Male      9     0    10
# 2 1     Male      4     0    10
# 3 1     Male      6     0    10

    # Note that run is a factor

# Examine GLM w binomial distribution
m <- glm(cbind(y, n - y) ~ dose*sex + dose*Run + sex*Run,
         family = binomial,
         data = df_flight)
summary(m)
# 
# Call:
#   glm(formula = cbind(y, n - y) ~ dose * sex + dose * Run + sex * 
#         Run, family = binomial, data = df_flight)
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -0.5173895  0.2179101  -2.374   0.0176 *  
# dose         -0.0004715  0.0009367  -0.503   0.6148    
# sexMale       1.3360689  0.2775512   4.814 1.48e-06 ***
# Run2          0.2682830  0.2838491   0.945   0.3446    
# dose:sexMale  0.0003283  0.0009953   0.330   0.7415    
# dose:Run2     0.0004882  0.0010205   0.478   0.6324    
# sexMale:Run2 -0.0582098  0.3222142  -0.181   0.8566    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 212.07  on 71  degrees of freedom
# Residual deviance: 130.20  on 65  degrees of freedom
# AIC: 315.24
# 
# Number of Fisher Scoring iterations: 4

# Look at model fit
disp <- sum(residuals(m, type="pearson")^2) / m$df.residual
disp 
#[1] 1.773443  Mild but acceptable

sim <- simulateResiduals(m)
plot(sim)
# QQ plot w significant deviation, resids vs. preds shows deviation

# Based on significant lack of fit for binomial model, try GLM w quasibonomial
# m_quasi <- glm(cbind(y, n - y) ~ dose*sex + dose*Run + sex*Run, 
#                family = quasibinomial,
#                data = df_flight)
# 
# summary(m_quasi)
# 
# Call:
#   glm(formula = cbind(y, n - y) ~ dose * sex + dose * Run + sex * 
#         Run, family = quasibinomial, data = df_flight)
# 
# Coefficients:
#                Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -0.5173895  0.2902010  -1.783 0.079277 .  
# dose         -0.0004715  0.0012475  -0.378 0.706719    
# sexMale       1.3360689  0.3696278   3.615 0.000587 ***
#   Run2          0.2682830  0.3780149   0.710 0.480419    
# dose:sexMale  0.0003283  0.0013255   0.248 0.805185    
# dose:Run2     0.0004882  0.0013590   0.359 0.720597    
# sexMale:Run2 -0.0582098  0.4291075  -0.136 0.892515    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for quasibinomial family taken to be 1.773548)
# 
# Null deviance: 212.07  on 71  degrees of freedom
# Residual deviance: 130.20  on 65  degrees of freedom
# AIC: NA
# 
# Number of Fisher Scoring iterations: 4

plot(m_quasi$fitted.values,m_quasi$residuals) # looks good

cooks <- cooks.distance(m_quasi)
cooks[cooks > 0.5]
# named numeric(0)

vif(m_quasi)
#     dose      sex      Run dose:sex dose:Run  sex:Run 
# 3.642205 3.018268 3.178758 3.127102 4.337054 2.753274 

# Use car::Anova() to get a type II Analysis of Deviance table



car::Anova(m_quasi,  type = 2)
# Analysis of Deviance Table (Type II tests)
# 
# Response: cbind(y, n - y)
#          LR Chisq Df Pr(>Chisq)    
# dose        0.001  1     0.9755    
# Run         2.379  1     0.1230    
# sex        43.830  1  3.583e-11 ***
# dose:Run    0.129  1     0.7193    
# dose:sex    0.061  1     0.8043    
# Run:sex     0.018  1     0.8921    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Examine plot 
  
