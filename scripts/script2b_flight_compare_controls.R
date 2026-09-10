# Check whether performance of sham controls (went to Fremont eBeam facility)
# and stay-at home controls were different

# Outcome--similar means, see lines 42-43. Treatment (0 vs 6 Gy) and run ns,
# but sex < 0.001.

library(tidyverse)
library(DHARMa)


###### Get Flight data and identify controls ############################## #

df_flight <- readRDS("./data/dat_expt1_flt_cylinder_cleaned.Rds")

# Get controls
df_cntrl <- df_flight %>% 
  filter(dose < 10)

df_cntrl %>% 
  group_by(Run,dose) %>% 
  summarise(nObs = n())
# # A tibble: 4 × 3
# # Groups:   Run [2]
#   Run    dose  nObs
#   <fct> <dbl> <int>
# 1 1         0     6
# 2 1         6     6
# 3 2         0     6
# 4 2         6     6

# Get character labels on controls (factor)
df_cntrl <- df_cntrl %>% 
  mutate(ctrl_type = (if_else(dose == 6,"travel","stay")))

# Display summary comparison of controls
df_cntrl %>% 
  group_by(ctrl_type) %>% 
  summarise(nObs = sum(!is.na(y)),
            mn = mean(y, na.rm = T),
            sem = FSA::se(y)
            )
# # A tibble: 2 × 4
#   ctrl_type  nObs    mn   sem
#   <chr>     <int> <dbl> <dbl>
# 1 stay         12  5.67 0.772
# 2 travel       12  5.33 0.689
  
# Determined y rather than y/n for convenience. Looks close. Confirm w GLM

# Test binomial distribution
m1 <- glm(cbind(y,n - y) ~ ctrl_type,
               family = binomial(),
               data = df_cntrl
)

# Dispersion test
dispersion <- sum(residuals(m1, type = "pearson")^2) / m1$df.residual
dispersion # 2.6--too high

# Residuals plot
sim <- simulateResiduals(m1, n = 1000)
plot(sim)  # Significant deviation in QQ plo

# Add additional predictors
m_binomial <- glm(cbind(y,n - y) ~ ctrl_type + Run + sex,
          family = binomial,
          data = df_cntrl
         )

# Dispersion test
dispersion <- sum(residuals(m_binomial, type = "pearson")^2) / m_binomial$df.residual
dispersion # 1.75 w run and sex included

# Plot residuals
sim <- simulateResiduals(m_binomial, n = 1000)
plot(sim)  # No significant deviation, no resid problems detected

# Examine model summary
summary(m_binomial)
# 
# Call:
#   glm(formula = cbind(y, n - y) ~ ctrl_type + Run + sex, family = binomial, 
#       data = df_cntrl)
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)      -0.5184     0.2710  -1.913   0.0557 .  
# ctrl_typetravel  -0.1507     0.2747  -0.549   0.5832    
# Run2              0.3008     0.2749   1.094   0.2738    
# sexMale           1.3366     0.2757   4.848 1.24e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 65.234  on 23  degrees of freedom
# Residual deviance: 38.975  on 20  degrees of freedom
# AIC: 105.46
# 
# Number of Fisher Scoring iterations: 4



