# Uses DescTools for KW analysis by Run * Sex, messes w
# parametric models to separate trts w/i Run * Sex

library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(car)
library(emmeans)
library(multcomp)
library(multcompView)                 # Used with Dunn test
library(DescTools)                    # nicer kruskal-wallis output
library(FSA)                          # for dunnTest()

# Reconstruct the data set
df_dev <- readRDS("./data/dat_exp2_pupae.Rds")

# # Log transform for Y axis of ggplot
# df_dev <- df_dev %>% 
#   mutate(logtrans = log2(pupae_larvae + 1))
    # That's crap!

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
# # A tibble: 349 × 6
#   Run      gy sex    pupae_larvae logtrans dose2
#   <fct> <dbl> <chr>         <dbl>    <dbl> <fct>
# 1 1         0 Female            1     1    0    
# 2 1         0 Female           57     5.86 0    
# 3 1         0 Female            2     1.58 0    
# 4 1         0 Female           33     5.09 0    
# 5 1         0 Female           15     4    0    

# Use DescTools and KW to test within Run*Sex.

############## Run 1. Test 0 and 5 in females ####################### #

df_run1fem <- df_dev %>% 
  filter(Run == 1) %>% 
  filter(sex == "Female" & dose2 %in% c(0,5))

df_run1fem %>% 
  group_by(dose2) %>% 
  summarise(n_gt0 = sum(pupae_larvae > 0),
            nObs = n(),
            pct_fer = 100*n_gt0/nObs)
# # A tibble: 2 × 4
#   dose2 n_gt0  nObs pct_fer
#   <fct> <int> <int>   <dbl>
# 1 0        13    15    86.7
# 2 5         2     8    25

Desc(pupae_larvae ~ dose2, data = df_run1fem)
# 
# pupae_larvae ~ dose2 (df_run1fem)
# 
# Summary: 
#   n pairs: 23, valid: 23 (100.0%), missings: 0 (0.0%), groups: 2
# 
# 
#               0        5
# mean     18.533    0.250
# median    8.000    0.000
# sd       19.508    0.463
# IQR      32.000    0.250
# n            15        8
# np      65.217%  34.783%
# NAs           0        0
# 0s            2        6
# 
# Kruskal-Wallis rank sum test:
#   Kruskal-Wallis chi-squared = 10.464, df = 1, p-value = 0.001217

# See remarks in Run2 fems concerning kruskal.wallis vs. parametric models


############## Run 1. Test 0, 5, 175,and 210 in males ####################### #

# Select data
df_run1male <- df_dev %>% 
  filter(Run == 1) %>% 
  filter(sex == "Male" & dose2 %in% c(0,5,175,210))

# Look at data distribution
df_run1male %>% 
  group_by(dose2) %>% 
  summarise(n_gt0 = sum(pupae_larvae > 0),
            nObs = n(),
            pct_fer = 100*n_gt0/nObs)
# # A tibble: 4 × 4
#   dose2 n_gt0  nObs pct_fer
#   <fct> <int> <int>   <dbl>
# 1 0        16    23   69.6 
# 2 5         8    16   50   
# 3 175       3    14   21.4 
# 4 210       2    23    8.70

# Non-parametric
Desc(pupae_larvae ~ dose2, data = df_run1male)
# pupae_larvae ~ dose2 (df_run1male)
# 
# Summary: 
#   n pairs: 76, valid: 76 (100.0%), missings: 0 (0.0%), groups: 4
# 
# 
#              0        5      175      210
# mean     26.739   12.625    0.929    0.130
# median   23.000    0.500    0.000    0.000
# sd       28.604   22.259    2.921    0.458
# IQR      49.500   14.250    0.000    0.000
# n            23       16       14       23
# np      30.263%  21.053%  18.421%  30.263%
# NAs           0        0        0        0
# 0s            7        8       11       21
# 
# Kruskal-Wallis rank sum test:
#   Kruskal-Wallis chi-squared = 24.462, df = 3, p-value = 2e-05

# Dunn test
out <- dunnTest(pupae_larvae ~ dose2, data = df_run1male, method = "holm")

# Convert result into a named matrix of p-values
pvals <- out$res$P.adj
names(pvals) <- paste(out$res$Comparison)
pvals
#      0 - 175      0 - 210    175 - 210        0 - 5      175 - 5      210 - 5 
# 3.118181e-03 2.650668e-05 5.670486e-01 2.139387e-01 2.482526e-01 4.353834e-02 




############## Run 2. Test 5 and 175 in females ####################### #

df_run2fem <- df_dev %>% 
  filter(Run == 2) %>% 
  filter(sex == "Female" & dose2 %in% c(5,175))

df_run2fem %>% 
  group_by(dose2) %>% 
  summarise(n_gt0 = sum(pupae_larvae > 0),
            nObs = n(),
            pct_fer = 100*n_gt0/nObs
  )
# A tibble: 2 × 4
# dose2 n_gt0  nObs pct_fer
#   <fct> <int> <int>   <dbl>
# 1 5        11    18    61.1
# 2 175       1     8    12.5


Desc(pupae_larvae ~ dose2, data = df_run2fem)
# pupae_larvae ~ dose2 (df_run2fem)
# 
# Summary: 
#   n pairs: 26, valid: 26 (100.0%), missings: 0 (0.0%), groups: 2
# 
# 
#               5      175
# mean     30.833    0.500
# median    3.000    0.000
# sd       36.378    1.414
# IQR      61.500    0.000
# n            18        8
# np      69.231%  30.769%
# NAs           0        0
# 0s            7        7
# 
# Kruskal-Wallis rank sum test:
#   Kruskal-Wallis chi-squared = 5.1419, df = 1, p-value = 0.02336


############## Run 2. Test 5, 175 and 210 in males ####################### #

df_run2male <- df_dev %>% 
  filter(Run == 2) %>% 
  filter(sex == "Male" & dose2 %in% c(5,175,210))

df_run2male %>% 
  group_by(dose2) %>% 
  summarise(n_gt0 = sum(pupae_larvae > 0),
            nObs = n(),
            pct_fert = 100*n_gt0/nObs)
# # A tibble: 3 × 4
#   dose2 n_gt0  nObs pct_fert
#   <fct> <int> <int>    <dbl>
# 1 5         5    14    35.7 
# 2 175       1    15     6.67
# 3 210       2    14    14.3 

Desc(pupae_larvae ~ dose2, data = df_run2male)
# pupae_larvae ~ dose2 (df_run2fem)
# 
# Summary: 
#   n pairs: 26, valid: 26 (100.0%), missings: 0 (0.0%), groups: 2
# 
# 
#               5      175
# mean     30.833    0.500
# median    3.000    0.000
# sd       36.378    1.414
# IQR      61.500    0.000
# n            18        8
# np      69.231%  30.769%
# NAs           0        0
# 0s            7        7
# 
# Kruskal-Wallis rank sum test:
#   Kruskal-Wallis chi-squared = 5.1419, df = 1, p-value = 0.02336

