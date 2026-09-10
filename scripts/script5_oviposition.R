# Examine oviposition data from experiment 2. Examines for and finds over-
# dispersion, and determines that GLM w negative binomial is suitable.
# Initially examines oviposition vs. dose with dose as a continuous variable,
# but then grouped similar doses as categorical predictors because that 
# approach is more informative. Output includes:
#  - A Type II Analysis of Deviance for dose, sex, run, dose*sex & dose*run
#    (all predictors categorical, GLM w nb) (contained in lines 318-322)
#  - A vertical bar chart with SE for eggs per female vs. dose, arranged
#    as facet_grid(sex ~ cohort) (output lines 249-256)
#  - emmeans/cld tukey groupings for all 18 levels of dose*sex*run (line 385)

# Load libraries
library(tidyverse)
library(DHARMa)
library(MASS)     # for overdispersion and nb
library(car)
library(FSA)      # for se()
library(DescTools)
library(emmeans)
library(estimability)
library(multcomp)

source("./scripts/themes.R")

###### 1. Load data set and obtain an overview ######################## #

df_ovip <- readRDS("./data/dat_expt2_oviposition.Rds")
df_ovip
# # A tibble: 349 × 4
#     run  dose sex    total_eggs
#   <dbl> <dbl> <chr>       <dbl>
# 1     1     0 Female         11
# 2     1     0 Female        162
# 3     1     0 Female         21

df_ovip %>%
  group_by(run,dose) %>% 
  summarise(nObs = n())
# # A tibble: 11 × 3
# # Groups:   run [2]
#      run  dose  nObs
#    <dbl> <dbl> <int>
#  1     1     0    38
#  2     1     6    24
#  3     1   177    29
#  4     1   211    40
#  5     1   405    30
#  6     1   445    37
#  7     2     6    32
#  8     2   176    23
#  9     2   208    34
# 10     2   309    39
# 11     2   310    23

# A half an egg is an obvious entry error, and throws off count-based statistics
df_ovip$total_eggs[df_ovip$total_eggs == 40.5] <- 40

## nObs < 60 because unmated moths are excluded

##### 2. Plot total fecundity vs dose by sex ########## #

# Make both run and dose into factors
df_ovip$run <- factor(df_ovip$run, levels = c(1,2))

# Get dose into discreet groups
df_ovip <- df_ovip %>% 
  mutate(dose2 = case_when(
    dose < 5 ~ 0,              # 0 only
    dose < 100 ~ 5,            # 6 only
    dose < 200 ~ 175,          # 176, 177
    dose < 250 ~ 210,          # 208, 211
    dose < 315 ~ 310,          # 310
    .default = 420             # 405, 445
  ))

# Make dose2 a factor
df_ovip$dose2 <- factor(df_ovip$dose2)

df_ovip <- df_ovip |> 
  mutate(cohort_name =ifelse(run == 1,"Run 1","Run 2")) 

# Box plots, sex by run
p1 <- ggplot(df_ovip, aes(x = dose2, y = total_eggs)) +
  geom_boxplot()  +
  xlab("Dose (Gy)") +
  ylab("Eggs per female") +
  facet_grid(sex ~ cohort_name) + 
  theme_powerpoint()

p1

ggsave(
  "fig_ovipos_vbar_error_sex_by_run.jpg",
  device = "jpg",
  plot = p1,
  path = "results",
  width = 5.83,
  height = 5.83,
  units = "in",
  dpi = 300
)


####### Examine data distribution and neg. bin. vs Poisson ############# #

## Fit a Poison model to test for overdispersion
m_pois <- glm(total_eggs ~ dose2*sex + dose2*run + sex*run, family = poisson, data = df_ovip )
summary(m_pois)
# 
# Call:
#   glm(formula = total_eggs ~ dose2 * sex + dose2 * run + sex * 
#         run, family = poisson, data = df_ovip)
# 
# Coefficients: (3 not defined because of singularities)
# Estimate Std. Error z value             Pr(>|z|)    
# (Intercept)       4.21262    0.03142 134.078 < 0.0000000000000002 ***
# dose25           -0.88864    0.06481 -13.712 < 0.0000000000000002 ***
# dose2175         -1.32404    0.06226 -21.265 < 0.0000000000000002 ***
# dose2210         -1.40636    0.06064 -23.192 < 0.0000000000000002 ***
# dose2310         -1.91314    0.07266 -26.331 < 0.0000000000000002 ***
# dose2420         -0.79639    0.04452 -17.888 < 0.0000000000000002 ***
# sexMale          -0.02297    0.04057  -0.566             0.571328    
# run2              1.18057    0.05574  21.178 < 0.0000000000000002 ***
# dose25:sexMale    0.26646    0.07359   3.621             0.000293 ***
# dose2175:sexMale  0.84705    0.07299  11.606 < 0.0000000000000002 ***
# dose2210:sexMale  0.76041    0.06981  10.893 < 0.0000000000000002 ***
# dose2310:sexMale  0.95595    0.08955  10.675 < 0.0000000000000002 ***
# dose2420:sexMale  0.16544    0.05898   2.805             0.005031 ** 
# dose25:run2      -0.13101    0.07000  -1.872             0.061252 .  
# dose2175:run2    -0.31530    0.07212  -4.372            0.0000123 ***
# dose2210:run2          NA         NA      NA                   NA    
# dose2310:run2          NA         NA      NA                   NA    
# dose2420:run2          NA         NA      NA                   NA    
# sexMale:run2     -1.63872    0.05878 -27.877 < 0.0000000000000002 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 15577  on 348  degrees of freedom
# Residual deviance: 12605  on 333  degrees of freedom
# AIC: 14195
# 
# Number of Fisher Scoring iterations: 5

overdispersion_ratio <- sum(residuals(m_pois, type = "pearson")^2) / m_pois$df.residual
overdispersion_ratio
#[1] 49.19061 -- Overdispersed if > 1.5, Poisson inappropriate if > 2!

# Run a negative binomial model and compare to poisson

m_nb <- glm.nb(total_eggs ~ dose2*sex + dose2*run + sex*run, data = df_ovip)
summary(m_nb)
# 
# Call:
#   glm.nb(formula = total_eggs ~ dose2 * sex + dose2 * run + sex * 
#            run, data = df_ovip, init.theta = 0.751339568, link = log)
# 
# Coefficients: (3 not defined because of singularities)
# Estimate Std. Error z value             Pr(>|z|)    
# (Intercept)       4.21262    0.29953  14.064 < 0.0000000000000002 ***
# dose25           -0.89026    0.46137  -1.930             0.053655 .  
# dose2175         -1.45056    0.41106  -3.529             0.000417 ***
# dose2210         -1.34788    0.39838  -3.383             0.000716 ***
# dose2310         -1.92579    0.49896  -3.860             0.000114 ***
# dose2420         -0.79639    0.36200  -2.200             0.027810 *  
# sexMale          -0.02297    0.38502  -0.060             0.952435    
# run2              1.19322    0.32690   3.650             0.000262 ***
# dose25:sexMale    0.28245    0.54504   0.518             0.604304    
# dose2175:sexMale  0.89750    0.53008   1.693             0.090432 .  
# dose2210:sexMale  0.73481    0.50136   1.466             0.142751    
# dose2310:sexMale  0.98331    0.60617   1.622             0.104765    
# dose2420:sexMale  0.16544    0.47912   0.345             0.729876    
# dose25:run2      -0.10599    0.42721  -0.248             0.804060    
# dose2175:run2    -0.32854    0.43335  -0.758             0.448366    
# dose2210:run2          NA         NA      NA                   NA    
# dose2310:run2          NA         NA      NA                   NA    
# dose2420:run2          NA         NA      NA                   NA    
# sexMale:run2     -1.66607    0.35793  -4.655           0.00000324 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for Negative Binomial(0.7513) family taken to be 1)
# 
# Null deviance: 475.9  on 348  degrees of freedom
# Residual deviance: 413.5  on 333  degrees of freedom
# AIC: 3144.7
# 
# Number of Fisher Scoring iterations: 1
# 
# 
# Theta:  0.7513 
# Std. Err.:  0.0549 
# 
# 2 x log-likelihood:  -3110.7450 


AIC(m_pois, m_nb)
#       df       AIC
# m_pois 16 14194.526
# m_nb   17  3144.745     -- makes the same point
  # Model changed with script iterations. This was the same simple model for
  # both distributions

sim_nb <- simulateResiduals(m_nb, n = 1000)

# Residual diagnostics
plot(sim_nb)
testDispersion(sim_nb)
testZeroInflation(sim_nb)
testOutliers(sim_nb)

# Residuals
plotResiduals(sim_nb, form = df_ovip$dose)
plotResiduals(sim_nb, form = df_ovip$run)
plotResiduals(sim_nb, form = df_ovip$sex)
    # The negative binomial model passes all the tests


# Analysis of Deviance Table
car::Anova(m_nb, type = 2)
# Analysis of Deviance Table (Type II tests)
# 
# Response: total_eggs
# LR Chisq Df  Pr(>Chisq)    
# dose2      23.5143  5   0.0002691 ***
# sex         2.5662  1   0.1091701    
# run         0.9206  1   0.3373168    
# dose2:sex   5.7253  5   0.3338758    
# dose2:run   0.5928  2   0.7435012    
# sex:run    21.2655  1 0.000003999 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

   # Dose and the sex:run interaction were the only two significant effects
   # Examine slices of sex:run separately

####### Examine each level of Run*Sex separate ############# #

# Create subset data frames
df_f1 <- df_ovip %>% filter(sex == "Female", run == "1")
df_f2 <- df_ovip %>% filter(sex == "Female", run == "2")
df_m1 <- df_ovip %>% filter(sex == "Male",   run == "1")
df_m2 <- df_ovip %>% filter(sex == "Male",   run == "2")

# Use Desc to explore each
Desc(total_eggs ~ dose2, df_f1)
# 
# Summary: 
#   n pairs: 88, valid: 88 (100.0%), missings: 0 (0.0%), groups: 5
# 
# 
#               0        5      175      210      420
# mean     67.533   31.125   12.133   20.118   30.455
# median   53.000    9.000    8.000   14.000   22.000
# sd       56.639   44.147   14.086   19.580   27.049
# IQR      72.000   33.500    8.500   26.000   35.000
# n            15        8       15       17       33
# np      17.045%   9.091%  17.045%  19.318%  37.500%
# NAs           0        0        0        0        0
# 0s            0        1        0        1        2
# 
# Kruskal-Wallis rank sum test:
#   Kruskal-Wallis chi-squared = 17.189, df = 4, p-value = 0.001776

x <- dunnTest(total_eggs ~ dose2, data=df_f1, method = "holm")
x  # A list

x <- x$res %>% 
  dplyr::select(Comparison, P.adj)
  
Desc(total_eggs ~ dose2, df_m1)
# 
# Summary: 
#   n pairs: 110, valid: 110 (100.0%), missings: 0 (0.0%), groups: 5
# 
# 
#               0        5      175      210      420
# mean     66.000   33.750   47.214   31.957   35.118
# median   52.000   18.000   33.500   17.000   21.500
# sd       58.938   41.258   56.096   34.563   40.001
# IQR      76.500   34.000   27.000   37.000   36.500
# n            23       16       14       23       34
# np      20.909%  14.545%  12.727%  20.909%  30.909%
# NAs           0        0        0        0        0
# 0s            0        0        1        1        2
# 
# Kruskal-Wallis rank sum test:
#   Kruskal-Wallis chi-squared = 6.2428, df = 4, p-value = 0.1817

Desc(total_eggs ~ dose2, df_f2)
# total_eggs ~ dose2 (df_f2)
# 
# Summary: 
#   n pairs: 72, valid: 72 (100.0%), missings: 0 (0.0%), groups: 4
# 
# 
#               5      175      210      310
# mean     77.833   53.625   50.850   32.462
# median   25.500   42.000   43.500   12.000
# sd       76.048   41.214   44.152   43.210
# IQR     146.750   58.750   67.000   52.500
# n            18        8       20       26
# np      25.000%  11.111%  27.778%  36.111%
# NAs           0        0        0        0
# 0s            2        0        3        3
# 
# Kruskal-Wallis rank sum test:
#   Kruskal-Wallis chi-squared = 5.4584, df = 3, p-value = 0.1411

Desc(total_eggs ~ dose2, df_m2)
# 
# Summary: 
#   n pairs: 79, valid: 79 (100.0%), missings: 0 (0.0%), groups: 4
# 
# 
#               5      175      210      310
# mean     21.571   13.067   26.214   16.028
# median    5.500    9.000   13.000   10.000
# sd       29.589   15.102   38.120   23.069
# IQR      24.750    8.500   19.000   14.500
# n            14       15       14       36
# np      17.722%  18.987%  17.722%  45.570%
# NAs           0        0        0        0
# 0s            3        1        2        4
# 
# Kruskal-Wallis rank sum test:
#   Kruskal-Wallis chi-squared = 1.0143, df = 3, p-value = 0.7978
