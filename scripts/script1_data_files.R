library(tidyverse)

# For each of the 5 data sets in /data, provide a statement to load into
# Global Enivornment (active memory), and provide a short example and
# destription


##### Open and describe Experiment 1 flight cylinder data ################# #

df_flight <- readRDS("./data/dat_expt1_flt_cylinder.Rds")
df_flight
# # A tibble: 72 × 5
#   Run   sex       y  dose     n
#   <fct> <chr> <dbl> <dbl> <dbl>
# 1 1     Male      9     0    10
# 2 1     Male      4     0    10
# 3 1     Male      6     0    10
# 4 2     Male      8     0    10


# Descriptions
# Run  -- Cohort or rep in time, levels = 1 and 2
# sex  -- Male, Female
# y    -- Number of moths leaving the flight cylinder
# dose -- Amount in Gy of the e-beam irradiation received
# n    -- denominator; each flight cylinder had 10 moths

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

# Descriptions
# run -- as above
# sex -- as above (m,f)
# dose -- as above
# cage -- physical cage for flt cyl & longevity test, 10 moths/cage
# individual_moth -- distinguishes moths withi the same cage
# survival_time -- when moth was found dead, between 1 and 14
# status -- Completed if found dead, Censored if still alive on day 14

# NB df_flight and df_surv descripe the same group of moths

##### Open and describe Experiment 2 mating ################# #

# Different moths and different experiment from last two data sets

df_mating <- readRDS("./data/dat_expt2_mating.Rds")
df_mating
# # A tibble: 660 × 4
#     run  Dose Sex   Mated
#   <dbl> <dbl> <fct> <fct>
# 1     1     6 Male  No   
# 2     1     6 Male  Yes  
# 3     1     6 Male  No   
# 4     1     6 Male  Yes

# Labels are as above

##### Open and describe Experiment 2 oviposition ################# #

df_ovip <- read_rds("./data/dat_expt2_oviposition.Rds")

df_ovip
# # A tibble: 349 × 4
#     run  dose sex    total_eggs
#   <dbl> <dbl> <chr>       <dbl>
# 1     2     0 Female         11
# 2     2     0 Female        162
# 3     2     0 Female         21
# 4     2     0 Female         93
# 5     2     0 Female         35

# Same as moths in df_mating
# Fewer moths because unmated moths dropped
# Sex refers to which parent was irradiated (obviously all eggs from females!)


##### Open and describe Experiment 2 development ################# #

df_dev <- readRDS("./data/dat_expt2_development.Rds")
df_dev
# # A tibble: 349 × 4
#   Run    dose sex    pupae_larvae
#   <fct> <dbl> <chr>         <dbl>
# 1 1         0 Female            1
# 2 1         0 Female           57
# 3 1         0 Female            2
# 4 1         0 Female           33
# 5 1         0 Female           15
# 6 1         0 Female           35
# 7 1         0 Female           31

# These are progeny from the same moths as the other expt2 data sets