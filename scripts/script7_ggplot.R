library(tidyverse)
library(scales)

source("./scripts/themes.R")

# p1 -- Flight cylinder -- saved line 38
# p2 -- Longevity -- saved line 91
# p3 -- Proportion mated -- saved line 137
# p4 -- Total feduncity -- saved line 189
# p5 -- Pupal count -- saved line 

########### Expt 1 Flight cylinder figure ################################# #

# Retrieve data
df_flight <- readRDS("./data/dat_expt1_flt_cylinder.Rds")
df_flight
# # A tibble: 72 × 5
#   Run   sex       y  dose     n
#   <fct> <chr> <dbl> <dbl> <dbl>
# 1 1     Male      9     0    10
# 2 1     Male      4     0    10
# 3 1     Male      6     0    10

df_esc_x_sex <- df_flight %>% 
  group_by(sex,Run,n) %>% 
  mutate(pct_esc = 100*y/n)

df_esc_x_sex <- df_esc_x_sex %>% 
  rename(Sex = sex)

p1 <- ggplot(df_esc_x_sex, aes(x = dose, y = pct_esc, group = Sex, color = Sex)) +
  geom_point(aes(shape = Sex), stroke = 1.5, size = 2) +
  geom_smooth(aes(linetype = Sex), method = "lm", se = TRUE) +
  scale_shape_manual(values = c("Female" = 1, "Male" = 4)) +
  scale_linetype_manual(values = c("Female" = "solid", "Male" = "dashed")) + 
  labs(x = "Dose in Gy",
       y = "Percent of adults leaving the cylinder") +
  theme_fullwidth()

p1

ggsave(filename = "./results/fig2_surv.jpg",
       plot = p2, device = "jpg",
       dpi = 300, width = 5.83, height = 3.5, units = "in")

########### Expt 1 Longevity figure ####################################### #

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

df_medians <- df_surv %>% 
  group_by(run, sex, cage, dose) %>% 
  summarise(nObs = n(),
            med = median(survival_time),
            n_surv = sum(status == "Censored"))

df_medians <- df_medians %>% 
  mutate(pct_14dsurv = 100*n_surv/nObs)

df_survmeans <- df_medians %>% 
  group_by(sex,dose) %>% 
  summarise(nObs = n(),
            mn = mean(pct_14dsurv),
            sem = FSA::se(pct_14dsurv))

df_survmeans
## Prettier variable names for plotting
df_survmeans$sex[df_survmeans$sex == "f"] <- "Female"
df_survmeans$sex[df_survmeans$sex == "m"] <- "Male"



p2 <- ggplot(df_survmeans,aes(x = dose, y = mn)) +
  geom_point() +
  geom_errorbar(aes(ymin = mn - sem, ymax = mn + sem)) +
  geom_smooth(method = "lm") +
  facet_grid(sex ~ .) +
  # ylim(c(0,15)) +
  xlab("Dose (Gy)") +
  ylab("Survival at 14 days (%)") +
  facet_grid(sex ~ .)+
  theme_powerpoint()

p2
 
ggsave("./results/fig2_surv_by_sex_run.jpg",
       plot = p2,
       device = "jpg",
       width = 5.83,
       height = 5.83,
       units = "in",
       dpi = 300)

########### Expt 2 Mating figure ####################################### #

df_mating <- readRDS("./data/dat_expt2_mating.Rds")

df_pmated <- df_mating %>% 
  group_by(run,Dose,Sex,Mated) %>%
  summarise(nObs = n()) %>% 
  pivot_wider(names_from = Mated, values_from = nObs) %>% 
  rename(mated = Yes) %>% 
  mutate(trials = No + mated, prop = mated/trials) 


#%>% 
#  select(run,Dose,Sex,mated,trials,prop)

# Make run into a factor
df_pmated <- df_pmated %>% 
  rename(Run = run)

df_pmated$Run <- factor(df_pmated$Run, levels = c(1,2))


# Plot proportion mated by vs dose by sex, with run color-coded

p3 <- ggplot(df_pmated2, aes(x = Dose, y = prop, color = Run)) +
  geom_point(aes(shape = Run), stroke = 1.5, size = 2) +
  geom_smooth(aes(group = Run, linetype = Run), method = "lm", se = TRUE) +
  scale_shape_manual(values = c("1" = 1, "2" = 4)) +
  scale_linetype_manual(values = c("1" = "solid", "2" = "dashed")) + 
  labs(x = "Dose in Gy",
       y = "Proportion mated") +
  ylim(c(0,1)) +
  facet_grid(Sex ~ .) +
  theme_powerpoint()  

p3

ggsave("./results/fig3_mating_by_sex_run.jpg",
       plot = p3,
       device = "jpg",
       width = 5.83,
       height = 5.83,
       units = "in",
       dpi = 300)

########### Expt 2 Oviposition figure ##################################### #

df_ovip <- readRDS("./data/dat_expt2_oviposition.Rds")

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
p4 <- ggplot(df_ovip, aes(x = dose2, y = total_eggs)) +
  geom_boxplot()  +
  geom_text(
    data = sig_letters2,
    aes(x = dose2, y = ypos, label = label),
    size = 5
  ) +
  coord_cartesian(ylim = c(0, 250)) +
  xlab("Dose (Gy)") +
  ylab("Eggs per female") +
  facet_grid(sex ~ cohort_name) + 
  theme_powerpoint()

p4

ggsave(
  "fig_ovipos_vbar_error_sex_by_run.jpg",
  device = "jpg",
  plot = p4,
  path = "results",
  width = 5.83,
  height = 5.83,
  units = "in",
  dpi = 300
)

########### Pupal devopment figure ############################### #

df_dev <- readRDS("./data/dat_expt2_pupae.Rds")

# Apply dose categories
df_dev <- df_dev %>% 
  mutate(dose2 = case_when(
    gy < 5 ~ 0,              # 0 only
    gy < 100 ~ 5,            # 6 only
    gy < 200 ~ 175,          # 176, 177
    gy < 250 ~ 210,          # 208, 211
    gy < 315 ~ 310,          # 310
    .default = 420             # 405, 445
  ))

# Set dose2 data type as factor
df_dev$dose2 <- factor(df_dev$dose2)

df_dev <- df_dev %>% 
  mutate(cohort = if_else(Run == 1,"Run 1","Run 2"))

p5 <- ggplot(df_dev, aes(x = dose2, y = pupae_larvae)) +
  geom_boxplot()   +
  geom_text(
    data = sig_labels,
    aes(x = dose2, y = ypos, label = label),
    size = 5
  )  +
  coord_cartesian(ylim = c(0, 100)) +
  xlab("Dose (Gy)") +
  ylab("Pupae per mated pair") +
  facet_grid(sex ~ cohort) +
  theme_powerpoint()

p5

ggsave("./results/fig5_ebeam_fertility.png",
       plot = p5,
       device = "png",
       width = 6,
       heigh = 4,
       units = "in",
       dpi = 300)
