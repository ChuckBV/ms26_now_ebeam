

#### Create 2nd tibble for ovip plot to create sig. labels ################ #

sig_letters  <- read_csv("dose_labels.csv")

glimpse(sig_letters)

sig_letters$dose2 <- factor(sig_letters$dose2)

sig_letters <- sig_letters %>% 
  rename(cohort_name = cohort)

sig_letters$cohort_name[sig_letters$cohort_name == "Run1"] <- "Run 1"
sig_letters$cohort_name[sig_letters$cohort_name == "Run2"] <- "Run 2"

df_ovip

sig_letters <- sig_letters %>% 
  rename(sex = Sex)

sig_letters2 <- sig_letters %>% 
  mutate(ypos = 210)

saveRDS(sig_letters2,"./data/dat_expt2_ovip_sig_labels.Rds")

#### Create 2nd tibble for pupa plot to create sig. labels ################ #

# Get data set from script1

# Create scaffold tibble--first get proper var names
colnames(df_dev)
# [1] "Run"          "gy"           "sex"          "pupae_larvae"
# [5] "dose2"        "cohort"  

# assemble components
sexes <- c(rep("Female",6),rep("Male",6))

dose_levels <- as.character(unique(df_dev$dose2))
dose_levels
# [1] "0"   "5"   "175" "210" "310" "420"
levels(dose_levels)
# NULL

labels_r1fem <- c("a","b","","","","")
labels_r1male <-  c("a","ab","bc","c"," "," ")
labels_r2fem <- c(" ","a","b"," "," "," ")
labels_r2male <- c(" ","a","a","a"," "," ")

# Assemble tibble
sig_labels <- tibble(
  cohort = c(rep("Run 1",12),rep("Run 2",12)),
  sex = c(rep(sexes,2)),
  dose2 = rep(dose_levels,4),
  label = c(labels_r1fem,labels_r1male,labels_r2fem,labels_r2male),
  ypos = rep(95,24)
)
# [1] 0   5   175 210 310 420
# Levels: 0 5 175 210 310 420