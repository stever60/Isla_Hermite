# Code for compositing different sites into pne csv file

# -------------------------------------------------------------------------

# SECTION 1: SET UP 

# -------------------------------------------------------------------------

#clear previous console
remove (list = ls())
#clear plot window
dev.off()

# setup workspace ----
library(itraxR)
library(tidyverse) # all Section tidyverse packages
library(tidypaleo) # Dewey Dunnington's ggplot extensions for palaeo-style plots
library(dplyr)
library(readr)
library(ggpubr) # for ggarrange plotting
library(GGally) # for correlation and Prob density matrix plotting
library(PeriodicTable)
library(errors)
library(chemometrics)
library(psych) # for generating summary stats
library(compositions) # for generating clr 

# -------------------------------------------------------------------------

# SECTION 2: IMPORT DATA & COMBINE FILES 

# -------------------------------------------------------------------------

# set working directory
setwd("/Users/Steve/Dropbox/BAS/Data/R/")


# MSCL Overlaps dataset for ALL sites from original files ------------------------------------------

# import site data
HER14L <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HER14L_OVERLAPS.csv") 
HER24L <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HER24L_OVERLAPS.csv") 
HER34L <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HER34L_OVERLAPS.csv")
HER42L <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HER42L_OVERLAPS.csv") #NEEDS FIXING
HER42PB <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HER42PB_OVERLAPS.csv")
HER44L <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HER44L_OVERLAPS.csv") #NEEDS FIXING
HER49L <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HER49L_OVERLAPS.csv") #NEEDS FIXING
HERPT1S1 <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HERPT1S1_OVERLAPS.csv") #NEEDS FIXING
HERPT1S4 <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HERPT1S4_OVERLAPS.csv") #NEEDS FIXING
HERPT1S8 <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HERPT1S8_OVERLAPS.csv") #NEEDS FIXING
HERPT2S1 <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HERPT2S1_OVERLAPS.csv") #NEEDS FIXING
# bind all data and rename to match itraxR column names
HER_MSCL <-  bind_rows(HER14L, HER24L, HER34L, HER42L, HER42PB, HER44L, HER49L, HERPT1S1, HERPT1S4, HERPT1S8, HERPT2S1) %>% 
  select(SITE:SECTION_ID, Strat_depth:WMAR_Mean) %>% 
  rename(Site = SITE, depth = Strat_depth, label = SECTION_ID,
         min_age = SH20_min_age_95CI, max_age = SH20_max_age_95CI,
         median_age = SH20_median_age, mean_age = SH20_mean_age, accum_rate = Accum_rate)
HER14L_0ka <- HER_MSCL %>% 
  filter(Site=='HER14L' & mean_age <10000) %>% 
  mutate(Site = recode(Site, 'HER14L' = 'HER14L_0ka'))
HER14L_10ka <- HER_MSCL %>% 
  filter(Site=='HER14L' & mean_age >10000) %>% 
  mutate(Site = recode(Site, 'HER14L' = 'HER14L_10ka'))
HER_MSCL_overlap <- HER_MSCL %>% 
  filter(!(Site=='HER14L')) %>%
  bind_rows(HER14L_0ka, HER14L_10ka) %>% 
  arrange(Site)
head(HER_MSCL_overlaps)
tail(HER_MSCL_overlaps)
# remove site records from Global Environment
rm(HER14L, HER14L_0ka, HER14L_10ka, HER24L, HER34L, HER42L, HER42PB, HER44L, HER49L, HERPT1S1, HERPT1S4, HERPT1S8, HERPT2S1, HER_MSCL)
# write to file
write.csv(HER_MSCL_overlaps,"Papers/Hodgson_Hermite_2022/Data/MSCL/Outputs/HER_MSCL_overlaps.csv", row.names = FALSE)

# MSCL Composite dataset for ALL sites from original files ------------------------------------------

# import site data
HER14L <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HER14L_COMPOSITE.csv") 
HER24L <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HER24L_COMPOSITE.csv") 
HER34L <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HER34L_COMPOSITE.csv")
HER42L <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HER42L_COMPOSITE.csv")
HER42PB <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HER42PB_COMPOSITE.csv")
HER44L <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HER44L_COMPOSITE.csv")
HER49L <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HER49L_COMPOSITE.csv")
HERPT1S1 <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HERPT1S1_COMPOSITE.csv")
HERPT1S4 <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HERPT1S4_COMPOSITE.csv")
HERPT1S8 <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HERPT1S8_COMPOSITE.csv")
HERPT2S1 <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HERPT2S1_COMPOSITE.csv")
# bind all data and rename to match itraxR column names
HER_MSCL <-  bind_rows(HER14L, HER24L, HER34L, HER42L, HER42PB, HER44L, HER49L, HERPT1S1, HERPT1S4, HERPT1S8, HERPT2S1) %>% 
  select(SITE:SECTION_ID, Strat_depth:WMAR_Mean) %>% 
  rename(Site = SITE, depth = Strat_depth, label = SECTION_ID,
         min_age = SH20_min_age_95CI, max_age = SH20_max_age_95CI,
         median_age = SH20_median_age, mean_age = SH20_mean_age, accum_rate = Accum_rate)
HER14L_0ka <- HER_MSCL %>% 
  filter(Site=='HER14L' & mean_age <10000) %>% 
  mutate(Site = recode(Site, 'HER14L' = 'HER14L_0ka'))
HER14L_10ka <- HER_MSCL %>% 
  filter(Site=='HER14L' & mean_age >10000) %>% 
  mutate(Site = recode(Site, 'HER14L' = 'HER14L_10ka'))
HER_MSCL_comp <- HER_MSCL %>% 
  filter(!(Site=='HER14L')) %>%
  bind_rows(HER14L_0ka, HER14L_10ka) %>% 
  arrange(Site)
head(HER_MSCL_comp)
tail(HER_MSCL_comp)
# remove site records from Global Environment
rm(HER14L, HER14L_0ka, HER14L_10ka, HER24L, HER34L, HER42L, HER42PB, HER44L, HER49L, HERPT1S1, HERPT1S4, HERPT1S8, HERPT2S1, HER_MSCL)
# write to file
write.csv(HER_MSCL_comp,"Papers/Hodgson_Hermite_2022/Data/MSCL/Outputs/HER_MSCL_comp.csv", row.names = FALSE)