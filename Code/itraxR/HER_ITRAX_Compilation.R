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


# HER_ITRAX_COMPOSITE_cps  ------------------------------------------------

# COMPOSITE files are pre-existing downcore profiles with no overlaps
# Section join points determined from MSCL orexisting Excel ITRAX analysis
# Imports all raw data COMPOSITE cps files to use as inputs to HER-ITRAX site by site processing

# import COMPOSITE site datasets using a previously made cps composite file
HER14L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER14L/HER14L_ITRAX_COMPOSITE_cps.csv") 
HER24L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER24L/HER24L_ITRAX_COMPOSITE_cps.csv") 
HER34L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER34L/HER34L_ITRAX_COMPOSITE_cps.csv")
HER42L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42L/HER42L_ITRAX_COMPOSITE_cps.csv")
HER42PB <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB_ITRAX_COMPOSITE_cps.csv")
HER44L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER44L/HER44L_ITRAX_COMPOSITE_cps.csv")
HER49L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER49L/HER49L_ITRAX_COMPOSITE_cps.csv")
# Assign elements from PeriodicTable package to 'elementsList' - to select all elements present across all sites
data(periodicTable)
elementsList <- c(periodicTable$symb, "Mo_inc", "Mo-coh")
# bind all data and rename to match itraxR column names - 
HER_ITRAX_bind <-  bind_rows(HER14L, HER24L, HER34L, HER42L, HER42PB, HER44L, HER49L) %>% 
  rename(
    label = Section,
    position = position_mm,
    field_depth = Field_depth_cm,
    section_depth = Section_depth_cm,
    depth = Strat_depth_cm,
    min_age = SH20_min_age_95CI,
    max_age = SH20_max_age_95CI,
    median_age = SH20_median_age,
    mean_age = SH20_mean_age,
    `Fe a*2` = D1
         ) %>% 
  select(Site:ID, depth, min_age, max_age, median_age, mean_age, surface, 
         validity, kcps, MSE, any_of(elementsList), Mo_inc:last_col()) %>% 
  relocate(label, .after = Mo_coh) %>% 
  relocate(surface, .after = last_col()) %>% 
  relocate(validity, .after = last_col()) %>%
  relocate(ID, .after = last_col()) %>% 
  mutate(cps = rowSums(across(Mg:Mo_coh))) %>% #mutate(cps = rowSums(.[df1_rowsums])) %>%
  relocate(cps, .before = MSE) %>% 
  filter(validity=='1')

# split HER14L into before and after 10 ka hiatus for plotting
HER14L_0ka <- HER_ITRAX_bind %>% 
  filter(Site=='HER14L' & mean_age <10000) %>% 
  mutate(Site = recode(Site, 'HER14L' = 'HER14L_0ka'))
HER14L_10ka <- HER_ITRAX_bind %>% 
  filter(Site=='HER14L' & mean_age >10000) %>% 
  mutate(Site = recode(Site, 'HER14L' = 'HER14L_10ka'))
HER_ITRAX_COMPOSITE_cps <- HER_ITRAX_bind %>% 
  filter(!(Site=='HER14L')) %>%
  bind_rows(HER14L_0ka, HER14L_10ka) %>% 
  arrange(Site)
HER_ITRAX_COMPOSITE_cps
tail(HER_ITRAX_COMPOSITE_cps)
# remove site records from Global Environment
rm(HER14L, HER14L_0ka, HER14L_10ka, HER24L, HER34L, HER42L, HER42PB, HER44L, HER49L, HER_ITRAX_bind)
# write to file
write.csv(HER_ITRAX_COMPOSITE_cps,"Papers/Hodgson_Hermite_2022/Data/ITRAX/Outputs/HER_ITRAX_COMPOSITE_cps.csv", row.names = FALSE)


# -------------------------------------------------------------------------

# Post-processed ITRAX COMPOSITE site files 

# -------------------------------------------------------------------------

# Import xrf = all elements
# Import xrf1 = selected elements and ratios

# HER_ITRAX_cps_comp ----------------------------------------------------

# import processed composite site files xrf or xrf1
HER14L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER14L/Output/Composite/HER14L_xrf1.csv") 
HER24L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER24L/Output/Composite/HER24L_xrf1.csv") 
HER34L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER34L/Output/Composite/HER34L_xrf1.csv")
HER42L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42L/Output/Composite/HER42L_xrf1.csv")
HER42PB <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/Composite/HER42PB_xrf1.csv")
HER44L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER44L/Output/Composite/HER44L_xrf1.csv")
HER49L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER49L/Output/Composite/HER49L_xrf1.csv")
# Assign elements from PeriodicTable package to 'elementsList' - to select all elements present across all sites
data(periodicTable)
elementsList <- c(periodicTable$symb, "Mo_inc", "Mo-coh")
# bind all data and rename to match itraxR column names - 
# all NA's in a column means that element failed filtering test at that site
HER_ITRAX_bind <-  bind_rows(HER14L, HER24L, HER34L, HER42L, HER42PB, HER44L, HER49L) %>% 
  rename(mean_age = SH20_age) %>% 
  select(Site, depth, mean_age, kcps:MSE, any_of(elementsList), Mo_inc: last_col())
# split HER14L into before and after 10 ka hiatus for plotting
HER14L_0ka <- HER_ITRAX_bind %>% 
  filter(Site=='HER14L' & mean_age <10000) %>% 
  mutate(Site = recode(Site, 'HER14L' = 'HER14L_0ka'))
HER14L_10ka <- HER_ITRAX_bind %>% 
  filter(Site=='HER14L' & mean_age >10000) %>% 
  mutate(Site = recode(Site, 'HER14L' = 'HER14L_10ka'))
HER_ITRAX_xrf1_cps_comp <- HER_ITRAX_bind %>% 
  filter(!(Site=='HER14L')) %>%
  bind_rows(HER14L_0ka, HER14L_10ka) %>% 
  arrange(Site)
HER_ITRAX_xrf1_cps_comp
tail(HER_ITRAX_xrf1_cps_comp)
# remove site records from Global Environment
rm(HER14L, HER14L_0ka, HER14L_10ka, HER24L, HER34L, HER42L, HER42PB, HER44L, HER49L, HER_ITRAX_bind)
# write to file
write.csv(HER_ITRAX_xrf1_cps_comp,"Papers/Hodgson_Hermite_2022/Data/ITRAX/Outputs/HER_ITRAX_xrf1_cps_comp.csv", row.names = FALSE)

# HER_ITRAX_inc_comp ----------------------------------------------------

# import processed composite site files xrf or xrf1
HER14L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER14L/Output/Composite/HER14L_xrf1_inc.csv") 
HER24L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER24L/Output/Composite/HER24L_xrf1_inc.csv") 
HER34L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER34L/Output/Composite/HER34L_xrf1_inc.csv")
HER42L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42L/Output/Composite/HER42L_xrf1_inc.csv")
HER42PB <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/Composite/HER42PB_xrf1_inc.csv")
HER44L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER44L/Output/Composite/HER44L_xrf1_inc.csv")
HER49L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER49L/Output/Composite/HER49L_xrf1_inc.csv")
# Assign elements from PeriodicTable package to 'elementsList' - to select all elements present across all sites
data(periodicTable)
elementsList <- c(periodicTable$symb, "Mo_inc", "Mo-coh")
# bind all data and rename to match itraxR column names - 
# all NA's in a column means that element failed filtering test at that site
HER_ITRAX_bind <-  bind_rows(HER14L, HER24L, HER34L, HER42L, HER42PB, HER44L, HER49L) %>% 
  rename(mean_age = SH20_age) %>% 
  select(Site, depth, mean_age, kcps:MSE, any_of(elementsList), Mo_inc: last_col())
# split HER14L into before and after 10 ka hiatus for plotting
HER14L_0ka <- HER_ITRAX_bind %>% 
  filter(Site=='HER14L' & mean_age <10000) %>% 
  mutate(Site = recode(Site, 'HER14L' = 'HER14L_0ka'))
HER14L_10ka <- HER_ITRAX_bind %>% 
  filter(Site=='HER14L' & mean_age >10000) %>% 
  mutate(Site = recode(Site, 'HER14L' = 'HER14L_10ka'))
HER_ITRAX_xrf1_inc_comp <- HER_ITRAX_bind %>% 
  filter(!(Site=='HER14L')) %>%
  bind_rows(HER14L_0ka, HER14L_10ka) %>% 
  arrange(Site)
HER_ITRAX_xrf1_inc_comp
tail(HER_ITRAX_xrf1_inc_comp)
# remove site records from Global Environment
rm(HER14L, HER14L_0ka, HER14L_10ka, HER24L, HER34L, HER42L, HER42PB, HER44L, HER49L, HER_ITRAX_bind)
# write to file
write.csv(HER_ITRAX_xrf1_inc_comp,"Papers/Hodgson_Hermite_2022/Data/ITRAX/Outputs/HER_ITRAX_xrf1_inc_comp.csv", row.names = FALSE)


# HER_ITRAX_%cps_sum_comp ----------------------------------------------------

# import processed composite site files xrf or xrf1
HER14L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER14L/Output/Composite/HER14L_xrf1_pc_cps.csv") 
HER24L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER24L/Output/Composite/HER24L_xrf1_pc_cps.csv") 
HER34L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER34L/Output/Composite/HER34L_xrf1_pc_cps.csv")
HER42L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42L/Output/Composite/HER42L_xrf1_pc_cps.csv")
HER42PB <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/Composite/HER42PB_xrf1_pc_cps.csv")
HER44L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER44L/Output/Composite/HER44L_xrf1_pc_cps.csv")
HER49L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER49L/Output/Composite/HER49L_xrf1_pc_cps.csv")
# Assign elements from PeriodicTable package to 'elementsList' - to select all elements present across all sites
data(periodicTable)
elementsList <- c(periodicTable$symb, "Mo_inc", "Mo-coh")
# bind all data and rename to match itraxR column names - 
# all NA's in a column means that element failed filtering test at that site
HER_ITRAX_bind <-  bind_rows(HER14L, HER24L, HER34L, HER42L, HER42PB, HER44L, HER49L) %>% 
  rename(mean_age = SH20_age) %>% 
  select(Site, depth, mean_age, kcps:MSE, any_of(elementsList), Mo_inc: last_col())
# split HER14L into before and after 10 ka hiatus for plotting
HER14L_0ka <- HER_ITRAX_bind %>% 
  filter(Site=='HER14L' & mean_age <10000) %>% 
  mutate(Site = recode(Site, 'HER14L' = 'HER14L_0ka'))
HER14L_10ka <- HER_ITRAX_bind %>% 
  filter(Site=='HER14L' & mean_age >10000) %>% 
  mutate(Site = recode(Site, 'HER14L' = 'HER14L_10ka'))
HER_ITRAX_xrf1_pc_cps_comp <- HER_ITRAX_bind %>% 
  filter(!(Site=='HER14L')) %>%
  bind_rows(HER14L_0ka, HER14L_10ka) %>% 
  arrange(Site)
HER_ITRAX_xrf1_pc_cps_comp
tail(HER_ITRAX_xrf1_pc_cps_comp)
# remove site records from Global Environment
rm(HER14L, HER14L_0ka, HER14L_10ka, HER24L, HER34L, HER42L, HER42PB, HER44L, HER49L, HER_ITRAX_bind)
# write to file
write.csv(HER_ITRAX_xrf1_pc_cps_comp,"Papers/Hodgson_Hermite_2022/Data/ITRAX/Outputs/HER_ITRAX_xrf1_pc_cps_comp.csv", row.names = FALSE)

# HER_ITRAX_clr_comp ----------------------------------------------------

# import processed composite site files xrf or xrf1
HER14L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER14L/Output/Composite/HER14L_xrf1_clr.csv") 
HER24L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER24L/Output/Composite/HER24L_xrf1_clr.csv") 
HER34L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER34L/Output/Composite/HER34L_xrf1_clr.csv")
HER42L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42L/Output/Composite/HER42L_xrf1_clr.csv")
HER42PB <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/Composite/HER42PB_xrf1_clr.csv")
HER44L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER44L/Output/Composite/HER44L_xrf1_clr.csv")
HER49L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER49L/Output/Composite/HER49L_xrf1_clr.csv")
# Assign elements from PeriodicTable package to 'elementsList' - to select all elements present across all sites
data(periodicTable)
elementsList <- c(periodicTable$symb, "Mo_inc", "Mo-coh")
# bind all data and rename to match itraxR column names - 
# all NA's in a column means that element failed filtering test at that site
HER_ITRAX_bind <-  bind_rows(HER14L, HER24L, HER34L, HER42L, HER42PB, HER44L, HER49L) %>% 
  rename(mean_age = SH20_age) %>% 
  select(Site, depth, mean_age, kcps:MSE, any_of(elementsList), Mo_inc: last_col())
# split HER14L into before and after 10 ka hiatus for plotting
HER14L_0ka <- HER_ITRAX_bind %>% 
  filter(Site=='HER14L' & mean_age <10000) %>% 
  mutate(Site = recode(Site, 'HER14L' = 'HER14L_0ka'))
HER14L_10ka <- HER_ITRAX_bind %>% 
  filter(Site=='HER14L' & mean_age >10000) %>% 
  mutate(Site = recode(Site, 'HER14L' = 'HER14L_10ka'))
HER_ITRAX_xrf1_clr_comp <- HER_ITRAX_bind %>% 
  filter(!(Site=='HER14L')) %>%
  bind_rows(HER14L_0ka, HER14L_10ka) %>% 
  arrange(Site)
HER_ITRAX_xrf1_clr_comp
tail(HER_ITRAX_xrf1_clr_comp)
# remove site records from Global Environment
rm(HER14L, HER14L_0ka, HER14L_10ka, HER24L, HER34L, HER42L, HER42PB, HER44L, HER49L, HER_ITRAX_bind)
# write to file
write.csv(HER_ITRAX_xrf1_clr_comp,"Papers/Hodgson_Hermite_2022/Data/ITRAX/Outputs/HER_ITRAX_xrf1_clr_comp.csv", row.names = FALSE)


# HER_ITRAX_cpsZ_comp ----------------------------------------------------

# import processed composite site files xrf or xrf1
HER14L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER14L/Output/Composite/HER14L_xrf1_Z.csv") 
HER24L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER24L/Output/Composite/HER24L_xrf1_Z.csv") 
HER34L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER34L/Output/Composite/HER34L_xrf1_Z.csv")
HER42L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42L/Output/Composite/HER42L_xrf1_Z.csv")
HER42PB <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/Composite/HER42PB_xrf1_Z.csv")
HER44L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER44L/Output/Composite/HER44L_xrf1_Z.csv")
HER49L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER49L/Output/Composite/HER49L_xrf1_Z.csv")
# Assign elements from PeriodicTable package to 'elementsList' - to select all elements present across all sites
data(periodicTable)
elementsList <- c(periodicTable$symb, "Mo_inc", "Mo-coh")
# bind all data and rename to match itraxR column names - 
# all NA's in a column means that element failed filtering test at that site
HER_ITRAX_bind <-  bind_rows(HER14L, HER24L, HER34L, HER42L, HER42PB, HER44L, HER49L) %>% 
  rename(mean_age = SH20_age) %>% 
  select(Site, depth, mean_age, kcps:MSE, any_of(elementsList), Mo_inc: last_col())
# split HER14L into before and after 10 ka hiatus for plotting
HER14L_0ka <- HER_ITRAX_bind %>% 
  filter(Site=='HER14L' & mean_age <10000) %>% 
  mutate(Site = recode(Site, 'HER14L' = 'HER14L_0ka'))
HER14L_10ka <- HER_ITRAX_bind %>% 
  filter(Site=='HER14L' & mean_age >10000) %>% 
  mutate(Site = recode(Site, 'HER14L' = 'HER14L_10ka'))
HER_ITRAX_xrf1_cpsZ_comp <- HER_ITRAX_bind %>% 
  filter(!(Site=='HER14L')) %>%
  bind_rows(HER14L_0ka, HER14L_10ka) %>% 
  arrange(Site)
HER_ITRAX_xrf1_cpsZ_comp
tail(HER_ITRAX_xrf1_cpsZ_comp)
# remove site records from Global Environment
rm(HER14L, HER14L_0ka, HER14L_10ka, HER24L, HER34L, HER42L, HER42PB, HER44L, HER49L, HER_ITRAX_bind)
# write to file
write.csv(HER_ITRAX_xrf1_cpsZ_comp,"Papers/Hodgson_Hermite_2022/Data/ITRAX/Outputs/HER_ITRAX_xrf1_cpsZ_comp.csv", row.names = FALSE)

# HER_ITRAX_Ln_Ti_Z_comp ----------------------------------------------------

# import processed composite site files xrf or xrf1
HER14L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER14L/Output/Composite/HER14L_xrf1_Ln_Ti_Z.csv") 
HER24L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER24L/Output/Composite/HER24L_xrf1_Ln_Ti_Z.csv") 
HER34L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER34L/Output/Composite/HER34L_xrf1_Ln_Ti_Z.csv")
HER42L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42L/Output/Composite/HER42L_xrf1_Ln_Ti_Z.csv")
HER42PB <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/Composite/HER42PB_xrf1_Ln_Ti_Z.csv")
HER44L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER44L/Output/Composite/HER44L_xrf1_Ln_Ti_Z.csv")
HER49L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER49L/Output/Composite/HER49L_xrf1_Ln_Ti_Z.csv")
# Assign elements from PeriodicTable package to 'elementsList' - to select all elements present across all sites
data(periodicTable)
elementsList <- c(periodicTable$symb, "Mo_inc", "Mo-coh")
# bind all data and rename to match itraxR column names - 
# all NA's in a column means that element failed filtering test at that site
HER_ITRAX_bind <-  bind_rows(HER14L, HER24L, HER34L, HER42L, HER42PB, HER44L, HER49L) %>% 
  rename(mean_age = SH20_age) %>% 
  select(Site, depth, mean_age, kcps:MSE, any_of(elementsList), Mo_inc: last_col())
# split HER14L into before and after 10 ka hiatus for plotting
HER14L_0ka <- HER_ITRAX_bind %>% 
  filter(Site=='HER14L' & mean_age <10000) %>% 
  mutate(Site = recode(Site, 'HER14L' = 'HER14L_0ka'))
HER14L_10ka <- HER_ITRAX_bind %>% 
  filter(Site=='HER14L' & mean_age >10000) %>% 
  mutate(Site = recode(Site, 'HER14L' = 'HER14L_10ka'))
HER_ITRAX_xrf1_Ln_Ti_Z_comp <- HER_ITRAX_bind %>% 
  filter(!(Site=='HER14L')) %>%
  bind_rows(HER14L_0ka, HER14L_10ka) %>% 
  arrange(Site)
HER_ITRAX_xrf1_Ln_Ti_Z_comp
tail(HER_ITRAX_xrf1_Ln_Ti_Z_comp)
# remove site records from Global Environment
rm(HER14L, HER14L_0ka, HER14L_10ka, HER24L, HER34L, HER42L, HER42PB, HER44L, HER49L, HER_ITRAX_bind)
# write to file
write.csv(HER_ITRAX_xrf1_Ln_Ti_Z_comp,"Papers/Hodgson_Hermite_2022/Data/ITRAX/Outputs/HER_ITRAX_xrf1_Ln_Ti_Z_comp.csv", row.names = FALSE)

