# -------------------------------------------------------------------------

# HERMITE: Site by site data import for MSCL facet plotting

# -------------------------------------------------------------------------

# Set up & clear ------------------------------------------------------------------

#clear previous console
remove (list = ls())
#clear plot window
dev.off()

# Load libraries & colours ----------------------------------------------------------

##load libraries
library(tidyverse)
library(tidypaleo)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(GGally)
library(ggExtra)
library(gridGraphics)
library(dplyr)
library(broom)
library(rioja)
library(itraxR)
library(PeriodicTable)
library(errors)
library(chemometrics)
library(psych) # for generating summary stats
library(compositions) # for generating clr 
#colour palettes
library(ggsci) #for npg etc
library(wesanderson) 
library(viridis)        
library(RColorBrewer)
library(scales)
#RColorBrewer display test
display.brewer.all()
display.brewer.pal(5,"BrBG")
display.brewer.all(colorblindFriendly = TRUE)

# Colours! keep minimized  --------
display.brewer.pal(11,"BrBG")
nb.cols <- 11
BrBG1 <- colorRampPalette(brewer.pal(11, "BrBG"))(nb.cols)
BrBG1
BrBG2 <- c("#A6611A", "#DFC27D", "#018571")
BrBG3 <- c("#8C510A", "#D8B365", "#F6E8C3", "#C7EAE5", "#5AB4AC", "#01665E")
BrBG4 <- c("#543005", "#8C510A", "#BF812D", "#DFC27D", "#F6E8C3", "#F5F5F5", "#C7EAE5", "#80CDC1", "#35978F", "#01665E", "#003C30", "black")

Paired1 <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
Paired1

Paired1 <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
Paired1

RdYlBu1 <-colorRampPalette(brewer.pal(11, "RdYlBu"))(nb.cols)
RdYlBu1
RdYlBu2 <- c("black","#313695","#4575B4","#74ADD1","#ABD9E9","#E0F3F8","#FFFFBF","#FEE090","#FDAE61","#F46D43","#D73027","#A50026")
RdYlBu3 <- c("#313695","#313695","#4575B4","#74ADD1","#ABD9E9","#C7EAE5", "#E0F3F8","#F5F5F5", "#F6E8C3", "#DFC27D", "#BF812D", "#8C510A", "#543005")  

Set2 <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)
Set2

Dark2.1 <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)
Dark2.1
Dark2.2 <- c("#7D54A5", "#A66753")
Dark2.3 <- c("#7D54A5", "#A66753", "#868686FF", "#A73030FF")

show_col(pal_jco("default")(10))
show_col(pal_jco("default")(10))
show_col(pal_npg("nrc")(10))
show_col(pal_jama("default")(7))

# -------------------------------------------------------------------------

# Site by site data import

# -------------------------------------------------------------------------

# set working directory
setwd("/Users/Steve/Dropbox/BAS/Data/R/")

# HER14L - SHW2015
HER14L_over <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HER14L_OVERLAPS.csv") %>% 
  mutate(SITE = 'HER14L', .before = SECTION_ID) %>% 
  select(SITE:MSCL_ID, Strat_depth:RES_SAT) %>% 
  filter(!if_all(PWAmp_SAT:RES_SAT, is.na))
HER14L_over

HER14L_comp <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HER14L_COMPOSITE.csv") %>% 
  mutate(SITE = 'HER14L', .before = SECTION_ID) %>% 
  select(SITE:MSCL_ID, Strat_depth:RES_SAT, WMAR_SAT) %>% 
  filter(!if_all(PWAmp_SAT:RES_SAT, is.na))
HER14L_comp

# HER24L - SHW2015
HER24L_over <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HER24L_OVERLAPS.csv") %>% 
  mutate(SITE = 'HER24L', .before = SECTION_ID) %>% 
  filter(!if_all(PWAmp_SAT:RES_SAT, is.na))
HER24L_over

HER24L_comp <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HER24L_COMPOSITE.csv") %>% 
  mutate(SITE = 'HER24L', .before = SECTION_ID) %>% 
  select(SITE:MSCL_ID, Strat_depth:RES_SAT, WMAR_SAT)
HER24L_comp

# HER34L - SHW2015
HER34L_over <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HER34L_OVERLAPS.csv") %>% 
  mutate(SITE = 'HER34L', .before = SECTION_ID) %>% 
  select(SITE:MSCL_ID, Strat_depth:RES_SAT) %>% 
  filter(!if_all(PWAmp_SAT:RES_SAT, is.na))
HER34L_over

HER34L_comp <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HER34L_COMPOSITE.csv") %>% 
  mutate(SITE = 'HER34L', .before = SECTION_ID) %>% 
  select(SITE:MSCL_ID, Strat_depth:RES_SAT, WMAR_SAT) %>% 
  filter(!if_all(PWAmp_SAT:RES_SAT, is.na))
HER34L_comp

# HER42L - 2015
HER42L_over <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HER42L_OVERLAPS.csv") %>% 
  mutate(SITE = 'HER42L', .before = SECTION_ID) %>% 
  select(SITE:MSCL_ID, Strat_depth:RES_SAT) %>% 
  filter(!if_all(PWAmp_SAT:RES_SAT, is.na))
HER42L_over

HER42L_comp <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HER42L_COMPOSITE.csv") %>% 
  mutate(SITE = 'HER42L', .before = SECTION_ID) %>% 
  select(SITE:MSCL_ID, Strat_depth:RES_SAT, WMAR_SAT) %>% 
  filter(!if_all(PWAmp_SAT:RES_SAT, is.na))
HER42L_comp

# HER42PB - ACE2018
HER42PB_over <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HER42PB_OVERLAPS.csv") %>% 
  mutate(SITE = 'HER42PB', .before = SECTION_ID) %>% 
  select(SITE:MSCL_ID, Strat_depth:RES_SAT) %>% 
  filter(!if_all(PWAmp_SAT:RES_SAT, is.na))
HER42PB_over

HER42PB_comp <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HER42PB_COMPOSITE.csv") %>% 
  mutate(SITE = 'HER42PB', .before = SECTION_ID) %>% 
  select(SITE:MSCL_ID, Strat_depth:RES_SAT, WMAR_SAT) %>% 
  filter(!if_all(PWAmp_SAT:RES_SAT, is.na))
HER42PB_comp

# HER44L - SHW2015 
HER44L_over <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HER44L_OVERLAPS.csv") %>% 
  mutate(SITE = 'HER44L', .before = SECTION_ID) %>% 
  select(SITE:MSCL_ID, Strat_depth:RES_SAT) %>% 
  filter(!if_all(PWAmp_SAT:RES_SAT, is.na))
HER44L_over

HER44L_comp <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HER44L_COMPOSITE.csv") %>% 
  mutate(SITE = 'HER44L', .before = SECTION_ID) %>% 
  select(SITE:MSCL_ID, Strat_depth:RES_SAT, WMAR_SAT) %>% 
  filter(!if_all(PWAmp_SAT:RES_SAT, is.na))
HER44L_comp

# HER49L -2015
HER49L_over <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HER49L_OVERLAPS.csv") %>% 
  mutate(SITE = 'HER49L', .before = SECTION_ID) %>% 
  select(SITE:MSCL_ID, Strat_depth:RES_SAT) %>% 
  filter(!if_all(PWAmp_SAT:RES_SAT, is.na))
HER49L_over

HER49L_comp <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HER49L_COMPOSITE.csv") %>% 
  mutate(SITE = 'HER49L', .before = SECTION_ID) %>% 
  select(SITE:MSCL_ID, Strat_depth:RES_SAT, WMAR_SAT) %>% 
  filter(!if_all(PWAmp_SAT:RES_SAT, is.na))
HER49L_comp

# HERPT1S1 - single section - no OVERLAP - ACE2018
HERPT1S1_comp <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HERPT1S1_COMPOSITE.csv") %>% 
  mutate(SITE = 'HERPT1S1', .before = SECTION_ID) %>% 
  select(SITE:MSCL_ID, Strat_depth:RES_SAT, WMAR_SAT) %>% 
  filter(!if_all(PWAmp_SAT:RES_SAT, is.na))
HERPT1S1_comp

# HERPT1S4 - ACE2018
HERPT1S4_over <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HERPT1S4_OVERLAPS.csv") %>% 
  mutate(SITE = 'HERPT1S4', .before = SECTION_ID) %>% 
  select(SITE:MSCL_ID, Strat_depth:RES_SAT) %>% 
  filter(!if_all(PWAmp_SAT:RES_SAT, is.na))
HERPT1S4_over

HERPT1S4_comp <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HERPT1S4_COMPOSITE.csv") %>% 
  mutate(SITE = 'HERPT1S4', .before = SECTION_ID) %>% 
  select(SITE:MSCL_ID, Strat_depth:RES_SAT) %>% 
  filter(!if_all(PWAmp_SAT:RES_SAT, is.na))
HERPT1S4_comp

# HERPT1S8 - ACE2018
HERPT1S8_over <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HERPT1S8_OVERLAPS.csv") %>% 
  mutate(SITE = 'HERPT1S8', .before = SECTION_ID) %>% 
  select(SITE:MSCL_ID, Strat_depth:RES_SAT)
HERPT1S8_over

HERPT1S8_comp <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HERPT1S8_COMPOSITE.csv") %>% 
  mutate(SITE = 'HERPT1S8', .before = SECTION_ID) %>% 
  select(SITE:MSCL_ID, Strat_depth:RES_SAT) %>% 
  filter(!if_all(PWAmp_SAT:RES_SAT, is.na))
HERPT1S8_comp

# HERPT2S1_SS1_SEC42 - ACE2018
HERPT2S1_over <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HERPT2S1_OVERLAPS.csv") %>% 
  mutate(SITE = 'HERPT2S1', .before = SECTION_ID) %>% 
  select(SITE:MSCL_ID, Strat_depth:RES_SAT) %>% 
  filter(!if_all(PWAmp_SAT:RES_SAT, is.na))
HERPT2S1_over

HERPT2S1_comp <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Inputs/HERPT2S1_COMPOSITE.csv") %>% 
  mutate(SITE = 'HERPT2S1', .before = SECTION_ID) %>% 
  select(SITE:MSCL_ID, Strat_depth:RES_SAT, WMAR_SAT) %>% 
  filter(!if_all(PWAmp_SAT:RES_SAT, is.na))
HERPT2S1_comp

# -------------------------------------------------------------------------

# LONG FORMAT DATA CONVERSION

# -------------------------------------------------------------------------

# for ggplot and tidypaleo facet plotting

# HER14L - long
# all parameters measured - overlaps
HER14L_over_long <- select(HER14L_over, SITE:MSCL_ID, Strat_depth:RES_SAT) %>%
  pivot_longer(c(PWVel_SAT,	Den1_SAT,	MS1_SAT,	DCMS1_SAT,	Imp_SAT,	FP_SAT,	RES_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HER14L_over_long

# all parameters measured - composite
HER14L_comp_long1 <- select(HER14L_comp, SITE:MSCL_ID,  Strat_depth:RES_SAT) %>%
  pivot_longer(c(PWVel_SAT,	Den1_SAT,	MS1_SAT,	DCMS1_SAT,	Imp_SAT,	FP_SAT,	RES_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HER14L_comp_long1

# for key plots in paper - composite 
HER14L_comp_long2 <- select(HER14L_comp, SITE:Accum_rate,  Den1_SAT, MS1_SAT, DCMS1_SAT, WMAR_SAT) %>%
  pivot_longer(c(Den1_SAT, MS1_SAT, DCMS1_SAT, WMAR_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HER14L_comp_long2

# for CONISS only - composite 
HER14L_comp_long3 <- select(HER14L_comp, SITE:Accum_rate,  Den1_SAT, DCMS1_SAT) %>%
  pivot_longer(c(Den1_SAT, DCMS1_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HER14L_comp_long3

# HER24L - long
HER24L_over_long <- select(HER24L_over, SITE:MSCL_ID, Strat_depth:RES_SAT) %>%
  pivot_longer(c(PWVel_SAT,	Den1_SAT,	MS1_SAT,	DCMS1_SAT,	Imp_SAT,	FP_SAT,	RES_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HER24L_over_long

HER24L_comp_long1 <- select(HER24L_comp, SITE:MSCL_ID,  Strat_depth:RES_SAT) %>%
  pivot_longer(c(PWVel_SAT,	Den1_SAT,	MS1_SAT,	DCMS1_SAT,	Imp_SAT,	FP_SAT,	RES_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HER24L_comp_long1

HER24L_comp_long2 <- select(HER24L_comp, SITE:Accum_rate,  Den1_SAT, MS1_SAT, DCMS1_SAT, WMAR_SAT) %>%
  pivot_longer(c(Den1_SAT, MS1_SAT, DCMS1_SAT, WMAR_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HER24L_comp_long2

# for CONISS only - composite 
HER24L_comp_long3 <- select(HER24L_comp, SITE:Accum_rate,  Den1_SAT, DCMS1_SAT) %>%
  pivot_longer(c(Den1_SAT, DCMS1_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HER24L_comp_long3

# HER34L - long
HER34L_over_long <- select(HER34L_over, SITE:MSCL_ID, Strat_depth:RES_SAT) %>%
  pivot_longer(c(PWVel_SAT,	Den1_SAT,	MS1_SAT,	DCMS1_SAT,	Imp_SAT,	FP_SAT,	RES_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HER34L_over_long

HER34L_comp_long1 <- select(HER34L_comp, SITE:MSCL_ID,  Strat_depth:RES_SAT) %>%
  pivot_longer(c(PWVel_SAT,	Den1_SAT,	MS1_SAT,	DCMS1_SAT,	Imp_SAT,	FP_SAT,	RES_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HER34L_comp_long1

HER34L_comp_long2 <- select(HER34L_comp, SITE:Accum_rate, Den1_SAT, MS1_SAT, DCMS1_SAT, WMAR_SAT) %>%
  pivot_longer(c(Den1_SAT, MS1_SAT, DCMS1_SAT, WMAR_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HER34L_comp_long2

# for CONISS only - composite 
HER34L_comp_long3 <- select(HER34L_comp, SITE:Accum_rate,  Den1_SAT, DCMS1_SAT) %>%
  pivot_longer(c(Den1_SAT, DCMS1_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HER34L_comp_long3

# HER42L - long
HER42L_over_long <- select(HER42L_over, SITE:MSCL_ID, Strat_depth:RES_SAT) %>%
  pivot_longer(c(PWVel_SAT,	Den1_SAT,	MS1_SAT,	DCMS1_SAT,	Imp_SAT,	FP_SAT,	RES_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HER42L_over_long

HER42L_comp_long1 <- select(HER42L_comp, SITE:MSCL_ID,  Strat_depth:RES_SAT) %>%
  pivot_longer(c(PWVel_SAT,	Den1_SAT,	MS1_SAT,	DCMS1_SAT,	Imp_SAT,	FP_SAT,	RES_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HER42L_comp_long1

HER42L_comp_long2 <- select(HER42L_comp, SITE:Accum_rate,  Den1_SAT, MS1_SAT, DCMS1_SAT, WMAR_SAT) %>%
  pivot_longer(c(Den1_SAT, MS1_SAT, DCMS1_SAT, WMAR_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HER42L_comp_long2

# for CONISS only - composite 
HER42L_comp_long3 <- select(HER42L_comp, SITE:Accum_rate,  Den1_SAT, DCMS1_SAT) %>%
  pivot_longer(c(Den1_SAT, DCMS1_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HER42L_comp_long3

# HER42PB - long
HER42PB_over_long <- select(HER42PB_over, SITE:MSCL_ID, Strat_depth:RES_SAT) %>%
  pivot_longer(c(PWVel_SAT,	Den1_SAT,	MS1_SAT,	DCMS1_SAT,	Imp_SAT,	FP_SAT,	RES_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HER42PB_over_long

HER42PB_comp_long1 <- select(HER42PB_comp, SITE:MSCL_ID,  Strat_depth:RES_SAT) %>%
  pivot_longer(c(PWVel_SAT,	Den1_SAT,	MS1_SAT,	DCMS1_SAT,	Imp_SAT,	FP_SAT,	RES_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HER42PB_comp_long1

HER42PB_comp_long2 <- select(HER42PB_comp, SITE:Accum_rate,  Den1_SAT, MS1_SAT, DCMS1_SAT, WMAR_SAT) %>%
  pivot_longer(c(Den1_SAT, MS1_SAT, DCMS1_SAT, WMAR_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HER42PB_comp_long2

# for CONISS only - composite 
HER42PB_comp_long3 <- select(HER42PB_comp, SITE:Accum_rate,  Den1_SAT, DCMS1_SAT) %>%
  pivot_longer(c(Den1_SAT, DCMS1_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HER42PB_comp_long3

# HER44L - long
HER44L_over_long <- select(HER44L_over, SITE:MSCL_ID, Strat_depth:RES_SAT) %>%
  pivot_longer(c(PWVel_SAT,	Den1_SAT,	MS1_SAT,	DCMS1_SAT,	Imp_SAT,	FP_SAT,	RES_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HER44L_over_long

HER44L_comp_long1 <- select(HER44L_comp, SITE:MSCL_ID,  Strat_depth:RES_SAT) %>%
  pivot_longer(c(PWVel_SAT,	Den1_SAT,	MS1_SAT,	DCMS1_SAT,	Imp_SAT,	FP_SAT,	RES_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HER44L_comp_long1

HER44L_comp_long2 <- select(HER44L_comp, SITE:Accum_rate,  Den1_SAT, MS1_SAT, DCMS1_SAT, WMAR_SAT) %>%
  pivot_longer(c(Den1_SAT, MS1_SAT, DCMS1_SAT, WMAR_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HER44L_comp_long2

# for CONISS only - composite 
HER44L_comp_long3 <- select(HER44L_comp, SITE:Accum_rate,  Den1_SAT, DCMS1_SAT) %>%
  pivot_longer(c(Den1_SAT, DCMS1_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HER44L_comp_long3

# HER49L - long
HER49L_over_long <- select(HER49L_over, SITE:MSCL_ID, Strat_depth:RES_SAT) %>%
  pivot_longer(c(PWVel_SAT,	Den1_SAT,	MS1_SAT,	DCMS1_SAT,	Imp_SAT,	FP_SAT,	RES_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HER49L_over_long

HER49L_comp_long1 <- select(HER49L_comp, SITE:MSCL_ID,  Strat_depth:RES_SAT) %>%
  pivot_longer(c(PWVel_SAT,	Den1_SAT,	MS1_SAT,	DCMS1_SAT,	Imp_SAT,	FP_SAT,	RES_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HER49L_comp_long1

HER49L_comp_long2 <- select(HER49L_comp, SITE:Accum_rate,  Den1_SAT, MS1_SAT, DCMS1_SAT, WMAR_SAT) %>%
  pivot_longer(c(Den1_SAT, MS1_SAT, DCMS1_SAT, WMAR_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HER49L_comp_long2

# for CONISS only - composite 
HER49L_comp_long3 <- select(HER49L_comp, SITE:Accum_rate,  Den1_SAT, DCMS1_SAT) %>%
  pivot_longer(c(Den1_SAT, DCMS1_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HER49L_comp_long3

# HERPT1S1 - long
HERPT1S1_comp_long1 <- select(HERPT1S1_comp, SITE:MSCL_ID, Strat_depth:RES_SAT) %>%
  pivot_longer(c(PWVel_SAT,	Den1_SAT,	MS1_SAT,	DCMS1_SAT,	Imp_SAT,	FP_SAT,	RES_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HERPT1S1_comp_long1

HERPT1S1_comp_long2 <- select(HERPT1S1_comp, SITE:Accum_rate,  Den1_SAT, MS1_SAT, DCMS1_SAT, WMAR_SAT) %>%
  pivot_longer(c(Den1_SAT, MS1_SAT, DCMS1_SAT, WMAR_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
tail(HERPT1S1_comp_long2)

# for CONISS only - composite 
HERPT1S1_comp_long3 <- select(HERPT1S1_comp, SITE:Accum_rate,  Den1_SAT, DCMS1_SAT) %>%
  pivot_longer(c(Den1_SAT, DCMS1_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HERPT1S1_comp_long3

# HERPT1S4 - long
HERPT1S4_over_long <- select(HERPT1S4_over, SITE:MSCL_ID, Strat_depth:RES_SAT) %>%
  pivot_longer(c(PWVel_SAT,	Den1_SAT,	MS1_SAT,	DCMS1_SAT,	Imp_SAT,	FP_SAT,	RES_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HERPT1S4_over_long

HERPT1S4_comp_long1 <- select(HERPT1S4_comp, SITE:MSCL_ID, Strat_depth:RES_SAT) %>%
  pivot_longer(c(PWVel_SAT,	Den1_SAT,	MS1_SAT,	DCMS1_SAT,	Imp_SAT,	FP_SAT,	RES_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HERPT1S4_comp_long1

HERPT1S4_comp_long2 <- select(HERPT1S4_comp, SITE:Strat_depth,  Den1_SAT, MS1_SAT, DCMS1_SAT) %>%
  pivot_longer(c(Den1_SAT, MS1_SAT, DCMS1_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HERPT1S4_comp_long2

# for CONISS only - composite 
HERPT1S4_comp_long3 <- select(HERPT1S4_comp, SITE:Strat_depth,  Den1_SAT, DCMS1_SAT) %>%
  pivot_longer(c(Den1_SAT, DCMS1_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HERPT1S4_comp_long3

# HERPT1S8 - long
HERPT1S8_over_long <- select(HERPT1S8_over, SITE:MSCL_ID, Strat_depth:RES_SAT) %>%
  pivot_longer(c(PWVel_SAT,	Den1_SAT,	MS1_SAT,	DCMS1_SAT,	Imp_SAT,	FP_SAT,	RES_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HERPT1S8_over_long

HERPT1S8_comp_long1 <- select(HERPT1S8_comp, SITE:MSCL_ID, Strat_depth:RES_SAT) %>%
  pivot_longer(c(PWVel_SAT,	Den1_SAT,	MS1_SAT,	DCMS1_SAT,	Imp_SAT,	FP_SAT,	RES_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HERPT1S8_comp_long1

HERPT1S8_comp_long2 <- select(HERPT1S8_comp, SITE:Strat_depth,  Den1_SAT, MS1_SAT, DCMS1_SAT) %>%
  pivot_longer(c(Den1_SAT, MS1_SAT, DCMS1_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HERPT1S8_comp_long2

# for CONISS only - composite 
HERPT1S8_comp_long3 <- select(HERPT1S8_comp, SITE:Strat_depth,  Den1_SAT, DCMS1_SAT) %>%
  pivot_longer(c(Den1_SAT, DCMS1_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HERPT1S8_comp_long3

# HERPT2S1_SS2_XRF6_SEC42 - long
HERPT2S1_over_long <- select(HERPT2S1_over, SITE:MSCL_ID, Strat_depth:RES_SAT) %>%
  pivot_longer(c(PWVel_SAT,	Den1_SAT,	MS1_SAT,	DCMS1_SAT,	Imp_SAT,	FP_SAT,	RES_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HERPT2S1_over_long

HERPT2S1_comp_long1 <- select(HERPT2S1_comp, SITE:MSCL_ID, Strat_depth:RES_SAT) %>%
  pivot_longer(c(PWVel_SAT,	Den1_SAT,	MS1_SAT,	DCMS1_SAT,	Imp_SAT,	FP_SAT,	RES_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HERPT2S1_comp_long1

HERPT2S1_comp_long2 <- select(HERPT2S1_comp, SITE:Accum_rate,  Den1_SAT, MS1_SAT, DCMS1_SAT, WMAR_SAT) %>%
  pivot_longer(c(Den1_SAT, MS1_SAT, DCMS1_SAT, WMAR_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HERPT2S1_comp_long2

# for CONISS only - composite 
HERPT2S1_comp_long3 <- select(HERPT2S1_comp, SITE:Accum_rate,  Den1_SAT, DCMS1_SAT) %>%
  pivot_longer(c(Den1_SAT, DCMS1_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HERPT2S1_comp_long3


# -------------------------------------------------------------------------

# PLOTTING   ----------------------------------------------

# set theme for Tidypalaeo - 7 pt to fit - scales up when saving to A4 pdf
theme_set(theme_paleo(base_size=7) + theme(
  plot.title = element_text(color="black", size=7, face="bold"),
  axis.title.x = element_text(color="black", size=7),
  axis.title.y = element_text(color="black", size=7),
  axis.text.x = element_text(color="black", size = 7),
  axis.text.y = element_text(color="black", size = 7)
))

## HER14L Plots ----------------------------------------------------------

# Overlap and Composite plots for XRAD mapping ---------------------

# All data - showing overlaps
HER14L_p1 <- HER14L_over_long %>% 
  mutate(param = fct_relevel(param, "PWVel_SAT", "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"Imp_SAT",	"FP_SAT",	"RES_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "Imp_SAT" = "Imp.", "FP_SAT" = "FP", "RES_SAT" = "Res.")) %>% 
  ggplot(aes(x = value, y = Strat_depth, colour = SECTION_ID, 
             fill = SECTION_ID, shape = SECTION_ID), size = 0.3) +
  geom_lineh(lwd = 0.3) +
  #geom_point(size = 0.2) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HER14L [SAT]: Overlaps") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("PWave Vel." = "ms-1", "GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "Res." = "Ohm.m"))
HER14L_p1
ggsave("Outputs/Figures/HER14L/Fig1.1_ALL_over.pdf",
       height = c(10), width = c(20), dpi = 600, units = "cm")

# All data - composite colour by core for XRAD single core mapping
HER14L_p2 <- HER14L_comp_long1 %>% 
  mutate(param = fct_relevel(param, "PWVel_SAT", "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"Imp_SAT",	"FP_SAT",	"RES_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "Imp_SAT" = "Imp.", "FP_SAT" = "FP", "RES_SAT" = "Res.")) %>% 
  ggplot(aes(x = value, y = Strat_depth, colour = SECTION_ID, 
             fill = SECTION_ID, shape = SECTION_ID), size = 0.3) +
  geom_lineh(lwd = 0.3) +
  #geom_point(size = 0.2) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HER14L [SAT]: Composite") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("PWave Vel." = "ms-1", "GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "10^-8 m3 kg-1", "Res." = "Ohm.m"))
HER14L_p2
ggsave("Outputs/Figures/HER14L/Fig1.2_ALL_comp.pdf",
       height = c(10), width = c(20), dpi = 600, units = "cm")

# All data - composite plotted as black
HER14L_p3 <- HER14L_comp_long1 %>% 
  mutate(param = fct_relevel(param, "PWVel_SAT", "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"Imp_SAT",	"FP_SAT",	"RES_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "Imp_SAT" = "Imp.", "FP_SAT" = "FP", "RES_SAT" = "Res.")) %>% 
  ggplot(aes(x = value, y = Strat_depth), size = 0.3) +
  geom_lineh(lwd = 0.3) +
  #geom_point(size = 0.2) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HER14L [SAT]: Composite") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("PWave Vel." = "ms-1", "GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "Res." = "Ohm.m"))
HER14L_p3
ggsave("Outputs/Figures/HER14L/Fig1.3_ALL_comp_black.pdf",
       height = c(10), width = c(16), dpi = 600, units = "cm")

# plot on A4
ggarrange(HER14L_p1, HER14L_p2, nrow = 2)
ggsave("Outputs/Figures/HER14L/Fig1.4_over_comp_coresx2.pdf",
       height = c(20), width = c(20), dpi = 600, units = "cm")

ggarrange(HER14L_p1, HER14L_p2, HER14L_p3, nrow = 3)
ggsave("Outputs/Figures/HER14L/Fig1.5_over_comp_coresx3.pdf",
       height = c(30), width = c(20), dpi = 600, units = "cm")

# Plot the age depth model to check ---------------------------------------

HER14L_adm_plot <- age_depth_model(
  HER14L_comp,
  depth = Strat_depth, age = SH20_mean_age,
  age_max = SH20_max_age_95CI,
  age_min = SH20_min_age_95CI
) 

plot(HER14L_adm_plot)

# save plot to file 
pdf(file = "Outputs/Figures/HER14L/Fig1.6_age_depth_model.pdf",
    width = 6, # The width of the plot in inches
    height = 6) # The height of the plot in inches
HER14L_adm_plot <- age_depth_model(
  HER14L_comp,
  depth = Strat_depth, age = SH20_mean_age,
  age_max = SH20_max_age_95CI,
  age_min = SH20_min_age_95CI
) 
plot(HER14L_adm_plot)
dev.off()

# Composite plots for paper ------------------------------------------

HER14L_plot <- HER14L_comp_long2 %>% 
  mutate(param = fct_relevel(param, "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"WMAR_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "WMAR_SAT" = "WMAR")) %>%
  ggplot(aes(x = value, y = Strat_depth), size = 0.01) +
  geom_lineh(lwd = 0.3) +
  geom_point(size = 0.01) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HER14L [SAT]") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "WMAR" = "g cm-2  a-1"))
HER14L_plot

# OPTIONAL
# add hiatus (red) horizontal line and other points of interest (green) - or comment out

HER14L_plot1 <- HER14L_plot

# OR

HER14L_plot1 <- HER14L_plot +
  geom_hline(yintercept = c(76), col = "red", lty = 1, alpha = 1, lwd = 0.5)
HER14L_plot1

# Add age on secondary age y axis
HER14L_adm <- age_depth_model(
  HER14L_comp, 
  depth = Strat_depth, age = SH20_mean_age)

HER14L_plot_age <- HER14L_plot1 + 
  scale_y_depth_age(HER14L_adm, age_name = "Age (cal a BP)")
HER14L_plot_age

# CONISS
# Uses Den1_SAT, DCMS_SAT only in _com_long3, but plotted with MS1_SAT and WMAR
# CONISS using nested_chclust_coniss() and the rioja package for constrained cluster analysis 
# By default, this calculates the number of significant zones based on a broken stick simulation
# TO DO <- add number of groups and their depth/ages

HER14L_coniss <- HER14L_comp_long3 %>%
  nested_data(qualifiers = c(SH20_mean_age, Strat_depth), 
              key = param, value = value, trans = scale) %>%
  nested_chclust_coniss()

HER14L_plot_coniss <- HER14L_plot_age +
  layer_dendrogram(HER14L_coniss, aes(y = Strat_depth), 
                   param = "CONISS", size = 0.2, lwd = 0.2) +
  layer_zone_boundaries(HER14L_coniss, aes(y = Strat_depth), 
                        size = 0.3, lwd = 0.3, colour = "red")

ggarrange(HER14L_plot_age, HER14L_plot_coniss, nrow = 2)
ggsave("Outputs/Figures/HER14L/Fig2.1_comp_depth_coniss.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")

# Plot versus age
HER14L_age <- HER14L_comp_long2 %>% 
  mutate(param = fct_relevel(param, "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"WMAR_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "WMAR_SAT" = "WMAR")) %>%
  ggplot(aes(x = SH20_mean_age, y = value), size = 0.01) +
  geom_line(lwd = 0.3) +
  geom_point(size = 0.01) +
  #scale_y_reverse() +
  # add units to the graph headers
  facet_geochem_grid(
    vars(param),
    units = c("GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "WMAR" = "g cm-2  a-1")) +
  #facet_geochem_grid(vars(param)) +
  labs(x = "Age (cal a BP)", y = NULL)
HER14L_age

# OPTIONAL
# shading and lines for HER14 - added hiatus as a shaded area - min-max = rounded to the nearest 10 years - not needed for other plots
zone_data <- tibble(ymin = -Inf, ymax = Inf , xmin = 1800, xmax = 10790)

HER14L_age_lines <- HER14L_age +
  geom_rect(
    mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
    data = zone_data, 
    #alpha = 0.2,
    fill = "white",
    inherit.aes = FALSE
  ) +
  # add a red dotted lines of interest at 1 ka, 2 ka, 11.75 ka and 13.8 ka - need to label these later on
  geom_vline(xintercept = c(1800, 10790), col = "red", lty = 1, alpha = 1, lwd = 0.5) +
  geom_vline(xintercept = c(0, 1000, 2000, 3000, 10000, 13800), col = "blue", lty = 3, alpha = 0.7, lwd = 0.3) +
  geom_vline(xintercept = c(4200, 8200, 11750, 16000), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  labs(title = "HER14L [SAT]")
HER14L_age_lines

# Save final figure
ggarrange(HER14L_plot_coniss, HER14L_age_lines, nrow = 2)
ggsave("Outputs/Figures/HER14L/Fig2.2_depth_age_summary.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")

# -------------------------------------------------------------------------

# HER24L plots ------------------------------------------------------------------
# Overlap and Composite plots for XRAD mapping ---------------------

# All data - showing overlaps
HER24L_p1 <- HER24L_over_long %>% 
  mutate(param = fct_relevel(param, "PWVel_SAT", "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"Imp_SAT",	"FP_SAT",	"RES_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "Imp_SAT" = "Imp.", "FP_SAT" = "FP", "RES_SAT" = "Res.")) %>% 
  ggplot(aes(x = value, y = Strat_depth, colour = SECTION_ID, 
             fill = SECTION_ID, shape = SECTION_ID), size = 0.3) +
  geom_lineh(lwd = 0.3) +
  #geom_point(size = 0.2) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HER24L [SAT]: Overlaps") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("PWave Vel." = "ms-1", "GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "Res." = "Ohm.m"))
HER24L_p1
ggsave("Outputs/Figures/HER24L/Fig1.1_ALL_over.pdf",
       height = c(10), width = c(20), dpi = 600, units = "cm")

# All data - composite colour by core for XRAD single core mapping
HER24L_p2 <- HER24L_comp_long1 %>% 
  mutate(param = fct_relevel(param, "PWVel_SAT", "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"Imp_SAT",	"FP_SAT",	"RES_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "Imp_SAT" = "Imp.", "FP_SAT" = "FP", "RES_SAT" = "Res.")) %>% 
  ggplot(aes(x = value, y = Strat_depth, colour = SECTION_ID, 
             fill = SECTION_ID, shape = SECTION_ID), size = 0.3) +
  geom_lineh(lwd = 0.3) +
  #geom_point(size = 0.2) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HER24L [SAT]: Composite") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("PWave Vel." = "ms-1", "GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "10^-8 m3 kg-1", "Res." = "Ohm.m"))
HER24L_p2
ggsave("Outputs/Figures/HER24L/Fig1.2_ALL_comp.pdf",
       height = c(10), width = c(20), dpi = 600, units = "cm")

# All data - composite plotted as black
HER24L_p3 <- HER24L_comp_long1 %>% 
  mutate(param = fct_relevel(param, "PWVel_SAT", "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"Imp_SAT",	"FP_SAT",	"RES_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "Imp_SAT" = "Imp.", "FP_SAT" = "FP", "RES_SAT" = "Res.")) %>% 
  ggplot(aes(x = value, y = Strat_depth), size = 0.3) +
  geom_lineh(lwd = 0.3) +
  #geom_point(size = 0.2) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HER24L [SAT]: Composite") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("PWave Vel." = "ms-1", "GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "Res." = "Ohm.m"))
HER24L_p3
ggsave("Outputs/Figures/HER24L/Fig1.3_ALL_comp_black.pdf",
       height = c(10), width = c(16), dpi = 600, units = "cm")

# Multi-plots per page
ggarrange(HER24L_p1, HER24L_p2, nrow = 2)
ggsave("Outputs/Figures/HER24L/Fig1.4_over_comp_coresx2.pdf",
       height = c(20), width = c(20), dpi = 600, units = "cm")

ggarrange(HER24L_p1, HER24L_p2, HER24L_p3, nrow = 3)
ggsave("Outputs/Figures/HER24L/Fig1.5_over_comp_coresx3.pdf",
       height = c(30), width = c(20), dpi = 600, units = "cm")

# Age depth model - plot to check ---------------------------------------

#TO DO - add adm plot to core plots
HER24L_adm_plot <- age_depth_model(
  HER24L_comp,
  depth = Strat_depth, age = SH20_mean_age,
  age_max = SH20_max_age_95CI,
  age_min = SH20_min_age_95CI
) 

plot(HER24L_adm_plot)

# save plot to file 
pdf(file = "Outputs/Figures/HER24L/Fig1.6_age_depth_model.pdf",
    width = 6, # The width of the plot in inches
    height = 6) # The height of the plot in inches
HER24L_adm_plot <- age_depth_model(
  HER24L_comp,
  depth = Strat_depth, age = SH20_mean_age,
  age_max = SH20_max_age_95CI,
  age_min = SH20_min_age_95CI
) 
plot(HER24L_adm_plot)
dev.off()

# Composite plots for paper ------------------------------------------

HER24L_plot <- HER24L_comp_long2 %>% 
  mutate(param = fct_relevel(param, "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"WMAR_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "WMAR_SAT" = "WMAR")) %>%
  ggplot(aes(x = value, y = Strat_depth), size = 0.01) +
  geom_lineh(lwd = 0.3) +
  geom_point(size = 0.01) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HER24L [SAT]") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "WMAR" = "g cm-2  a-1"))
HER24L_plot

# OPTIONAL
# add hiatus (red) horizontal line and other points of interest (green) - or comment out

HER24L_plot1 <- HER24L_plot

# OR

#HER24L_plot1 <- HER24L_plot +
#geom_hline(yintercept = c(76), col = "red", lty = 1, alpha = 1, lwd = 0.5)
#HER24L_plot1

# Add age on secondary age y axis
HER24L_adm <- age_depth_model(
  HER24L_comp, 
  depth = Strat_depth, age = SH20_mean_age)

HER24L_plot_age <- HER24L_plot1 + 
  scale_y_depth_age(HER24L_adm, age_name = "Age (cal a BP)")
HER24L_plot_age

# CONISS
# Uses Den1_SAT, DCMS_SAT only in _com_long3, but plotted with MS1_SAT and WMAR
# CONISS using nested_chclust_coniss() and the rioja package for constrained cluster analysis 
# By default, this calculates the number of significant zones based on a broken stick simulation
# TO DO <- add number of groups and their depth/ages

HER24L_coniss <- HER24L_comp_long3 %>%
  nested_data(qualifiers = c(SH20_mean_age, Strat_depth), 
              key = param, value = value, trans = scale) %>%
  nested_chclust_coniss()

HER24L_plot_coniss <- HER24L_plot_age +
  layer_dendrogram(HER24L_coniss, aes(y = Strat_depth), 
                   param = "CONISS", size = 0.2, lwd = 0.2) +
  layer_zone_boundaries(HER24L_coniss, aes(y = Strat_depth), 
                        size = 0.3, lwd = 0.3, colour = "red")

ggarrange(HER24L_plot_age, HER24L_plot_coniss, nrow = 2)
ggsave("Outputs/Figures/HER24L/Fig2.1_comp_depth_coniss.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")

# Plot versus age
HER24L_age <- HER24L_comp_long2 %>% 
  mutate(param = fct_relevel(param, "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"WMAR_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "WMAR_SAT" = "WMAR")) %>%
  ggplot(aes(x = SH20_mean_age, y = value), size = 0.01) +
  geom_line(lwd = 0.3) +
  geom_point(size = 0.01) +
  #scale_y_reverse() +
  # add units to the graph headers
  facet_geochem_grid(
    vars(param),
    units = c("GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "WMAR" = "g cm-2  a-1")) +
  #facet_geochem_grid(vars(param)) +
  labs(x = "Age (cal a BP)", y = NULL)
HER24L_age

# OPTIONAL
# shading and lines for HER14 - added hiatus as a shaded area - min-max = rounded to the nearest 10 years - not needed for other plots
#zone_data <- tibble(ymin = -Inf, ymax = Inf , xmin = 1800, xmax = 10790)

HER24L_age_lines1 <- HER24L_age +
  #geom_rect(
  #mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
  #data = zone_data, 
  #alpha = 0.2,
  #fill = "white",
  #inherit.aes = FALSE
  #) +
  # add lines of interest at 1 ka, 2 ka, 11.75 ka and 13.8 ka etc. - need to label these later on
  #geom_vline(xintercept = c(0, 1000, 2000, 3000, 7500, 10000, 13800), col = "blue", lty = 3, alpha = 0.7, lwd = 0.3) +
  geom_vline(xintercept = c(4200), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  #geom_vline(xintercept = c(4200), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  labs(title = "HER24L [SAT]")
HER24L_age_lines1

# Save final figure
ggarrange(HER24L_plot_coniss, HER24L_age_lines1, nrow = 2)
ggsave("Outputs/Figures/HER24L/Fig2.2_depth_age_summary.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")

HER24L_age_lines2 <- HER24L_age +
  #geom_rect(
  #mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
  #data = zone_data, 
  #alpha = 0.2,
  #fill = "white",
  #inherit.aes = FALSE
  #) +
  # add lines of interest at 1 ka, 2 ka, 11.75 ka and 13.8 ka etc. - need to label these later on
  geom_vline(xintercept = c(0, 1000, 2000, 3000, 7500, 10000, 13800), col = "blue", lty = 3, alpha = 0.7, lwd = 0.3) +
  geom_vline(xintercept = c(4200, 8200, 11750, 16000), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  #geom_vline(xintercept = c(4200), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  labs(title = "HER24L [SAT]")
HER24L_age_lines2

# Save final figure
ggarrange(HER24L_plot_coniss, HER24L_age_lines2, nrow = 2)
ggsave("Outputs/Figures/HER24L/Fig2.3_depth_age_summary.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")

# -------------------------------------------------------------------------

# HERPT1S1 plots ----------------------------------------------

#  Composite plots for XRAD mapping ---------------------

# All data - composite colour by core for XRAD single core mapping
HERPT1S1_p2 <- HERPT1S1_comp_long1 %>% 
  mutate(param = fct_relevel(param, "PWVel_SAT", "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"Imp_SAT",	"FP_SAT",	"RES_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "Imp_SAT" = "Imp.", "FP_SAT" = "FP", "RES_SAT" = "Res.")) %>% 
  ggplot(aes(x = value, y = Strat_depth, colour = SECTION_ID, 
             fill = SECTION_ID, shape = SECTION_ID), size = 0.3) +
  geom_lineh(lwd = 0.3) +
  #geom_point(size = 0.2) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HERPT1S1 [SAT]: Composite") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("PWave Vel." = "ms-1", "GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "10^-8 m3 kg-1", "Res." = "Ohm.m"))
HERPT1S1_p2
ggsave("Outputs/Figures/HERPT1S1/Fig1.2_ALL_comp.pdf",
       height = c(10), width = c(20), dpi = 600, units = "cm")

# All data - composite plotted as black
HERPT1S1_p3 <- HERPT1S1_comp_long1 %>% 
  mutate(param = fct_relevel(param, "PWVel_SAT", "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"Imp_SAT",	"FP_SAT",	"RES_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "Imp_SAT" = "Imp.", "FP_SAT" = "FP", "RES_SAT" = "Res.")) %>% 
  ggplot(aes(x = value, y = Strat_depth), size = 0.3) +
  geom_lineh(lwd = 0.3) +
  #geom_point(size = 0.2) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HERPT1S1 [SAT]: Composite") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("PWave Vel." = "ms-1", "GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "Res." = "Ohm.m"))
HERPT1S1_p3
ggsave("Outputs/Figures/HERPT1S1/Fig1.3_ALL_comp_black.pdf",
       height = c(10), width = c(16), dpi = 600, units = "cm")

# Multi-plots per page
ggarrange(HERPT1S1_p2, HERPT1S1_p2, nrow = 2)
ggsave("Outputs/Figures/HERPT1S1/Fig1.4_over_comp_coresx2.pdf",
       height = c(20), width = c(20), dpi = 600, units = "cm")

ggarrange(HERPT1S1_p2, HERPT1S1_p2, HERPT1S1_p3, nrow = 3)
ggsave("Outputs/Figures/HERPT1S1/Fig1.5_over_comp_coresx3.pdf",
       height = c(30), width = c(20), dpi = 600, units = "cm")

# Age depth model ---------------------------------------

#TO DO - add adm plot to core plots
HERPT1S1_adm_plot <- age_depth_model(
  HERPT1S1_comp,
  depth = Strat_depth, age = SH20_mean_age,
  age_max = SH20_max_age_95CI,
  age_min = SH20_min_age_95CI
) 

plot(HERPT1S1_adm_plot)

# save plot to file 
pdf(file = "Outputs/Figures/HERPT1S1/Fig1.6_age_depth_model.pdf",
    width = 6, # The width of the plot in inches
    height = 6) # The height of the plot in inches
HERPT1S1_adm_plot <- age_depth_model(
  HERPT1S1_comp,
  depth = Strat_depth, age = SH20_mean_age,
  age_max = SH20_max_age_95CI,
  age_min = SH20_min_age_95CI
) 
plot(HERPT1S1_adm_plot)
dev.off()

# Composite plots for paper ------------------------------------------

HERPT1S1_plot <- HERPT1S1_comp_long2 %>% 
  mutate(param = fct_relevel(param, "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"WMAR_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "WMAR_SAT" = "WMAR")) %>%
  ggplot(aes(x = value, y = Strat_depth), size = 0.01) +
  geom_lineh(lwd = 0.3) +
  geom_point(size = 0.01) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HERPT1S1 [SAT]") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "WMAR" = "g cm-2  a-1"))
HERPT1S1_plot

# OPTIONAL
# add hiatus (red) horizontal line and other points of interest (green) - or comment out

HERPT1S1_plot1 <- HERPT1S1_plot

# OR

#HERPT1S1_plot1 <- HERPT1S1_plot +
#geom_hline(yintercept = c(76), col = "red", lty = 1, alpha = 1, lwd = 0.5)
#HERPT1S1_plot1

# Add age on secondary age y axis
HERPT1S1_adm <- age_depth_model(
  HERPT1S1_comp, 
  depth = Strat_depth, age = SH20_mean_age)

HERPT1S1_plot_age <- HERPT1S1_plot1 + 
  scale_y_depth_age(HERPT1S1_adm, age_name = "Age (cal a BP)")
HERPT1S1_plot_age

# CONISS
# Uses Den1_SAT, DCMS_SAT only in _com_long3, but plotted with MS1_SAT and WMAR
# CONISS using nested_chclust_coniss() and the rioja package for constrained cluster analysis 
# By default, this calculates the number of significant zones based on a broken stick simulation
# TO DO <- add number of groups and their depth/ages

HERPT1S1_coniss <- HERPT1S1_comp_long3 %>%
  nested_data(qualifiers = c(SH20_mean_age, Strat_depth), 
              key = param, value = value, trans = scale) %>%
  nested_chclust_coniss()

HERPT1S1_plot_coniss <- HERPT1S1_plot_age +
  layer_dendrogram(HERPT1S1_coniss, aes(y = Strat_depth), 
                   param = "CONISS", size = 0.2, lwd = 0.2) +
  layer_zone_boundaries(HERPT1S1_coniss, aes(y = Strat_depth), 
                        size = 0.3, lwd = 0.3, colour = "red")
HERPT1S1_plot_coniss

ggarrange(HERPT1S1_plot_age, HERPT1S1_plot_coniss, nrow = 2)
ggsave("Outputs/Figures/HERPT1S1/Fig2.1_comp_depth_coniss.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")

# Plot versus age
HERPT1S1_age <- HERPT1S1_comp_long2 %>% 
  mutate(param = fct_relevel(param, "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"WMAR_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "WMAR_SAT" = "WMAR")) %>%
  ggplot(aes(x = SH20_mean_age, y = value), size = 0.01) +
  geom_line(lwd = 0.3) +
  geom_point(size = 0.01) +
  #scale_y_reverse() +
  # add units to the graph headers
  facet_geochem_grid(
    vars(param),
    units = c("GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "WMAR" = "g cm-2  a-1")) +
  #facet_geochem_grid(vars(param)) +
  labs(x = "Age (cal a BP)", y = NULL)
HERPT1S1_age

# OPTIONAL
# shading and lines for HER14 - added hiatus as a shaded area - min-max = rounded to the nearest 10 years - not needed for other plots
#zone_data <- tibble(ymin = -Inf, ymax = Inf , xmin = 1800, xmax = 10790)

HERPT1S1_age_lines1 <- HERPT1S1_age +
  #geom_rect(
  #mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
  #data = zone_data, 
  #alpha = 0.2,
  #fill = "white",
  #inherit.aes = FALSE
  #) +
  # add lines of interest at 1 ka, 2 ka, 11.75 ka and 13.8 ka etc. - need to label these later on
  #geom_vline(xintercept = c(0, 1000, 2000, 3000, 7500, 10000, 13800), col = "blue", lty = 3, alpha = 0.7, lwd = 0.3) +
  geom_vline(xintercept = c(4200, 8200), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  #geom_vline(xintercept = c(4200), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  labs(title = "HERPT1S1 [SAT]")
HERPT1S1_age_lines1

# Save final figure
ggarrange(HERPT1S1_plot_coniss, HERPT1S1_age_lines1, nrow = 2)
ggsave("Outputs/Figures/HERPT1S1/Fig2.2_depth_age_summary.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")

HERPT1S1_age_lines2 <- HERPT1S1_age +
  #geom_rect(
  #mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
  #data = zone_data, 
  #alpha = 0.2,
  #fill = "white",
  #inherit.aes = FALSE
  #) +
  # add lines of interest at 1 ka, 2 ka, 11.75 ka and 13.8 ka etc. - need to label these later on
  geom_vline(xintercept = c(0, 1000, 2000, 3000, 7500, 10000, 13800), col = "blue", lty = 3, alpha = 0.7, lwd = 0.3) +
  geom_vline(xintercept = c(4200, 8200, 11750, 16000), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  #geom_vline(xintercept = c(4200), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  labs(title = "HERPT1S1 [SAT]")
HERPT1S1_age_lines2

# Save final figure
ggarrange(HERPT1S1_plot_coniss, HERPT1S1_age_lines2, nrow = 2)
ggsave("Outputs/Figures/HER42L/Fig2.3_depth_age_summary.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")

# -------------------------------------------------------------------------

# HERPT1S4 ----------------------------------------------

# Overlap and Composite plots for XRAD mapping ---------------------

# All data - showing overlaps
HERPT1S4_p1 <- HERPT1S4_over_long %>% 
  mutate(param = fct_relevel(param, "PWVel_SAT", "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"Imp_SAT",	"FP_SAT",	"RES_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "Imp_SAT" = "Imp.", "FP_SAT" = "FP", "RES_SAT" = "Res.")) %>% 
  ggplot(aes(x = value, y = Strat_depth, colour = SECTION_ID, 
             fill = SECTION_ID, shape = SECTION_ID), size = 0.3) +
  geom_lineh(lwd = 0.3) +
  #geom_point(size = 0.2) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HERPT1S4 [SAT]: Overlaps") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("PWave Vel." = "ms-1", "GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "Res." = "Ohm.m"))
HERPT1S4_p1
ggsave("Outputs/Figures/HERPT1S4/Fig1.1_ALL_over.pdf",
       height = c(10), width = c(20), dpi = 600, units = "cm")

# All data - composite colour by core for XRAD single core mapping
HERPT1S4_p2 <- HERPT1S4_comp_long1 %>% 
  mutate(param = fct_relevel(param, "PWVel_SAT", "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"Imp_SAT",	"FP_SAT",	"RES_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "Imp_SAT" = "Imp.", "FP_SAT" = "FP", "RES_SAT" = "Res.")) %>% 
  ggplot(aes(x = value, y = Strat_depth, colour = SECTION_ID, 
             fill = SECTION_ID, shape = SECTION_ID), size = 0.3) +
  geom_lineh(lwd = 0.3) +
  #geom_point(size = 0.2) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HERPT1S4 [SAT]: Composite") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("PWave Vel." = "ms-1", "GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "10^-8 m3 kg-1", "Res." = "Ohm.m"))
HERPT1S4_p2
ggsave("Outputs/Figures/HERPT1S4/Fig1.2_ALL_comp.pdf",
       height = c(10), width = c(20), dpi = 600, units = "cm")

# All data - composite plotted as black
HERPT1S4_p3 <- HERPT1S4_comp_long1 %>% 
  mutate(param = fct_relevel(param, "PWVel_SAT", "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"Imp_SAT",	"FP_SAT",	"RES_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "Imp_SAT" = "Imp.", "FP_SAT" = "FP", "RES_SAT" = "Res.")) %>% 
  ggplot(aes(x = value, y = Strat_depth), size = 0.3) +
  geom_lineh(lwd = 0.3) +
  #geom_point(size = 0.2) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HERPT1S4 [SAT]: Composite") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("PWave Vel." = "ms-1", "GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "Res." = "Ohm.m"))
HERPT1S4_p3
ggsave("Outputs/Figures/HERPT1S4/Fig1.3_ALL_comp_black.pdf",
       height = c(10), width = c(16), dpi = 600, units = "cm")

# Multi-plots per page
ggarrange(HERPT1S4_p1, HERPT1S4_p2, nrow = 2)
ggsave("Outputs/Figures/HERPT1S4/Fig1.4_over_comp_coresx2.pdf",
       height = c(20), width = c(20), dpi = 600, units = "cm")

ggarrange(HERPT1S4_p1, HERPT1S4_p2, HERPT1S4_p3, nrow = 3)
ggsave("Outputs/Figures/HERPT1S4/Fig1.5_over_comp_coresx3.pdf",
       height = c(30), width = c(20), dpi = 600, units = "cm")

# Code to end doesnt work without an age depth model ----------------------

# Age depth model - plot to check ---------------------------------------

#TO DO - add adm plot to core plots
HERPT1S4_adm_plot <- age_depth_model(
  HERPT1S4_comp,
  depth = Strat_depth, age = SH20_mean_age,
  age_max = SH20_max_age_95CI,
  age_min = SH20_min_age_95CI
) 

plot(HERPT1S4_adm_plot)

# save plot to file 
pdf(file = "Outputs/Figures/HERPT1S4/Fig1.6_age_depth_model.pdf",
    width = 6, # The width of the plot in inches
    height = 6) # The height of the plot in inches
HERPT1S4_adm_plot <- age_depth_model(
  HERPT1S4_comp,
  depth = Strat_depth, age = SH20_mean_age,
  age_max = SH20_max_age_95CI,
  age_min = SH20_min_age_95CI
) 
plot(HERPT1S4_adm_plot)
dev.off()

# Composite plots for paper ------------------------------------------

HERPT1S4_plot <- HERPT1S4_comp_long2 %>% 
  mutate(param = fct_relevel(param, "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"WMAR_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "WMAR_SAT" = "WMAR")) %>%
  ggplot(aes(x = value, y = Strat_depth), size = 0.01) +
  geom_lineh(lwd = 0.3) +
  geom_point(size = 0.01) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HERPT1S4 [SAT]") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "WMAR" = "g cm-2  a-1"))
HERPT1S4_plot

# OPTIONAL
# add hiatus (red) horizontal line and other points of interest (green) - or comment out

HERPT1S4_plot1 <- HERPT1S4_plot

# OR

#HERPT1S4_plot1 <- HERPT1S4_plot +
#geom_hline(yintercept = c(76), col = "red", lty = 1, alpha = 1, lwd = 0.5)
#HERPT1S4_plot1

# Add age on secondary age y axis
HERPT1S4_adm <- age_depth_model(
  HERPT1S4_comp, 
  depth = Strat_depth, age = SH20_mean_age)

HERPT1S4_plot_age <- HERPT1S4_plot1 + 
  scale_y_depth_age(HERPT1S4_adm, age_name = "Age (cal a BP)")
HERPT1S4_plot_age

# CONISS
# Uses Den1_SAT, DCMS_SAT only in _com_long3, but plotted with MS1_SAT and WMAR
# CONISS using nested_chclust_coniss() and the rioja package for constrained cluster analysis 
# By default, this calculates the number of significant zones based on a broken stick simulation
# TO DO <- add number of groups and their depth/ages

HERPT1S4_coniss <- HERPT1S4_comp_long3 %>%
  nested_data(qualifiers = c(SH20_mean_age, Strat_depth), 
              key = param, value = value, trans = scale) %>%
  nested_chclust_coniss()

HERPT1S4_plot_coniss <- HERPT1S4_plot_age +
  layer_dendrogram(HERPT1S4_coniss, aes(y = Strat_depth), 
                   param = "CONISS", size = 0.2, lwd = 0.2) +
  layer_zone_boundaries(HERPT1S4_coniss, aes(y = Strat_depth), 
                        size = 0.3, lwd = 0.3, colour = "red")

ggarrange(HERPT1S4_plot_age, HERPT1S4_plot_coniss, nrow = 2)
ggsave("Outputs/Figures/HERPT1S4/Fig2.1_comp_depth_coniss.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")

# Plot versus age
HERPT1S4_age <- HERPT1S4_comp_long2 %>% 
  mutate(param = fct_relevel(param, "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"WMAR_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "WMAR_SAT" = "WMAR")) %>%
  ggplot(aes(x = SH20_mean_age, y = value), size = 0.01) +
  geom_line(lwd = 0.3) +
  geom_point(size = 0.01) +
  #scale_y_reverse() +
  # add units to the graph headers
  facet_geochem_grid(
    vars(param),
    units = c("GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "WMAR" = "g cm-2  a-1")) +
  #facet_geochem_grid(vars(param)) +
  labs(x = "Age (cal a BP)", y = NULL)
HERPT1S4_age

# OPTIONAL
# shading and lines for HER14 - added hiatus as a shaded area - min-max = rounded to the nearest 10 years - not needed for other plots
#zone_data <- tibble(ymin = -Inf, ymax = Inf , xmin = 1800, xmax = 10790)

HERPT1S4_age_lines1 <- HERPT1S4_age +
  #geom_rect(
  #mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
  #data = zone_data, 
  #alpha = 0.2,
  #fill = "white",
  #inherit.aes = FALSE
  #) +
  # add lines of interest at 1 ka, 2 ka, 11.75 ka and 13.8 ka etc. - need to label these later on
  #geom_vline(xintercept = c(0, 1000, 2000, 3000, 7500, 10000, 13800), col = "blue", lty = 3, alpha = 0.7, lwd = 0.3) +
  geom_vline(xintercept = c(4200), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  #geom_vline(xintercept = c(4200), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  labs(title = "HERPT1S4 [SAT]")
HERPT1S4_age_lines1

# Save final figure
ggarrange(HERPT1S4_plot_coniss, HERPT1S4_age_lines1, nrow = 2)
ggsave("Outputs/Figures/HERPT1S4/Fig2.2_depth_age_summary.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")

HERPT1S4_age_lines2 <- HERPT1S4_age +
  #geom_rect(
  #mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
  #data = zone_data, 
  #alpha = 0.2,
  #fill = "white",
  #inherit.aes = FALSE
  #) +
  # add lines of interest at 1 ka, 2 ka, 11.75 ka and 13.8 ka etc. - need to label these later on
  geom_vline(xintercept = c(0, 1000, 2000, 3000, 7500, 10000, 13800), col = "blue", lty = 3, alpha = 0.7, lwd = 0.3) +
  geom_vline(xintercept = c(4200, 8200, 11750, 16000), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  #geom_vline(xintercept = c(4200), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  labs(title = "HERPT1S4 [SAT]")
HERPT1S4_age_lines2

# Save final figure
ggarrange(HERPT1S4_plot_coniss, HERPT1S4_age_lines2, nrow = 2)
ggsave("Outputs/Figures/HERPT1S4/Fig2.3_depth_age_summary.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")






# HERPT1S8 - plots  ----------------------------------------------

# Overlap and Composite plots for XRAD mapping ---------------------

# All data - showing overlaps
HERPT1S8_p1 <- HERPT1S8_over_long %>% 
  mutate(param = fct_relevel(param, "PWVel_SAT", "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"Imp_SAT",	"FP_SAT",	"RES_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "Imp_SAT" = "Imp.", "FP_SAT" = "FP", "RES_SAT" = "Res.")) %>% 
  ggplot(aes(x = value, y = Strat_depth, colour = SECTION_ID, 
             fill = SECTION_ID, shape = SECTION_ID), size = 0.3) +
  geom_lineh(lwd = 0.3) +
  #geom_point(size = 0.2) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HERPT1S8 [SAT]: Overlaps") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("PWave Vel." = "ms-1", "GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "Res." = "Ohm.m"))
HERPT1S8_p1
ggsave("Outputs/Figures/HERPT1S8/Fig1.1_ALL_over.pdf",
       height = c(10), width = c(20), dpi = 600, units = "cm")

# All data - composite colour by core for XRAD single core mapping
HERPT1S8_p2 <- HERPT1S8_comp_long1 %>% 
  mutate(param = fct_relevel(param, "PWVel_SAT", "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"Imp_SAT",	"FP_SAT",	"RES_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "Imp_SAT" = "Imp.", "FP_SAT" = "FP", "RES_SAT" = "Res.")) %>% 
  ggplot(aes(x = value, y = Strat_depth, colour = SECTION_ID, 
             fill = SECTION_ID, shape = SECTION_ID), size = 0.3) +
  geom_lineh(lwd = 0.3) +
  #geom_point(size = 0.2) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HERPT1S8 [SAT]: Composite") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("PWave Vel." = "ms-1", "GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "10^-8 m3 kg-1", "Res." = "Ohm.m"))
HERPT1S8_p2
ggsave("Outputs/Figures/HERPT1S8/Fig1.2_ALL_comp.pdf",
       height = c(10), width = c(20), dpi = 600, units = "cm")

# All data - composite plotted as black
HERPT1S8_p3 <- HERPT1S8_comp_long1 %>% 
  mutate(param = fct_relevel(param, "PWVel_SAT", "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"Imp_SAT",	"FP_SAT",	"RES_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "Imp_SAT" = "Imp.", "FP_SAT" = "FP", "RES_SAT" = "Res.")) %>% 
  ggplot(aes(x = value, y = Strat_depth), size = 0.3) +
  geom_lineh(lwd = 0.3) +
  #geom_point(size = 0.2) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HERPT1S8 [SAT]: Composite") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("PWave Vel." = "ms-1", "GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "Res." = "Ohm.m"))
HERPT1S8_p3
ggsave("Outputs/Figures/HERPT1S8/Fig1.3_ALL_comp_black.pdf",
       height = c(10), width = c(16), dpi = 600, units = "cm")

# Multi-plots per page
ggarrange(HERPT1S8_p1, HERPT1S8_p2, nrow = 2)
ggsave("Outputs/Figures/HERPT1S8/Fig1.4_over_comp_coresx2.pdf",
       height = c(20), width = c(20), dpi = 600, units = "cm")

ggarrange(HERPT1S8_p1, HERPT1S8_p2, HERPT1S8_p3, nrow = 3)
ggsave("Outputs/Figures/HERPT1S8/Fig1.5_over_comp_coresx3.pdf",
       height = c(30), width = c(20), dpi = 600, units = "cm")

# Code to end doesnt work without an age depth model ----------------------

# Age depth model - plot to check ---------------------------------------

#TO DO - add adm plot to core plots
HERPT1S8_adm_plot <- age_depth_model(
  HERPT1S8_comp,
  depth = Strat_depth, age = SH20_mean_age,
  age_max = SH20_max_age_95CI,
  age_min = SH20_min_age_95CI
) 

plot(HERPT1S8_adm_plot)

# save plot to file 
pdf(file = "Outputs/Figures/HERPT1S8/Fig1.6_age_depth_model.pdf",
    width = 6, # The width of the plot in inches
    height = 6) # The height of the plot in inches
HERPT1S8_adm_plot <- age_depth_model(
  HERPT1S8_comp,
  depth = Strat_depth, age = SH20_mean_age,
  age_max = SH20_max_age_95CI,
  age_min = SH20_min_age_95CI
) 
plot(HERPT1S8_adm_plot)
dev.off()

# Composite plots for paper ------------------------------------------

HERPT1S8_plot <- HERPT1S8_comp_long2 %>% 
  mutate(param = fct_relevel(param, "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"WMAR_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "WMAR_SAT" = "WMAR")) %>%
  ggplot(aes(x = value, y = Strat_depth), size = 0.01) +
  geom_lineh(lwd = 0.3) +
  geom_point(size = 0.01) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HERPT1S8 [SAT]") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "WMAR" = "g cm-2  a-1"))
HERPT1S8_plot

# OPTIONAL
# add hiatus (red) horizontal line and other points of interest (green) - or comment out

HERPT1S8_plot1 <- HERPT1S8_plot

# OR

#HERPT1S8_plot1 <- HERPT1S8_plot +
#geom_hline(yintercept = c(76), col = "red", lty = 1, alpha = 1, lwd = 0.5)
#HERPT1S8_plot1

# Add age on secondary age y axis
HERPT1S8_adm <- age_depth_model(
  HERPT1S8_comp, 
  depth = Strat_depth, age = SH20_mean_age)

HERPT1S8_plot_age <- HERPT1S8_plot1 + 
  scale_y_depth_age(HERPT1S8_adm, age_name = "Age (cal a BP)")
HERPT1S8_plot_age

# CONISS
# Uses Den1_SAT, DCMS_SAT only in _com_long3, but plotted with MS1_SAT and WMAR
# CONISS using nested_chclust_coniss() and the rioja package for constrained cluster analysis 
# By default, this calculates the number of significant zones based on a broken stick simulation
# TO DO <- add number of groups and their depth/ages

HERPT1S8_coniss <- HERPT1S8_comp_long3 %>%
  nested_data(qualifiers = c(SH20_mean_age, Strat_depth), 
              key = param, value = value, trans = scale) %>%
  nested_chclust_coniss()

HERPT1S8_plot_coniss <- HERPT1S8_plot_age +
  layer_dendrogram(HERPT1S8_coniss, aes(y = Strat_depth), 
                   param = "CONISS", size = 0.2, lwd = 0.2) +
  layer_zone_boundaries(HERPT1S8_coniss, aes(y = Strat_depth), 
                        size = 0.3, lwd = 0.3, colour = "red")

ggarrange(HERPT1S8_plot_age, HERPT1S8_plot_coniss, nrow = 2)
ggsave("Outputs/Figures/HERPT1S8/Fig2.1_comp_depth_coniss.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")

# Plot versus age
HERPT1S8_age <- HERPT1S8_comp_long2 %>% 
  mutate(param = fct_relevel(param, "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"WMAR_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "WMAR_SAT" = "WMAR")) %>%
  ggplot(aes(x = SH20_mean_age, y = value), size = 0.01) +
  geom_line(lwd = 0.3) +
  geom_point(size = 0.01) +
  #scale_y_reverse() +
  # add units to the graph headers
  facet_geochem_grid(
    vars(param),
    units = c("GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "WMAR" = "g cm-2  a-1")) +
  #facet_geochem_grid(vars(param)) +
  labs(x = "Age (cal a BP)", y = NULL)
HERPT1S8_age

# OPTIONAL
# shading and lines for HER14 - added hiatus as a shaded area - min-max = rounded to the nearest 10 years - not needed for other plots
#zone_data <- tibble(ymin = -Inf, ymax = Inf , xmin = 1800, xmax = 10790)

HERPT1S8_age_lines1 <- HERPT1S8_age +
  #geom_rect(
  #mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
  #data = zone_data, 
  #alpha = 0.2,
  #fill = "white",
  #inherit.aes = FALSE
  #) +
  # add lines of interest at 1 ka, 2 ka, 11.75 ka and 13.8 ka etc. - need to label these later on
  #geom_vline(xintercept = c(0, 1000, 2000, 3000, 7500, 10000, 13800), col = "blue", lty = 3, alpha = 0.7, lwd = 0.3) +
  geom_vline(xintercept = c(4200), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  #geom_vline(xintercept = c(4200), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  labs(title = "HERPT1S8 [SAT]")
HERPT1S8_age_lines1

# Save final figure
ggarrange(HERPT1S8_plot_coniss, HERPT1S8_age_lines1, nrow = 2)
ggsave("Outputs/Figures/HERPT1S8/Fig2.2_depth_age_summary.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")

HERPT1S8_age_lines2 <- HERPT1S8_age +
  #geom_rect(
  #mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
  #data = zone_data, 
  #alpha = 0.2,
  #fill = "white",
  #inherit.aes = FALSE
  #) +
  # add lines of interest at 1 ka, 2 ka, 11.75 ka and 13.8 ka etc. - need to label these later on
  geom_vline(xintercept = c(0, 1000, 2000, 3000, 7500, 10000, 13800), col = "blue", lty = 3, alpha = 0.7, lwd = 0.3) +
  geom_vline(xintercept = c(4200, 8200, 11750, 16000), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  #geom_vline(xintercept = c(4200), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  labs(title = "HERPT1S8 [SAT]")
HERPT1S8_age_lines2

# Save final figure
ggarrange(HERPT1S8_plot_coniss, HERPT1S8_age_lines2, nrow = 2)
ggsave("Outputs/Figures/HERPT1S8/Fig2.3_depth_age_summary.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")




# -------------------------------------------------------------------------

# # HER34L - plots to do --------------------------------------------------

# Overlap and Composite plots for XRAD mapping ---------------------

# All data - showing overlaps
HER34L_p1 <- HER34L_over_long %>% 
  mutate(param = fct_relevel(param, "PWVel_SAT", "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"Imp_SAT",	"FP_SAT",	"RES_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "Imp_SAT" = "Imp.", "FP_SAT" = "FP", "RES_SAT" = "Res.")) %>% 
  ggplot(aes(x = value, y = Strat_depth, colour = SECTION_ID, 
             fill = SECTION_ID, shape = SECTION_ID), size = 0.3
  ) +
  geom_lineh(lwd = 0.3) +
  #geom_point(size = 0.2) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HER34L [SAT]: Overlaps") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("PWave Vel." = "ms-1", "GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "Res." = "Ohm.m"))
HER34L_p1
ggsave("Outputs/Figures/HER34L/Fig1.1_ALL_over.pdf",
       height = c(10), width = c(20), dpi = 600, units = "cm")

# All data - composite colour by core for XRAD single core mapping
HER34L_p2 <- HER34L_comp_long1 %>% 
  mutate(param = fct_relevel(param, "PWVel_SAT", "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"Imp_SAT",	"FP_SAT",	"RES_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "Imp_SAT" = "Imp.", "FP_SAT" = "FP", "RES_SAT" = "Res.")) %>% 
  ggplot(aes(x = value, y = Strat_depth, colour = SECTION_ID, 
             fill = SECTION_ID, shape = SECTION_ID), size = 0.3) +
  geom_lineh(lwd = 0.3) +
  #geom_point(size = 0.2) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HER34L [SAT]: Composite") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("PWave Vel." = "ms-1", "GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "10^-8 m3 kg-1", "Res." = "Ohm.m"))
HER34L_p2
ggsave("Outputs/Figures/HER34L/Fig1.2_ALL_comp.pdf",
       height = c(10), width = c(20), dpi = 600, units = "cm")

# All data - composite plotted as black
HER34L_p3 <- HER34L_comp_long1 %>% 
  mutate(param = fct_relevel(param, "PWVel_SAT", "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"Imp_SAT",	"FP_SAT",	"RES_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "Imp_SAT" = "Imp.", "FP_SAT" = "FP", "RES_SAT" = "Res.")) %>% 
  ggplot(aes(x = value, y = Strat_depth), size = 0.5) +
  geom_lineh(lwd = 0.5) +
  #geom_point(size = 0.2) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HER34L [SAT]: Composite") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("PWave Vel." = "ms-1", "GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "Res." = "Ohm.m"))
HER34L_p3
ggsave("Outputs/Figures/HER34L/Fig1.3_ALL_comp_black.pdf",
       height = c(10), width = c(16), dpi = 600, units = "cm")

# Multi-plots per page
ggarrange(HER34L_p1, HER34L_p2, nrow = 2)
ggsave("Outputs/Figures/HER34L/Fig1.4_over_comp_coresx2.pdf",
       height = c(20), width = c(20), dpi = 600, units = "cm")

ggarrange(HER34L_p1, HER34L_p2, HER34L_p3, nrow = 3)
ggsave("Outputs/Figures/HER34L/Fig1.5_over_comp_coresx3.pdf",
       height = c(30), width = c(20), dpi = 600, units = "cm")

# Age depth model - plot to check ---------------------------------------

#TO DO - add adm plot to core plots
HER34L_adm_plot <- age_depth_model(
  HER34L_comp,
  depth = Strat_depth, age = SH20_mean_age,
  age_max = SH20_max_age_95CI,
  age_min = SH20_min_age_95CI
) 

plot(HER34L_adm_plot)

# save plot to file 
pdf(file = "Outputs/Figures/HER34L/Fig1.6_age_depth_model.pdf",
    width = 6, # The width of the plot in inches
    height = 6) # The height of the plot in inches
HER34L_adm_plot <- age_depth_model(
  HER34L_comp,
  depth = Strat_depth, age = SH20_mean_age,
  age_max = SH20_max_age_95CI,
  age_min = SH20_min_age_95CI
) 
plot(HER34L_adm_plot)
dev.off()

# Composite plots for paper ------------------------------------------

HER34L_plot <- HER34L_comp_long2 %>% 
  mutate(param = fct_relevel(param, "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"WMAR_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "WMAR_SAT" = "WMAR")) %>%
  ggplot(aes(x = value, y = Strat_depth), size = 0.01) +
  geom_lineh(lwd = 0.3) +
  geom_point(size = 0.01) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HER34L [SAT]") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "WMAR" = "g cm-2  a-1"))
HER34L_plot

# OPTIONAL
# add hiatus (red) horizontal line and other points of interest (green) - or comment out

HER34L_plot1 <- HER34L_plot

# OR

#HER34L_plot1 <- HER34L_plot +
#geom_hline(yintercept = c(76), col = "red", lty = 1, alpha = 1, lwd = 0.5)
#HER34L_plot1

# Add age on secondary age y axis
HER34L_adm <- age_depth_model(
  HER34L_comp, 
  depth = Strat_depth, age = SH20_mean_age)

HER34L_plot_age <- HER34L_plot1 + 
  scale_y_depth_age(HER34L_adm, age_name = "Age (cal a BP)")
HER34L_plot_age

# CONISS
# Uses Den1_SAT, DCMS_SAT only in _com_long3, but plotted with MS1_SAT and WMAR
# CONISS using nested_chclust_coniss() and the rioja package for constrained cluster analysis 
# By default, this calculates the number of significant zones based on a broken stick simulation
# TO DO <- add number of groups and their depth/ages

HER34L_coniss <- HER34L_comp_long3 %>%
  nested_data(qualifiers = c(SH20_mean_age, Strat_depth), 
              key = param, value = value, trans = scale) %>%
  nested_chclust_coniss()

HER34L_plot_coniss <- HER34L_plot_age +
  layer_dendrogram(HER34L_coniss, aes(y = Strat_depth), 
                   param = "CONISS", size = 0.2, lwd = 0.2) +
  layer_zone_boundaries(HER34L_coniss, aes(y = Strat_depth), 
                        size = 0.3, lwd = 0.3, colour = "red")
HER34L_plot_coniss

ggarrange(HER34L_plot_age, HER34L_plot_coniss, nrow = 2)
ggsave("Outputs/Figures/HER34L/Fig2.1_comp_depth_coniss.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")

# Plot versus age
HER34L_age <- HER34L_comp_long2 %>% 
  mutate(param = fct_relevel(param, "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"WMAR_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "WMAR_SAT" = "WMAR")) %>%
  ggplot(aes(x = SH20_mean_age, y = value), size = 0.01) +
  geom_line(lwd = 0.3) +
  geom_point(size = 0.01) +
  #scale_y_reverse() +
  # add units to the graph headers
  facet_geochem_grid(
    vars(param),
    units = c("GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "WMAR" = "g cm-2  a-1")) +
  #facet_geochem_grid(vars(param)) +
  labs(x = "Age (cal a BP)", y = NULL)
HER34L_age

# OPTIONAL
# shading and lines for HER14 - added hiatus as a shaded area - min-max = rounded to the nearest 10 years - not needed for other plots
#zone_data <- tibble(ymin = -Inf, ymax = Inf , xmin = 1800, xmax = 10790)

HER34L_age_lines1 <- HER34L_age +
  #geom_rect(
  #mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
  #data = zone_data, 
  #alpha = 0.2,
  #fill = "white",
  #inherit.aes = FALSE
  #) +
  # add lines of interest at 1 ka, 2 ka, 11.75 ka and 13.8 ka etc. - need to label these later on
  #geom_vline(xintercept = c(0, 1000, 2000, 3000, 7500, 10000, 13800), col = "blue", lty = 3, alpha = 0.7, lwd = 0.3) +
  geom_vline(xintercept = c(4200), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  labs(title = "HER34L [SAT]")
HER34L_age_lines1

# Save final figure
ggarrange(HER34L_plot_coniss, HER34L_age_lines1, nrow = 2)
ggsave("Outputs/Figures/HER34L/Fig2.2_depth_age_summary.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")

HER34L_age_lines2 <- HER34L_age +
  #geom_rect(
  #mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
  #data = zone_data, 
  #alpha = 0.2,
  #fill = "white",
  #inherit.aes = FALSE
  #) +
  # add lines of interest at 1 ka, 2 ka, 11.75 ka and 13.8 ka etc. - need to label these later on
  geom_vline(xintercept = c(0, 1000, 2000, 3000, 7500, 10000, 13800), col = "blue", lty = 3, alpha = 0.7, lwd = 0.3) +
  geom_vline(xintercept = c(4200, 8200, 11750, 16000), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  #geom_vline(xintercept = c(4200), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  labs(title = "HER34L [SAT]")
HER34L_age_lines2

# Save final figure
ggarrange(HER34L_plot_coniss, HER34L_age_lines2, nrow = 2)
ggsave("Outputs/Figures/HER34L/Fig2.3_depth_age_summary.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")

# -------------------------------------------------------------------------

# # HER42L --------------------------------------------------

# Overlap and Composite plots for XRAD mapping ---------------------

# All data - showing overlaps
HER42L_p1 <- HER42L_over_long %>% 
  mutate(param = fct_relevel(param, "PWVel_SAT", "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"Imp_SAT",	"FP_SAT",	"RES_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "Imp_SAT" = "Imp.", "FP_SAT" = "FP", "RES_SAT" = "Res.")) %>% 
  ggplot(aes(x = value, y = Strat_depth, colour = SECTION_ID, 
             fill = SECTION_ID, shape = SECTION_ID), size = 0.3
  ) +
  geom_lineh(lwd = 0.3) +
  #geom_point(size = 0.2) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HER42PB [SAT]: Overlaps") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("PWave Vel." = "ms-1", "GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "Res." = "Ohm.m"))
HER42L_p1
ggsave("Outputs/Figures/HER42L/Fig1.1_ALL_over.pdf",
       height = c(10), width = c(20), dpi = 600, units = "cm")

# All data - composite colour by core for XRAD single core mapping
HER42L_p2 <- HER42L_comp_long1 %>% 
  mutate(param = fct_relevel(param, "PWVel_SAT", "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"Imp_SAT",	"FP_SAT",	"RES_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "Imp_SAT" = "Imp.", "FP_SAT" = "FP", "RES_SAT" = "Res.")) %>% 
  ggplot(aes(x = value, y = Strat_depth, colour = SECTION_ID, 
             fill = SECTION_ID, shape = SECTION_ID), size = 0.3
  ) +
  geom_lineh(lwd = 0.3) +
  #geom_point(size = 0.2) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HER42L [SAT]: Composite") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("PWave Vel." = "ms-1", "GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "10^-8 m3 kg-1", "Res." = "Ohm.m"))
HER42L_p2
ggsave("Outputs/Figures/HER42L/Fig1.2_ALL_comp.pdf",
       height = c(10), width = c(20), dpi = 600, units = "cm")

# All data - composite plotted as black
HER42L_p3 <- HER42L_comp_long1 %>% 
  mutate(param = fct_relevel(param, "PWVel_SAT", "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"Imp_SAT",	"FP_SAT",	"RES_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "Imp_SAT" = "Imp.", "FP_SAT" = "FP", "RES_SAT" = "Res.")) %>% 
  ggplot(aes(x = value, y = Strat_depth), size = 0.3) +
  geom_lineh(lwd = 0.3) +
  #geom_point(size = 0.2) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HER42L [SAT]: Composite") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("PWave Vel." = "ms-1", "GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "Res." = "Ohm.m"))
HER42L_p3
ggsave("Outputs/Figures/HER42L/Fig1.3_ALL_comp_black.pdf",
       height = c(10), width = c(16), dpi = 600, units = "cm")

# Multi-plots per page
ggarrange(HER42L_p1, HER42L_p2, nrow = 2)
ggsave("Outputs/Figures/HER42L/Fig1.4_over_comp_coresx2.pdf",
       height = c(20), width = c(20), dpi = 600, units = "cm")

ggarrange(HER42L_p1, HER42L_p2, HER42L_p3, nrow = 3)
ggsave("Outputs/Figures/HER42L/Fig1.5_over_comp_coresx3.pdf",
       height = c(30), width = c(20), dpi = 600, units = "cm")

# Age depth model - plot to check ---------------------------------------

#TO DO - add adm plot to core plots
HER42L_adm_plot <- age_depth_model(
  HER42L_comp,
  depth = Strat_depth, age = SH20_mean_age,
  age_max = SH20_max_age_95CI,
  age_min = SH20_min_age_95CI
) 

plot(HER42L_adm_plot)

# save plot to file 
pdf(file = "Outputs/Figures/HER42L/Fig1.6_age_depth_model.pdf",
    width = 6, # The width of the plot in inches
    height = 6) # The height of the plot in inches
HER42L_adm_plot <- age_depth_model(
  HER42L_comp,
  depth = Strat_depth, age = SH20_mean_age,
  age_max = SH20_max_age_95CI,
  age_min = SH20_min_age_95CI
) 
plot(HER42L_adm_plot)
dev.off()

# Composite plots for paper ------------------------------------------

HER42L_plot <- HER42L_comp_long2 %>% 
  mutate(param = fct_relevel(param, "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"WMAR_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "WMAR_SAT" = "WMAR")) %>%
  ggplot(aes(x = value, y = Strat_depth), size = 0.01) +
  geom_lineh(lwd = 0.3) +
  geom_point(size = 0.01) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HER42L [SAT]") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "WMAR" = "g cm-2  a-1"))
HER42L_plot

# OPTIONAL
# add hiatus (red) horizontal line and other points of interest (green) - or comment out

HER42L_plot1 <- HER42L_plot

# OR

#HER42L_plot1 <- HER42L_plot +
#geom_hline(yintercept = c(76), col = "red", lty = 1, alpha = 1, lwd = 0.5)
#HER42L_plot1

# Add age on secondary age y axis
HER42L_adm <- age_depth_model(
  HER42L_comp, 
  depth = Strat_depth, age = SH20_mean_age)

HER42L_plot_age <- HER42L_plot1 + 
  scale_y_depth_age(HER42L_adm, age_name = "Age (cal a BP)")
HER42L_plot_age

# CONISS
# Uses Den1_SAT, DCMS_SAT only in _com_long3, but plotted with MS1_SAT and WMAR
# CONISS using nested_chclust_coniss() and the rioja package for constrained cluster analysis 
# By default, this calculates the number of significant zones based on a broken stick simulation
# TO DO <- add number of groups and their depth/ages

HER42L_coniss <- HER42L_comp_long3 %>%
  nested_data(qualifiers = c(SH20_mean_age, Strat_depth), 
              key = param, value = value, trans = scale) %>%
  nested_chclust_coniss()

HER42L_plot_coniss <- HER42L_plot_age +
  layer_dendrogram(HER42L_coniss, aes(y = Strat_depth), 
                   param = "CONISS", size = 0.2, lwd = 0.2) +
  layer_zone_boundaries(HER42L_coniss, aes(y = Strat_depth), 
                        size = 0.3, lwd = 0.3, colour = "red")
HER42L_plot_coniss

ggarrange(HER42L_plot_age, HER42L_plot_coniss, nrow = 2)
ggsave("Outputs/Figures/HER42L/Fig2.1_comp_depth_coniss.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")

# Plot versus age
HER42L_age <- HER42L_comp_long2 %>% 
  mutate(param = fct_relevel(param, "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"WMAR_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "WMAR_SAT" = "WMAR")) %>%
  ggplot(aes(x = SH20_mean_age, y = value), size = 0.01) +
  geom_line(lwd = 0.3) +
  geom_point(size = 0.01) +
  #scale_y_reverse() +
  # add units to the graph headers
  facet_geochem_grid(
    vars(param),
    units = c("GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "WMAR" = "g cm-2  a-1")) +
  #facet_geochem_grid(vars(param)) +
  labs(x = "Age (cal a BP)", y = NULL)
HER42L_age

# OPTIONAL
# shading and lines for HER14 - added hiatus as a shaded area - min-max = rounded to the nearest 10 years - not needed for other plots
#zone_data <- tibble(ymin = -Inf, ymax = Inf , xmin = 1800, xmax = 10790)

HER42L_age_lines1 <- HER42L_age +
  #geom_rect(
  #mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
  #data = zone_data, 
  #alpha = 0.2,
  #fill = "white",
  #inherit.aes = FALSE
  #) +
  # add lines of interest at 1 ka, 2 ka, 11.75 ka and 13.8 ka etc. - need to label these later on
  #geom_vline(xintercept = c(0, 1000, 2000, 3000, 7500, 10000, 13800), col = "blue", lty = 3, alpha = 0.7, lwd = 0.3) +
  geom_vline(xintercept = c(4200), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  #geom_vline(xintercept = c(4200), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  labs(title = "HER42L [SAT]")
HER42L_age_lines1

# Save final figure
ggarrange(HER42L_plot_coniss, HER42L_age_lines1, nrow = 2)
ggsave("Outputs/Figures/HER42L/Fig2.2_depth_age_summary.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")

HER42L_age_lines2 <- HER42L_age +
  #geom_rect(
  #mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
  #data = zone_data, 
  #alpha = 0.2,
  #fill = "white",
  #inherit.aes = FALSE
  #) +
  # add lines of interest at 1 ka, 2 ka, 11.75 ka and 13.8 ka etc. - need to label these later on
  geom_vline(xintercept = c(0, 1000, 2000, 3000, 7500, 10000, 13800), col = "blue", lty = 3, alpha = 0.7, lwd = 0.3) +
  geom_vline(xintercept = c(4200, 8200, 11750, 16000), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  #geom_vline(xintercept = c(4200), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  labs(title = "HER42L [SAT]")
HER42L_age_lines2

# Save final figure
ggarrange(HER42L_plot_coniss, HER42L_age_lines2, nrow = 2)
ggsave("Outputs/Figures/HER42L/Fig2.3_depth_age_summary.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")

# -------------------------------------------------------------------------

# # HER42PB -------------------------------------------------

# Overlap and Composite plots for XRAD mapping ---------------------

# All data - showing overlaps
HER42PB_p1 <- HER42PB_over_long %>% 
  mutate(param = fct_relevel(param, "PWVel_SAT", "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"Imp_SAT",	"FP_SAT",	"RES_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "Imp_SAT" = "Imp.", "FP_SAT" = "FP", "RES_SAT" = "Res.")) %>% 
  ggplot(aes(x = value, y = Strat_depth, colour = SECTION_ID, 
             fill = SECTION_ID, shape = SECTION_ID), size = 0.3
  ) +
  geom_lineh(lwd = 0.3) +
  #geom_point(size = 0.2) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HER42PB [SAT]: Overlaps") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("PWave Vel." = "ms-1", "GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "Res." = "Ohm.m"))
HER42PB_p1
ggsave("Outputs/Figures/HER42PB/Fig1.1_ALL_over.pdf",
       height = c(10), width = c(20), dpi = 600, units = "cm")

# All data - composite colour by core for XRAD single core mapping
HER42PB_p2 <- HER42PB_comp_long1 %>% 
  mutate(param = fct_relevel(param, "PWVel_SAT", "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"Imp_SAT",	"FP_SAT",	"RES_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "Imp_SAT" = "Imp.", "FP_SAT" = "FP", "RES_SAT" = "Res.")) %>% 
  ggplot(aes(x = value, y = Strat_depth, colour = SECTION_ID, 
             fill = SECTION_ID, shape = SECTION_ID), size = 0.3
  ) +
  geom_lineh(lwd = 0.3) +
  #geom_point(size = 0.2) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HER42PB [SAT]: Composite") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("PWave Vel." = "ms-1", "GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "10^-8 m3 kg-1", "Res." = "Ohm.m"))
HER42PB_p2
ggsave("Outputs/Figures/HER42PB/Fig1.2_ALL_comp.pdf",
       height = c(10), width = c(20), dpi = 600, units = "cm")

# All data - composite plotted as black
HER42PB_p3 <- HER42PB_comp_long1 %>% 
  mutate(param = fct_relevel(param, "PWVel_SAT", "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"Imp_SAT",	"FP_SAT",	"RES_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "Imp_SAT" = "Imp.", "FP_SAT" = "FP", "RES_SAT" = "Res.")) %>% 
  ggplot(aes(x = value, y = Strat_depth), size = 0.3) +
  geom_lineh(lwd = 0.3) +
  #geom_point(size = 0.2) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HER42PB [SAT]: Composite") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("PWave Vel." = "ms-1", "GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "Res." = "Ohm.m"))
HER42PB_p3
ggsave("Outputs/Figures/HER42PB/Fig1.3_ALL_comp_black.pdf",
       height = c(10), width = c(16), dpi = 600, units = "cm")

# Multi-plots per page
ggarrange(HER42PB_p1, HER42PB_p2, nrow = 2)
ggsave("Outputs/Figures/HER42PB/Fig1.4_over_comp_coresx2.pdf",
       height = c(20), width = c(20), dpi = 600, units = "cm")

ggarrange(HER42PB_p1, HER42PB_p2, HER42PB_p3, nrow = 3)
ggsave("Outputs/Figures/HER42PB/Fig1.5_over_comp_coresx3.pdf",
       height = c(30), width = c(20), dpi = 600, units = "cm")

# Age depth model - plot to check ---------------------------------------

#TO DO - add adm plot to core plots
HER42PB_adm_plot <- age_depth_model(
  HER42PB_comp,
  depth = Strat_depth, age = SH20_mean_age,
  age_max = SH20_max_age_95CI,
  age_min = SH20_min_age_95CI
) 

plot(HER42PB_adm_plot)

# save plot to file 
pdf(file = "Outputs/Figures/HER42PB/Fig1.6_age_depth_model.pdf",
    width = 6, # The width of the plot in inches
    height = 6) # The height of the plot in inches
HER42PB_adm_plot <- age_depth_model(
  HER42PB_comp,
  depth = Strat_depth, age = SH20_mean_age,
  age_max = SH20_max_age_95CI,
  age_min = SH20_min_age_95CI
) 
plot(HER42PB_adm_plot)
dev.off()

# Composite plots for paper ------------------------------------------

HER42PB_plot <- HER42PB_comp_long2 %>% 
  mutate(param = fct_relevel(param, "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"WMAR_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "WMAR_SAT" = "WMAR")) %>%
  ggplot(aes(x = value, y = Strat_depth), size = 0.01) +
  geom_lineh(lwd = 0.3) +
  geom_point(size = 0.01) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HER42PB [SAT]") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "WMAR" = "g cm-2  a-1"))
HER42PB_plot

# OPTIONAL
# add hiatus (red) horizontal line and other points of interest (green) - or comment out

HER42PB_plot1 <- HER42PB_plot

# OR

#HER42PB_plot1 <- HER42PB_plot +
#geom_hline(yintercept = c(76), col = "red", lty = 1, alpha = 1, lwd = 0.5)
#HER42PB_plot1

# Add age on secondary age y axis
HER42PB_adm <- age_depth_model(
  HER42PB_comp, 
  depth = Strat_depth, age = SH20_mean_age)

HER42PB_plot_age <- HER42PB_plot1 + 
  scale_y_depth_age(HER42PB_adm, age_name = "Age (cal a BP)")
HER42PB_plot_age

# CONISS
# Uses Den1_SAT, DCMS_SAT only in _com_long3, but plotted with MS1_SAT and WMAR
# CONISS using nested_chclust_coniss() and the rioja package for constrained cluster analysis 
# By default, this calculates the number of significant zones based on a broken stick simulation
# TO DO <- add number of groups and their depth/ages

HER42PB_coniss <- HER42PB_comp_long3 %>%
  nested_data(qualifiers = c(SH20_mean_age, Strat_depth), 
              key = param, value = value, trans = scale) %>%
  nested_chclust_coniss()

HER42PB_plot_coniss <- HER42PB_plot_age +
  layer_dendrogram(HER42PB_coniss, aes(y = Strat_depth), 
                   param = "CONISS", size = 0.2, lwd = 0.2) +
  layer_zone_boundaries(HER42PB_coniss, aes(y = Strat_depth), 
                        size = 0.3, lwd = 0.3, colour = "red")
HER42PB_plot_coniss

ggarrange(HER42PB_plot_age, HER42PB_plot_coniss, nrow = 2)
ggsave("Outputs/Figures/HER42PB/Fig2.1_comp_depth_coniss.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")

# Plot versus age
HER42PB_age <- HER42PB_comp_long2 %>% 
  mutate(param = fct_relevel(param, "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"WMAR_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "WMAR_SAT" = "WMAR")) %>%
  ggplot(aes(x = SH20_mean_age, y = value), size = 0.01) +
  geom_line(lwd = 0.3) +
  geom_point(size = 0.01) +
  #scale_y_reverse() +
  # add units to the graph headers
  facet_geochem_grid(
    vars(param),
    units = c("GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "WMAR" = "g cm-2  a-1")) +
  #facet_geochem_grid(vars(param)) +
  labs(x = "Age (cal a BP)", y = NULL)
HER42PB_age

# OPTIONAL
# shading and lines for HER14 - added hiatus as a shaded area - min-max = rounded to the nearest 10 years - not needed for other plots
#zone_data <- tibble(ymin = -Inf, ymax = Inf , xmin = 1800, xmax = 10790)

HER42PB_age_lines1 <- HER42PB_age +
  #geom_rect(
  #mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
  #data = zone_data, 
  #alpha = 0.2,
  #fill = "white",
  #inherit.aes = FALSE
  #) +
  # add lines of interest at 1 ka, 2 ka, 11.75 ka and 13.8 ka etc. - need to label these later on
  #geom_vline(xintercept = c(0, 1000, 2000, 3000, 7500, 10000, 13800), col = "blue", lty = 3, alpha = 0.7, lwd = 0.3) +
  geom_vline(xintercept = c(4200), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  #geom_vline(xintercept = c(4200), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  labs(title = "HER42PB [SAT]")
HER42PB_age_lines1

# Save final figure
ggarrange(HER42PB_plot_coniss, HER42PB_age_lines1, nrow = 2)
ggsave("Outputs/Figures/HER42PB/Fig2.2_depth_age_summary.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")

HER42PB_age_lines2 <- HER42PB_age +
  #geom_rect(
  #mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
  #data = zone_data, 
  #alpha = 0.2,
  #fill = "white",
  #inherit.aes = FALSE
  #) +
  # add lines of interest at 1 ka, 2 ka, 11.75 ka and 13.8 ka etc. - need to label these later on
  geom_vline(xintercept = c(0, 1000, 2000, 3000, 7500, 10000, 13800), col = "blue", lty = 3, alpha = 0.7, lwd = 0.3) +
  geom_vline(xintercept = c(4200, 8200, 11750, 16000), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  #geom_vline(xintercept = c(4200), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  labs(title = "HER42PB [SAT]")
HER42PB_age_lines2

# Save final figure
ggarrange(HER42PB_plot_coniss, HER42PB_age_lines2, nrow = 2)
ggsave("Outputs/Figures/HER42L/Fig2.3_depth_age_summary.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")

# -------------------------------------------------------------------------

# # HERPT2S1_SS2_XRF6_SEC42 plots  ----------------------------------------------

# Overlap and Composite plots for XRAD mapping ---------------------

# All data - showing overlaps
HERPT2S1_p1 <- HERPT2S1_over_long %>% 
  mutate(param = fct_relevel(param, "PWVel_SAT", "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"Imp_SAT",	"FP_SAT",	"RES_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "Imp_SAT" = "Imp.", "FP_SAT" = "FP", "RES_SAT" = "Res.")) %>% 
  ggplot(aes(x = value, y = Strat_depth, colour = SECTION_ID, 
             fill = SECTION_ID, shape = SECTION_ID), size = 0.3
  ) +
  geom_lineh(lwd = 0.3) +
  #geom_point(size = 0.2) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HERPT2S1 [SAT]: Overlaps") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("PWave Vel." = "ms-1", "GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "Res." = "Ohm.m"))
HERPT2S1_p1
ggsave("Outputs/Figures/HERPT2S1/Fig1.1_ALL_over.pdf",
       height = c(10), width = c(20), dpi = 600, units = "cm")

# All data - composite colour by core for XRAD single core mapping
HERPT2S1_p2 <- HERPT2S1_comp_long1 %>% 
  mutate(param = fct_relevel(param, "PWVel_SAT", "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"Imp_SAT",	"FP_SAT",	"RES_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "Imp_SAT" = "Imp.", "FP_SAT" = "FP", "RES_SAT" = "Res.")) %>% 
  ggplot(aes(x = value, y = Strat_depth, colour = SECTION_ID, 
             fill = SECTION_ID, shape = SECTION_ID), size = 0.3
  ) +
  geom_lineh(lwd = 0.3) +
  #geom_point(size = 0.2) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HERPT2S1 [SAT]: Composite") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("PWave Vel." = "ms-1", "GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "10^-8 m3 kg-1", "Res." = "Ohm.m"))
HERPT2S1_p2
ggsave("Outputs/Figures/HERPT2S1/Fig1.2_ALL_comp.pdf",
       height = c(10), width = c(20), dpi = 600, units = "cm")

# All data - composite plotted as black
HERPT2S1_p3 <- HERPT2S1_comp_long1 %>% 
  mutate(param = fct_relevel(param, "PWVel_SAT", "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"Imp_SAT",	"FP_SAT",	"RES_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "Imp_SAT" = "Imp.", "FP_SAT" = "FP", "RES_SAT" = "Res.")) %>% 
  ggplot(aes(x = value, y = Strat_depth), size = 0.5) +
  geom_lineh(lwd = 0.5) +
  #geom_point(size = 0.2) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HERPT2S1 [SAT]: Composite") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("PWave Vel." = "ms-1", "GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "Res." = "Ohm.m"))
HERPT2S1_p3
ggsave("Outputs/Figures/HERPT2S1/Fig1.3_ALL_comp_black.pdf",
       height = c(10), width = c(16), dpi = 600, units = "cm")

# Multi-plots per page
ggarrange(HERPT2S1_p1, HERPT2S1_p2, nrow = 2)
ggsave("Outputs/Figures/HERPT2S1/Fig1.4_over_comp_coresx2.pdf",
       height = c(20), width = c(20), dpi = 600, units = "cm")

ggarrange(HERPT2S1_p1, HERPT2S1_p2, HERPT2S1_p3, nrow = 3)
ggsave("Outputs/Figures/HERPT2S1/Fig1.5_over_comp_coresx3.pdf",
       height = c(30), width = c(20), dpi = 600, units = "cm")

# Age depth model - plot to check ---------------------------------------

#TO DO - add adm plot to core plots
HERPT2S1_adm_plot <- age_depth_model(
  HERPT2S1_comp,
  depth = Strat_depth, age = SH20_mean_age,
  age_max = SH20_max_age_95CI,
  age_min = SH20_min_age_95CI
) 

plot(HERPT2S1_adm_plot)

# save plot to file 
pdf(file = "Outputs/Figures/HERPT2S1/Fig1.6_age_depth_model.pdf",
    width = 6, # The width of the plot in inches
    height = 6) # The height of the plot in inches
HERPT2S1_adm_plot <- age_depth_model(
  HERPT2S1_comp,
  depth = Strat_depth, age = SH20_mean_age,
  age_max = SH20_max_age_95CI,
  age_min = SH20_min_age_95CI
) 
plot(HERPT2S1_adm_plot)
dev.off()

# Composite plots for paper ------------------------------------------

HERPT2S1_plot <- HERPT2S1_comp_long2 %>% 
  mutate(param = fct_relevel(param, "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"WMAR_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "WMAR_SAT" = "WMAR")) %>%
  ggplot(aes(x = value, y = Strat_depth), size = 0.01) +
  geom_lineh(lwd = 0.3) +
  geom_point(size = 0.01) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HERPT2S1 [SAT]") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "WMAR" = "g cm-2  a-1"))
HERPT2S1_plot

# OPTIONAL
# add hiatus (red) horizontal line and other points of interest (green) - or comment out

HERPT2S1_plot1 <- HERPT2S1_plot

# OR

#HERPT2S1_plot1 <- HERPT2S1_plot +
#geom_hline(yintercept = c(76), col = "red", lty = 1, alpha = 1, lwd = 0.5)
#HERPT2S1_plot1

# Add age on secondary age y axis
HERPT2S1_adm <- age_depth_model(
  HERPT2S1_comp, 
  depth = Strat_depth, age = SH20_mean_age)

HERPT2S1_plot_age <- HERPT2S1_plot1 + 
  scale_y_depth_age(HERPT2S1_adm, age_name = "Age (cal a BP)")
HERPT2S1_plot_age

# CONISS
# Uses Den1_SAT, DCMS_SAT only in _com_long3, but plotted with MS1_SAT and WMAR
# CONISS using nested_chclust_coniss() and the rioja package for constrained cluster analysis 
# By default, this calculates the number of significant zones based on a broken stick simulation
# TO DO <- add number of groups and their depth/ages

HERPT2S1_coniss <- HERPT2S1_comp_long3 %>%
  nested_data(qualifiers = c(SH20_mean_age, Strat_depth), 
              key = param, value = value, trans = scale) %>%
  nested_chclust_coniss()

HERPT2S1_plot_coniss <- HERPT2S1_plot_age +
  layer_dendrogram(HERPT2S1_coniss, aes(y = Strat_depth), 
                   param = "CONISS", size = 0.2, lwd = 0.2) +
  layer_zone_boundaries(HERPT2S1_coniss, aes(y = Strat_depth), 
                        size = 0.3, lwd = 0.3, colour = "red")
HERPT2S1_plot_coniss

ggarrange(HERPT2S1_plot_age, HERPT2S1_plot_coniss, nrow = 2)
ggsave("Outputs/Figures/HERPT2S1/Fig2.1_comp_depth_coniss.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")

# Plot versus age
HERPT2S1_age <- HERPT2S1_comp_long2 %>% 
  mutate(param = fct_relevel(param, "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"WMAR_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "WMAR_SAT" = "WMAR")) %>%
  ggplot(aes(x = SH20_mean_age, y = value), size = 0.01) +
  geom_line(lwd = 0.3) +
  geom_point(size = 0.01) +
  #scale_y_reverse() +
  # add units to the graph headers
  facet_geochem_grid(
    vars(param),
    units = c("GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "WMAR" = "g cm-2  a-1")) +
  #facet_geochem_grid(vars(param)) +
  labs(x = "Age (cal a BP)", y = NULL)
HERPT2S1_age

# OPTIONAL
# shading and lines for HER14 - added hiatus as a shaded area - min-max = rounded to the nearest 10 years - not needed for other plots
#zone_data <- tibble(ymin = -Inf, ymax = Inf , xmin = 1800, xmax = 10790)

HERPT2S1_age_lines1 <- HERPT2S1_age +
  #geom_rect(
  #mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
  #data = zone_data, 
  #alpha = 0.2,
  #fill = "white",
  #inherit.aes = FALSE
  #) +
  # add lines of interest at 1 ka, 2 ka, 11.75 ka and 13.8 ka etc. - need to label these later on
  #geom_vline(xintercept = c(0, 1000, 2000, 3000, 7500, 10000, 13800), col = "blue", lty = 3, alpha = 0.7, lwd = 0.3) +
  geom_vline(xintercept = c(4200), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  #geom_vline(xintercept = c(4200), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  labs(title = "HERPT2S1 [SAT]")
HERPT2S1_age_lines1

# Save final figure
ggarrange(HERPT2S1_plot_coniss, HERPT2S1_age_lines1, nrow = 2)
ggsave("Outputs/Figures/HERPT2S1/Fig2.2_depth_age_summary.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")

HERPT2S1_age_lines2 <- HERPT2S1_age +
  #geom_rect(
  #mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
  #data = zone_data, 
  #alpha = 0.2,
  #fill = "white",
  #inherit.aes = FALSE
  #) +
  # add lines of interest at 1 ka, 2 ka, 11.75 ka and 13.8 ka etc. - need to label these later on
  geom_vline(xintercept = c(0, 1000, 2000, 3000, 7500, 10000, 13800), col = "blue", lty = 3, alpha = 0.7, lwd = 0.3) +
  geom_vline(xintercept = c(4200, 8200, 11750, 16000), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  #geom_vline(xintercept = c(4200), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  labs(title = "HERPT2S1 [SAT]")
HERPT2S1_age_lines2

# Save final figure
ggarrange(HERPT2S1_plot_coniss, HERPT2S1_age_lines2, nrow = 2)
ggsave("Outputs/Figures/HERPT2S1/Fig2.3_depth_age_summary.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")

# -------------------------------------------------------------------------

# # HER44L plots --------------------------------------------------

# Overlap and Composite plots for XRAD mapping ---------------------

# All data - showing overlaps
HER44L_p1 <- HER44L_over_long %>% 
  mutate(param = fct_relevel(param, "PWVel_SAT", "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"Imp_SAT",	"FP_SAT",	"RES_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "Imp_SAT" = "Imp.", "FP_SAT" = "FP", "RES_SAT" = "Res.")) %>% 
  ggplot(aes(x = value, y = Strat_depth, colour = SECTION_ID, 
             fill = SECTION_ID, shape = SECTION_ID), size = 0.3
  ) +
  geom_lineh(lwd = 0.3) +
  #geom_point(size = 0.2) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HER44L [SAT]: Overlaps") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("PWave Vel." = "ms-1", "GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "Res." = "Ohm.m"))
HER44L_p1
ggsave("Outputs/Figures/HER44L/Fig1.1_ALL_over.pdf",
       height = c(10), width = c(20), dpi = 600, units = "cm")

# All data - composite colour by core for XRAD single core mapping
HER44L_p2 <- HER44L_comp_long1 %>% 
  mutate(param = fct_relevel(param, "PWVel_SAT", "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"Imp_SAT",	"FP_SAT",	"RES_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "Imp_SAT" = "Imp.", "FP_SAT" = "FP", "RES_SAT" = "Res.")) %>% 
  ggplot(aes(x = value, y = Strat_depth, colour = SECTION_ID, 
             fill = SECTION_ID, shape = SECTION_ID), size = 0.3
  ) +
  geom_lineh(lwd = 0.3) +
  #geom_point(size = 0.2) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HER44L [SAT]: Composite") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("PWave Vel." = "ms-1", "GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "10^-8 m3 kg-1", "Res." = "Ohm.m"))
HER44L_p2
ggsave("Outputs/Figures/HER44L/Fig1.2_ALL_comp.pdf",
       height = c(10), width = c(20), dpi = 600, units = "cm")

# All data - composite plotted as black
HER44L_p3 <- HER44L_comp_long1 %>% 
  mutate(param = fct_relevel(param, "PWVel_SAT", "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"Imp_SAT",	"FP_SAT",	"RES_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "Imp_SAT" = "Imp.", "FP_SAT" = "FP", "RES_SAT" = "Res.")) %>% 
  ggplot(aes(x = value, y = Strat_depth), size = 0.5) +
  geom_lineh(lwd = 0.5) +
  #geom_point(size = 0.2) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HER44L [SAT]: Composite") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("PWave Vel." = "ms-1", "GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "Res." = "Ohm.m"))
HER44L_p3
ggsave("Outputs/Figures/HER44L/Fig1.3_ALL_comp_black.pdf",
       height = c(10), width = c(16), dpi = 600, units = "cm")

# Multi-plots per page
ggarrange(HER44L_p1, HER44L_p2, nrow = 2)
ggsave("Outputs/Figures/HER44L/Fig1.4_over_comp_coresx2.pdf",
       height = c(20), width = c(20), dpi = 600, units = "cm")

ggarrange(HER44L_p1, HER44L_p2, HER44L_p3, nrow = 3)
ggsave("Outputs/Figures/HER44L/Fig1.5_over_comp_coresx3.pdf",
       height = c(30), width = c(20), dpi = 600, units = "cm")

# Age depth model - plot to check ---------------------------------------

#TO DO - add adm plot to core plots
HER44L_adm_plot <- age_depth_model(
  HER44L_comp,
  depth = Strat_depth, age = SH20_mean_age,
  age_max = SH20_max_age_95CI,
  age_min = SH20_min_age_95CI
) 

plot(HER44L_adm_plot)

# save plot to file 
pdf(file = "Outputs/Figures/HER44L/Fig1.6_age_depth_model.pdf",
    width = 6, # The width of the plot in inches
    height = 6) # The height of the plot in inches
HER44L_adm_plot <- age_depth_model(
  HER44L_comp,
  depth = Strat_depth, age = SH20_mean_age,
  age_max = SH20_max_age_95CI,
  age_min = SH20_min_age_95CI
) 
plot(HER44L_adm_plot)
dev.off()

# Composite plots for paper ------------------------------------------

HER44L_plot <- HER44L_comp_long2 %>% 
  mutate(param = fct_relevel(param, "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"WMAR_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "WMAR_SAT" = "WMAR")) %>%
  ggplot(aes(x = value, y = Strat_depth), size = 0.01) +
  geom_lineh(lwd = 0.3) +
  geom_point(size = 0.01) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HER44L [SAT]") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "WMAR" = "g cm-2  a-1"))
HER44L_plot

# OPTIONAL
# add hiatus (red) horizontal line and other points of interest (green) - or comment out

HER44L_plot1 <- HER44L_plot

# OR

#HER44L_plot1 <- HER44L_plot +
#geom_hline(yintercept = c(76), col = "red", lty = 1, alpha = 1, lwd = 0.5)
#HER44L_plot1

# Add age on secondary age y axis
HER44L_adm <- age_depth_model(
  HER44L_comp, 
  depth = Strat_depth, age = SH20_mean_age)

HER44L_plot_age <- HER44L_plot1 + 
  scale_y_depth_age(HER44L_adm, age_name = "Age (cal a BP)")
HER44L_plot_age

# CONISS
# Uses Den1_SAT, DCMS_SAT only in _com_long3, but plotted with MS1_SAT and WMAR
# CONISS using nested_chclust_coniss() and the rioja package for constrained cluster analysis 
# By default, this calculates the number of significant zones based on a broken stick simulation
# TO DO <- add number of groups and their depth/ages

HER44L_coniss <- HER44L_comp_long3 %>%
  nested_data(qualifiers = c(SH20_mean_age, Strat_depth), 
              key = param, value = value, trans = scale) %>%
  nested_chclust_coniss()

HER44L_plot_coniss <- HER44L_plot_age +
  layer_dendrogram(HER44L_coniss, aes(y = Strat_depth), 
                   param = "CONISS", size = 0.2, lwd = 0.2) +
  layer_zone_boundaries(HER44L_coniss, aes(y = Strat_depth), 
                        size = 0.3, lwd = 0.3, colour = "red")
HER44L_plot_coniss

ggarrange(HER44L_plot_age, HER44L_plot_coniss, nrow = 2)
ggsave("Outputs/Figures/HER44L/Fig2.1_comp_depth_coniss.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")

# Plot versus age
HER44L_age <- HER44L_comp_long2 %>% 
  mutate(param = fct_relevel(param, "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"WMAR_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "WMAR_SAT" = "WMAR")) %>%
  ggplot(aes(x = SH20_mean_age, y = value), size = 0.01) +
  geom_line(lwd = 0.3) +
  geom_point(size = 0.01) +
  #scale_y_reverse() +
  # add units to the graph headers
  facet_geochem_grid(
    vars(param),
    units = c("GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "WMAR" = "g cm-2  a-1")) +
  #facet_geochem_grid(vars(param)) +
  labs(x = "Age (cal a BP)", y = NULL)
HER44L_age

# OPTIONAL
# shading and lines for HER14 - added hiatus as a shaded area - min-max = rounded to the nearest 10 years - not needed for other plots
#zone_data <- tibble(ymin = -Inf, ymax = Inf , xmin = 1800, xmax = 10790)

HER44L_age_lines1 <- HER44L_age +
  #geom_rect(
  #mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
  #data = zone_data, 
  #alpha = 0.2,
  #fill = "white",
  #inherit.aes = FALSE
  #) +
  # add lines of interest at 1 ka, 2 ka, 11.75 ka and 13.8 ka etc. - need to label these later on
  #geom_vline(xintercept = c(0, 1000, 2000, 3000, 7500, 10000, 13800), col = "blue", lty = 3, alpha = 0.7, lwd = 0.3) +
  geom_vline(xintercept = c(4200, 8200), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  #geom_vline(xintercept = c(4200), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  labs(title = "HER44L [SAT]")
HER44L_age_lines1

# Save final figure
ggarrange(HER44L_plot_coniss, HER44L_age_lines1, nrow = 2)
ggsave("Outputs/Figures/HER44L/Fig2.2_depth_age_summary.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")

HER44L_age_lines2 <- HER44L_age +
  #geom_rect(
  #mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
  #data = zone_data, 
  #alpha = 0.2,
  #fill = "white",
  #inherit.aes = FALSE
  #) +
  # add lines of interest at 1 ka, 2 ka, 11.75 ka and 13.8 ka etc. - need to label these later on
  geom_vline(xintercept = c(0, 1000, 2000, 3000, 7500, 10000, 13800), col = "blue", lty = 3, alpha = 0.7, lwd = 0.3) +
  geom_vline(xintercept = c(4200, 8200, 11750, 16000), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  #geom_vline(xintercept = c(4200), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  labs(title = "HER44L [SAT]")
HER44L_age_lines2

# Save final figure
ggarrange(HER44L_plot_coniss, HER44L_age_lines2, nrow = 2)
ggsave("Outputs/Figures/HER44L/Fig2.3_depth_age_summary.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")


# -------------------------------------------------------------------------

## HER49L plots --------------------------------------------------

# Overlap and Composite plots for XRAD mapping ---------------------

# All data - showing overlaps
HER49L_p1 <- HER49L_over_long %>% 
  mutate(param = fct_relevel(param, "PWVel_SAT", "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"Imp_SAT",	"FP_SAT",	"RES_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "Imp_SAT" = "Imp.", "FP_SAT" = "FP", "RES_SAT" = "Res.")) %>% 
  ggplot(aes(x = value, y = Strat_depth, colour = SECTION_ID, 
             fill = SECTION_ID, shape = SECTION_ID), size = 0.3
  ) +
  geom_lineh(lwd = 0.3) +
  #geom_point(size = 0.2) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HER49L [SAT]: Overlaps") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("PWave Vel." = "ms-1", "GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "Res." = "Ohm.m"))
HER49L_p1
ggsave("Outputs/Figures/HER49L/Fig1.1_ALL_over.pdf",
       height = c(10), width = c(20), dpi = 600, units = "cm")

# All data - composite colour by core for XRAD single core mapping
HER49L_p2 <- HER49L_comp_long1 %>% 
  mutate(param = fct_relevel(param, "PWVel_SAT", "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"Imp_SAT",	"FP_SAT",	"RES_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "Imp_SAT" = "Imp.", "FP_SAT" = "FP", "RES_SAT" = "Res.")) %>% 
  ggplot(aes(x = value, y = Strat_depth, colour = SECTION_ID, 
             fill = SECTION_ID, shape = SECTION_ID), size = 0.3
  ) +
  geom_lineh(lwd = 0.3) +
  #geom_point(size = 0.2) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HER49L [SAT]: Composite") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("PWave Vel." = "ms-1", "GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "10^-8 m3 kg-1", "Res." = "Ohm.m"))
HER49L_p2
ggsave("Outputs/Figures/HER49L/Fig1.2_ALL_comp.pdf",
       height = c(10), width = c(20), dpi = 600, units = "cm")

# All data - composite plotted as black
HER49L_p3 <- HER49L_comp_long1 %>% 
  mutate(param = fct_relevel(param, "PWVel_SAT", "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"Imp_SAT",	"FP_SAT",	"RES_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "Imp_SAT" = "Imp.", "FP_SAT" = "FP", "RES_SAT" = "Res.")) %>% 
  ggplot(aes(x = value, y = Strat_depth), size = 0.3) +
  geom_lineh(lwd = 0.3) +
  #geom_point(size = 0.2) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HER49L [SAT]: Composite") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("PWave Vel." = "ms-1", "GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "Res." = "Ohm.m"))
HER49L_p3
ggsave("Outputs/Figures/HER49L/Fig1.3_ALL_comp_black.pdf",
       height = c(10), width = c(16), dpi = 600, units = "cm")

# Multi-plots per page
ggarrange(HER49L_p1, HER49L_p2, nrow = 2)
ggsave("Outputs/Figures/HER49L/Fig1.4_over_comp_coresx2.pdf",
       height = c(20), width = c(20), dpi = 600, units = "cm")

ggarrange(HER49L_p1, HER49L_p2, HER49L_p3, nrow = 3)
ggsave("Outputs/Figures/HER49L/Fig1.5_over_comp_coresx3.pdf",
       height = c(30), width = c(20), dpi = 600, units = "cm")

# Age depth model - plot to check ---------------------------------------

#TO DO - add adm plot to core plots
HER49L_adm_plot <- age_depth_model(
  HER49L_comp,
  depth = Strat_depth, age = SH20_mean_age,
  age_max = SH20_max_age_95CI,
  age_min = SH20_min_age_95CI
) 

plot(HER49L_adm_plot)

# save plot to file 
pdf(file = "Outputs/Figures/HER49L/Fig1.6_age_depth_model.pdf",
    width = 6, # The width of the plot in inches
    height = 6) # The height of the plot in inches
HER49L_adm_plot <- age_depth_model(
  HER49L_comp,
  depth = Strat_depth, age = SH20_mean_age,
  age_max = SH20_max_age_95CI,
  age_min = SH20_min_age_95CI
) 
plot(HER49L_adm_plot)
dev.off()

# Composite plots for paper ------------------------------------------

HER49L_plot <- HER49L_comp_long2 %>% 
  mutate(param = fct_relevel(param, "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"WMAR_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "WMAR_SAT" = "WMAR")) %>%
  ggplot(aes(x = value, y = Strat_depth), size = 0.01) +
  geom_lineh(lwd = 0.3) +
  geom_point(size = 0.01) +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HER49L [SAT]") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "WMAR" = "g cm-2  a-1"))
HER49L_plot

# OPTIONAL
# add hiatus (red) horizontal line and other points of interest (green) - or comment out

HER49L_plot1 <- HER49L_plot

# OR

#HER49L_plot1 <- HER49L_plot +
#geom_hline(yintercept = c(76), col = "red", lty = 1, alpha = 1, lwd = 0.5)
#HER49L_plot1

# Add age on secondary age y axis
HER49L_adm <- age_depth_model(
  HER49L_comp, 
  depth = Strat_depth, age = SH20_mean_age)

HER49L_plot_age <- HER49L_plot1 + 
  scale_y_depth_age(HER49L_adm, age_name = "Age (cal a BP)")
HER49L_plot_age

# CONISS
# Uses Den1_SAT, DCMS_SAT only in _com_long3, but plotted with MS1_SAT and WMAR
# CONISS using nested_chclust_coniss() and the rioja package for constrained cluster analysis 
# By default, this calculates the number of significant zones based on a broken stick simulation
# TO DO <- add number of groups and their depth/ages

HER49L_coniss <- HER49L_comp_long3 %>%
  nested_data(qualifiers = c(SH20_mean_age, Strat_depth), 
              key = param, value = value, trans = scale) %>%
  nested_chclust_coniss()

HER49L_plot_coniss <- HER49L_plot_age +
  layer_dendrogram(HER49L_coniss, aes(y = Strat_depth), 
                   param = "CONISS", size = 0.2, lwd = 0.2) +
  layer_zone_boundaries(HER49L_coniss, aes(y = Strat_depth), 
                        size = 0.3, lwd = 0.3, colour = "red")
HER49L_plot_coniss

ggarrange(HER49L_plot_age, HER49L_plot_coniss, nrow = 2)
ggsave("Outputs/Figures/HER49L/Fig2.1_comp_depth_coniss.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")

# Plot versus age
HER49L_age <- HER49L_comp_long2 %>% 
  mutate(param = fct_relevel(param, "Den1_SAT",	"MS1_SAT",	
                             "DCMS1_SAT",	"WMAR_SAT")) %>% 
  mutate(param = recode(param, "PWVel_SAT" = "PWave Vel.", "Den1_SAT" = "GRD", "MS1_SAT" = "MS", 
                        "DCMS1_SAT" = "DCMS", "WMAR_SAT" = "WMAR")) %>%
  ggplot(aes(x = SH20_mean_age, y = value), size = 0.01) +
  geom_line(lwd = 0.3) +
  geom_point(size = 0.01) +
  #scale_y_reverse() +
  # add units to the graph headers
  facet_geochem_grid(
    vars(param),
    units = c("GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^-8 m3 kg-1", "WMAR" = "g cm-2  a-1")) +
  #facet_geochem_grid(vars(param)) +
  labs(x = "Age (cal a BP)", y = NULL)
HER49L_age

# OPTIONAL
# shading and lines for HER14 - added hiatus as a shaded area - min-max = rounded to the nearest 10 years - not needed for other plots
#zone_data <- tibble(ymin = -Inf, ymax = Inf , xmin = 1800, xmax = 10790)

HER49L_age_lines1 <- HER49L_age +
  #geom_rect(
  #mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
  #data = zone_data, 
  #alpha = 0.2,
  #fill = "white",
  #inherit.aes = FALSE
  #) +
  # add lines of interest at 1 ka, 2 ka, 11.75 ka and 13.8 ka etc. - need to label these later on
  #geom_vline(xintercept = c(0, 1000, 2000, 3000, 7500, 10000, 13800), col = "blue", lty = 3, alpha = 0.7, lwd = 0.3) +
  geom_vline(xintercept = c(4200, 8200), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  #geom_vline(xintercept = c(4200), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  labs(title = "HER49L [SAT]")
HER49L_age_lines1

# Save final figure
ggarrange(HER49L_plot_coniss, HER49L_age_lines1, nrow = 2)
ggsave("Outputs/Figures/HER49L/Fig2.2_depth_age_summary.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")

HER49L_age_lines2 <- HER49L_age +
  #geom_rect(
  #mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
  #data = zone_data, 
  #alpha = 0.2,
  #fill = "white",
  #inherit.aes = FALSE
  #) +
  # add lines of interest at 1 ka, 2 ka, 11.75 ka and 13.8 ka etc. - need to label these later on
  geom_vline(xintercept = c(0, 1000, 2000, 3000, 7500, 10000, 13800), col = "blue", lty = 3, alpha = 0.7, lwd = 0.3) +
  geom_vline(xintercept = c(4200, 8200, 11750, 16000), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  #geom_vline(xintercept = c(4200), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  labs(title = "HER49L [SAT]")
HER49L_age_lines2

# Save final figure
ggarrange(HER49L_plot_coniss, HER49L_age_lines2, nrow = 2)
ggsave("Outputs/Figures/HER42L/Fig2.3_depth_age_summary.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")

# -------------------------------------------------------------------------

# COMBINED SITE PLOTS -------------------------------------------------------------------------

# Age ---------------------------------------------------------------------

# 10ka+ sites
ggarrange(HER14L_age_lines, HER24L_age_Lines2, HER44L_age_Lines2, HER49L_age_Lines2, nrow = 4)
ggsave("Outputs/Figures/Combined_sites/Fig_10ka.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")


# 1.1 NORTHERN_LAKES: HER14L + HER24L + HER34L ------------------------------

theme_set(theme_classic(base_size=12) + theme(
  plot.title = element_text(color="black", size=12, face="bold"),
  axis.title.x = element_text(color="black", size=12),
  axis.title.y = element_text(color="black", size=12)
))

# GRD plot
HER14L_GRD <- ggplot() + 
  geom_line(data = HER14L_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "black", lwd = 0.5) +
  geom_point(data = HER14L_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "black", size = 0.5) +
  #scale_color_manual(values = colours) +
  labs(title = "Northern Coastal Lakes: GRD",
       subtitle = "HER14L (black), HER24L (red), HER34L (blue) GRD") + 
  xlab("Age (cal a BP") +
  ylab(as.expression(expression(paste("GRDensity [SAT] " (kg/m ^ 3 )) ))) +
  #coord_flip() +
  ylim(0.5,2.5)
HER14L_GRD

# add hiatus as white shading and plot HER24L and HER34L on top
zone_data <- tibble(ymin = -Inf, ymax = Inf , xmin = 1800, xmax = 10790)

HER14L_24L_34L_GRD  <- HER14L_GRD +
  geom_rect(
    mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
    data = zone_data, 
    #alpha = 0.2,
    fill = "white",
    inherit.aes = FALSE
  ) +
  geom_line(data = HER24L_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "red", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER24L_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "red", size = 0.5, alpha = 0.5) +
  geom_line(data = HER34L_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "blue", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER34L_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "blue", size = 0.5, alpha = 0.5) +
  # add a red dotted line at point of interests (mid HER24L-3B core jump to higher MS and GRD at ~7.4 ka)
  geom_vline(xintercept = c(7400, 13800), col = "grey", lty = 2, alpha = 0.7, lwd = 0.75)
HER14L_24L_34L_GRD

# MS plot
HER14L_MS <- ggplot() + 
  geom_line(data = HER14L_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "black", lwd = 0.5) +
  geom_point(data = HER14L_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "black", size = 0.5) +
  #scale_color_manual(values = colours) +
  labs(title = "Northern Coastal Lakes: MS",
       subtitle = "HER14L (black), HER24L (red), HER34L (blue)") + 
  xlab("Age (cal a BP") +
  ylab(as.expression(expression(paste("MS [SAT] SI " (x10 ^ -5)) ))) +
  #coord_flip() +
  ylim(-5,15)
HER14L_MS

# add hiatus as white shading and plot HER24L and HER34L on top
zone_data <- tibble(ymin = -Inf, ymax = Inf , xmin = 1800, xmax = 10790)

HER14L_24L_34L_MS  <- HER14L_MS +
  geom_rect(
    mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
    data = zone_data, 
    #alpha = 0.2,
    fill = "white",
    inherit.aes = FALSE
  ) +
  geom_line(data = HER24L_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "red", lwd = 0.5, alpha = 1) +
  geom_point(data = HER24L_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "red", size = 0.5, alpha = 1) +
  geom_line(data = HER34L_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "blue", lwd = 0.5, alpha = 1) +
  geom_point(data = HER34L_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "blue", size = 0.5, alpha = 15) +
  # add a red dotted line at point of interests (mid HER24L-3B core jump to higher MS and GRD at ~7.5 ka)
  geom_vline(xintercept = c(7500, 13800), col = "grey", lty = 2, alpha = 0.7, lwd = 0.75)
HER14L_24L_34L_MS

# DCMS plot
HER14L_DCMS <- ggplot() + 
  geom_line(data = HER14L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "black", lwd = 0.5) +
  geom_point(data = HER14L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "black", size = 0.5) +
  #scale_color_manual(values = colours) +
  labs(title = "Northern Coastal Lakes: DCMS",
       subtitle = "HER14L (black), HER24L (red), HER34L (blue)") + 
  xlab("Age (cal a BP") +
  ylab(as.expression(expression(paste("DCMS [SAT]" (x10 ^ -8 ~m ^ 3/kg)) ))) +
  #coord_flip() +
  ylim(-5,15)
HER14L_DCMS

# add hiatus as white shading and plot HER24L and HER34L on top
zone_data <- tibble(ymin = -Inf, ymax = Inf , xmin = 1800, xmax = 10790)

HER14L_24L_34L_DCMS  <- HER14L_DCMS +
  geom_rect(
    mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
    data = zone_data, 
    #alpha = 0.2,
    fill = "white",
    inherit.aes = FALSE
  ) +
  geom_line(data = HER24L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "red", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER24L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "red", size = 0.5, alpha = 0.5) +
  geom_line(data = HER34L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "blue", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER34L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "blue", size = 0.5, alpha = 0.5) +
  # add a red dotted line at point of interests (mid HER24L-3B core jump to higher MS and GRD at ~7.5 ka)
  geom_vline(xintercept = c(7500, 13800), col = "grey", lty = 2, alpha = 0.7, lwd = 0.75)
HER14L_24L_34L_DCMS

ggarrange(HER14L_24L_34L_GRD, HER14L_24L_34L_MS, HER14L_24L_34L_DCMS, nrow = 3)
ggsave("Outputs/Figures/Combined_sites/Fig1.1_NorthernLakes_HER14_24_34.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")

# 1.2 NORTHERN + SOUTHERN 0-10 ka SITES: HER14L + HER24L + HER44L >7.5 ka ------------------------------

theme_set(theme_classic(base_size=12) + theme(
  plot.title = element_text(color="black", size=12, face="bold"),
  axis.title.x = element_text(color="black", size=12),
  axis.title.y = element_text(color="black", size=12)
))

#filter HER44L to remove data <6ka - as there's no age data less then 6ka - pot HER14/24/44L>6ka
HER44L_comp_6ka <- filter(HER44L_comp, SH20_mean_age >7500)

# GRD plot
HER14L_GRD <- ggplot() + 
  geom_line(data = HER14L_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "black", lwd = 0.5) +
  geom_point(data = HER14L_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "black", size = 0.5) +
  #scale_color_manual(values = colours) +
  labs(title = "Northern & Southern Coastal Lakes 0-10 ka records: GRD",
       subtitle = "HER14L (black), HER24L (red), HER44L >7.5 ka (blue) GRD") + 
  xlab("Age (cal a BP") +
  ylab(as.expression(expression(paste("GRDensity [SAT] " (kg/m ^ 3 )) ))) +
  #coord_flip() +
  ylim(0.5,2.5)
HER14L_GRD

# add hiatus as white shading and plot HER24L and HER44L >7.5 ka on top
zone_data <- tibble(ymin = -Inf, ymax = Inf , xmin = 1800, xmax = 10790)

HER14L_24L_44L_GRD  <- HER14L_GRD +
  geom_rect(
    mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
    data = zone_data, 
    #alpha = 0.2,
    fill = "white",
    inherit.aes = FALSE
  ) +
  geom_line(data = HER24L_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "red", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER24L_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "red", size = 0.5, alpha = 0.5) +
  geom_line(data = HER44L_comp_6ka, aes(x=SH20_mean_age, y=Den1_SAT), colour = "blue", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER44L_comp_6ka, aes(x=SH20_mean_age, y=Den1_SAT),colour = "blue", size = 0.5, alpha = 0.5) +
  # add a red dotted line at point of interests (mid HER24L-3B core jump to higher MS and GRD at ~7.4 ka)
  geom_vline(xintercept = c(7400, 13800), col = "grey", lty = 2, alpha = 0.7, lwd = 0.75)
HER14L_24L_44L_GRD

# MS plot
HER14L_MS <- ggplot() + 
  geom_line(data = HER14L_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "black", lwd = 0.5) +
  geom_point(data = HER14L_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "black", size = 0.5) +
  #scale_color_manual(values = colours) +
  labs(title = "Northern & Southern Coastal Lakes 0-10 ka records: MS",
       subtitle = "HER14L (black), HER24L (red), HER44L >7.5 ka (blue)") + 
  xlab("Age (cal a BP") +
  ylab(as.expression(expression(paste("MS [SAT] SI " (x10 ^ -5)) ))) +
  #coord_flip() +
  ylim(-5,15)
HER14L_MS

# add hiatus as white shading and plot HER24L and HER44L on top
zone_data <- tibble(ymin = -Inf, ymax = Inf , xmin = 1800, xmax = 10790)

HER14L_24L_44L_MS  <- HER14L_MS +
  geom_rect(
    mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
    data = zone_data, 
    #alpha = 0.2,
    fill = "white",
    inherit.aes = FALSE
  ) +
  geom_line(data = HER24L_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "red", lwd = 0.5, alpha = 1) +
  geom_point(data = HER24L_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "red", size = 0.5, alpha = 1) +
  geom_line(data = HER44L_comp_6ka, aes(x=SH20_mean_age, y=MS1_SAT), colour = "blue", lwd = 0.5, alpha = 1) +
  geom_point(data = HER44L_comp_6ka, aes(x=SH20_mean_age, y=MS1_SAT),colour = "blue", size = 0.5, alpha = 15) +
  # add a red dotted line at point of interests (mid HER24L-3B core jump to higher MS and GRD at ~7.5 ka)
  geom_vline(xintercept = c(7500, 13800), col = "grey", lty = 2, alpha = 0.7, lwd = 0.75)
HER14L_24L_44L_MS

# DCMS plot
HER14L_DCMS <- ggplot() + 
  geom_line(data = HER14L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "black", lwd = 0.5) +
  geom_point(data = HER14L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "black", size = 0.5) +
  #scale_color_manual(values = colours) +
  labs(title = "Northern & Southern Coastal Lakes 0-10 ka records: DCMS",
       subtitle = "HER14L (black), HER24L (red), HER44L >7.5 (blue)") + 
  xlab("Age (cal a BP") +
  ylab(as.expression(expression(paste("DCMS [SAT]" (x10 ^ -8 ~m ^ 3/kg)) ))) +
  #coord_flip() +
  ylim(-5,15)
HER14L_DCMS

# add hiatus as white shading and plot HER24L and HER44L on top
zone_data <- tibble(ymin = -Inf, ymax = Inf , xmin = 1800, xmax = 10790)

HER14L_24L_44L_DCMS  <- HER14L_DCMS +
  geom_rect(
    mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
    data = zone_data, 
    #alpha = 0.2,
    fill = "white",
    inherit.aes = FALSE
  ) +
  geom_line(data = HER24L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "red", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER24L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "red", size = 0.5, alpha = 0.5) +
  geom_line(data = HER44L_comp_6ka, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "blue", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER44L_comp_6ka, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "blue", size = 0.5, alpha = 0.5) +
  # add a red dotted line at point of interests (mid HER24L-3B core jump to higher MS and GRD at ~7.5 ka)
  geom_vline(xintercept = c(7500, 13800), col = "grey", lty = 2, alpha = 0.7, lwd = 0.75)
HER14L_24L_44L_DCMS

ggarrange(HER14L_24L_44L_GRD, HER14L_24L_44L_MS, HER14L_24L_44L_DCMS, nrow = 3)
ggsave("Outputs/Figures/Combined_sites/Fig1.2_0-10ka_NorthernLakes_HER14_24_44.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")


# 1.3 NORTHERN + SOUTHERN 5-10 ka LAKES: HER14L + HER24L + HER42PB + HER44L + HER49L------------------------------

theme_set(theme_classic(base_size=12) + theme(
  plot.title = element_text(color="black", size=12, face="bold"),
  axis.title.x = element_text(color="black", size=12),
  axis.title.y = element_text(color="black", size=12)
))

#filter HER44L to remove data <6ka - as there's no age data less then 6ka - pot HER14/24/44L>6ka
HER44L_comp_6ka <- filter(HER44L_comp, SH20_mean_age >7500)

# GRD plot
HER14L_GRD <- ggplot() + 
  geom_line(data = HER14L_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "black", lwd = 0.5) +
  geom_point(data = HER14L_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "black", size = 0.5) +
  #scale_color_manual(values = colours) +
  labs(title = "Northern & Southern Coastal Lakes 0-10 ka records + 0-15 ka HER42L: GRD",
       subtitle = "HER14L (black), HER24L (red), HER42L (orange), HER44L >7.5 ka (blue), HER49L (grey)") + 
  xlab("Age (cal a BP") +
  ylab(as.expression(expression(paste("GRDensity [SAT] " (kg/m ^ 3 )) ))) +
  #coord_flip() +
  ylim(0.5,2.5)
HER14L_GRD

# add hiatus as white shading and plot HER24L and HER44L >9.4 ka on top
zone_data <- tibble(ymin = -Inf, ymax = Inf , xmin = 1800, xmax = 10790)

HER14L_24L_44L_42L_49L_GRD  <- HER14L_GRD +
  geom_rect(
    mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
    data = zone_data, 
    #alpha = 0.2,
    fill = "white",
    inherit.aes = FALSE
  ) +
  geom_line(data = HER24L_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "red", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER24L_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "red", size = 0.5, alpha = 0.5) +
  geom_line(data = HER44L_comp_6ka, aes(x=SH20_mean_age, y=Den1_SAT), colour = "blue", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER44L_comp_6ka, aes(x=SH20_mean_age, y=Den1_SAT),colour = "blue", size = 0.5, alpha = 0.5) +
  geom_line(data = HER42L_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "orange", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER42L_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "orange", size = 0.5, alpha = 0.5) +
  geom_line(data = HER49L_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "grey", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER49L_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "grey", size = 0.5, alpha = 0.5) +
  # add a red dotted line at point of interests (mid HER24L-3B core jump to higher MS and GRD at ~7.4 ka)
  geom_vline(xintercept = c(7400, 13800), col = "grey", lty = 2, alpha = 0.7, lwd = 0.75)
HER14L_24L_44L_42L_49L_GRD

# MS plot
HER14L_MS <- ggplot() + 
  geom_line(data = HER14L_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "black", lwd = 0.5) +
  geom_point(data = HER14L_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "black", size = 0.5) +
  #scale_color_manual(values = colours) +
  labs(title = "Northern & Southern Coastal Lakes 0-10 ka records + 0-15 ka HER42L: MS",
       subtitle = "HER14L (black), HER24L (red), HER42L (orange), HER44L >7.5 ka (blue), HER49L (grey)") + 
  xlab("Age (cal a BP") +
  ylab(as.expression(expression(paste("MS [SAT] SI " (x10 ^ -5)) ))) +
  #coord_flip() +
  ylim(-5,15)
HER14L_MS

# add hiatus as white shading and plot HER24L and HER44L on top
zone_data <- tibble(ymin = -Inf, ymax = Inf , xmin = 1800, xmax = 10790)

HER14L_24L_44L_42L_49L_MS  <- HER14L_MS +
  geom_rect(
    mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
    data = zone_data, 
    #alpha = 0.2,
    fill = "white",
    inherit.aes = FALSE
  ) +
  geom_line(data = HER24L_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "red", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER24L_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "red", size = 0.5, alpha = 0.5) +
  geom_line(data = HER44L_comp_6ka, aes(x=SH20_mean_age, y=MS1_SAT), colour = "blue", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER44L_comp_6ka, aes(x=SH20_mean_age, y=MS1_SAT),colour = "blue", size = 0.5, alpha = 0.5) +
  geom_line(data = HER42L_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "orange", lwd = 0.5, alpha = 0.3) +
  geom_point(data = HER42L_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "orange", size = 0.5, alpha = 0.3) +
  geom_line(data = HER49L_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "grey", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER49L_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "grey", size = 0.5, alpha = 0.5) +
  # add a red dotted line at point of interests (mid HER24L-3B core jump to higher MS and GRD at ~7.5 ka)
  geom_vline(xintercept = c(7500, 13800), col = "grey", lty = 2, alpha = 0.7, lwd = 0.75)
HER14L_24L_44L_42L_49L_MS

# DCMS plot
HER14L_DCMS <- ggplot() + 
  geom_line(data = HER14L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "black", lwd = 0.5) +
  geom_point(data = HER14L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "black", size = 0.5) +
  #scale_color_manual(values = colours) +
  labs(title = "Northern & Southern Coastal Lakes 0-10 ka records + 0-15 ka HER42L: DCMS",
       subtitle = "HER14L (black), HER24L (red), HER42L (orange), HER44L >7.5 (blue), HER49L (grey)") + 
  xlab("Age (cal a BP") +
  ylab(as.expression(expression(paste("DCMS [SAT]" (x10 ^ -8 ~m ^ 3/kg)) ))) +
  #coord_flip() +
  ylim(-5,15)
HER14L_DCMS

# add hiatus as white shading and plot HER24L and HER44L on top
zone_data <- tibble(ymin = -Inf, ymax = Inf , xmin = 1800, xmax = 10790)

HER14L_24L_44L_42L_49L_DCMS  <- HER14L_DCMS +
  geom_rect(
    mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
    data = zone_data, 
    #alpha = 0.2,
    fill = "white",
    inherit.aes = FALSE
  ) +
  geom_line(data = HER24L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "red", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER24L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "red", size = 0.5, alpha = 0.5) +
  geom_line(data = HER44L_comp_6ka, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "blue", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER44L_comp_6ka, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "blue", size = 0.5, alpha = 0.5) +
  geom_line(data = HER42L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "orange", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER42L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "orange", size = 0.5, alpha = 0.5) +
  geom_line(data = HER49L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "grey", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER49L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "grey", size = 0.5, alpha = 0.5) +
  # add a red dotted line at point of interests (mid HER24L-3B core jump to higher MS and GRD at ~7.5 ka)
  geom_vline(xintercept = c(7500, 13800), col = "grey", lty = 2, alpha = 0.7, lwd = 0.75)
HER14L_24L_44L_42L_49L_DCMS

ggarrange(HER14L_24L_44L_42L_49L_GRD, HER14L_24L_44L_42L_49L_MS, HER14L_24L_44L_42L_49L_DCMS, nrow = 3)
ggsave("Outputs/Figures/Combined_sites/Fig1.3_0-10ka_Lakes_24L_42L_44L_49L.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")


# 1.4 SOUTHERN LAKES: HER42L + HER44L + HER49L  ------------------------------

theme_set(theme_classic(base_size=12) + theme(
  plot.title = element_text(color="black", size=12, face="bold"),
  axis.title.x = element_text(color="black", size=12),
  axis.title.y = element_text(color="black", size=12)
))

#filter HER44L to remove data <6ka - as there's no age data less then 6ka - pot HER14/24/44L>6ka
HER44L_comp_9ka <- filter(HER44L_comp, SH20_mean_age >7500)

# GRD plot
HER14L_GRD <- ggplot() + 
  geom_line(data = HER14L_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "white", lwd = 0.5) +
  geom_point(data = HER14L_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "white", size = 0.5) +
  #scale_color_manual(values = colours) +
  labs(title = "Southern Coastal Lakes: GRD",
       subtitle = "HER42L (orange), HER44L >7.5 ka (blue), HER49L (grey)") + 
  xlab("Age (cal a BP") +
  ylab(as.expression(expression(paste("GRDensity [SAT] " (kg/m ^ 3 )) ))) +
  #coord_flip() +
  ylim(0.5,2.5)
HER14L_GRD

# add hiatus as white shading and plot HER24L and HER44L >9.4 ka on top
zone_data <- tibble(ymin = -Inf, ymax = Inf , xmin = 1800, xmax = 10790)

HER14L_42L_44L_42L_49L_GRD  <- HER14L_GRD +
  geom_rect(
    mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
    data = zone_data, 
    #alpha = 0.2,
    fill = "white",
    inherit.aes = FALSE
  ) +
  #geom_line(data = HER24L_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "red", lwd = 0.5, alpha = 0.5) +
  #geom_point(data = HER24L_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "red", size = 0.5, alpha = 0.5) +
  geom_line(data = HER44L_comp_9ka, aes(x=SH20_mean_age, y=Den1_SAT), colour = "blue", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER44L_comp_9ka, aes(x=SH20_mean_age, y=Den1_SAT),colour = "blue", size = 0.5, alpha = 0.5) +
  geom_line(data = HER42L_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "orange", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER42L_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "orange", size = 0.5, alpha = 0.5) +
  geom_line(data = HER49L_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "grey", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER49L_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "grey", size = 0.5, alpha = 0.5) +
  # add a red dotted line at point of interests (mid HER24L-3B core jump to higher MS and GRD at ~7.4 ka)
  geom_vline(xintercept = c(7400, 13800), col = "grey", lty = 2, alpha = 0.7, lwd = 0.75)
HER14L_42L_44L_42L_49L_GRD

# MS plot
HER14L_MS <- ggplot() + 
  geom_line(data = HER14L_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "white", lwd = 0.5) +
  geom_point(data = HER14L_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "white", size = 0.5) +
  #scale_color_manual(values = colours) +
  labs(title = "Southern Coastal Lakes: MS",
       subtitle = "HER42L (orange), HER44L >7.5 ka (blue), HER49L (grey)") + 
  xlab("Age (cal a BP") +
  ylab(as.expression(expression(paste("MS [SAT] SI " (x10 ^ -5)) ))) +
  #coord_flip() +
  ylim(-5,15)
HER14L_MS

# add hiatus as white shading and plot HER24L and HER44L on top
zone_data <- tibble(ymin = -Inf, ymax = Inf , xmin = 1800, xmax = 10790)

HER14L_42L_44L_42L_49L_MS  <- HER14L_MS +
  geom_rect(
    mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
    data = zone_data, 
    #alpha = 0.2,
    fill = "white",
    inherit.aes = FALSE
  ) +
  #geom_line(data = HER24L_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "red", lwd = 0.5, alpha = 0.5) +
  #geom_point(data = HER24L_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "red", size = 0.5, alpha = 0.5) +
  geom_line(data = HER44L_comp_9ka, aes(x=SH20_mean_age, y=MS1_SAT), colour = "blue", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER44L_comp_9ka, aes(x=SH20_mean_age, y=MS1_SAT),colour = "blue", size = 0.5, alpha = 0.5) +
  geom_line(data = HER42L_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "orange", lwd = 0.5, alpha = 0.3) +
  geom_point(data = HER42L_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "orange", size = 0.5, alpha = 0.3) +
  geom_line(data = HER49L_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "grey", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER49L_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "grey", size = 0.5, alpha = 0.5) +
  # add a red dotted line at point of interests (mid HER24L-3B core jump to higher MS and GRD at ~7.5 ka)
  geom_vline(xintercept = c(7500, 13800), col = "grey", lty = 2, alpha = 0.7, lwd = 0.75)
HER14L_42L_44L_42L_49L_MS

# DCMS plot
HER14L_DCMS <- ggplot() + 
  geom_line(data = HER14L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "white", lwd = 0.5) +
  geom_point(data = HER14L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "white", size = 0.5) +
  #scale_color_manual(values = colours) +
  labs(title = "Southern Coastal Lakes: DCMS",
       subtitle = "HER42L (orange), HER44L >7.5 ka (blue), HER49L (grey)") + 
  xlab("Age (cal a BP") +
  ylab(as.expression(expression(paste("DCMS [SAT]" (x10 ^ -8 ~m ^ 3/kg)) ))) +
  #coord_flip() +
  ylim(-5,15)
HER14L_DCMS

# add hiatus as white shading and plot HER24L and HER44L on top
zone_data <- tibble(ymin = -Inf, ymax = Inf , xmin = 1800, xmax = 10790)

HER14L_42L_44L_42L_49L_DCMS  <- HER14L_DCMS +
  geom_rect(
    mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
    data = zone_data, 
    #alpha = 0.2,
    fill = "white",
    inherit.aes = FALSE
  ) +
  #geom_line(data = HER24L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "red", lwd = 0.5, alpha = 0.5) +
  #geom_point(data = HER24L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "red", size = 0.5, alpha = 0.5) +
  geom_line(data = HER44L_comp_9ka, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "blue", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER44L_comp_9ka, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "blue", size = 0.5, alpha = 0.5) +
  geom_line(data = HER42L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "orange", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER42L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "orange", size = 0.5, alpha = 0.5) +
  geom_line(data = HER49L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "grey", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER49L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "grey", size = 0.5, alpha = 0.5) +
  # add a red dotted line at point of interests (mid HER24L-3B core jump to higher MS and GRD at ~7.5 ka)
  geom_vline(xintercept = c(7500, 13800), col = "grey", lty = 2, alpha = 0.7, lwd = 0.75)
HER14L_42L_44L_42L_49L_DCMS

ggarrange(HER14L_24L_44L_42L_49L_GRD, HER14L_24L_44L_42L_49L_MS, HER14L_24L_44L_42L_49L_DCMS, nrow = 3)
ggsave("Outputs/Figures/Combined_sites/Fig1.4_SouthernLakes_HER42L_44L_49L.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")


# 1.5 SOUTHERN LAKES & PEAT BOGS: HER42L + HER44L + HER49L + HER42PB + PT2S1_SS1_XRF6_SEC42  ------------------------------ 

theme_set(theme_classic(base_size=12) + theme(
  plot.title = element_text(color="black", size=12, face="bold"),
  axis.title.x = element_text(color="black", size=12),
  axis.title.y = element_text(color="black", size=12)
))

#filter HER44L to remove data <6ka - as there's no age data less then 6ka - pot HER14/24/44L>6ka
HER44L_comp_9ka <- filter(HER44L_comp, SH20_mean_age >7500)

# GRD plot
HER14L_GRD <- ggplot() + 
  geom_line(data = HER14L_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "white", lwd = 0.5) +
  geom_point(data = HER14L_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "white", size = 0.5) +
  #scale_color_manual(values = colours) +
  labs(title = "Southern Coastal Lakes & Peat Bogs: GRD",
       subtitle = "HER42L (orange), HER42PB(black), HERPT2S1_SS1_SEC42(brown), HER44L >7.5 ka (blue), HER49L (grey)") + 
  xlab("Age (cal a BP") +
  ylab(as.expression(expression(paste("GRDensity [SAT] " (kg/m ^ 3 )) ))) +
  #coord_flip() +
  ylim(0.5,2.5)
HER14L_GRD

# add hiatus as white shading and plot HER24L and HER44L >9.4 ka on top
zone_data <- tibble(ymin = -Inf, ymax = Inf , xmin = 1800, xmax = 10790)

HER14L_42L_44L_42L_49L_42PB_PT2S1_GRD  <- HER14L_GRD +
  geom_rect(
    mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
    data = zone_data, 
    #alpha = 0.2,
    fill = "white",
    inherit.aes = FALSE
  ) +
  #geom_line(data = HER24L_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "red", lwd = 0.5, alpha = 0.5) +
  #geom_point(data = HER24L_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "red", size = 0.5, alpha = 0.5) +
  geom_line(data = HER44L_comp_9ka, aes(x=SH20_mean_age, y=Den1_SAT), colour = "blue", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER44L_comp_9ka, aes(x=SH20_mean_age, y=Den1_SAT),colour = "blue", size = 0.5, alpha = 0.5) +
  geom_line(data = HER42L_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "orange", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER42L_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "orange", size = 0.5, alpha = 0.5) +
  geom_line(data = HER49L_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "grey", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER49L_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "grey", size = 0.5, alpha = 0.5) +
  geom_line(data = HER42PB_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "black", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER42PB_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "black", size = 0.5, alpha = 0.5) +
  geom_line(data = HERPT2S1_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "brown", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HERPT2S1_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "brown", size = 0.5, alpha = 0.5) +
  # add a red dotted line at point of interests (mid HER24L-3B core jump to higher MS and GRD at ~7.4 ka)
  geom_vline(xintercept = c(7400, 13800), col = "grey", lty = 2, alpha = 0.7, lwd = 0.75)
HER14L_42L_44L_42L_49L_42PB_PT2S1_GRD

# MS plot
HER14L_MS <- ggplot() + 
  geom_line(data = HER14L_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "white", lwd = 0.5) +
  geom_point(data = HER14L_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "white", size = 0.5) +
  #scale_color_manual(values = colours) +
  labs(title = "Southern Coastal Lakes & Peat Bogs: MS",
       subtitle = "HER42L (orange), HER42PB(black), HERPT2S1_SS1_SEC42(brown), HER44L >7.5 ka (blue), HER49L (grey)") + 
  xlab("Age (cal a BP") +
  ylab(as.expression(expression(paste("MS [SAT] SI " (x10 ^ -5)) ))) +
  #coord_flip() +
  ylim(-5,15)
HER14L_MS

# add hiatus as white shading and plot HER24L and HER44L on top
zone_data <- tibble(ymin = -Inf, ymax = Inf , xmin = 1800, xmax = 10790)

HER14L_42L_44L_42L_49L_42PB_PT2S1_MS  <- HER14L_MS +
  geom_rect(
    mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
    data = zone_data, 
    #alpha = 0.2,
    fill = "white",
    inherit.aes = FALSE
  ) +
  #geom_line(data = HER24L_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "red", lwd = 0.5, alpha = 0.5) +
  #geom_point(data = HER24L_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "red", size = 0.5, alpha = 0.5) +
  geom_line(data = HER44L_comp_9ka, aes(x=SH20_mean_age, y=MS1_SAT), colour = "blue", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER44L_comp_9ka, aes(x=SH20_mean_age, y=MS1_SAT),colour = "blue", size = 0.5, alpha = 0.5) +
  geom_line(data = HER42L_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "orange", lwd = 0.5, alpha = 0.3) +
  geom_point(data = HER42L_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "orange", size = 0.5, alpha = 0.3) +
  geom_line(data = HER49L_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "grey", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER49L_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "grey", size = 0.5, alpha = 0.5) +
  geom_line(data = HER42PB_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "black", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER42PB_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "black", size = 0.5, alpha = 0.5) +
  geom_line(data = HERPT2S1_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "brown", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HERPT2S1_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "brown", size = 0.5, alpha = 0.5) +
  # add a red dotted line at point of interests (mid HER24L-3B core jump to higher MS and GRD at ~7.5 ka)
  geom_vline(xintercept = c(7500, 13800), col = "grey", lty = 2, alpha = 0.7, lwd = 0.75)
HER14L_42L_44L_42L_49L_42PB_PT2S1_MS

# DCMS plot
HER14L_DCMS <- ggplot() + 
  geom_line(data = HER14L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "white", lwd = 0.5) +
  geom_point(data = HER14L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "white", size = 0.5) +
  #scale_color_manual(values = colours) +
  labs(title = "Southern Coastal Lakes & Peat Bogs: DCMS",
       subtitle = "HER42L (orange), HER42PB(black), HERPT2S1_SS1_SEC42(brown), HER44L >7.5 ka (blue), HER49L (grey)") + 
  xlab("Age (cal a BP") +
  ylab(as.expression(expression(paste("DCMS [SAT]" (x10 ^ -8 ~m ^ 3/kg)) ))) +
  #coord_flip() +
  ylim(-5,15)
HER14L_DCMS

# add hiatus as white shading and plot HER24L and HER44L on top
zone_data <- tibble(ymin = -Inf, ymax = Inf , xmin = 1800, xmax = 10790)

HER14L_42L_44L_42L_49L_42PB_PT2S1_DCMS  <- HER14L_DCMS +
  geom_rect(
    mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
    data = zone_data, 
    #alpha = 0.2,
    fill = "white",
    inherit.aes = FALSE
  ) +
  #geom_line(data = HER24L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "red", lwd = 0.5, alpha = 0.5) +
  #geom_point(data = HER24L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "red", size = 0.5, alpha = 0.5) +
  geom_line(data = HER44L_comp_9ka, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "blue", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER44L_comp_9ka, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "blue", size = 0.5, alpha = 0.5) +
  geom_line(data = HER42L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "orange", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER42L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "orange", size = 0.5, alpha = 0.5) +
  geom_line(data = HER49L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "grey", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER49L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "grey", size = 0.5, alpha = 0.5) +
  geom_line(data = HER42PB_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "black", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER42PB_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "black", size = 0.5, alpha = 0.5) +
  geom_line(data = HERPT2S1_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "brown", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HERPT2S1_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "brown", size = 0.5, alpha = 0.5) +
  # add a red dotted line at point of interests (mid HER24L-3B core jump to higher MS and GRD at ~7.5 ka)
  geom_vline(xintercept = c(7500, 13800), col = "grey", lty = 2, alpha = 0.7, lwd = 0.75)
HER14L_42L_44L_42L_49L_42PB_PT2S1_DCMS

ggarrange(HER14L_42L_44L_42L_49L_42PB_PT2S1_GRD, 
          HER14L_42L_44L_42L_49L_42PB_PT2S1_MS, 
          HER14L_42L_44L_42L_49L_42PB_PT2S1_DCMS, nrow = 3)
ggsave("Outputs/Figures/Combined_sites/Fig1.5_SouthernLakes&Peat_HER14L_42L_44L_42L_49L_42PB_PT2S1.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")


# 1.6 SOUTHERN PEAT BOGS: HER42PB + PT2S1_SS1_XRF6_SEC42  ------------------------------ 

theme_set(theme_classic(base_size=12) + theme(
  plot.title = element_text(color="black", size=12, face="bold"),
  axis.title.x = element_text(color="black", size=12),
  axis.title.y = element_text(color="black", size=12)
))

#filter HER44L to remove data <6ka - as there's no age data less then 6ka - pot HER14/24/44L>6ka
HER44L_comp_9ka <- filter(HER44L_comp, SH20_mean_age >7500)
#HERPT2S1_comp_5ka <- filter(HERPT2S1_comp, SH20_mean_age >4800)

# GRD plot
HER14L_GRD <- ggplot() + 
  geom_line(data = HER14L_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "white", lwd = 0.5) +
  geom_point(data = HER14L_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "white", size = 0.5) +
  #scale_color_manual(values = colours) +
  labs(title = "Southern Peat Bogs: GRD",
       subtitle = "HER42PB(black), HERPT2S1_SS1_SEC42(brown) >5ka") + 
  xlab("Age (cal a BP") +
  ylab(as.expression(expression(paste("GRDensity [SAT] " (kg/m ^ 3 )) ))) +
  #coord_flip() +
  ylim(0.5,2.5)
HER14L_GRD

# add hiatus as white shading and plot HER24L and HER44L >9.4 ka on top
zone_data <- tibble(ymin = -Inf, ymax = Inf , xmin = 1800, xmax = 10790)

HER42PB_PT2S15ka_GRD  <- HER14L_GRD +
  geom_rect(
    mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
    data = zone_data, 
    #alpha = 0.2,
    fill = "white",
    inherit.aes = FALSE
  ) +
  #geom_line(data = HER24L_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "red", lwd = 0.5, alpha = 0.5) +
  #geom_point(data = HER24L_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "red", size = 0.5, alpha = 0.5) +
  #geom_line(data = HER44L_comp_9ka, aes(x=SH20_mean_age, y=Den1_SAT), colour = "blue", lwd = 0.5, alpha = 0.5) +
  #geom_point(data = HER44L_comp_9ka, aes(x=SH20_mean_age, y=Den1_SAT),colour = "blue", size = 0.5, alpha = 0.5) +
  #geom_line(data = HER42L_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "orange", lwd = 0.5, alpha = 0.5) +
  #geom_point(data = HER42L_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "orange", size = 0.5, alpha = 0.5) +
  #geom_line(data = HER49L_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "grey", lwd = 0.5, alpha = 0.5) +
  #geom_point(data = HER49L_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "grey", size = 0.5, alpha = 0.5) +
  geom_line(data = HER42PB_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "black", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER42PB_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "black", size = 0.5, alpha = 0.5) +
  geom_line(data = HERPT2S1_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "brown", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HERPT2S1_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "brown", size = 0.5, alpha = 0.5) +
  # add a red dotted line at point of interests (mid HER24L-3B core jump to higher MS and GRD at ~7.4 ka)
  geom_vline(xintercept = c(7400, 13800), col = "grey", lty = 2, alpha = 0.7, lwd = 0.75)
HER42PB_PT2S15ka_GRD

# MS plot
HER14L_MS <- ggplot() + 
  geom_line(data = HER14L_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "white", lwd = 0.5) +
  geom_point(data = HER14L_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "white", size = 0.5) +
  #scale_color_manual(values = colours) +
  labs(title = "Southern Peat Bogs: MS",
       subtitle = "HER42PB(black), HERPT2S1_SS1_SEC42(brown) >5ka") + 
  xlab("Age (cal a BP") +
  ylab(as.expression(expression(paste("MS [SAT] SI " (x10 ^ -5)) ))) +
  #coord_flip() +
  ylim(-5,15)
HER14L_MS

# add hiatus as white shading and plot HER24L and HER44L on top
zone_data <- tibble(ymin = -Inf, ymax = Inf , xmin = 1800, xmax = 10790)

HER42PB_PT2S15ka_MS  <- HER14L_MS +
  geom_rect(
    mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
    data = zone_data, 
    #alpha = 0.2,
    fill = "white",
    inherit.aes = FALSE
  ) +
  #geom_line(data = HER24L_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "red", lwd = 0.5, alpha = 0.5) +
  #geom_point(data = HER24L_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "red", size = 0.5, alpha = 0.5) +
  #geom_line(data = HER44L_comp_9ka, aes(x=SH20_mean_age, y=MS1_SAT), colour = "blue", lwd = 0.5, alpha = 0.5) +
  #geom_point(data = HER44L_comp_9ka, aes(x=SH20_mean_age, y=MS1_SAT),colour = "blue", size = 0.5, alpha = 0.5) +
  #geom_line(data = HER42L_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "orange", lwd = 0.5, alpha = 0.3) +
  #geom_point(data = HER42L_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "orange", size = 0.5, alpha = 0.3) +
  #geom_line(data = HER49L_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "grey", lwd = 0.5, alpha = 0.5) +
  #geom_point(data = HER49L_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "grey", size = 0.5, alpha = 0.5) +
  geom_line(data = HER42PB_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "black", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER42PB_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "black", size = 0.5, alpha = 0.5) +
  geom_line(data = HERPT2S1_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "brown", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HERPT2S1_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "brown", size = 0.5, alpha = 0.5) +
  # add a red dotted line at point of interests (mid HER24L-3B core jump to higher MS and GRD at ~7.5 ka)
  geom_vline(xintercept = c(7500, 13800), col = "grey", lty = 2, alpha = 0.7, lwd = 0.75)
HER42PB_PT2S15ka_MS

# DCMS plot
HER14L_DCMS <- ggplot() + 
  geom_line(data = HER14L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "white", lwd = 0.5) +
  geom_point(data = HER14L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "white", size = 0.5) +
  #scale_color_manual(values = colours) +
  labs(title = "Southern Peat Bogs: DCMS",
       subtitle = "HER42PB(black), HERPT2S1_SS1_SEC42(brown) >5ka") + 
  xlab("Age (cal a BP") +
  ylab(as.expression(expression(paste("DCMS [SAT]" (x10 ^ -8 ~m ^ 3/kg)) ))) +
  #coord_flip() +
  ylim(-5,15)
HER14L_DCMS

# add hiatus as white shading and plot HER24L and HER44L on top
zone_data <- tibble(ymin = -Inf, ymax = Inf , xmin = 1800, xmax = 10790)

HER42PB_PT2S15ka_DCMS  <- HER14L_DCMS +
  geom_rect(
    mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
    data = zone_data, 
    #alpha = 0.2,
    fill = "white",
    inherit.aes = FALSE
  ) +
  #geom_line(data = HER24L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "red", lwd = 0.5, alpha = 0.5) +
  #geom_point(data = HER24L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "red", size = 0.5, alpha = 0.5) +
  #geom_line(data = HER44L_comp_9ka, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "blue", lwd = 0.5, alpha = 0.5) +
  #geom_point(data = HER44L_comp_9ka, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "blue", size = 0.5, alpha = 0.5) +
  #geom_line(data = HER42L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "orange", lwd = 0.5, alpha = 0.5) +
  #geom_point(data = HER42L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "orange", size = 0.5, alpha = 0.5) +
  #geom_line(data = HER49L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "grey", lwd = 0.5, alpha = 0.5) +
  #geom_point(data = HER49L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "grey", size = 0.5, alpha = 0.5) +
  geom_line(data = HER42PB_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "black", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER42PB_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "black", size = 0.5, alpha = 0.5) +
  geom_line(data = HERPT2S1_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "brown", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HERPT2S1_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "brown", size = 0.5, alpha = 0.5) +
  # add a red dotted line at point of interests (mid HER24L-3B core jump to higher MS and GRD at ~7.5 ka)
  geom_vline(xintercept = c(7500, 13800), col = "grey", lty = 2, alpha = 0.7, lwd = 0.75)
HER42PB_PT2S15ka_DCMS

ggarrange(HER42PB_PT2S15ka_GRD, 
          HER42PB_PT2S15ka_MS, 
          HER42PB_PT2S15ka_DCMS, nrow = 3)
ggsave("Outputs/Figures/Combined_sites/Fig1.6_Southern Coastal Peat_HER42PB_PT2S1.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")


# 1.7 SOUTHERN PEAT BOGS + >10 ka Lakes: PT2S1 + HER44L >10 ka + HER14L >10 ka ------------------------------ 

#clear plot window
dev.off()

theme_set(theme_classic(base_size=12) + theme(
  plot.title = element_text(color="black", size=12, face="bold"),
  axis.title.x = element_text(color="black", size=12),
  axis.title.y = element_text(color="black", size=12)
))

#filter HER44L to remove data <6ka - as there's no age data less then 6ka - pot HER14/24/44L>6ka
HER44L_comp_9ka <- filter(HER44L_comp, SH20_mean_age >9000)
#HERPT2S1_comp_5ka <- filter(HERPT2S1_comp, SH20_mean_age >4800)

# GRD plot
HER14L_GRD <- ggplot() + 
  geom_line(data = HER14L_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "blue", lwd = 0.5) +
  geom_point(data = HER14L_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "blue", size = 0.5) +
  #scale_color_manual(values = colours) +
  labs(title = "Southern Peat Bogs + >10 ka Northern Lakes: GRD",
       subtitle = "HER42PB(black), HERPT2S1_SS1_SEC42(brown), HER44L>9ka(red), HER49L(grey), HER14L>10 ka(blue)") +
  xlab("Age (cal a BP") +
  ylab(as.expression(expression(paste("GRDensity [SAT] " (kg/m ^ 3 )) ))) +
  #coord_flip() +
  ylim(0.5,2.5)
HER14L_GRD

# add hiatus as white shading and plot HER24L and HER44L >9.4 ka on top
zone_data <- tibble(ymin = -Inf, ymax = Inf , xmin = 0, xmax = 10790)

HER42PB_PT2S15ka_HER44L_HER14L_10K_GRD  <- HER14L_GRD +
  geom_rect(
    mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
    data = zone_data, 
    #alpha = 0.2,
    fill = "white",
    inherit.aes = FALSE
  ) +
  #geom_line(data = HER24L_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "red", lwd = 0.5, alpha = 0.5) +
  #geom_point(data = HER24L_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "red", size = 0.5, alpha = 0.5) +
  geom_line(data = HER44L_comp_9ka, aes(x=SH20_mean_age, y=Den1_SAT), colour = "red", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER44L_comp_9ka, aes(x=SH20_mean_age, y=Den1_SAT),colour = "red", size = 0.5, alpha = 0.5) +
  #geom_line(data = HER42L_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "orange", lwd = 0.5, alpha = 0.5) +
  #geom_point(data = HER42L_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "orange", size = 0.5, alpha = 0.5) +
  geom_line(data = HER49L_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "grey", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER49L_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "grey", size = 0.5, alpha = 0.5) +
  geom_line(data = HER42PB_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "black", lwd = 0.5, alpha = 1) +
  geom_point(data = HER42PB_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "black", size = 0.5, alpha = 1) +
  geom_line(data = HERPT2S1_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "brown", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HERPT2S1_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "brown", size = 0.5, alpha = 0.5) +
  # add a red dotted line at point of interests (mid HER24L-3B core jump to higher MS and GRD at ~7.4 ka)
  geom_vline(xintercept = c(7400, 13800), col = "grey", lty = 2, alpha = 0.7, lwd = 0.75)
HER42PB_PT2S15ka_HER44L_HER14L_10K_GRD

# MS plot
HER14L_MS <- ggplot() + 
  geom_line(data = HER14L_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "blue", lwd = 0.5) +
  geom_point(data = HER14L_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "blue", size = 0.5) +
  #scale_color_manual(values = colours) +
  labs(title = "Southern Peat Bogs + >10 ka Northern Lakes: MS",
       subtitle = "HER42PB(black), HERPT2S1_SS1_SEC42(brown), HER44L>9ka(red), HER49L(grey), HER14L>10 ka(blue)") +
  xlab("Age (cal a BP") +
  ylab(as.expression(expression(paste("MS [SAT] SI " (x10 ^ -5)) ))) +
  #coord_flip() +
  ylim(-5,15)
HER14L_MS

# add hiatus as white shading and plot HER24L and HER44L on top
zone_data <- tibble(ymin = -Inf, ymax = Inf , xmin = 0, xmax = 10790)

HER42PB_PT2S15ka_HER44L_HER14L_10K_MS  <- HER14L_MS +
  geom_rect(
    mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
    data = zone_data, 
    #alpha = 0.2,
    fill = "white",
    inherit.aes = FALSE
  ) +
  #geom_line(data = HER24L_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "red", lwd = 0.5, alpha = 0.5) +
  #geom_point(data = HER24L_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "red", size = 0.5, alpha = 0.5) +
  geom_line(data = HER44L_comp_9ka, aes(x=SH20_mean_age, y=MS1_SAT), colour = "red", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER44L_comp_9ka, aes(x=SH20_mean_age, y=MS1_SAT),colour = "red", size = 0.5, alpha = 0.5) +
  #geom_line(data = HER42L_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "orange", lwd = 0.5, alpha = 0.3) +
  #geom_point(data = HER42L_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "orange", size = 0.5, alpha = 0.3) +
  geom_line(data = HER49L_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "grey", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER49L_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "grey", size = 0.5, alpha = 0.5) +
  geom_line(data = HER42PB_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "black", lwd = 0.5, alpha = 1) +
  geom_point(data = HER42PB_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "black", size = 0.5, alpha = 1) +
  geom_line(data = HERPT2S1_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "brown", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HERPT2S1_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "brown", size = 0.5, alpha = 0.5) +
  # add a red dotted line at point of interests (mid HER24L-3B core jump to higher MS and GRD at ~7.5 ka)
  geom_vline(xintercept = c(7500, 13800), col = "grey", lty = 2, alpha = 0.7, lwd = 0.75)
HER42PB_PT2S15ka_HER44L_HER14L_10K_MS

# DCMS plot
HER14L_DCMS <- ggplot() + 
  geom_line(data = HER14L_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "blue", lwd = 0.5) +
  geom_point(data = HER14L_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "blue", size = 0.5) +
  #scale_color_manual(values = colours) +
  labs(title = "Southern Peat Bogs + >10 ka Northern Lakes: DCMS",
       subtitle = "HER42PB(black), HERPT2S1_SS1_SEC42(brown), HER44L>9ka(blue), HER49L(grey), HER14L>10 ka(blue)") +
  xlab("Age (cal a BP") +
  ylab(as.expression(expression(paste("DCMS [SAT]" (x10 ^ -8 ~m ^ 3/kg)) ))) +
  #coord_flip() +
  ylim(-5,15)
HER14L_DCMS

# add hiatus as white shading and plot HER24L and HER44L on top
zone_data <- tibble(ymin = -Inf, ymax = Inf , xmin = 0, xmax = 10790)

HER42PB_PT2S15ka_HER44L_HER14L_10K_DCMS  <- HER14L_DCMS +
  geom_rect(
    mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
    data = zone_data, 
    #alpha = 0.2,
    fill = "white",
    inherit.aes = FALSE
  ) +
  #geom_line(data = HER24L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "red", lwd = 0.5, alpha = 0.5) +
  #geom_point(data = HER24L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "red", size = 0.5, alpha = 0.5) +
  geom_line(data = HER44L_comp_9ka, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "red", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER44L_comp_9ka, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "red", size = 0.5, alpha = 0.5) +
  #geom_line(data = HER42L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "orange", lwd = 0.5, alpha = 0.5) +
  #geom_point(data = HER42L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "orange", size = 0.5, alpha = 0.5) +
  geom_line(data = HER49L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "grey", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER49L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "grey", size = 0.5, alpha = 0.5) +
  geom_line(data = HER42PB_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "black", lwd = 0.5, alpha = 1) +
  geom_point(data = HER42PB_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "black", size = 0.5, alpha = 1) +
  geom_line(data = HERPT2S1_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "brown", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HERPT2S1_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "brown", size = 0.5, alpha = 0.5) +
  # add a red dotted line at point of interests (mid HER24L-3B core jump to higher MS and GRD at ~7.5 ka)
  geom_vline(xintercept = c(7500, 13800), col = "grey", lty = 2, alpha = 0.7, lwd = 0.75)
HER42PB_PT2S15ka_HER44L_HER14L_10K_DCMS

ggarrange(HER42PB_PT2S15ka_HER44L_HER14L_10K_GRD, 
          HER42PB_PT2S15ka_HER44L_HER14L_10K_MS, 
          HER42PB_PT2S15ka_HER44L_HER14L_10K_DCMS, nrow = 3)
ggsave("Outputs/Figures/Combined_sites/Fig1.7_SouthernPeat_Lakes10K_HER42PB_PT2S1_5K_HER44L_HER14L_10K.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")


# 1.8 SOUTHERN PEAT BOGS: PT2S1_SS1_XRF6_SEC42 +HER49L + HER44L>9ka ------------------------------ 

theme_set(theme_classic(base_size=12) + theme(
  plot.title = element_text(color="black", size=12, face="bold"),
  axis.title.x = element_text(color="black", size=12),
  axis.title.y = element_text(color="black", size=12)
))

#filter HER44L to remove data <6ka - as there's no age data less then 6ka - pot HER14/24/44L>6ka
HER44L_comp_9ka <- filter(HER44L_comp, SH20_mean_age >9000)
#HERPT2S1_comp_5ka <- filter(HERPT2S1_comp, SH20_mean_age >4800)

# GRD plot
HER14L_GRD <- ggplot() + 
  geom_line(data = HER14L_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "blue", lwd = 0.5) +
  geom_point(data = HER14L_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "blue", size = 0.5) +
  #scale_color_manual(values = colours) +
  labs(title = "Southern Peat Bog + HER49L + HER44L + HER14L:GRD",
       subtitle = "HERPT2S1_SS1_SEC42(brown), HER49L(grey), HER44L>8ka(red), HER14L>10ka(blue)") + 
  xlab("Age (cal a BP") +
  ylab(as.expression(expression(paste("GRDensity [SAT] " (kg/m ^ 3 )) ))) +
  #coord_flip() +
  ylim(0.5,2.5)
HER14L_GRD

# add hiatus as white shading and plot HER24L and HER44L >9.4 ka on top
zone_data <- tibble(ymin = -Inf, ymax = Inf , xmin = 0, xmax = 10790)

HER49L_PT2S1_GRD  <- HER14L_GRD +
  geom_rect(
    mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
    data = zone_data, 
    #alpha = 0.2,
    fill = "white",
    inherit.aes = FALSE
  ) +
  #geom_line(data = HER24L_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "red", lwd = 0.5, alpha = 0.5) +
  #geom_point(data = HER24L_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "red", size = 0.5, alpha = 0.5) +
  geom_line(data = HER44L_comp_9ka, aes(x=SH20_mean_age, y=Den1_SAT), colour = "blue", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER44L_comp_9ka, aes(x=SH20_mean_age, y=Den1_SAT),colour = "blue", size = 0.5, alpha = 0.5) +
  #geom_line(data = HER42L_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "orange", lwd = 0.5, alpha = 0.5) +
  #geom_point(data = HER42L_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "orange", size = 0.5, alpha = 0.5) +
  geom_line(data = HER49L_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "grey", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER49L_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "grey", size = 0.5, alpha = 0.5) +
  #geom_line(data = HER42PB_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "black", lwd = 0.5, alpha = 0.5) +
  #geom_point(data = HER42PB_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "black", size = 0.5, alpha = 0.5) +
  geom_line(data = HERPT2S1_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "brown", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HERPT2S1_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "brown", size = 0.5, alpha = 0.5) +
  # add a red dotted line at point of interests (mid HER24L-3B core jump to higher MS and GRD at ~7.4 ka)
  geom_vline(xintercept = c(7400, 13800), col = "grey", lty = 2, alpha = 0.7, lwd = 0.75)
HER49L_PT2S1_GRD

# MS plot
HER14L_MS <- ggplot() + 
  geom_line(data = HER14L_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "blue", lwd = 0.5) +
  geom_point(data = HER14L_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "blue", size = 0.5) +
  #scale_color_manual(values = colours) +
  labs(title = "Southern Peat Bog + HER49L + HER44L + HER14L: MS",
       subtitle = "HERPT2S1_SS1_SEC42(brown), HER49L(grey), HER44L>8ka(red), HER14L>10ka(blue)") +   xlab("Age (cal a BP") +
  ylab(as.expression(expression(paste("MS [SAT] SI " (x10 ^ -5)) ))) +
  #coord_flip() +
  ylim(-5,15)
HER14L_MS

# add hiatus as white shading and plot HER24L and HER44L on top
zone_data <- tibble(ymin = -Inf, ymax = Inf , xmin = 0, xmax = 10790)

HER49L_PT2S1_MS  <- HER14L_MS +
  geom_rect(
    mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
    data = zone_data, 
    #alpha = 0.2,
    fill = "white",
    inherit.aes = FALSE
  ) +
  #geom_line(data = HER24L_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "red", lwd = 0.5, alpha = 0.5) +
  #geom_point(data = HER24L_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "red", size = 0.5, alpha = 0.5) +
  geom_line(data = HER44L_comp_9ka, aes(x=SH20_mean_age, y=MS1_SAT), colour = "blue", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER44L_comp_9ka, aes(x=SH20_mean_age, y=MS1_SAT),colour = "blue", size = 0.5, alpha = 0.5) +
  #geom_line(data = HER42L_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "orange", lwd = 0.5, alpha = 0.3) +
  #geom_point(data = HER42L_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "orange", size = 0.5, alpha = 0.3) +
  geom_line(data = HER49L_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "grey", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER49L_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "grey", size = 0.5, alpha = 0.5) +
  #geom_line(data = HER42PB_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "black", lwd = 0.5, alpha = 0.5) +
  #geom_point(data = HER42PB_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "black", size = 0.5, alpha = 0.5) +
  geom_line(data = HERPT2S1_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "brown", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HERPT2S1_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "brown", size = 0.5, alpha = 0.5) +
  # add a red dotted line at point of interests (mid HER24L-3B core jump to higher MS and GRD at ~7.5 ka)
  geom_vline(xintercept = c(7500, 13800), col = "grey", lty = 2, alpha = 0.7, lwd = 0.75)
HER49L_PT2S1_MS

# DCMS plot
HER14L_DCMS <- ggplot() + 
  geom_line(data = HER14L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "blue", lwd = 0.5) +
  geom_point(data = HER14L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "blue", size = 0.5) +
  #scale_color_manual(values = colours) +
  labs(title = "Southern Peat Bog + HER49L + HER44L + HER14L:DCMS",
       subtitle = "HERPT2S1_SS1_SEC42(brown), HER49L(grey), HER44L>8ka(red), HER14L>10ka(blue)") +   xlab("Age (cal a BP") +
  ylab(as.expression(expression(paste("DCMS [SAT]" (x10 ^ -8 ~m ^ 3/kg)) ))) +
  #coord_flip() +
  ylim(-5,15)
HER14L_DCMS

# add hiatus as white shading and plot HER24L and HER44L on top
zone_data <- tibble(ymin = -Inf, ymax = Inf , xmin = 0, xmax = 10790)

HER49L_PT2S1_DCMS  <- HER14L_DCMS +
  geom_rect(
    mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
    data = zone_data, 
    #alpha = 0.2,
    fill = "white",
    inherit.aes = FALSE
  ) +
  #geom_line(data = HER24L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "red", lwd = 0.5, alpha = 0.5) +
  #geom_point(data = HER24L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "red", size = 0.5, alpha = 0.5) +
  geom_line(data = HER44L_comp_9ka, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "blue", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER44L_comp_9ka, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "blue", size = 0.5, alpha = 0.5) +
  #geom_line(data = HER42L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "orange", lwd = 0.5, alpha = 0.5) +
  #geom_point(data = HER42L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "orange", size = 0.5, alpha = 0.5) +
  geom_line(data = HER49L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "grey", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER49L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "grey", size = 0.5, alpha = 0.5) +
  #geom_line(data = HER42PB_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "black", lwd = 0.5, alpha = 0.5) +
  #geom_point(data = HER42PB_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "black", size = 0.5, alpha = 0.5) +
  geom_line(data = HERPT2S1_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "brown", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HERPT2S1_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "brown", size = 0.5, alpha = 0.5) +
  # add a red dotted line at point of interests (mid HER24L-3B core jump to higher MS and GRD at ~7.5 ka)
  geom_vline(xintercept = c(7500, 13800), col = "grey", lty = 2, alpha = 0.7, lwd = 0.75)
HER49L_PT2S1_DCMS

ggarrange(HER49L_PT2S1_GRD, 
          HER49L_PT2S1_MS, 
          HER49L_PT2S1_DCMS, nrow = 3)
ggsave("Outputs/Figures/Combined_sites/Fig1.8_Southern Coastal Peat_HERPT2S1+49L+44L+14L.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")


# 1.9 SOUTHERN PEAT BOGS: HER42PB +HER49L + HER44L>9ka ------------------------------ 

theme_set(theme_classic(base_size=12) + theme(
  plot.title = element_text(color="black", size=12, face="bold"),
  axis.title.x = element_text(color="black", size=12),
  axis.title.y = element_text(color="black", size=12)
))

#filter HER44L to remove data <6ka - as there's no age data less then 6ka - pot HER14/24/44L>6ka
HER44L_comp_9ka <- filter(HER44L_comp, SH20_mean_age >9000)
#HERPT2S1_comp_5ka <- filter(HERPT2S1_comp, SH20_mean_age >4800)

# GRD plot
HER14L_GRD <- ggplot() + 
  geom_line(data = HER14L_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "blue", lwd = 0.5) +
  geom_point(data = HER14L_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "blue", size = 0.5) +
  #scale_color_manual(values = colours) +
  labs(title = "Southern Peat Bog HER42PB + HER49L + HER44L + HER14L: GRD",
       subtitle = "HER42PB(brown), HER49L(grey), HER44L>8ka(blue), HER14L>10ka(blue)") +   
  xlab("Age (cal a BP") +
  ylab(as.expression(expression(paste("GRDensity [SAT] " (kg/m ^ 3 )) ))) +
  #coord_flip() +
  ylim(0.5,2.5)
HER14L_GRD

# add hiatus as white shading and plot HER24L and HER44L >9.4 ka on top
zone_data <- tibble(ymin = -Inf, ymax = Inf , xmin = 0, xmax = 10790)

HER49L_HER42PB_GRD  <- HER14L_GRD +
  geom_rect(
    mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
    data = zone_data, 
    #alpha = 0.2,
    fill = "white",
    inherit.aes = FALSE
  ) +
  #geom_line(data = HER24L_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "red", lwd = 0.5, alpha = 0.5) +
  #geom_point(data = HER24L_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "red", size = 0.5, alpha = 0.5) +
  geom_line(data = HER44L_comp_9ka, aes(x=SH20_mean_age, y=Den1_SAT), colour = "blue", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER44L_comp_9ka, aes(x=SH20_mean_age, y=Den1_SAT),colour = "blue", size = 0.5, alpha = 0.5) +
  #geom_line(data = HER42L_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "orange", lwd = 0.5, alpha = 0.5) +
  #geom_point(data = HER42L_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "orange", size = 0.5, alpha = 0.5) +
  geom_line(data = HER49L_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "grey", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER49L_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "grey", size = 0.5, alpha = 0.5) +
  geom_line(data = HER42PB_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "black", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER42PB_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "black", size = 0.5, alpha = 0.5) +
  #geom_line(data = HERPT2S1_comp, aes(x=SH20_mean_age, y=Den1_SAT), colour = "brown", lwd = 0.5, alpha = 0.5) +
  #geom_point(data = HERPT2S1_comp, aes(x=SH20_mean_age, y=Den1_SAT),colour = "brown", size = 0.5, alpha = 0.5) +
  # add a red dotted line at point of interests (mid HER24L-3B core jump to higher MS and GRD at ~7.4 ka)
  geom_vline(xintercept = c(7400, 13800), col = "grey", lty = 2, alpha = 0.7, lwd = 0.75)
HER49L_HER42PB_GRD

# MS plot
HER14L_MS <- ggplot() + 
  geom_line(data = HER14L_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "blue", lwd = 0.5) +
  geom_point(data = HER14L_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "blue", size = 0.5) +
  #scale_color_manual(values = colours) +
  labs(title = "Southern Peat Bog HER42PB + HER49L + HER44L + HER14L: MS",
       subtitle = "HER42PB(brown), HER49L(grey), HER44L>8ka(blue), HER14L>10ka(blue)") +   
  xlab("Age (cal a BP") +
  ylab(as.expression(expression(paste("MS [SAT] SI " (x10 ^ -5)) ))) +
  #coord_flip() +
  ylim(-5,15)
HER14L_MS

# add hiatus as white shading and plot HER24L and HER44L on top
zone_data <- tibble(ymin = -Inf, ymax = Inf , xmin = 0, xmax = 10790)

HER49L_HER42PB_MS  <- HER14L_MS +
  geom_rect(
    mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
    data = zone_data, 
    #alpha = 0.2,
    fill = "white",
    inherit.aes = FALSE
  ) +
  #geom_line(data = HER24L_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "red", lwd = 0.5, alpha = 0.5) +
  #geom_point(data = HER24L_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "red", size = 0.5, alpha = 0.5) +
  geom_line(data = HER44L_comp_9ka, aes(x=SH20_mean_age, y=MS1_SAT), colour = "blue", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER44L_comp_9ka, aes(x=SH20_mean_age, y=MS1_SAT),colour = "blue", size = 0.5, alpha = 0.5) +
  #geom_line(data = HER42L_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "orange", lwd = 0.5, alpha = 0.3) +
  #geom_point(data = HER42L_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "orange", size = 0.5, alpha = 0.3) +
  geom_line(data = HER49L_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "grey", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER49L_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "grey", size = 0.5, alpha = 0.5) +
  geom_line(data = HER42PB_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "black", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER42PB_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "black", size = 0.5, alpha = 0.5) +
  #geom_line(data = HERPT2S1_comp, aes(x=SH20_mean_age, y=MS1_SAT), colour = "brown", lwd = 0.5, alpha = 0.5) +
  #geom_point(data = HERPT2S1_comp, aes(x=SH20_mean_age, y=MS1_SAT),colour = "brown", size = 0.5, alpha = 0.5) +
  # add a red dotted line at point of interests (mid HER24L-3B core jump to higher MS and GRD at ~7.5 ka)
  geom_vline(xintercept = c(7500, 13800), col = "grey", lty = 2, alpha = 0.7, lwd = 0.75)
HER49L_HER42PB_MS

# DCMS plot
HER14L_DCMS <- ggplot() + 
  geom_line(data = HER14L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "blue", lwd = 0.5) +
  geom_point(data = HER14L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "blue", size = 0.5) +
  #scale_color_manual(values = colours) +
  labs(title = "Southern Peat Bog HER42PB + HER49L + HER44L + HER14L: DCMS",
       subtitle = "HER42PB(brown), HER49L(grey), HER44L>8ka(blue), HER14L>10ka(blue)") +   
  xlab("Age (cal a BP") +
  ylab(as.expression(expression(paste("DCMS [SAT]" (x10 ^ -8 ~m ^ 3/kg)) ))) +
  #coord_flip() +
  ylim(-5,15)
HER14L_DCMS

# add hiatus as white shading and plot HER24L and HER44L on top
zone_data <- tibble(ymin = -Inf, ymax = Inf , xmin = 0, xmax = 10790)

HER49L_HER42PB_DCMS  <- HER14L_DCMS +
  geom_rect(
    mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
    data = zone_data, 
    #alpha = 0.2,
    fill = "white",
    inherit.aes = FALSE
  ) +
  #geom_line(data = HER24L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "red", lwd = 0.5, alpha = 0.5) +
  #geom_point(data = HER24L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "red", size = 0.5, alpha = 0.5) +
  geom_line(data = HER44L_comp_9ka, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "blue", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER44L_comp_9ka, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "blue", size = 0.5, alpha = 0.5) +
  #geom_line(data = HER42L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "orange", lwd = 0.5, alpha = 0.5) +
  #geom_point(data = HER42L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "orange", size = 0.5, alpha = 0.5) +
  geom_line(data = HER49L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "grey", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER49L_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "grey", size = 0.5, alpha = 0.5) +
  geom_line(data = HER42PB_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "black", lwd = 0.5, alpha = 0.5) +
  geom_point(data = HER42PB_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "black", size = 0.5, alpha = 0.5) +
  #geom_line(data = HERPT2S1_comp, aes(x=SH20_mean_age, y=DCMS1_SAT), colour = "brown", lwd = 0.5, alpha = 0.5) +
  #geom_point(data = HERPT2S1_comp, aes(x=SH20_mean_age, y=DCMS1_SAT),colour = "brown", size = 0.5, alpha = 0.5) +
  # add a red dotted line at point of interests (mid HER24L-3B core jump to higher MS and GRD at ~7.5 ka)
  geom_vline(xintercept = c(7500, 13800), col = "grey", lty = 2, alpha = 0.7, lwd = 0.75)
HER49L_HER42PB_DCMS

ggarrange(HER49L_HER42PB_GRD, 
          HER49L_HER42PB_MS, 
          HER49L_HER42PB_DCMS, nrow = 3)
ggsave("Outputs/Figures/Combined_sites/Fig1.9_Southern Coastal Peat_HER42PB+49L+44L+14L.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")
