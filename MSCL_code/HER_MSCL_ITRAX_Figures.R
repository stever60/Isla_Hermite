# HERMITE PAPER

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

# -------------------------------------------------------------------------

# Data input

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



# ITRAX Composite dataset for ALL sites from original files ----------------------------

# import clr site datasets
HER14L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER14L/Output/Composite/HER14L_xrf1_clr.csv") 
HER24L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER24L/Output/Composite/HER24L_xrf1_clr.csv") 
HER34L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER34L/Output/Composite/HER34L_xrf1_clr.csv")
HER42L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42L/Output/Composite/HER42L_xrf1_clr.csv")
HER42PB <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/Composite/HER42PB_xrf1_clr.csv")
HER44L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER44L/Output/Composite/HER44L_xrf1_clr.csv")
HER49L <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER49L/Output/Composite/HER49L_xrf1_clr.csv")
# Assign elements from PeriodicTable package to 'elementsList' - to select all elements present across all sites
data(periodicTable)
elementsList <- periodicTable$symb
# bind all data and rename to match itraxR column names - 
# all NA's in a column means that element failed filtering test at that site
HER_ITRAX_bind <-  bind_rows(HER14L, HER24L, HER34L, HER42L, HER42PB, HER44L, HER49L) %>% 
  rename(mean_age = SH20_age) %>% 
  select(Site, depth, mean_age, any_of(elementsList), Mo_inc: LnTi_Ca)
HER14L_0ka <- HER_ITRAX_bind %>% 
  filter(Site=='HER14L' & mean_age <10000) %>% 
  mutate(Site = recode(Site, 'HER14L' = 'HER14L_0ka'))
HER14L_10ka <- HER_ITRAX_clr %>% 
  filter(Site=='HER14L' & mean_age >10000) %>% 
  mutate(Site = recode(Site, 'HER14L' = 'HER14L_10ka'))
HER_ITRAX_xrf1_clr_comp <- HER_ITRAX_bind %>% 
  filter(!(Site=='HER14L')) %>%
  bind_rows(HER14L_0ka, HER14L_10ka) %>% 
  arrange(Site)
HER_ITRAX_xrf1_clr_comp
tail(HER_ITRAX_xrf1_clr_comp)
# remove site records from Global Environment
rm(HER14L, HER14L_0ka, HER14L_10ka, HER24L, HER34L, HER42L, HER42PB, HER44L, HER49L, HER_ITRAX_clr)
# write to file
write.csv(HER_ITRAX_xrf1_clr_comp,"Papers/Hodgson_Hermite_2022/Data/ITRAX/Outputs/HER_ITRAX_xrf1_clr_comp.csv", row.names = FALSE)

HER_ITRAX_xrf1_clr_comp <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/Outputs/HER_ITRAX_xrf1_clr_comp.csv") 
HER_MSCL_comp <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Outputs/HER_MSCL_comp.csv") 

# ACE09 Subsample density dataset  -----------------------------------------------------

ACE09_density <- read_csv("Papers/Hodgson_Hermite_2022/Data/Subsample/ACE09_density.csv") %>% 
  select(-c(MSCL_ID, depth_min, depth_max)) %>%
  rename(Density = Wet_den) %>% 
  filter(Site == 'HER42PB')
ACE09_density

# Add age-depth model - run this once - minimise and move to next section  -----------------------------

# set up depth file for RBacon
ACE09_density_depth <- ACE09_density %>% 
  select(depth) 
ACE09_density_depth
write.table(ACE09_density_depth, "RBacon/Bacon_runs/HER42PB/HER42PB_depths.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE)

# switch dir to RBacon
setwd("/Users/Steve/Dropbox/BAS/Data/R/RBacon")
library(rbacon)
#clear plot window
dev.off()
Bacon("HER42PB",depths.file=TRUE,d.max=500,thick=10, cc=3,
      postbomb=4,rotate.axes=TRUE,mem.mean=0.4,rounded=2,
      mem.strength=20,acc.mean=10, yr.max=15000,
      C14.border=rgb(0, 0, 1, 1), C14.col=rgb(0, 0, 1, 0.5), 
      cal.col=rgb(0, 0.5, 0.5, 0.5), cal.border=rgb(0, 0.5, 0.5, 1),
      dark=0.5, height=15, plot.pdf = TRUE, model.only=TRUE,
      #mean.col="red", mean.lty="1", mean.lwd ="2",
      #range.col=grey(0.5), range.lty="1", range.lwd = "1",
      title.location='topright', mgp=c(1.5, 0.7, 0))
file.copy("RBacon/Bacon_runs/HER42PB/HER42PB_52_ages.txt",
          to = "/Users/Steve/Dropbox/BAS/Data/R/Papers/Hodgson_Hermite_2022/Data/HER42PB_52_ages.txt", recursive = TRUE,
          overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)
file.copy("RBacon/Bacon_runs/HER42PB/HER42PB_52.pdf",
          to = "/Users/Steve/Dropbox/BAS/Data/R/Papers/Hodgson_Hermite_2022/Data/HER42PB_52.pdf", recursive = TRUE,
          overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)
#clear plot window
dev.off()
# switch dir back from RBacon
setwd("/Users/Steve/Dropbox/BAS/Data/R/")


# Add ages data to dataframe and rename to match MSCL Wet density --------------------------------------------
# Make sure switched back from RBacon directory
setwd("/Users/Steve/Dropbox/BAS/Data/R/")
ACE09_HER42_density <- read_tsv("Papers/Hodgson_Hermite_2022/Data/Subsample/HER42PB_52_ages.txt") %>% 
left_join(ACE09_density, HER42PB_ages, by = c("depth" = "depth")) %>% 
  rename(label = Section, min_age = min, max_age = max, median_age = median, mean_age = mean) %>% 
  relocate(min_age:mean_age, .after = depth) %>% 
  relocate(c(Site, label), .before = depth) %>%
  mutate(Site = recode(Site, 'HER42PB' = 'HER42PB_WD'))
  
ACE09_WD <- ACE09_HER42_density %>% 
  select(Site, label, depth, mean_age, Density)

# Write  to file
write.csv(ACE09_HER42_density,"Papers/Hodgson_Hermite_2022/Data/Subsample/Outputs/ACE09_HER42_density.csv", row.names = FALSE)
write.csv(ACE09_WD,"Papers/Hodgson_Hermite_2022/Data/Subsample/Outputs/ACE09_WD.csv", row.names = FALSE)

HER_ITRAX_xrf1_clr_comp <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/Outputs/HER_ITRAX_xrf1_clr_comp.csv") 
HER_MSCL_comp <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Outputs/HER_MSCL_comp.csv") 
ACE09_HER42_density<- read_csv("Papers/Hodgson_Hermite_2022/Data/Subsample/Outputs/ACE09_HER42_density.csv") 
ACE09_WD <- read_csv("Papers/Hodgson_Hermite_2022/Data/Subsample/Outputs/ACE09_WD.csv") 

# -------------------------------------------------------------------------

# START FROM HERE IF MSCL & TRAX COMPOSITE FILES HAVE AREADY BEEN MADE

# -------------------------------------------------------------------------

# Import files needed below

# set working directory
setwd("/Users/Steve/Dropbox/BAS/Data/R/")

# MSCL
HER_MSCL_comp <- read_csv("Papers/Hodgson_Hermite_2022/Data/MSCL/Outputs/HER_MSCL_comp.csv") 
ACE09_HER42_density<- read_csv("Papers/Hodgson_Hermite_2022/Data/Subsample/Outputs/ACE09_HER42_density.csv") 
ACE09_WD <- read_csv("Papers/Hodgson_Hermite_2022/Data/Subsample/Outputs/ACE09_WD.csv") 
# ITRAX - choose one of these - or others from folder
HER_ITRAX_xrf1_cps_comp <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/Outputs/HER_ITRAX_xrf1_cps_comp.csv") 
HER_ITRAX_xrf1_inc_comp <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/Outputs/HER_ITRAX_xrf1_inc_comp.csv") 
HER_ITRAX_xrf1_pc_cps_comp <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/Outputs/HER_ITRAX_xrf1_pc_cps_comp.csv") 
HER_ITRAX_xrf1_Ln_Ti_Z_comp <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/Outputs/HER_ITRAX_xrf1_Ln_Ti_Z_comp.csv")
HER_ITRAX_xrf1_clr_comp <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/Outputs/HER_ITRAX_xrf1_clr_comp.csv") 

# Plotting set up--------------------------------------------------------------

library(tidypaleo)

# Define plot order
Site_plot_order <- c('HER14L_0ka', 'HER14L_10ka','HER24L','HER34L',
                     'HER42L', 'HER42PB', 'HER44L','HER49L', 
                     'HERPT1S1','HERPT1S4', 'HERPT1S8', 'HERPT2S1')
jco_violin0 <- c("#0073C2FF", "#0073C2FF", "#003C67FF","#313695", 
                "#7AA6DCFF","#B24745FF", "#003C67FF", "#374E55FF",
                "#DFC27D", "#BF812D", "#8C510A", "#543005")

Site_rename_GRD <- c('HER14L_0ka' = 'HER14L_0ka_GRD', 'HER14L_10ka' = 'HER14L_10ka_GRD','HER24L' = 'HER24L_GRD', 
                 'HER34L' = 'HER34L_GRD','HER42L' = 'HER42L_GRD', 
                 'HER44L' = 'HER44L_GRD', 'HER49L' = 'HER49L_GRD', 
                 'HER42PB' = 'HER42PB_GRD','HER42PB_WD' = 'HER42PB_WDSS', 
                 'HERPT1S1' = 'HERPT1S1_GRD', 'HERPT1S4' = 'HERPT1S4_GRD',  
                 'HERPT1S8' = 'HERPT1S8_GRD', 'HERPT2S1' = 'HERPT2S1_GRD')

Site_plot_order_violin <- c('HER14L_0ka_GRD', 'HER14L_10ka_GRD','HER24L_GRD', 'HER34L_GRD', 
                            'HER42L_GRD', 'HER42PB_GRD', 'HER44L_GRD', 'HER49L_GRD', 
                            'HERPT1S1_GRD', 'HERPT1S4_GRD', 'HERPT1S8_GRD','HERPT2S1_GRD', 
                            'HER42PB_WDSS')
jco_violin <- c("#0073C2FF", "#0073C2FF", "#003C67FF","#313695", 
                "#7AA6DCFF","#B24745FF", "#003C67FF", "#374E55FF", 
                "#DFC27D", "#BF812D", "#8C510A", 
                "#543005", "#8F7700FF")
show_col(jco_violin)

Site_plot_order_GRD  <- c('HER14L_0ka_GRD', 'HER14L_10ka_GRD','HER24L_GRD', 'HER34L_GRD', 
                     'HER42L_GRD', 'HER42PB_WDSS','HER42PB_GRD', 'HER44L_GRD', 
                     'HER49L_GRD', 'HERPT1S1_GRD', 'HERPT1S4_GRD', 'HERPT1S8_GRD',
                     'HERPT2S1_GRD')
jco1 <- c("#0073C2FF", "#0073C2FF", "#003C67FF","#313695", 
          "#7AA6DCFF","#8F7700FF","#B24745FF", "#003C67FF", 
          "#374E55FF", "#DFC27D", "#BF812D", "#8C510A", 
          "#543005")
show_col(jco1)

Site_plot_order_GRD1 <- c('HER14L_0ka_GRD', 'HER14L_10ka_GRD', 'HER24L_GRD', 
                          'HER34L_GRD','HER42L_GRD', 'HER44L_GRD', 'HER49L_GRD', 
                          'HER42PB_GRD', 'HER42PB_WDSS', 'HERPT1S1_GRD', 'HERPT1S4_GRD',  
                          'HERPT1S8_GRD', 'HERPT2S1_GRD')

northern <- c('HER14L_0ka_GRD', 'HER14L_10ka_GRD','HER24L_GRD')
jcoN <- c("#0073C2FF", "#0073C2FF", "#003C67FF")
show_col(jcoN)

southern <- c('HER42L_GRD', 'HER49L_GRD')
jamaS <- c("#7AA6DCFF", "#003C67FF")
show_col(jamaS)

HER42LP <- c('HER42L_GRD', 'HER42PB_WDSS', 'HER42PB_GRD')
jamaPB <- c("#7AA6DCFF", "#8F7700FF","#B24745FF")
show_col(jamaPB)

jamaHER42PB <- c("#8F7700FF","#B24745FF")
show_col(jamaHER42PB)

peat_bog_S <- c('HER42L_GRD', 'HER42PB_WDSS', 'HER42PB_GRD', 'HER49L_GRD')
jamaPBS <- c("#7AA6DCFF", "#8F7700FF","#B24745FF", "#003C67FF")
show_col(jamaPBS)

ITRAX_plot_order <- c('HER14L_0ka', 'HER14L_10ka', 'HER24L', 'HER34L',
                      'HER42PB', 'HER42L',  'HER44L', 'HER49L')
jco_violin1 <- c("#0073C2FF", "#0073C2FF", "#003C67FF", "#313695",
                 "#B24745FF", "#7AA6DCFF", "#003C67FF","#374E55FF")
show_col(jco_violin1)

northern1 <- c('HER14L_0ka', 'HER14L_10ka', 'HER24L')
jcoN1 <- c("#0073C2FF", "#0073C2FF", "#003C67FF")
show_col(jcoN1)

southern1 <- c('HER42L', 'HER49L')
jamaS1 <- c("#7AA6DCFF", "#003C67FF")
show_col(jamaS1)

HER42LP1 <- c('HER42PB', 'HER42L')
jamaPB1 <- c("#B24745FF", "#7AA6DCFF")
show_col(jamaPB1)

peat_bog_S1 <- c('HER42PB', 'HER42L', 'HER49L')
jamaPBS1 <- c("#B24745FF", "#7AA6DCFF", "#003C67FF")
show_col(jamaPBS1)

# define themes
theme_set(theme_classic(base_size=8))
HER_view_theme8 <- theme(plot.title = element_text(size = 10, face = "bold"),
                   axis.title.x=element_text(size = 8),
                   axis.text.x=element_text(size = 8, color="black",angle=0, hjust=0),
                   axis.text.y=element_text(size = 8, color="black",angle=0, hjust=0),
                   axis.ticks.length = unit(0.25, "cm"),
                   panel.grid.major.x = element_line(color = "grey", size = 0.1, linetype = 1),
                   legend.position = "right"
                   )
HER_view_theme7 <- theme(plot.title = element_text(size = 9, face = "bold"),
                        axis.title.x=element_text(size = 7, color="black"),
                        axis.text.x=element_text(size = 7, color="black", angle=0, hjust=0),
                        axis.text.y=element_text(size = 7, color="black", angle=0, hjust=0),
                        axis.ticks.length = unit(0.12, "cm"),
                        panel.grid.major.x = element_line(color = "grey", size = 0.1, linetype = 1),
                        legend.position = "right"
)

# define patchwork theme items
remove_x <- theme(
  axis.line.x =element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.x = element_blank()
)

remove_y <- theme(
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.title.y = element_blank()
)

# Density plotting -----------------------------------------------------------------

# Set up for plots
HER_MSCL_comp_density <- HER_MSCL_comp %>% 
  filter(!if_any(everything(), is.na)) %>% #remove lines with NAs
  select(Site, label, depth, mean_age, accum_rate, PWAmp_SAT:RES_SAT) %>%
  rename(PWAMp = PWAmp_SAT, Density = Den1_SAT, MS = MS1_SAT, DCMS = DCMS1_SAT, Impedance = Imp_SAT, FP = FP_SAT, Resistivity = RES_SAT) %>% 
  select(Site: mean_age, Density) %>% 
  bind_rows(ACE09_WD) %>%
  mutate(mean_age = mean_age / 1000) %>%
  mutate(Site = recode(Site,'HER14L_0ka' = 'HER14L_0ka_GRD', 'HER14L_10ka' = 'HER14L_10ka_GRD', 'HER24L' = 'HER24L_GRD', 
                       'HER34L' = 'HER34L_GRD','HER42L' = 'HER42L_GRD', 
                       'HER44L' = 'HER44L_GRD', 'HER49L' = 'HER49L_GRD', 
                       'HER42PB' = 'HER42PB_GRD','HER42PB_WD' = 'HER42PB_WDSS', 
                       'HERPT1S1' = 'HERPT1S1_GRD', 'HERPT1S4' = 'HERPT1S4_GRD',  
                       'HERPT1S8' = 'HERPT1S8_GRD', 'HERPT2S1' = 'HERPT2S1_GRD')) %>% 
  mutate(Site = fct_relevel(Site, Site_plot_order_GRD))  #reorder based on site rename list
HER_MSCL_comp_density
write.csv(HER_MSCL_comp_density,"Papers/Hodgson_Hermite_2022/Data/MSCL/Outputs/HER_MSCL_comp_density.csv", row.names = FALSE)
# reorder based on median values
#HER_MSCL_comp1$Site = with(HER_MSCL_comp1 , reorder(Site, Density, median))

# Violin density plots
HER_density_violin <- HER_MSCL_comp_density %>%
  ggplot(aes(x=Site, y=Density, fill=Site)) + 
  geom_violin() +
  xlab("Site") +
  ylab(as.expression(expression(paste("Wet Density " (kg/m ^ 3 )) ))) +
  theme(legend.position="none", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("") +
  scale_fill_manual(values = jco_violin) +
  scale_color_manual(values = jco_violin)
HER_density_violin
ggsave("Papers/Hodgson_Hermite_2022/Figures/HER_density_violin.pdf",
       height = c(6), width = c(12), dpi = 600, units = "cm")

# Density vs age 
HER_Density_N  <- HER_MSCL_comp_density %>% 
  filter(Site==northern) %>% 
  mutate(Site = fct_relevel(Site, northern)) %>% #sets plot order as northern 
  pivot_longer(c(Density), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth) %>% 
  arrange(param) %>% 
  ggplot(aes(x = mean_age, y = value)) +
  geom_line(aes(group = Site, colour = Site), size = 0.5, alpha = 0.3) +
  #geom_point(aes(group = Site, colour = Site), size = 0.5, alpha = 0.2) +
  geom_smooth(aes(group = Site, colour = Site), size = 0.75 , span = 0.1, method = "loess") +
  #scale_fill_jco() +
  #scale_color_jco() +
  scale_fill_manual(values = jcoN) +
  scale_color_manual(values = jcoN) +
  xlab("Age (cal ka BP)") +
  ylab(as.expression(expression(paste("Wet density " (kg/m ^ 3 )) ))) +
  HER_view_theme8 +   
  scale_x_continuous(limits= c(-0.1, 16), breaks=seq(0, 16, 2)) + 
  scale_y_continuous(limits= c(0.5, 2.5), breaks=seq(0.5, 2.5, 0.5)) +
  #geom_vline(xintercept = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16), col = "grey", lty = 3, alpha = 0.5, lwd = 0.3) +
  geom_vline(xintercept = c(4.2, 8.2, 11.75), col = "blue", lty = 3, alpha = 0.5, lwd = 0.5) +
  removeGrid(x = TRUE, y = TRUE)
HER_Density_N

HER_Density_S <- HER_MSCL_comp_density %>% 
  filter(Site == southern) %>%
  mutate(Site = fct_relevel(Site, southern)) %>% 
  pivot_longer(c(Density), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth) %>% 
  arrange(param) %>% 
  ggplot(aes(x = mean_age, y = value)) +
  geom_line(aes(group = Site, colour = Site), size = 0.5, alpha = 0.3) +
  #geom_point(aes(group = Site, colour = Site), size = 0.5, alpha = 0.2) +
  geom_smooth(aes(group = Site, colour = Site), size = 0.75 , span = 0.1, method = "loess") +
  #scale_fill_jama() +
  #scale_color_jama() +
  scale_fill_manual(values = jamaS) +
  scale_color_manual(values = jamaS) +
  xlab("Age (cal a BP)") +
  ylab(as.expression(expression(paste("Wet density " (kg/m ^ 3 )) ))) +
  HER_view_theme8 +   
  scale_x_continuous(limits= c(-0.1, 16), breaks=seq(0, 16, 2)) + 
  scale_y_continuous(limits= c(0.5, 2.5), breaks=seq(0.5, 2.5, 0.5)) +
  #geom_vline(xintercept = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16), col = "grey", lty = 3, alpha = 0.5, lwd = 0.3) +
  geom_vline(xintercept = c(4.2, 8.2, 11.75), col = "blue", lty = 3, alpha = 0.5, lwd = 0.5) +
  removeGrid(x = TRUE, y = TRUE)
HER_Density_S

HER42_Density_S <- HER_MSCL_comp_density %>% 
  bind_rows(ACE09_WD)%>%
  filter(Site==peat_bog_S) %>% #peat_bog
  mutate(Site = fct_relevel(Site, peat_bog_S)) %>%
  pivot_longer(c(Density), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth) %>% 
  arrange(param) %>% 
  ggplot(aes(x = mean_age, y = value)) +
  geom_line(aes(group = Site, colour = Site), size = 0.5, alpha = 0.3) +
  #geom_point(aes(group = Site, colour = Site), size = 0.5, alpha = 0.3) +
  geom_smooth(aes(group = Site, colour = Site), size = 0.75 , span = 0.1, method = "loess") +
  #scale_fill_jama() +
  #scale_color_jama() +
  scale_fill_manual(values = jamaPBS) +
  scale_color_manual(values = jamaPBS) +
  xlab("Age (cal ka BP)") +
  ylab(as.expression(expression(paste("Wet Density " (kg/m ^ 3 )) ))) +
  HER_view_theme8 +   
  scale_x_continuous(limits= c(-0.1, 16), breaks=seq(0, 16, 2)) + 
  scale_y_continuous(limits= c(0.5, 2), breaks=seq(0.5, 2.5, 0.5)) +
  #geom_vline(xintercept = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16), col = "grey", lty = 3, alpha = 0.5, lwd = 0.3) +
  geom_vline(xintercept = c(4.2, 8.2, 11.75), col = "blue", lty = 3, alpha = 0.5, lwd = 0.5) +
  removeGrid(x = TRUE, y = TRUE)
HER42_Density_S

# plot data with patchwork and align y-axis using wrap_plots
library(patchwork)

# Make plot list & plot
HER_density_p1 <- list(HER_Density_N + remove_x, HER_Density_S + remove_x, HER42_Density_S)

HER_density <- wrap_plots(HER_density_p1) +
  plot_layout(
    nrow = 3, ncol = 1) + #widths = c(0.5, 8)) + widths controls the relative width of two or more plots
  plot_annotation(
    title = "HER Density",
    caption = NULL,
    theme = HER_view_theme8
  )
HER_density
ggsave("Papers/Hodgson_Hermite_2022/Figures/HER_density.pdf",
       height = c(12), width = c(12), dpi = 600, units = "cm")


# DCMS plotting -----------------------------------------------------------------

# Set up for plots
HER_MSCL_comp_DCMS <- HER_MSCL_comp %>% 
  select(Site, label, depth, mean_age, accum_rate, PWAmp_SAT:RES_SAT) %>%
  rename(PWAMp = PWAmp_SAT, Density = Den1_SAT, MS = MS1_SAT, DCMS = DCMS1_SAT, Impedance = Imp_SAT, FP = FP_SAT, Resistivity = RES_SAT) %>% 
  select(Site: mean_age, DCMS) %>% 
  #bind_rows(ACE09_WD) %>%
  mutate(mean_age = mean_age / 1000) %>%
  mutate(Site = fct_relevel(Site, Site_plot_order))  #reorder based on site rename list
HER_MSCL_comp_DCMS
write.csv(HER_MSCL_comp_DCMS,"Papers/Hodgson_Hermite_2022/Data/MSCL/Outputs/HER_MSCL_comp_DCMS.csv", row.names = FALSE)
# reorder based on median values
#HER_MSCL_comp1$Site = with(HER_MSCL_comp1 , reorder(Site, DCMS, median))

# Violin DCMS plots
HER_DCMS_violin <- HER_MSCL_comp_DCMS %>%
  ggplot(aes(x=Site, y=DCMS, fill=Site)) + 
  geom_violin() +
  xlab("Site") +
  ylab(as.expression(expression(paste("log10 DCMS " (x10 ^ -8 ~m ^ 3/kg)) ))) +
  theme(legend.position="none", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("") +
  scale_fill_manual(values = jco_violin0) +
  scale_color_manual(values = jco_violin0)
HER_DCMS_violin
ggsave("Papers/Hodgson_Hermite_2022/Figures/HER_DCMS_violin.pdf",
       height = c(6), width = c(12), dpi = 600, units = "cm")


# DCMS vs age 

HER_DCMS_N  <- HER_MSCL_comp_DCMS %>% 
  filter(Site==northern1) %>% 
  mutate(Site = fct_relevel(Site, northern1)) %>% #sets plot order as northern 
  pivot_longer(c(DCMS), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth) %>% 
  arrange(param) %>% 
  ggplot(aes(x = mean_age, y = value)) +
  geom_line(aes(group = Site, colour = Site), size = 0.5, alpha = 0.3) +
  #geom_point(aes(group = Site, colour = Site), size = 0.5, alpha = 0.2) +
  geom_smooth(aes(group = Site, colour = Site), size = 0.75 , span = 0.1, method = "loess") +
  #scale_fill_jco() +
  #scale_color_jco() +
  scale_fill_manual(values = jcoN1) +
  scale_color_manual(values = jcoN1) +
  xlab("Age (cal ka BP)") +
  ylab(as.expression(expression(paste("DCMS " (x10 ^ -8 ~m ^ 3/kg)) ))) +
  HER_view_theme8 +   
  scale_x_continuous(limits= c(-0.1, 16), breaks=seq(0, 16, 2)) + 
  #scale_y_continuous(limits= c(0.5, 2.5), breaks=seq(0.5, 2.5, 0.5)) +
  #geom_vline(xintercept = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16), col = "grey", lty = 3, alpha = 0.5, lwd = 0.3) +
  geom_vline(xintercept = c(4.2, 8.2, 11.75), col = "blue", lty = 3, alpha = 0.5, lwd = 0.5) +
  removeGrid(x = TRUE, y = TRUE)
HER_DCMS_N

HER_DCMS_S <- HER_MSCL_comp_DCMS %>% 
  filter(Site == southern1) %>%
  mutate(Site = fct_relevel(Site, southern1)) %>% 
  pivot_longer(c(DCMS), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth) %>% 
  arrange(param) %>% 
  ggplot(aes(x = mean_age, y = value)) +
  geom_line(aes(group = Site, colour = Site), size = 0.5, alpha = 0.3) +
  #geom_point(aes(group = Site, colour = Site), size = 0.5, alpha = 0.2) +
  geom_smooth(aes(group = Site, colour = Site), size = 0.75 , span = 0.1, method = "loess") +
  #scale_fill_jama() +
  #scale_color_jama() +
  scale_fill_manual(values = jamaS1) +
  scale_color_manual(values = jamaS1) +
  xlab("Age (cal ka BP)") +
  ylab(as.expression(expression(paste("DCMS " (x10 ^ -8 ~m ^ 3/kg)) ))) +
  HER_view_theme8 +   
  scale_x_continuous(limits= c(-0.1, 16), breaks=seq(0, 16, 2)) + 
  #scale_y_continuous(limits= c(0.5, 2.5), breaks=seq(0.5, 2.5, 0.5)) +
  #geom_vline(xintercept = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16), col = "grey", lty = 3, alpha = 0.5, lwd = 0.3) +
  geom_vline(xintercept = c(4.2, 8.2, 11.75), col = "blue", lty = 3, alpha = 0.5, lwd = 0.5) +
  removeGrid(x = TRUE, y = TRUE)
HER_DCMS_S

HER42_DCMS_S <- HER_MSCL_comp_DCMS %>%
  filter(Site==peat_bog_S1) %>% #peat_bog
  mutate(Site = fct_relevel(Site, peat_bog_S1)) %>%
  pivot_longer(c(DCMS), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth) %>% 
  arrange(param) %>% 
  ggplot(aes(x = mean_age, y = value)) +
  geom_line(aes(group = Site, colour = Site), size = 0.5, alpha = 0.3) +
  #geom_point(aes(group = Site, colour = Site), size = 0.5, alpha = 0.3) +
  geom_smooth(aes(group = Site, colour = Site), size = 0.75 , span = 0.1, method = "loess") +
  #scale_fill_jama() +
  #scale_color_jama() +
  scale_fill_manual(values = jamaPBS1) +
  scale_color_manual(values = jamaPBS1) +
  xlab("Age (cal ka BP)") +
  ylab(as.expression(expression(paste("DCMS " (x10 ^ -8 ~m ^ 3/kg)) ))) +
  HER_view_theme8 +   
  scale_x_continuous(limits= c(-0.1, 16), breaks=seq(0, 16, 2)) + 
  #scale_y_continuous(limits= c(0.5, 2.5), breaks=seq(0.5, 2.5, 0.5)) +
  #geom_vline(xintercept = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16), col = "grey", lty = 3, alpha = 0.5, lwd = 0.3) +
  geom_vline(xintercept = c(4.2, 8.2, 11.75), col = "blue", lty = 3, alpha = 0.5, lwd = 0.5) +
  removeGrid(x = TRUE, y = TRUE)
HER42_DCMS_S

# plot data with patchwork and align y-axis using wrap_plots
library(patchwork)

# Make plot list & plot
HER_DCMS_p1 <- list(HER_DCMS_N + remove_x, HER_DCMS_S + remove_x, HER42_DCMS_S)

HER_DCMS <- wrap_plots(HER_DCMS_p1) +
  plot_layout(
    nrow = 3, ncol = 1) + #widths = c(0.5, 8)) + widths controls the relative width of two or more plots
  plot_annotation(
    title = "HER MSCL DCMS",
    caption = NULL,
    theme = HER_view_theme8
  )
HER_DCMS
ggsave("Papers/Hodgson_Hermite_2022/Figures/HER_DCMS.pdf",
       height = c(12), width = c(12), dpi = 600, units = "cm")



# ITRAX Plotting --------------------------------------------------------------

library(tidypaleo)

# Assign elements from PeriodicTable package to 'elementsList' - to select all elements present across all sites
data(periodicTable)
elementsList <- c(periodicTable$symb, "Mo_inc", "Mo_inc")

# Set dataframe up for plots
HER_ITRAX_xrf1_clr_comp1 <- HER_ITRAX_xrf1_clr_comp %>% 
  mutate(mean_age = mean_age / 1000) %>%
  mutate(Site = fct_relevel(Site, ITRAX_plot_order)) %>% 
  mutate(Ln_inc_coh = log(inc_coh), Ln_coh_inc = log(coh_inc)) %>%  
  mutate(Site = fct_relevel(Site, ITRAX_plot_order)) %>%  #set up general plot order
  rename_with(~paste0(.x, "_clr"), Si:Mo_coh) # adds _clr to column names
HER_ITRAX_xrf1_clr_comp1
write.csv(HER_ITRAX_xrf1_clr_comp1,"Papers/Hodgson_Hermite_2022/Data/ITRAX/Outputs/HER_ITRAX_xrf1_clr_comp1.csv", row.names = FALSE)
# reorder based on median values
#HER_MSCL_comp1$Site = with(HER_MSCL_comp1 , reorder(Site, DCMS, median))

# Violin DCMS plots

# reorder based on median values
#HER_MSCL_comp1$Site = with(HER_MSCL_comp1 , reorder(Site, DCMS, median))

# Key elements
HER_ITRAX_violin_grid <- HER_ITRAX_xrf1_clr_comp1 %>%
pivot_longer(c(K_clr, Ca_clr, Ti_clr, Fe_clr, Br_clr, Sr_clr, Mo_inc_clr, Mo_coh_clr), 
             names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth) %>% 
  arrange(param) %>%
  ggplot(aes(x=Site, y=value, fill=Site)) + 
  geom_violin() +
  xlab("Site") +
  facet_geochem_grid(vars(param)) +
  #ylab(as.expression(expression(paste("Ln(K/Ca)")))) +
  theme(legend.position="none", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("") +
  #scale_y_continuous(breaks = seq(-10, 10, 0.5))
  scale_fill_manual(values = jco_violin1) +
  scale_color_manual(values = jco_violin1)
HER_ITRAX_violin_grid
ggsave("Papers/Hodgson_Hermite_2022/Figures/HER_ITRAX_violin_key.pdf",
       height = c(16), width = c(12), dpi = 600, units = "cm")
  
HER_ITRAX_violin_K_Ca <- HER_ITRAX_xrf1_clr_comp1 %>%
  ggplot(aes(x=Site, y=LnK_Ca, fill=Site)) + 
  geom_violin() +
  xlab("Site") +
  ylab(as.expression(expression(paste("Ln(K/Ca)")))) +
  theme(legend.position="none", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("") +
  scale_fill_manual(values = jco_violin1) +
  scale_color_manual(values = jco_violin1) +
  scale_y_continuous(breaks = seq(-10, 10, 0.5))
HER_ITRAX_violin_K_Ca 

HER_ITRAX_violin_coh_inc <- HER_ITRAX_xrf1_clr_comp1 %>%
  ggplot(aes(x=Site, y=Ln_coh_inc, fill=Site)) + 
  geom_violin() +
  xlab("Site") +
  ylab(as.expression(expression(paste("Ln(coh./inc.)")))) +
  theme(legend.position="none", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("") +
  scale_fill_manual(values = jco_violin1) +
  scale_color_manual(values = jco_violin1) +
  scale_y_continuous(breaks = seq(-10, 10, 0.25))
HER_ITRAX_violin_coh_inc

HER_ITRAX_violin_Ti <- HER_ITRAX_xrf1_clr_comp1 %>%
  ggplot(aes(x=Site, y=Ti_clr, fill=Site)) + 
  geom_violin() +
  xlab("Site") +
  ylab(as.expression(expression(paste("Ti [clr]")))) +
  theme(legend.position="none", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("") +
  scale_fill_manual(values = jco_violin1) +
  scale_color_manual(values = jco_violin1) +
  scale_y_continuous(breaks = seq(-10, 10, 1))
HER_ITRAX_violin_Ti 

# plot data with patchwork and align y-axis using wrap_plots
library(patchwork)

# Make plot list & plot
HER_ITRAX_violin <- list(HER_ITRAX_violin_K_Ca + remove_x, HER_ITRAX_violin_Ti + remove_x, HER_ITRAX_violin_coh_inc)

HER_ITRAX_violin3 <- wrap_plots(HER_ITRAX_violin) +
  plot_layout(
    nrow = 3, ncol = 1) + #widths = c(0.5, 8)) + widths controls the relative width of two or more plots
  plot_annotation(
    title = "HER ITRAX Violin plots",
    caption = NULL,
    theme = HER_view_theme8
  )
HER_ITRAX_violin3
ggsave("Papers/Hodgson_Hermite_2022/Figures/HER_ITRAX_violin3.pdf",
       height = c(16), width = c(12), dpi = 600, units = "cm")

# ITRAX age plots 

# northern sites - key elements & ratios 
HER_ITRAX_N  <- HER_ITRAX_xrf1_clr_comp1 %>% 
  filter(Site==northern1) %>% 
  mutate(Site = fct_relevel(Site, northern1)) %>% #sets plot order as northern 
  pivot_longer(c(Ti_clr, LnK_Ca, Ln_coh_inc), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth) %>% 
  arrange(param) %>% 
  ggplot(aes(x = mean_age, y = value)) +
  geom_line(aes(group = Site, colour = Site), size = 0.5, alpha = 0.3) +
  #geom_point(aes(group = Site, colour = Site), size = 0.5, alpha = 0.2) +
  geom_smooth(aes(group = Site, colour = Site), size = 0.75 , span = 0.1, method = "loess") +
  facet_geochem_grid(vars(param)) +
  #scale_fill_jco() +
  #scale_color_jco() +
  scale_fill_manual(values = jcoN1) +
  scale_color_manual(values = jcoN1) +
  xlab("Age (cal ka BP)") +
  #ylab(as.expression(expression(paste("coh./inc.")))) +
  HER_view_theme8 +   
  scale_x_continuous(limits= c(-0.1, 16), breaks=seq(0, 16, 2)) + 
  #scale_y_continuous(limits= c(0, 0.5), breaks=seq(0, 0.5, 0.1)) +
  #geom_vline(xintercept = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16), col = "grey", lty = 3, alpha = 0.5, lwd = 0.3) +
  geom_vline(xintercept = c(4.2, 8.2, 11.75), col = "blue", lty = 3, alpha = 0.5, lwd = 0.5) +
  removeGrid(x = TRUE, y = TRUE)
HER_ITRAX_N
ggsave("Papers/Hodgson_Hermite_2022/Figures/HER_ITRAX_N.pdf",
       height = c(16), width = c(12), dpi = 600, units = "cm")

# coh_inc only
HER_ITRAX_N_coh_inc <- HER_ITRAX_xrf1_clr_comp1 %>% 
  filter(Site == northern1) %>%
  mutate(Site = fct_relevel(Site, northern1)) %>% 
  pivot_longer(c(Ln_coh_inc), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth) %>% 
  arrange(param) %>% 
  ggplot(aes(x = mean_age, y = value)) +
  geom_line(aes(group = Site, colour = Site), size = 0.5, alpha = 0.3) +
  #geom_point(aes(group = Site, colour = Site), size = 0.5, alpha = 0.3) +
  geom_smooth(aes(group = Site, colour = Site), size = 0.75 , span = 0.1, method = "loess") +
  #scale_fill_jco() +
  #scale_color_jco() +
  scale_fill_manual(values = jcoN1) +
  scale_color_manual(values = jcoN1) +
  xlab("Age (cal ka BP)") +
  ylab(as.expression(expression(paste("Ln(coh./inc.)")))) +
  HER_view_theme8 +   
  scale_x_continuous(limits= c(-0.1,16), breaks=seq(0, 16, 2)) + 
  scale_y_continuous(limits= c(-2, -1), breaks=seq(-2, 1.5, 0.25)) +
  #geom_vline(xintercept = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16), col = "grey", lty = 3, alpha = 0.5, lwd = 0.3) +
  geom_vline(xintercept = c(4.2, 8.2, 11.75), col = "blue", lty = 3, alpha = 0.5, lwd = 0.5) +
  removeGrid(x = TRUE, y = TRUE)
HER_ITRAX_N_coh_inc

# coh_inc only
HER_ITRAX_N_LnK_Ca <- HER_ITRAX_xrf1_clr_comp1 %>% 
  filter(Site == northern1) %>%
  mutate(Site = fct_relevel(Site, northern1)) %>% 
  pivot_longer(c(LnK_Ca), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth) %>% 
  arrange(param) %>% 
  ggplot(aes(x = mean_age, y = value)) +
  geom_line(aes(group = Site, colour = Site), size = 0.5, alpha = 0.3) +
  #geom_point(aes(group = Site, colour = Site), size = 0.5, alpha = 0.3) +
  geom_smooth(aes(group = Site, colour = Site), size = 0.75 , span = 0.1, method = "loess") +
  #scale_fill_jco() +
  #scale_color_jco() +
  scale_fill_manual(values = jcoN1) +
  scale_color_manual(values = jcoN1) +
  xlab("Age (cal ka BP)") +
  ylab(as.expression(expression(paste("Ln K_Ca")))) +
  HER_view_theme8 +   
  scale_x_continuous(limits= c(-0.1,16), breaks=seq(0, 16, 2)) + 
  #scale_y_continuous(limits= c(-2, -1), breaks=seq(-2, 1.5, 0.25)) +
  #geom_vline(xintercept = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16), col = "grey", lty = 3, alpha = 0.5, lwd = 0.3) +
  geom_vline(xintercept = c(4.2, 8.2, 11.75), col = "blue", lty = 3, alpha = 0.5, lwd = 0.5) +
  removeGrid(x = TRUE, y = TRUE)
HER_ITRAX_N_LnK_Ca

# Southern lake sites - key elements & ratios 
HER_ITRAX_S  <- HER_ITRAX_xrf1_clr_comp1 %>% 
  filter(Site==southern1) %>% 
  mutate(Site = fct_relevel(Site, southern1)) %>% #sets plot order as northern 
  pivot_longer(c(Ti_clr, LnK_Ca, Ln_coh_inc), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth) %>% 
  arrange(param) %>% 
  ggplot(aes(x = mean_age, y = value)) +
  geom_line(aes(group = Site, colour = Site), size = 0.5, alpha = 0.3) +
  #geom_point(aes(group = Site, colour = Site), size = 0.5, alpha = 0.2) +
  geom_smooth(aes(group = Site, colour = Site), size = 0.75 , span = 0.1, method = "loess") +
  facet_geochem_grid(vars(param)) +
  #scale_fill_jama() +
  #scale_color_jama() +
  scale_fill_manual(values = jamaS1) +
  scale_color_manual(values = jamaS1) +
  xlab("Age (cal ka BP)") +
  #ylab(as.expression(expression(paste("Ln(coh./inc.)")))) +
  HER_view_theme8 +   
  scale_x_continuous(limits= c(-0.1, 16), breaks=seq(0, 16, 2)) + 
  #scale_y_continuous(limits= c(0, 0.5), breaks=seq(0, 0.5, 0.1)) +
  #geom_vline(xintercept = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16), col = "grey", lty = 3, alpha = 0.5, lwd = 0.3) +
  geom_vline(xintercept = c(4.2, 8.2, 11.75), col = "blue", lty = 3, alpha = 0.5, lwd = 0.5) +
  removeGrid(x = TRUE, y = TRUE)
HER_ITRAX_S
ggsave("Papers/Hodgson_Hermite_2022/Figures/HER_ITRAX_S.pdf",
       height = c(16), width = c(12), dpi = 600, units = "cm")

# coh_inc only
HER_ITRAX_S_coh_inc <- HER_ITRAX_xrf1_clr_comp1 %>% 
  filter(Site == southern1) %>%
  mutate(Site = fct_relevel(Site, southern1)) %>% 
  pivot_longer(c(Ln_coh_inc), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth) %>% 
  arrange(param) %>% 
  ggplot(aes(x = mean_age, y = value)) +
  geom_line(aes(group = Site, colour = Site), size = 0.5, alpha = 0.3) +
  #geom_point(aes(group = Site, colour = Site), size = 0.5, alpha = 0.3) +
  geom_smooth(aes(group = Site, colour = Site), size = 0.75 , span = 0.1, method = "loess") +
  #scale_fill_jama() +
  #scale_color_jama() +
  scale_fill_manual(values = jamaS1) +
  scale_color_manual(values = jamaS1) +
  xlab("Age (cal ka BP)") +
  ylab(as.expression(expression(paste("Ln(coh./inc.)")))) +
  HER_view_theme8 +   
  scale_x_continuous(limits= c(-0.1,16), breaks=seq(0, 16, 2)) + 
  #scale_y_continuous(limits= c(-10, 60), breaks=seq(-10, 60, 10)) +
  #geom_vline(xintercept = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16), col = "grey", lty = 3, alpha = 0.5, lwd = 0.3) +
  geom_vline(xintercept = c(4.2, 8.2, 11.75), col = "blue", lty = 3, alpha = 0.5, lwd = 0.5) +
  removeGrid(x = TRUE, y = TRUE)
HER_ITRAX_S_coh_inc

# HER42 Southern sites - key elements & ratios 
HER42_ITRAX_S <- HER_ITRAX_xrf1_clr_comp1 %>% 
  filter(Site==peat_bog_S1 ) %>% 
  mutate(Site = fct_relevel(Site, peat_bog_S1)) %>% #sets plot order as northern 
  pivot_longer(c(Ti_clr, LnK_Ca, Ln_coh_inc), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth) %>% 
  arrange(param) %>% 
  ggplot(aes(x = mean_age, y = value)) +
  geom_line(aes(group = Site, colour = Site), size = 0.5, alpha = 0.3) +
  #geom_point(aes(group = Site, colour = Site), size = 0.5, alpha = 0.2) +
  geom_smooth(aes(group = Site, colour = Site), size = 0.75 , span = 0.1, method = "loess") +
  facet_geochem_grid(vars(param)) +
  #scale_fill_jama() +
  #scale_color_jama() +
  scale_fill_manual(values = jamaPBS1) +
  scale_color_manual(values = jamaPBS1) +
  xlab("Age (cal ka BP)") +
  #ylab(as.expression(expression(paste("Ln(coh./inc.)")))) +
  HER_view_theme8 +   
  scale_x_continuous(limits= c(-0.1, 16), breaks=seq(0, 16, 2)) + 
  #scale_y_continuous(limits= c(0, 0.5), breaks=seq(0, 0.5, 0.1)) +
  #geom_vline(xintercept = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16), col = "grey", lty = 3, alpha = 0.5, lwd = 0.3) +
  geom_vline(xintercept = c(4.2, 8.2, 11.75), col = "blue", lty = 3, alpha = 0.5, lwd = 0.5) +
  removeGrid(x = TRUE, y = TRUE)
HER42_ITRAX_S
ggsave("Papers/Hodgson_Hermite_2022/Figures/HER42_ITRAX_S.pdf",
       height = c(16), width = c(12), dpi = 600, units = "cm")

# coh_inc only
HER42_ITRAX_S_coh_inc <- HER_ITRAX_xrf1_clr_comp1 %>% 
  filter(Site == peat_bog_S1 ) %>%
  mutate(Site = fct_relevel(Site, peat_bog_S1)) %>% 
  pivot_longer(c(Ln_coh_inc), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth) %>% 
  arrange(param) %>% 
  ggplot(aes(x = mean_age, y = value)) +
  geom_line(aes(group = Site, colour = Site), size = 0.5, alpha = 0.3) +
  #geom_point(aes(group = Site, colour = Site), size = 0.5, alpha = 0.3) +
  geom_smooth(aes(group = Site, colour = Site), size = 0.75 , span = 0.1, method = "loess") +
  #scale_fill_jama() +
  #scale_color_jama() +
  scale_fill_manual(values = jamaPBS1) +
  scale_color_manual(values = jamaPBS1) +
  xlab("Age (cal ka BP)") +
  ylab(as.expression(expression(paste("Ln(coh./inc.)")))) +
  HER_view_theme8 +   
  scale_x_continuous(limits= c(-0.1,16), breaks=seq(0, 16, 2)) + 
  #scale_y_continuous(limits= c(-10, 60), breaks=seq(-10, 60, 10)) +
  #geom_vline(xintercept = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16), col = "grey", lty = 3, alpha = 0.5, lwd = 0.3) +
  geom_vline(xintercept = c(4.2, 8.2, 11.75), col = "blue", lty = 3, alpha = 0.5, lwd = 0.5) +
  removeGrid(x = TRUE, y = TRUE)
HER42_ITRAX_S_coh_inc 

# LnK_Ca only
HER42_ITRAX_S_LnK_Ca <- HER_ITRAX_xrf1_clr_comp1 %>% 
  filter(Site == peat_bog_S1 ) %>%
  mutate(Site = fct_relevel(Site, peat_bog_S1)) %>% 
  pivot_longer(c(LnK_Ca), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = depth) %>% 
  arrange(param) %>% 
  ggplot(aes(x = mean_age, y = value)) +
  geom_line(aes(group = Site, colour = Site), size = 0.5, alpha = 0.3) +
  #geom_point(aes(group = Site, colour = Site), size = 0.5, alpha = 0.3) +
  geom_smooth(aes(group = Site, colour = Site), size = 0.75 , span = 0.1, method = "loess") +
  #scale_fill_jama() +
  #scale_color_jama() +
  scale_fill_manual(values = jamaPBS1) +
  scale_color_manual(values = jamaPBS1) +
  xlab("Age (cal ka BP)") +
  ylab(as.expression(expression(paste("Ln K/Ca")))) +
  HER_view_theme8 +   
  scale_x_continuous(limits= c(-0.1,16), breaks=seq(0, 16, 2)) + 
  #scale_y_continuous(limits= c(-10, 60), breaks=seq(-10, 60, 10)) +
  #geom_vline(xintercept = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16), col = "grey", lty = 3, alpha = 0.5, lwd = 0.3) +
  geom_vline(xintercept = c(4.2, 8.2, 11.75), col = "blue", lty = 3, alpha = 0.5, lwd = 0.5) +
  removeGrid(x = TRUE, y = TRUE)
HER42_ITRAX_S_LnK_Ca

# plot data with patchwork and align y-axis using wrap_plots
library(patchwork)

# Make plot list & plot
HER_ITRAX_p1 <- list(HER_ITRAX_N_coh_inc + remove_x, HER_ITRAX_S_coh_inc + remove_x, HER42_ITRAX_S_coh_inc)

HER_ITRAX_coh_inc <- wrap_plots(HER_ITRAX_p1) +
  plot_layout(
    nrow = 3, ncol = 1) + #widths = c(0.5, 8)) + widths controls the relative width of two or more plots
  plot_annotation(
    title = "HER Ln(coh./inc.) (proxy for dry mass)",
    caption = NULL,
    theme = HER_view_theme8
  )
HER_ITRAX_coh_inc
ggsave("Papers/Hodgson_Hermite_2022/Figures/HER_ITRAX_coh_inc.pdf",
       height = c(16), width = c(12), dpi = 600, units = "cm")


# Paper plots -------------------------------------------------------------


# Figure 3 Summary age plots ----------------------------------------------------------------

# Northern sites

# Make plot list & plot
HER_northern_p1 <- list(HER_ITRAX_N_LnK_Ca + remove_x, HER_Density_N + remove_x, HER_ITRAX_N_coh_inc + remove_x, HER_DCMS_N)

HER_northern_4plots <- wrap_plots(HER_northern_p1) +
  plot_layout(
    nrow = 4, ncol = 1) + #widths = c(0.5, 8)) + widths controls the relative width of two or more plots
  plot_annotation(
    title = "Punta Momberg HER14L & HER24L",
    caption = NULL,
    theme = HER_view_theme8
  )
HER_northern_4plots
ggsave("Papers/Hodgson_Hermite_2022/Figures/Fig3B_HER_northern.pdf",
       height = c(16), width = c(10), dpi = 600, units = "cm")


HER_southern_p1 <- list(HER42_ITRAX_S_LnK_Ca + remove_x, HER42_Density_S + remove_x, HER42_ITRAX_S_coh_inc + remove_x, HER42_DCMS_S)

HER_southern_4plots <- wrap_plots(HER_southern_p1) +
  plot_layout(
    nrow = 4, ncol = 1) + #widths = c(0.5, 8)) + widths controls the relative width of two or more plots
  plot_annotation(
    title = "Cabo West (South) HER42L, HER42PB & HER49L",
    caption = NULL,
    theme = HER_view_theme8
  )
HER_southern_4plots
ggsave("Papers/Hodgson_Hermite_2022/Figures/Fig3C_HER_southern.pdf",
       height = c(16), width = c(10), dpi = 600, units = "cm")


# Figure 6 Mini Strat plots ----------------------------------------------------
library(tidypaleo)

theme_set(theme_classic(base_size=7))
#theme_set(theme_paleo(base_size=7))

# MSCL component

# Set up GRD tolerance filter - need to use this carefully as GRD can indicate different lithologies and HER are only peat

HER_MSCL_SAT0 <-  HER_MSCL_comp %>% 
  select(Site, label, depth, mean_age, accum_rate, PWAmp_SAT:RES_SAT) %>%
  rename(PWAMp = PWAmp_SAT, GRD = Den1_SAT, MS = MS1_SAT, DCMS = DCMS1_SAT, Impedance = Imp_SAT, FP = FP_SAT, Resistivity = RES_SAT) %>% 
  select(Site: mean_age, GRD, MS, DCMS) %>% 
  filter(!if_all(GRD:DCMS, is.na))  %>%
  mutate(mean_age = mean_age / 1000)
HER_MSCL_SAT0

# Remove data beyond GRD and MS tolerance - 2 std dev is too strict when GRD values are very similar but all below 2
GRD.mean <- mean(HER_MSCL_SAT0$GRD, na.rm = TRUE)
GRD.sd <- 2*sd(HER_MSCL_SAT0$GRD, na.rm = TRUE)
GRD.min.thres <- 0.5 #GRD.mean - GRD.sd 
GRD.max.thres <- 3 # GRD.mean + GRD.sd 

MS.mean <- mean(HER_MSCL_SAT0$MS, na.rm = TRUE)
MS.sd <- 2*sd(HER_MSCL_SAT0$MS, na.rm = TRUE)
MS.min.thres <- -50 #MS.mean - MS.sd 
MS.max.thres <- 499 #MS.mean + MS.sd

HER_MSCL_SAT <- HER_MSCL_SAT0 %>% 
  mutate(GRD_tolerance = ifelse(GRD <=GRD.min.thres | GRD >=GRD.max.thres | is.na(GRD) == TRUE, FALSE, TRUE)) %>% 
  filter(GRD_tolerance == TRUE) %>% 
  mutate(MS_tolerance = ifelse(MS <=MS.min.thres | MS >=MS.max.thres | is.na(MS) == TRUE, TRUE, TRUE)) %>%
  filter(GRD_tolerance == TRUE) %>% 
  select(-MS, -GRD_tolerance, -MS_tolerance)
HER_MSCL_SAT

# ITRAX %cps sum component - replace with other xfr1 dataframes from above - tolerance filtering already done
plot_elements <- c("K", "Ca", "Lncoh_inc")
HER_ITRAX_pc_cps <-  HER_ITRAX_xrf1_pc_cps_comp  %>% 
  mutate(Lncoh_inc = log(coh_inc)) %>% 
  select(Site: mean_age, all_of(plot_elements)) %>% 
  filter(!if_all(plot_elements, is.na)) %>%
  mutate(mean_age = mean_age / 1000)
HER_ITRAX_pc_cps

# create combined MSCL & ITRAX file
HER_MSCL_ITRAX_pc_cps <- bind_rows(HER_MSCL_SAT, HER_ITRAX_pc_cps)
HER_MSCL_ITRAX_pc_cps

# Site HER14L
HER14L <- c("HER14L_0ka", "HER14L_10ka")
HER14L_MSCL_ITRAX_pc_cps <- HER_MSCL_ITRAX_pc_cps %>% 
  select(-label) %>% 
  filter(Site == HER14L) 
HER14L_MSCL_ITRAX_pc_cps

write.csv(HER14L_MSCL_ITRAX_pc_cps,"Papers/Hodgson_Hermite_2022/Data/HER14L_MSCL_ITRAX_pc_cps.csv", row.names = FALSE)

# plot - to fix omitting 
HER14L_plot <- HER14L_MSCL_ITRAX_pc_cps %>%
  mutate(Site = fct_relevel(Site, HER14L)) %>%
  pivot_longer(c(GRD, DCMS, K, Ca, Lncoh_inc), 
               names_to = "param", values_to = "value") %>%
  drop_na() %>% 
  mutate(param = fct_relevel(param, "GRD", "DCMS", "K", "Ca", "Lncoh_inc")) %>% 
  relocate(param, .before = depth) %>%
  arrange(param) %>% 
  ggplot(aes(x = value, y = depth), size = 0.01) +
  geom_lineh(lwd = 0.3) +
  geom_point(size = 0.01) +
  scale_y_reverse() +
  HER_view_theme7 + 
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  labs(title = "HER14L") +
  # add units to the graph headers
  facet_geochem_gridh(
    vars(param),
    units = c("GRD" = "g cm-3", "DCMS" = "x10^-8 m3 kg-1", "K" = "%cps sum", "Ca" = "%cps sum", "Lncoh_inc" = " "))
HER14L_plot

# OPTIONAL
# add hiatus (red) horizontal line and other points of interest (green) - or comment out

#HER14L_plot1 <- HER14L_plot

# OR

HER14L_plot1 <- HER14L_plot +
  geom_hline(yintercept = c(76), col = "red", lty = 1, alpha = 1, lwd = 0.5)
HER14L_plot1

# Add age on secondary age y axis - TO DO FROM HERE
HER14L_adm <- age_depth_model(
  HER14L_MSCL_ITRAX_pc_cps, 
  depth = depth, age = mean_age)

HER14L_plot_age <- HER14L_plot1 + 
  scale_y_depth_age(HER14L_adm, age_name = "Age (cal a BP)")
HER14L_plot_age

# CONISS
# Uses Den1_SAT, DCMS_SAT only in _com_long3, but plotted with MS1_SAT and WMAR
# CONISS using nested_chclust_coniss() and the rioja package for constrained cluster analysis 
# By default, this calculates the number of significant zones based on a broken stick simulation
# TO DO <- add number of groups and their depth/ages

HER14L_coniss <- HER_MSCL_SAT %>% 
  filter(Site == HER14L) %>%
  mutate(Site = fct_relevel(Site, HER14L)) %>%
  pivot_longer(c(GRD, DCMS), 
               names_to = "param", values_to = "value") %>% 
  mutate(param = fct_relevel(param, "GRD", "DCMS")) %>% 
  relocate(param, .before = depth) %>% 
  arrange(param) %>%
  nested_data(qualifiers = c(mean_age, depth), 
              key = param, value = value, trans = scale) %>%
  nested_chclust_coniss()

HER14L_plot_coniss <- HER14L_plot_age +
  layer_dendrogram(HER14L_coniss, aes(y = depth), 
                   param = "CONISS", size = 0.2, lwd = 0.2) +
  layer_zone_boundaries(HER14L_coniss, aes(y = depth), 
                        size = 0.3, lwd = 0.3, colour = "blue")

ggarrange(HER14L_plot_age, HER14L_plot_coniss, nrow = 2)
ggsave("Outputs/Figures/HER14L/Fig2.1_comp_depth_coniss.pdf",
       height = c(20), width = c(16), dpi = 600, units = "cm")

# Plot versus age
HER14L_age <- HER_MSCL_SAT %>% 
  filter(Site == HER14L) %>%
  mutate(Site = fct_relevel(Site, HER14L)) %>%
  pivot_longer(c(GRD, DCMS), 
               names_to = "param", values_to = "value") %>% 
  mutate(param = fct_relevel(param, "GRD", "DCMS")) %>% 
  relocate(param, .before = depth) %>% 
  arrange(param) %>%
  ggplot(aes(x = mean_age, y = value), size = 0.01) +
  geom_line(lwd = 0.3) +
  geom_point(size = 0.01) +
  #scale_y_reverse() +
  # add units to the graph headers
  facet_geochem_grid(
    vars(param),
    units = c("GRD" = "g cm-3", "DCMS" = "x10^-8 m3 kg-1")) +
  #facet_geochem_grid(vars(param)) +
  labs(x = "Age (cal ka BP)", y = NULL)
HER14L_age








