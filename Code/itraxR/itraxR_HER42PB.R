# HER42PB itrax.R

# install latest github version of ITRAX R
#remotes::install_github("tombishop1/itraxR")

# Set up & clear ------------------------------------------------------------------

#clear previous console
remove (list = ls())
#clear plot window
dev.off()

# setup workspace ----
library(itraxR)
library(tidyverse) # all core tidyverse packages
library(tidypaleo) # Dewey Dunnington's ggplot extensions for palaeo-style plots
library(readr)
library(ggpubr)
library(GGally) # for correlation and Prob density matrix plotting
library(PeriodicTable)
library(errors)
library(chemometrics)
library(patchwork)

# Raw HER data are located here:https://www.dropbox.com/sh/m2iguheb80k1x4m/AACX3Eplgel2LkDGPVyaZcpaa?dl=0

#set working directory
setwd("/Users/Steve/Dropbox/BAS/Data/R/")
#check working directory
getwd() 

# -------------------------------------------------------------------------

# SECTION 1: DATA WRANGLING - run once

# -------------------------------------------------------------------------

# Renaming columns in result.txt file to match itrax.R formats before using itraximport

# kcps renamed as cps (code not used here - see original)
#dfA<- read_tsv("HER42_1A_XRF/result.txt", col_names = TRUE, skip = 2) %>% 
#rename(cps = kcps, `Fe a*2` = D1) %>% 
#select(dfA,-X55) %>% 
#  write_tsv("HER42PB_3A/XRF/Results.tsv")
# Some reanalysis data contains extra column callled revalidity; validity is the original run validity

# Create Results1 file for itraxR -USE

# kcps & cps; cps calculated as element and scatter sum 
# import, calculate cps_sum and rename as cps, retain kcps column - save as Results1.tsv
# df1_rowsums <- c(11:49, 53, 54) # all elements + coh and inc scatter - old import now using mutate ... across

dfA <- read_tsv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1A_XRF/result.txt", col_names = TRUE, skip = 2) %>%
  relocate(`Mo inc`, .after = U) %>%
  relocate(`Mo coh`, .after = `Mo inc`) %>% 
  mutate(cps = rowSums(across(Mg:`Mo coh`))) %>% 
  rename(`Fe a*2` = D1) %>%
  relocate(cps, .before = MSE) %>% 
  select(-`...55`) 
write.table(dfA, "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1A_XRF/Results1.txt", sep = "\t", append = TRUE,
            row.names = FALSE, col.names = TRUE)

dfA2 <- read.table(file = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1A_XRF/result.txt", header = F, nrows = 2, sep = "\t") %>% 
  as_tibble()

dfB <- read_tsv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1B_XRF/result.txt", col_names = TRUE, skip = 2) %>%
  relocate(`Mo inc`, .after = U) %>%
  relocate(`Mo coh`, .after = `Mo inc`) %>% 
  mutate(cps = rowSums(across(Mg:`Mo coh`))) %>% 
  rename(`Fe a*2` = D1) %>%
  relocate(cps, .before = MSE) %>% 
  select(-`...55`) 
write.table(dfA, "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1B_XRF/Results1.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE)

dfC <- read_tsv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1C_XRF/result.txt", col_names = TRUE, skip = 2) %>%
  relocate(`Mo inc`, .after = U) %>%
  relocate(`Mo coh`, .after = `Mo inc`) %>% 
  mutate(cps = rowSums(across(Mg:`Mo coh`))) %>% 
  rename(`Fe a*2` = D1) %>%
  relocate(cps, .before = MSE) %>% 
  select(-`...55`) 
write.table(dfA, "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1C_XRF/Results1.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE, )

dfD <- read_tsv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1D_XRF/result.txt", col_names = TRUE, skip = 2) %>%
  relocate(`Mo inc`, .after = U) %>%
  relocate(`Mo coh`, .after = `Mo inc`) %>% 
  mutate(cps = rowSums(across(Mg:`Mo coh`))) %>% 
  rename(`Fe a*2` = D1) %>%
  relocate(cps, .before = MSE) %>% 
  select(-`...55`) 
write.table(dfA, "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1D_XRF/Results1.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE)

dfE <- read_tsv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1E_XRF/result.txt", col_names = TRUE, skip = 2) %>%
  relocate(`Mo inc`, .after = U) %>%
  relocate(`Mo coh`, .after = `Mo inc`) %>% 
  mutate(cps = rowSums(across(Mg:`Mo coh`))) %>% 
  rename(`Fe a*2` = D1) %>%
  relocate(cps, .before = MSE) %>% 
  select(-`...55`) 
write.table(dfA, "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1E_XRF/Results1.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE)

dfF <- read_tsv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1F_XRF/result.txt", col_names = TRUE, skip = 2) %>%
  relocate(`Mo inc`, .after = U) %>%
  relocate(`Mo coh`, .after = `Mo inc`) %>% 
  mutate(cps = rowSums(across(Mg:`Mo coh`))) %>% 
  rename(`Fe a*2` = D1) %>%
  relocate(cps, .before = MSE) %>% 
  select(-`...55`) 
write.table(dfA, "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1F_XRF/Results1.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE)

dfG <- read_tsv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1G_XRF/result.txt", col_names = TRUE, skip = 2) %>%
  relocate(`Mo inc`, .after = U) %>%
  relocate(`Mo coh`, .after = `Mo inc`) %>% 
  mutate(cps = rowSums(across(Mg:`Mo coh`))) %>% 
  rename(`Fe a*2` = D1) %>%
  relocate(cps, .before = MSE) %>% 
  select(-`...55`) 
write.table(dfA, "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1G_XRF/Results1.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE)

dfH <- read_tsv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1H_XRF/result.txt", col_names = TRUE, skip = 2) %>%
  relocate(`Mo inc`, .after = U) %>%
  relocate(`Mo coh`, .after = `Mo inc`) %>% 
  mutate(cps = rowSums(across(Mg:`Mo coh`))) %>% 
  rename(`Fe a*2` = D1) %>%
  relocate(cps, .before = MSE) %>% 
  select(-`...55`) 
write.table(dfA, "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1H_XRF/Results1.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE)

dfI <- read_tsv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1I_XRF/result.txt", col_names = TRUE, skip = 2) %>%
  relocate(`Mo inc`, .after = U) %>%
  relocate(`Mo coh`, .after = `Mo inc`) %>% 
  mutate(cps = rowSums(across(Mg:`Mo coh`))) %>% 
  rename(`Fe a*2` = D1) %>%
  relocate(cps, .before = MSE) %>% 
  select(-`...55`) 
write.table(dfA, "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1I_XRF/Results1.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE)

dfJ <- read_tsv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1J_XRF/result.txt", col_names = TRUE, skip = 2) %>%
  relocate(`Mo inc`, .after = U) %>%
  relocate(`Mo coh`, .after = `Mo inc`) %>% 
  mutate(cps = rowSums(across(Mg:`Mo coh`))) %>% 
  rename(`Fe a*2` = D1) %>%
  relocate(cps, .before = MSE) %>% 
  select(-`...55`) 
write.table(dfA, "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1J_XRF/Results1.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE)

# AT END
# open Results1.tsv in excel 
# add first two rows back into all new txt files using copy/paste
# rename .tsv files to .txt in Finder
# now ready to import into itrax.R and runs as normal
# choose whether to use Results.txt or Results1.txt
# below uses Results1.txt to match Manchester output from new ITRAX software

# remove from global environment
rm(dfA, dfB, dfC, dfD, dfE, dfF,  dfG, dfH, dfI, dfJ, df1_rowsums)

# -------------------------------------------------------------------------

# SECTION 2 - ITRAX data import

# -------------------------------------------------------------------------

# Pre-processing notes:
# - corrected field depth scale to match field report depths 
# - Strat depth scale based on ITRAX inc/coh and Ti data matches
# - document.txt files - start/end optical depths have been adjusted to match optical1 image
# - optical 1 image is brightness +150 then auto tone, contrast and colour corrected
# - optical 2 is cropped vertically to XRF start/end and start/end of the core and ~3 cm edge cropped horizontally 
# - optical 3 is optical 2 auto tone, auto contrast and auto colour corrected
# - field logs are as below
# - radiograph1 is radiograph autotone and autocontrast - positive - >density = darker 
# - radiograph1 is radiograph1 inverted - negative - >density = lighter

# Data import notes
# - Import ITRAX XRF data, image, radiographs for each core using itrax_import ##
# - ITRAX folders need to be named as e.g., HER42PB_1A_XRF, HER42PB_1A_XRAD etc
# - Use start depth from e.g., HER42PB_ITRAX_OVERLAPS.csv file
# - Depth is different from position on xrf scanner
# - itrax.R retains all data in overlaps and merges data - different to e.g., HER42PB_ITRAX_COMPOSITE.csv file

HER42PB1A <- list(metadata   = itrax_meta("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1A_XRF/document.txt"),
                xrf        = itrax_import("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1A_XRF/Results1.txt", 
                                          depth = 11, 
                                          parameters = "all"),
                image      = itrax_image(file = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1A_XRF/optical3.tif",
                                         meta = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1A_XRF/document.txt"),
                radiograph = itrax_radiograph(file = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1A_XRAD/radiograph1.tif",
                                              meta = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1A_XRAD/document.txt")
)

# trim = as.numeric(itrax_meta("CD166_19_S1/CD166_19_S1/document.txt")[6:7,2])))                                                 
# Core 1B
HER42PB1B <- list(metadata   = itrax_meta("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1B_XRF/document.txt"),
                xrf        = itrax_import("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1B_XRF/Results1.txt", 
                                          depth = 410, 
                                          parameters = "all"),
                image      = itrax_image(file = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1B_XRF/optical3.tif",
                                         meta = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1B_XRF/document.txt"),
                radiograph = itrax_radiograph(file = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1B_XRAD/radiograph1.tif",
                                              meta = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1B_XRAD/document.txt")
)

# core 1C
HER42PB1C <- list(metadata   = itrax_meta("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1C_XRF/document.txt"),
                xrf        = itrax_import("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1C_XRF/Results1.txt", 
                                          depth = 810, 
                                          parameters = "all"),
                image      = itrax_image(file = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1C_XRF/optical3.tif",
                                         meta = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1C_XRF/document.txt"),
                radiograph = itrax_radiograph(file = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1C_XRAD/radiograph1.tif",
                                              meta = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1C_XRAD/document.txt")
)

# core 1D
HER42PB1D <- list(metadata   = itrax_meta("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1D_XRF/document.txt"),
                xrf        = itrax_import("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1D_XRF/Results1.txt", 
                                          depth = 1210, 
                                          parameters = "all"),
                image      = itrax_image(file = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1D_XRF/optical3.tif",
                                         meta = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1D_XRF/document.txt"),
                radiograph = itrax_radiograph(file = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1D_XRAD/radiograph1.tif",
                                              meta = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1D_XRAD/document.txt")
)

# core 1E
HER42PB1E <- list(metadata   = itrax_meta("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1E_XRF/document.txt"),
                xrf        = itrax_import("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1E_XRF/Results1.txt", 
                                          depth = 1610, 
                                          parameters = "all"),
                image      = itrax_image(file = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1E_XRF/optical3.tif",
                                         meta = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1E_XRF/document.txt"),
                radiograph = itrax_radiograph(file = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1E_XRAD/radiograph1.tif",
                                              meta = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1E_XRAD/document.txt")
)

# core 1F
HER42PB1F <- list(metadata   = itrax_meta("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1F_XRF/document.txt"),
                xrf        = itrax_import("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1F_XRF/Results1.txt", 
                                          depth = 2010, 
                                          parameters = "all"),
                image      = itrax_image(file = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1F_XRF/optical3.tif",
                                         meta = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1F_XRF/document.txt"),
                radiograph = itrax_radiograph(file = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1F_XRAD/radiograph1.tif",
                                              meta = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1F_XRAD/document.txt")
)

# core 1G
HER42PB1G <- list(metadata   = itrax_meta("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1G_XRF/document.txt"),
                xrf        = itrax_import("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1G_XRF/Results1.txt", 
                                          depth =2410, 
                                          parameters = "all"),
                image      = itrax_image(file = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1G_XRF/optical3.tif",
                                         meta = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1G_XRF/document.txt"),
                radiograph = itrax_radiograph(file = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1G_XRAD/radiograph1.tif",
                                              meta = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1G_XRAD/document.txt")
)  

# core 1H
HER42PB1H <- list(metadata   = itrax_meta("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1H_XRF/document.txt"),
                xrf        = itrax_import("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1H_XRF/Results1.txt", 
                                          depth = 2810, 
                                          parameters = "all"),
                image      = itrax_image(file = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1H_XRF/optical3.tif",
                                         meta = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1H_XRF/document.txt"),
                radiograph = itrax_radiograph(file = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1H_XRAD/radiograph1.tif",
                                              meta = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1H_XRAD/document.txt")
)

# core 1I
HER42PB1I <- list(metadata   = itrax_meta("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1I_XRF/document.txt"),
                xrf        = itrax_import("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1I_XRF/Results1.txt", 
                                          depth = 3210, 
                                          parameters = "all"),
                image      = itrax_image(file = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1I_XRF/optical3.tif",
                                         meta = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1I_XRF/document.txt"),
                radiograph = itrax_radiograph(file = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1I_XRAD/radiograph1.tif",
                                              meta = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1I_XRAD/document.txt")
)

# core 1J
HER42PB1J <- list(metadata   = itrax_meta("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1J_XRF/document.txt"),
                xrf        = itrax_import("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1J_XRF/Results1.txt", 
                                          depth = 3560, 
                                          parameters = "all"),
                image      = itrax_image(file = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1J_XRF/optical3.tif",
                                         meta = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1J_XRF/document.txt"),
                radiograph = itrax_radiograph(file = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1J_XRAD/radiograph1.tif",
                                              meta = "Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB1J_XRAD/document.txt")
)

# combining just the xrf data 
HER42PB_xrf_AM <- itrax_join(list(A = HER42PB1A$xrf, 
                             B = HER42PB1B$xrf, 
                             C = HER42PB1C$xrf,
                             D = HER42PB1D$xrf,
                             E = HER42PB1E$xrf,
                             F = HER42PB1F$xrf,
                             G = HER42PB1G$xrf,
                             H = HER42PB1H$xrf,
                             I = HER42PB1I$xrf,
                             J = HER42PB1J$xrf))

# Intensity to cps conversion if needed
# Assign elements to/from PeriodicTable package to 'elementsList' (then unload package)
#data(periodicTable)
#elementsList <- periodicTable$symb
#current <- as.numeric(BIM9_Mo1$metadata[18,2])
#BIM9_Mo1 <- BIM9_Mo1 %>%
#  select(any_of(elementsList), `Mo inc`, `Mo coh`)
#BIM9_Mo1_xrf <-  mutate(round(BIM9_Mo1_xrf*current))
#rm(current)

# filter to remove validity = FALSE rows
HER42PB_xrf_all <- filter(HER42PB_xrf_AM, validity=='TRUE') %>% 
  relocate(label, .before = MSE) %>% 
  rename(Mo_inc = `Mo coh`) %>% 
  rename(Mo_coh = `Mo inc`) %>% 
  rename(Fe_a2 = `Fe a*2`) %>% 
  rename(sample_surface = `sample surface`)

HER42PB_xrf_all

# Save xrf dataset to file with no additional quality controls 
# i.e., if happy to use ITRAX validity = TRUE in as measured data
write.csv(HER42PB_xrf_all,"Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/itraxR/HER42PB_xrf_all_cps.csv", row.names = FALSE)

# -------------------------------------------------------------------------

# SECTION 3.1: Data QC filtering - run once

# -------------------------------------------------------------------------

# Setting up data filtering criteria using quality control measures in itrax.R (if required)
# can skip this section if happy to accept validity = TRUE from ITRAX measurements  

# Assign elements to/from PeriodicTable package to 'elementsList' (then unload package)
data(periodicTable)
elementsList <- periodicTable$symb

# list of all possible elements only
data(periodicTable)
allelements <- c(symb(1:117))
rm(periodicTable)

#  MSE and cps filtering

# cps filtering - Fe a*2 & adjust cps to between 10,000-20,000
Fig3.1 <- ggplot(data = HER42PB_xrf_all, mapping = aes(x = cps, y = Fe_a2)) + 
  geom_point(alpha = 0.1) + 
  theme_bw()
print(Fig3.1)

# cps - 2 std dev is too strict for HER42PB and other peat cores with a combination of low and high count matrices
cps.mean <- mean(HER42PB_xrf_all$cps)
cps.sd <- 6*sd(HER42PB_xrf_all$cps)
cps.min.thres <- cps.mean - cps.sd 
cps.max.thres <- cps.mean + cps.sd 

#  OR use this 

# cps tolerance filter 
# cps.min.thres <- 40000
# cps.max.thres <- 90000

Fig3.2 <- HER42PB_xrf_all  %>%
  mutate(in_cps_tolerance = ifelse(cps <=cps.min.thres | cps >=cps.max.thres | is.na(cps) == TRUE, FALSE, TRUE)) %>%
  ggplot(mapping = aes(x = depth, y = cps, col = in_cps_tolerance)) + 
  geom_line(aes(group = 1)) +
  scale_x_reverse() +
  geom_hline(yintercept = c(cps.min.thres, cps.max.thres)) +
  geom_rug(sides = "b", data = . %>% filter(in_cps_tolerance == FALSE)) + 
  theme_bw()
print(Fig3.2)

# MSE tolerance filter - 3 SD or 2 used here but need to use this carefully as MSE can indicate different lithologies 

# MSE tolerance - 2 std dev is too strict when MSE values are very similar but all below 2
MSE.mean <- mean(HER42PB_xrf_all$MSE)
MSE.sd <- 4*sd(HER42PB_xrf_all$MSE)
MSE.thres <- MSE.mean + MSE.sd 
MSE.thres

#  OR

#MSE.thres <- 2 # use this as an established general threshold

Fig3.3 <- HER42PB_xrf_all %>%
  mutate(in_mse_tolerance = ifelse(MSE >=MSE.thres, FALSE, TRUE)) %>% 
  ggplot(mapping = aes(x = depth,  y = MSE, col = in_mse_tolerance)) +
  geom_line(aes(group = 1)) +
  scale_x_reverse() +
  geom_hline(yintercept = MSE.thres) +
  geom_rug(sides = "b", data = . %>% filter(in_mse_tolerance == FALSE)) + 
  theme_bw()
print(Fig3.3)

# Surface slope tolerance filter

# use this with wide margins for peat cores eg 6 std dev equivalent to +/-0.5 
#slope.min.thres = -0.5
#slope.max.thres = 0.5

#  OR
slope1 <-  HER42PB_xrf_all$sample_surface - lag(HER42PB_xrf_all$sample_surface)
s1 <- as_tibble(slope1) %>% 
  filter(!if_any(everything(), is.na))
slope.mean <- mean(s1$value)
slope.sd <- 4*sd(s1$value)
slope.min.thres <- slope.mean - slope.sd 
slope.max.thres <- slope.mean + slope.sd 

Fig3.4 <- HER42PB_xrf_all %>%
  mutate(slope = sample_surface - dplyr::lag(sample_surface)) %>%
  mutate(in_slope_tolerance = ifelse(slope <=slope.min.thres | slope >=slope.max.thres | is.na(slope) == TRUE, FALSE, TRUE)) %>% 
  ggplot(mapping = aes(x = depth, y = slope, col = in_slope_tolerance)) +
  scale_y_continuous(limits = c(-1, 1), oob = scales::squish) +
  geom_line(aes(group = 1)) +
  geom_hline(yintercept = c(slope.min.thres, slope.max.thres)) +
  geom_rug(data = . %>% filter(validity == FALSE)) +
  scale_x_reverse() +
  theme_bw()
print(Fig3.4)

# show Figure 3.1 - 3.4 together:
ggarrange(Fig3.1, Fig3.2, Fig3.3, Fig3.4, nrow = 2, ncol = 2, labels = c('a', 'b', 'c', 'd'))
ggsave("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Figures/Fig3.1_tolerance_combined.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Combining 'validity' flags

# important to tailor to each record 
# commented out slope error b/c of overlapping cores and >3 cores (can only run 3 at once) causing high slope error
HER42PB_xrf_qc <- HER42PB_xrf_all %>%
  mutate(slope = sample_surface - dplyr::lag(sample_surface)) %>%
  mutate(in_slope_tolerance = ifelse(slope <=slope.min.thres | slope >=slope.max.thres | is.na(slope) == TRUE, FALSE, TRUE)) %>%
  select(-slope) %>%
  mutate(in_cps_tolerance = ifelse(cps <=cps.min.thres | cps >=cps.max.thres | is.na(cps) == TRUE, FALSE, TRUE)) %>%
  mutate(in_mse_tolerance = ifelse(MSE <=MSE.thres, TRUE, FALSE)) %>%
  rowwise() %>%
  mutate(qc = !any(c(validity, in_slope_tolerance, in_cps_tolerance, in_mse_tolerance) == FALSE)) %>%
  ungroup() %>%
  select(-c(in_slope_tolerance, in_cps_tolerance, in_mse_tolerance)) # %>% filter(qc == TRUE) #to remove from HER42PB_xrf_all rows that dont pass QC


# plot summary
theme_set(theme_bw(12))
Fig3.5 <- ggplot(data = HER42PB_xrf_qc, aes(y = depth, x = Ti, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.width = unit(0.5, 'cm'),
        legend.box.spacing = unit(0.1, 'cm')) +
  scale_color_discrete(name = "Pass QC") +
  scale_y_reverse(name = "Depth [cm]")

Fig3.6 <- ggplot(data = HER42PB_xrf_qc, aes(y = depth, x = Fe, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.width = unit(0.5, 'cm'),
        legend.box.spacing = unit(0.1, 'cm')) +
  scale_color_discrete(name = "Pass QC") +
  scale_y_reverse(name = "Depth [cm]")

Fig3.7 <- ggplot(data = HER42PB_xrf_qc, aes(y = depth, x = Br, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.width = unit(0.5, 'cm'),
        legend.box.spacing = unit(0.1, 'cm')) +
  scale_color_discrete(name = "Pass QC") +
  scale_y_reverse(name = "Depth [cm]")

Fig3.8 <- ggplot(data = HER42PB_xrf_qc, aes(y = depth, x = Mo_inc/Mo_coh, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.width = unit(0.5, 'cm'),
        legend.box.spacing = unit(0.1, 'cm')) +
  scale_color_discrete(name = "Pass QC") +
  scale_y_reverse(name = "Depth [cm]")

# Smoothing - eby running mean xample

HER42PB_xrf_qc_smooth <- HER42PB_xrf_qc %>%
  # uses a 10 point running mean (1 cm for 1 mm data); 5 before, 5 after
  mutate(across(any_of(elementsList), 
                function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
  )
  ) 

Fig3.9a <- ggplot(HER42PB_xrf_qc_smooth, mapping = aes(y = depth, x = Ca)) + 
  geom_line(data = ARD_xrf_qc, col = "grey80") + 
  geom_line() + 
  scale_y_reverse(name = "Depth [cm]")
Fig3.9b <- ggplot(HER42PB_xrf_qc_smooth, mapping = aes(y = depth, x = K)) + 
  geom_line(data = ARD_xrf_qc, col = "grey80") + 
  geom_line() + 
  scale_y_reverse(name = "Depth [cm]")

ggarrange(Fig3.5, Fig3.6, Fig3.7, Fig3.8,Fig3.9a, Fig3.9b, ncol = 2, nrow = 3, common.legend = TRUE)
ggsave("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Figures/itraxR/Fig3.2_tolerance_smoothing.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

ggsave("Figures/Fig14_Ca_smoothed.pdf", 
       height = c(15), width = c(30), dpi = 600, units = "cm")

# Save files
write.csv(HER42PB_xrf_qc,"Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB_xrf_qc_cps.csv", row.names = FALSE)
write.csv(HER42PB_xrf_qc_smooth,"Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB_xrf_qc_cps_smooth.csv", row.names = FALSE)

# Remove machine elements and known matrix effect elements from dataframe
machine_elements <- select(HER42PB_xrf_qc, c(Ar, Ta, W)) %>% 
  names()
machine_elements

# Create list of elements - can remove Zr (matrix effect in peat/organic seds) at this stage 
HER42PB_elements <- select(HER42PB_xrf_qc, any_of(elementsList), Mo_inc, Mo_coh, -c(all_of(machine_elements), Fe_a2)) %>% # , Zr
  names()
HER42PB_elements

# REE if in dataframe
#REE <- select(HER42PB_xrf_qc, c(La:Ho)) %>% 
#  names()
#REE

# Choose dataset to take forward to element filtering section -------------------------------------------------------------------------

# Simple QC filtering based on validity = TRUE only
HER42PB_xrf <- HER42PB_xrf_all

#  OR use below - use this most of the time

# Complex QC filtering based on multiple parameters above
HER42PB_xrf <- HER42PB_xrf_qc

# Summary stats for xrf
HER42PB_xrf_stats <- HER42PB_xrf %>%
  select(kcps:MSE, any_of(HER42PB_elements)) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()
tail(HER42PB_xrf_stats)

# Save QC xrf dataset and stats to file
write.csv(HER42PB_xrf,"Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/itraxR/HER42PB_xrf.csv", row.names = FALSE)
write.csv(HER42PB_xrf_stats,"Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/itraxR/HER42PB_xrf_stats.csv", row.names = FALSE)


## Correcting for dead time & Uncertainties 
# This doesn't work on cores scanned at Aber as there is no Dt (Dwell time) column in results.txt file
# ggplot(data = HER42PB_xrf, aes(x = depth, y = Dt)) +
#  scale_x_reverse() + 
#  scale_y_continuous(sec.axis = sec_axis( trans=~(.+(1-mean(HER42PB_xrf$Dt, na.rm = TRUE))), name="Correction Factor")) +
#  geom_line() +
#  geom_hline(yintercept = mean(HER42PB_xrf$Dt, na.rm = TRUE), linetype = "dotted")

## Correcting for water content and grain size
# use Mo inc/Mo coh as a proxy to correct for water content 
# in minerogenic seidments Zr can be used as a proxy for grain size changes
# for Mo tube analysis and in organic rich sediments Zr is a matrix affected element due to proximity of spectra peak to Mo

## Uncertainties
# Use repeat runs from other duplicate cores for this section to assess errors
# see code in itraxR book - to do when needed

# -------------------------------------------------------------------------

# SECTION 3.2: Measurable element filtering - run once

# -------------------------------------------------------------------------

#  Autocorrelation based filtering of elements

# Use autocorrelation function (acf) and plots to explore noise in a time-series
library(forecast)
library(ggrepel)
library(directlabels)

# Adjust for any element of interest by changing $
# split into two groups below to visualise the most common elements measured by ITRAX
theme_set(theme_bw(8))
Fig3.10 <- ggarrange(
  ggAcf(HER42PB_xrf$Al) + ylim(c(NA,1)), ggAcf(HER42PB_xrf$Si) + ylim(c(NA,1)), 
  ggAcf(HER42PB_xrf$P) + ylim(c(NA,1)), ggAcf(HER42PB_xrf$S) + ylim(c(NA,1)),
  ggAcf(HER42PB_xrf$Cl) + ylim(c(NA,1)), ggAcf(HER42PB_xrf$K) + ylim(c(NA,1)), 
  ggAcf(HER42PB_xrf$Ca) + ylim(c(NA,1)), ggAcf(HER42PB_xrf$Ti) + ylim(c(NA,1)),
  ggAcf(HER42PB_xrf$V) + ylim(c(NA,1)), ggAcf(HER42PB_xrf$Cr) + ylim(c(NA,1)),
  ggAcf(HER42PB_xrf$Mn) + ylim(c(NA,1)), ggAcf(HER42PB_xrf$Fe) + ylim(c(NA,1)),
  nrow = 4, ncol = 3)
print(Fig3.10)
ggsave("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Figures/itraxR/Fig3.3_ACF_pt1.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

Fig3.11 <- ggarrange(
  ggAcf(HER42PB_xrf$Co) + ylim(c(NA,1)), ggAcf(HER42PB_xrf$Ni) + ylim(c(NA,1)),
  ggAcf(HER42PB_xrf$Cu) + ylim(c(NA,1)), ggAcf(HER42PB_xrf$Zn) + ylim(c(NA,1)),
  ggAcf(HER42PB_xrf$Se) + ylim(c(NA,1)), ggAcf(HER42PB_xrf$Br) + ylim(c(NA,1)),
  ggAcf(HER42PB_xrf$Rb) + ylim(c(NA,1)), ggAcf(HER42PB_xrf$Sr) + ylim(c(NA,1)),
  ggAcf(HER42PB_xrf$Zr) + ylim(c(NA,1)), ggAcf(HER42PB_xrf$Ba) + ylim(c(NA,1)),
  ggAcf(HER42PB_xrf$Mo_inc) + ylim(c(NA,1)), ggAcf(HER42PB_xrf$Mo_coh) + ylim(c(NA,1)),
  nrow = 4, ncol = 3)
print(Fig3.11)
ggsave("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Figures/itraxR/Fig3.4_ACF_pt2.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Filter elements based on acf lag thresholds

# set lag threshold to an arbitary value of 5 or 20 for whole dataset
# for 1 mm dataset, 5 checks autocorrelation +/-5 mm around each measurement 
# can be performed on different units later
# 0.2 or 0.1 is usual minimum threshold

# define filter and lag thresholds
acf_thres_min <- 0.1
acf_thres_max <- 0.5
lag_thres <- 20

Fig3.12a <- apply(HER42PB_xrf %>% select(any_of(HER42PB_elements), Mo_inc, Mo_coh), 
                  2, FUN = function(x){round(Acf(x, plot = F)$acf, 3)}) %>%
  as_tibble(rownames = "lag") %>%
  pivot_longer(!c("lag"), names_to = "elements", values_to = "value") %>%
  mutate(lag = as.numeric(lag),
         elements = factor(elements, levels = filter(., lag == lag_thres) %>% 
                             arrange(desc(value)) %>% pull(elements))) %>%
  group_by(elements) %>%
  ggplot(aes(x = lag, y = value, col = elements)) +
  geom_line() +
  geom_hline(yintercept= c(acf_thres_min, acf_thres_max), color="red", size=0.5, linetype = 3) +
  xlim(0, 42) +
  geom_dl(aes(label = elements), method = list(dl.trans(x = x + 0.5), "last.points", cex = 0.8)) # adds element labels to end
print(Fig3.12a)
ggsave("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Figures/itraxR/Fig3.5a_ACF_all_elements.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# remove machine elements and apply acf based filtering to - leaving acf elements
apply(HER42PB_xrf %>% select(any_of(HER42PB_elements), Mo_inc, Mo_coh), 2, FUN = function(x){round(Acf(x, plot = F)$acf, 3)}) %>%
  as_tibble(rownames = "lag") %>% 
  pivot_longer(!c("lag"), names_to = "elements", values_to = "value") %>% 
  mutate(lag = as.numeric(lag),
         elements = factor(elements, levels = filter(., lag == lag_thres) %>% 
                             arrange(desc(value)) %>% 
                             pull(elements))) %>%
  group_by(elements) %>%
  filter(lag == lag_thres) %>%
  filter(value >= acf_thres_min) %>% 
  pull(elements) %>% 
  ordered() -> acfElements 
acfElements

acfElementsList <- select(HER42PB_xrf, any_of(acfElements)) %>% 
  names()
acfElementsList

# Replot with acf filtered elements only
Fig3.12b <- apply(HER42PB_xrf %>% select(any_of(acfElements)), 
                  2, FUN = function(x){round(Acf(x, plot = F)$acf, 3)}) %>%
  as_tibble(rownames = "lag") %>%
  pivot_longer(!c("lag"), names_to = "elements", values_to = "value") %>%
  mutate(lag = as.numeric(lag),
         elements = factor(elements, levels = filter(., lag == lag_thres) %>% 
                             arrange(desc(value)) %>% pull(elements))) %>%
  group_by(elements) %>%
  ggplot(aes(x = lag, y = value, col = elements)) +
  geom_line() +
  geom_hline(yintercept= c(acf_thres_min, acf_thres_max), color="red", size=0.5, linetype = 3) +
  xlim(0, 42) +
  geom_dl(aes(label = elements), method = list(dl.trans(x = x + 0.5), "last.points", cex = 0.8)) # adds element labels to end
print(Fig3.12b)
ggsave("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Figures/itraxR/Fig3.5b_ACF_elements.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Summary acf element plot vs depth
library(tidypaleo)
theme_set(theme_bw(base_size=8))

Fig3.13 <- HER42PB_xrf %>% 
  filter(qc == TRUE) %>% 
  select(any_of(acfElements), depth, label) %>%
  pivot_longer(any_of(acfElements), names_to = "param", values_to = "element") %>%
  filter(param %in% acfElements) %>%
  mutate(param = fct_relevel(param, acfElementsList)) %>%
  ggplot(aes(x = element, y = depth)) +
  geom_lineh(aes(color = label)) +
  #geom_point(size = 0.01) + #don't use for ITRAX - too many datapoints 
  #geom_lineh(size = 0.5) + #this will make a single black line plot
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(param)) +
  labs(x = "Peak area [cps]", y = "Depth [cm]") +
  ggtitle("HER42PB-ITRAX: cps (ACF-filtered elements)")
print(Fig3.13)
ggsave("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Figures/itraxR/Fig3.6_ACF_elements_depth.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# nested ggarrange to produce split level summary plot with different number of plots pe row
ggarrange(Fig3.13, # First row
          ggarrange(Fig3.12a, Fig3.12b, ncol = 2, labels = c("B", "C")), # Second row with two plots
          nrow = 2, 
          labels = "A", common.legend = TRUE)
ggsave("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Figures/itraxR/Fig3.7_ACF_elements_depth.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# ACF elements - further filtering based on mean and max cps and %cps sum values

# check acf element list
acfElementsList

# cps filtering based on mean >50 cps, max>100 cps

# mean cps > 50 function
HER42PB_mean50 <- function(x){
  is.numeric(x) && (mean(x, na.rm = TRUE) > 50)
}

# mean cps >50 dataframe 
HER42PBmean50 <- select(HER42PB_xrf, depth:kcps,
                       any_of(acfElementsList) &
                         where(HER42PB_mean50))
# mean cps >50 element list - remove Zr for peat as it's a matrix effect
HER42PBmean50List <- select(HER42PBmean50, -c(depth:kcps, Zr)) %>%
  names() %>% 
  print()

# max cps >100
HER42PB_max100 <- function(x){
  is.numeric(x) && (max(x, na.rm = TRUE) > 100)
}
HER42PBmax100 <- select(HER42PB_xrf, depth:kcps,
                       any_of(acfElementsList) &
                         where(HER42PB_max100))
HER42PBmax100List <- select(HER42PBmax100, -c(depth:kcps, Zr)) %>%
  names() %>% 
  print()

# mean cps >50 and max cps >100
HER42PBmean50_max100 <- select(HER42PBmean50, depth:kcps,
                              any_of(acfElementsList) &
                                where(HER42PB_max100))
HER42PBmean50_max100List <- select(HER42PBmean50_max100, -c(depth:kcps, Zr)) %>%
  names() %>% 
  print()

# Select cps filtered dataframe to take forward
HER42PB_xrf_filter <- HER42PBmean50_max100
HER42PB_xrf_filter

HER42PB_filterList <- HER42PBmean50_max100List
HER42PB_filterList

# Generate cps summary stats for filtered dataframe
HER42PB_xrf_filter_stats <- HER42PB_xrf_filter %>%
  select(kcps:cps, any_of(HER42PB_filterList)) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()

# Save qc element filtered dataset & stats to file
write.csv(HER42PB_xrf_filter,"Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/itraxR/HER42PB_xrf_filter.csv", row.names = FALSE)
write.csv(HER42PB_xrf_filter_stats,"Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/itraxR/HER42PB_xrf_filter_stats.csv", row.names = FALSE)

# Create a new dataframe called xrf1 of qc, acf and filtered elements and add  additional useful parameters
HER42PB_xrf[HER42PB_xrf == 0] <- NA #set 0 to NA log value calculations

HER42PB_xrf1 <- HER42PB_xrf %>%
  select(depth:kcps, any_of(HER42PB_filterList), Mo_inc, Mo_coh) %>% 
  mutate(TS_sum = Mo_inc + Mo_coh) %>%
  mutate(inc_coh = Mo_inc / Mo_coh) %>%
  mutate(coh_inc = Mo_coh / Mo_inc)%>%
  mutate(LnK_Ti = log(K / Ti)) %>% 
  mutate(LnK_Ca = log(K / Ca)) %>% 
  mutate(LnTi_Ca = log(Ti / Ca)) %>% 
  relocate(label, .after = depth) %>%
  relocate(K, .after = Cl) %>%
  relocate(Ti, .after = Ca) %>%
  relocate(TS_sum, .after = Mo_coh) %>%
  relocate(inc_coh, .after = TS_sum) %>%
  relocate(coh_inc, .after = inc_coh) %>% 
  relocate(LnK_Ti, .after = coh_inc) %>% 
  relocate(LnK_Ca, .after = LnK_Ti) %>% 
  relocate(LnTi_Ca, .after = LnK_Ca)
HER42PB_xrf1

#set NA back to 0 for _xrf dataframe
HER42PB_xrf <- HER42PB_xrf %>% 
  replace(is.na(.), 0) 
HER42PB_xrf

# new element and ratio lists
HER42PB_elements1 <- select(HER42PB_xrf1, any_of(elementsList), Mo_inc, Mo_coh) %>% 
  names()
HER42PB_elements1

HER42PB_ratioList1 <- select(HER42PB_xrf1, inc_coh, coh_inc, LnK_Ti, LnK_Ca, LnTi_Ca) %>% 
  names()
HER42PB_ratioList1

# Summary stats for xrf1
HER42PB_xrf1_stats <- HER42PB_xrf1 %>%
  select(any_of(HER42PB_elements1), any_of(HER42PB_ratioList1)) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()
tail(HER42PB_xrf1_stats)

write.csv(HER42PB_xrf1,"Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/itraxR/HER42PB_xrf1.csv", row.names = FALSE)
write.csv(HER42PB_xrf1_stats,"Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/itraxR/HER42PB_xrf1_stats.csv", row.names = FALSE)

# -------------------------------------------------------------------------

# SECTION 3.3: Data Transformation using filtered cps dataframe - run once

# -------------------------------------------------------------------------

# cps as % of cps_sum

# original cps qc dataframe
HER42PB_xrf_pc_cps0 <- HER42PB_xrf_qc %>%
  select(depth:kcps, Mg:Mo_coh, qc) %>%
  mutate(across(Mg:Mo_coh) /`cps`) %>%
  mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
  replace(is.na(.), 0) %>%
  mutate(across(Mg:Mo_coh) *100) %>% 
  print()
# check sum = 100%
HER42PB_xrf_pc_cps1 <- HER42PB_xrf_pc_cps0 %>% 
  mutate(cps_sum_all = rowSums(across(Mg:Mo_coh))) 
HER42PB_xrf_pc_cps2 <- HER42PB_xrf_pc_cps1 %>% 
  mutate(cps_sum_filtered = rowSums(across(all_of(HER42PB_elements1)))) 
head(HER42PB_xrf_pc_cps2$cps_sum_all)
head(HER42PB_xrf_pc_cps2$cps_sum_filtered)

# %cps_sum summary stats table for all ITRAX elements
HER42PB_xrf_pc_cps_stats <- HER42PB_xrf_pc_cps1 %>%
  select(any_of(HER42PB_elements), cps_sum_all) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()
tail(HER42PB_xrf_pc_cps_stats)
write.csv(HER42PB_xrf_pc_cps_stats,"Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/itraxR/HER42PB_xrf_pc_cps_stats.csv", row.names = FALSE)

# Filtered elements as % of cps_sum
HER42PB_xrf1_pc_cps_filtered <- select(HER42PB_xrf_pc_cps2, 
                                      depth:kcps,all_of(HER42PB_elements1), cps_sum_all, cps_sum_filtered)

HER42PB_xrf1_pc_cps <- select(HER42PB_xrf1, inc_coh:LnTi_Ca) %>% 
  bind_cols(HER42PB_xrf1_pc_cps_filtered) %>%
  relocate(inc_coh:LnTi_Ca, .after = Mo_coh)
HER42PB_xrf1_pc_cps

# %cps_sum summary stats table for filtered ITRAX elements
HER42PB_xrf1_pc_cps_stats <- HER42PB_xrf1_pc_cps %>%
  select(any_of(HER42PB_elements1), cps_sum_all, cps_sum_filtered) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()
write.csv(HER42PB_xrf1_pc_cps_stats,"Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/itraxR/HER42PB_xrf1_pc_cps_stats.csv", row.names = FALSE)


# Normalised datasets

# Normalised by inc scatter
HER42PB_xrf1_inc <- HER42PB_xrf1 %>%
  mutate(across(any_of(HER42PB_elements1)) /`Mo_inc`) %>%
  mutate_if(is.numeric, list(~na_if(., Inf)))
HER42PB_xrf1_inc

# Normalised by Ti
HER42PB_xrf1_Ti <- HER42PB_xrf1 %>% 
  mutate(across(any_of(HER42PB_elements1), .fns = ~./ Ti))

# Natural log of Ti normalised dataframe with -Inf replaced by NA
HER42PB_xrf1_Ln_Ti <- HER42PB_xrf1_Ti %>% 
  mutate(across(any_of(HER42PB_elements1),log)) 
is.na(HER42PB_xrf1_Ti)<-sapply(HER42PB_xrf1_Ti, is.infinite)
# Replace NA with 0 - if needed
# HER42PB_Ln_Ti_norm[is.na(HER42PB_Ln_Ti_norm)]<-0
# can replace Ti normalized with other normalisation parameters above

#Check = 1
head(HER42PB_xrf1_Ti$Ti)
#Check = 0
head(HER42PB_xrf1_Ln_Ti$Ti)

# Z-scores

HER42PB_elementsZ <- c(HER42PB_elements1, HER42PB_ratioList1)

# Standardised and centre cps dataframe
HER42PB_xrf1_Z <- HER42PB_xrf1
HER42PB_xrf1_Z[, HER42PB_elementsZ] <- scale(HER42PB_xrf1[, HER42PB_elementsZ], center = TRUE, scale = TRUE)
HER42PB_xrf1_Z

# Standardise and centre inc normalised dataframe
HER42PB_xrf1_inc_Z <- HER42PB_xrf1_inc
HER42PB_xrf1_inc_Z[, HER42PB_elementsZ] <- scale(HER42PB_xrf1_inc[, HER42PB_elementsZ], center = TRUE, scale = TRUE)
HER42PB_xrf1_inc_Z

# Standardize and center Ln Ti-normalized data
HER42PB_xrf1_Ln_Ti_Z <- HER42PB_xrf1_Ln_Ti
HER42PB_xrf1_Ln_Ti_Z[, HER42PB_elementsZ] <- scale(HER42PB_xrf1_Ln_Ti[, HER42PB_elementsZ], center = TRUE, scale = TRUE)
HER42PB_xrf1_Ln_Ti_Z

#set NA back to 0 for _xrf dataframe
HER42PB_xrf1 <- HER42PB_xrf1 %>% 
  replace(is.na(.), 0) 
HER42PB_xrf1

# Save original (xrf), acf filtered (xrf1) and transfromed cps dataframes to csv
write.csv(HER42PB_xrf,"Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/itraxR/HER42PB_xrf.csv", row.names = FALSE)
write.csv(HER42PB_xrf1,"Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/itraxR/HER42PB_xrf1.csv", row.names = FALSE)
write.csv(HER42PB_xrf1_pc_cps,"Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/itraxR/HER42PB_xrf1_pc_cps.csv", row.names = FALSE)
write.csv(HER42PB_xrf1_inc,"Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/itraxR/HER42PB_xrf1_inc.csv", row.names = FALSE)
write.csv(HER42PB_xrf1_Ti,"Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/itraxR/HER42PB_xrf1_Ti.csv", row.names = FALSE)
write.csv(HER42PB_xrf1_Ln_Ti,"Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/itraxR/HER42PB_xrf1_Ln_Ti.csv", row.names = FALSE)
write.csv(HER42PB_xrf1_Z,"Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/itraxR/HER42PB_xrf1_Z.csv", row.names = FALSE)
write.csv(HER42PB_xrf1_inc_Z,"Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/itraxR/HER42PB_xrf1_inc_Z.csv", row.names = FALSE)
write.csv(HER42PB_xrf1_Ln_Ti_Z,"Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/itraxR/HER42PB_xrf1_Ln_Ti_Z.csv", row.names = FALSE)


# %cps sum filtering

# mean %cps sum > 0.05%

HER42PB_mean1 <- function(x){
  is.numeric(x) && (mean(x, na.rm = TRUE) > 0.05)
}

HER42PBmean1 <- select(HER42PB_xrf1_pc_cps, depth:kcps,
                      any_of(acfElementsList) &
                        where(HER42PB_mean1))
HER42PBmean1List <- select(HER42PBmean1, -c(depth:kcps)) %>%
  names() %>% 
  print()

# max %cps sum >0.5%
HER42PB_max1 <- function(x){
  is.numeric(x) && (max(x, na.rm = TRUE) > 0.5)
}
HER42PBmax1 <- select(HER42PB_xrf1_pc_cps,depth:kcps,
                     any_of(acfElementsList) &
                       where(HER42PB_max1))
HER42PBmax1List <- select(HER42PBmax1, -c(depth:kcps)) %>%
  names() %>% 
  print()

# mean %cps sum > 0.05% and max %cps sum >0.5%
HER42PBmean1_max1 <- select(HER42PBmean1,depth:kcps,
                           any_of(acfElementsList) &
                             where(HER42PB_max1))
HER42PBmean1_max1List <- select(HER42PBmean1_max1, -c(depth:kcps)) %>%
  names() %>% 
  print()

# Select cps filtered dataframe to take forward
HER42PB_xrf_filter_pc <- HER42PBmean1_max1
HER42PB_xrf_filter_pc

HER42PB_filterList_pc <- HER42PBmean1_max1List
HER42PB_filterList_pc

# Filtered elements as % of cps_sum
HER42PB_xrf_pc_cps3 <- HER42PB_xrf_pc_cps2 %>% 
  mutate(cps_sum_filtered_pc = rowSums(across(all_of(HER42PB_filterList_pc)))) 
HER42PB_xrf2_pc_cps_filtered <- select(HER42PB_xrf_pc_cps3, 
                                      depth:kcps,all_of(HER42PB_filterList_pc), cps_sum_all, cps_sum_filtered, cps_sum_filtered_pc) %>% 
  relocate(K, .before = Ca)
HER42PB_xrf2_pc_cps_filtered

HER42PB_xrf2_pc_cps <- select(HER42PB_xrf1, inc_coh:LnTi_Ca) %>% 
  bind_cols(HER42PB_xrf2_pc_cps_filtered) %>%
  relocate(inc_coh:LnTi_Ca, .after = Mo_coh)
HER42PB_xrf2_pc_cps

# %cps_sum summary stats table for filtered ITRAX elements
HER42PB_xrf2_pc_cps_stats <- HER42PB_xrf2_pc_cps %>%
  select(any_of(HER42PB_filterList_pc), cps_sum_all, cps_sum_filtered, cps_sum_filtered_pc) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()
write.csv(HER42PB_xrf2_pc_cps,"Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/itraxR/HER42PB_xrf2_pc_cps.csv", row.names = FALSE)
write.csv(HER42PB_xrf2_pc_cps_stats,"Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/itraxR/HER42PB_xrf2_pc_cps_stats.csv", row.names = FALSE)


# clr (centered log ratios) dataset - filtered cps pc elements

library(compositions)

# Cannot run clr with zeroes or NAs
# removing rows with zeros creates too many gaps in the record for peat cores 
# so +1 to all cps
HER42PB_clr0 <- HER42PB_xrf1 %>%
  replace(is.na(.), 0) %>%
  select(K, any_of(HER42PB_filterList_pc))
HER42PB_clr1 <- HER42PB_clr0 + 1 #can run without +1 but this reduces the number of elements to 8
as_tibble(HER42PB_clr1)

# Or replace all 0 with NA and run clr
#HER42PB_clr1[HER42PB_clr1 == 0] <- NA
HER42PB_clr1 <- as_tibble(HER42PB_clr1) %>%
  clr() 
head(HER42PB_clr1)
tail(HER42PB_clr1)

# Add columns back into the clr dataframe
HER42PB_xrf1_clr <- select(HER42PB_xrf1, depth:kcps, inc_coh:LnTi_Ca) %>% 
  bind_cols(HER42PB_clr1) %>%
  relocate(inc_coh:LnTi_Ca, .after = Mo_coh) %>% 
  print()

write.csv(HER42PB_xrf1_clr,"Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/itraxR/HER42PB_xrf1_clr.csv", row.names = FALSE)

# -------------------------------------------------------------------------

# SECTION 4: START HERE ONCE DATASETS ABOVE HAVE BEEN ESTABLISHED 

# -------------------------------------------------------------------------

#set working directory
setwd("/Users/Steve/Dropbox/BAS/Data/R/Papers/Hodgson_Hermite_2022/Data/ITRAX/")
#check working directory
getwd() 

# used below
HER42PB_xrf <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/itraxR/HER42PB_xrf.csv")
HER42PB_xrf1 <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/itraxR/HER42PB_xrf1.csv")
HER42PB_xrf1_pc_cps <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/itraxR/HER42PB_xrf1_pc_cps.csv")
HER42PB_xrf2_pc_cps <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/itraxR/HER42PB_xrf2_pc_cps.csv")
HER42PB_xrf1_inc_Z <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/itraxR/HER42PB_xrf1_inc_Z.csv")
HER42PB_xrf1_clr <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/itraxR/HER42PB_xrf1_clr.csv")

# not used below
#HER42PB_xrf_all <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/itraxR/HER42PB_xrf_all.csv")
#HER42PB_xrf1_inc <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/itraxR/HER42PB_xrf1_inc.csv")
#HER42PB_xrf1_Ti <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/itraxR/HER42PB_xrf1_Ti.csv")
#HER42PB_xrf1_Ln_Ti <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/itraxR/HER42PB_xrf1_Ln_Ti.csv")
#HER42PB_xrf1_Z <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/itraxR/HER42PB_xrf1_Z.csv")
#HER42PB_xrf1_Ln_Ti_Z <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/itraxR/HER42PB_xrf1_Ln_Ti_Z.csv")

# Assign elements to/from PeriodicTable package to 'elementsList' (then unload package)
data(periodicTable)
elementsList <- periodicTable$symb

# list of all possible elements only
data(periodicTable)
allelements <- c(symb(1:117))
rm(periodicTable)

# Create lists of elements
HER42PB_elements <- select(HER42PB_xrf, c(Mg:Mo_coh)) %>% # , Zr
  names()
HER42PB_elements

# acf & cps filtered elements
HER42PB_elements1 <- select(HER42PB_xrf1, any_of(elementsList), Mo_inc, Mo_coh, inc_coh, coh_inc) %>% 
  names()
HER42PB_elements1

# acf, cps & %cps_sum filtered elements
HER42PB_elements2 <- select(HER42PB_xrf2_pc_cps, any_of(elementsList), Mo_inc, Mo_coh, inc_coh, coh_inc) %>% 
  names()
HER42PB_elements2

# additional ratio list
HER42PB_ratioList1 <- select(HER42PB_xrf1, inc_coh, coh_inc, LnK_Ti, LnK_Ca, LnTi_Ca) %>% 
  names()
HER42PB_ratioList1

#set NA back to 0 for _xrf1 dataframe
HER42PB_xrf1 <- HER42PB_xrf1 %>% 
  replace(is.na(.), 0) 
HER42PB_xrf1

HER42PB_xrf1_inc_Z <- HER42PB_xrf1_inc_Z %>% 
  replace(is.na(.), 0) 
HER42PB_xrf1_inc_Z

# Pivot transformed dataframes to long format for ggplot using filetered elements 2 (%cps sum filtered)
HER42PB_xrf_long <- select(HER42PB_xrf,  depth, label, any_of(HER42PB_elements2)) %>%
  pivot_longer(any_of(HER42PB_elements2), names_to = "param", values_to = "element")
HER42PB_xrf1_long <- select(HER42PB_xrf1,  depth, label, any_of(HER42PB_elements2)) %>%
  pivot_longer(any_of(HER42PB_elements2), names_to = "param", values_to = "element")
HER42PB_xrf1_inc_Z_long <- select(HER42PB_xrf1_inc_Z,  depth, label, MSE, kcps, cps, any_of(HER42PB_elements2)) %>%
  pivot_longer(any_of(HER42PB_elements2), names_to = "param", values_to = "element")
HER42PB_xrf1_pc_cps_long <- select(HER42PB_xrf1_pc_cps,  depth, label, any_of(HER42PB_elements2)) %>%
  pivot_longer(any_of(HER42PB_elements2), names_to = "param", values_to = "element")
HER42PB_xrf2_pc_cps_long <- select(HER42PB_xrf2_pc_cps,  depth, label, any_of(HER42PB_elements2)) %>%
  pivot_longer(any_of(HER42PB_elements2), names_to = "param", values_to = "element")
HER42PB_xrf1_clr_long <- select(HER42PB_xrf1_clr,  depth, label, any_of(HER42PB_elements2)) %>%
  pivot_longer(any_of(HER42PB_elements2), names_to = "param", values_to = "element")

# Not used  --------------------------------------------------------
#HER42PB_xrf_all <- select(HER42PB_xrf_all,  depth, label, any_of(HER42PB_elements)) %>%
#  pivot_longer(any_of(HER42PB_elements2), names_to = "param", values_to = "element")
#HER42PB_xrf1_inc_long <- select(HER42PB_xrf1_inc,  depth, label, any_of(HER42PB_elements2)) %>%
#  pivot_longer(any_of(HER42PB_elements2), names_to = "param", values_to = "element")
#HER42PB_xrf1_Z_long <- select(HER42PB_xrf1_Z,  depth, label, MSE, kcps, cps, any_of(HER42PB_elements2)) %>%
#  pivot_longer(any_of(HER42PB_elements2), names_to = "param", values_to = "element")
#HER42PB_xrf1_Ln_Ti_Z_long <- select(HER42PB_xrf1_Ln_Ti_Z,  depth, label, MSE, kcps, cps, any_of(HER42PB_elements2)) %>%
#  pivot_longer(any_of(HER42PB_elements2), names_to = "param", values_to = "element")

# Correlation plots - only run once -------------------------------------------------------------------------

library(GGally)
 
# Correlation matrices

# cps - all elements
theme_set(theme_bw(base_size=2))
ggcorr(HER42PB_xrf[,HER42PB_elements], method = c("everything", "pearson"), 
       size = 3, label = TRUE, label_alpha = TRUE, label_round=2, label_size= 2)
ggsave("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Figures/itraxR/Fig4.1a_Corr_matrix_xrf1_cps.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# cps - filtered elements2 
theme_set(theme_bw(base_size=2))
ggcorr(HER42PB_xrf1[,HER42PB_elements2], method = c("everything", "pearson"), 
       size = 3, label = TRUE, label_alpha = TRUE, label_round=2, label_size = 3)
ggsave("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Figures/itraxR/Fig4.1b_Corr_matrix_xrf1_cps.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# cps/inc as Z-scores - filtered elements2 
theme_set(theme_bw(base_size=2))
ggcorr(HER42PB_xrf1_inc_Z[,HER42PB_elements2], method = c("everything", "pearson"), 
       size = 3, label = TRUE, label_alpha = TRUE, label_round=2, label_size = 3)
ggsave("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Figures/itraxR/Fig4.3_Corr_matrix_incZ.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# %cps sum - filtered elements2 
theme_set(theme_bw(base_size=2))
ggcorr(HER42PB_xrf1_pc_cps[,HER42PB_elements2], method = c("everything", "pearson"), 
       size = 3, label = TRUE, label_alpha = TRUE, label_round=2, label_size = 3)
ggsave("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Figures/itraxR/Fig4.2_Corr_matrix_xrf1_pc_cps_sum.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# clr - filtered elements2 
ggcorr(HER42PB_xrf1_clr[,HER42PB_elements2], method = c("everything", "pearson"), 
       size = 3, label = TRUE, label_alpha = TRUE, label_round=2, label_size = 3) 
ggsave("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Figures/itraxR/Fig4.4_Corr_matrix_xrf1_clr.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Correlation density matrix plots
theme_set(theme_bw(base_size=8))
# cps - filtered elements2 
ggpairs(HER42PB_xrf1, columns = HER42PB_elements2, upper = list(continuous = wrap("cor", size = 3)),
        lower = list(continuous = wrap("points", alpha = 0.5, size=0.2)),
        title="Correlation-density plot: cps for ACF filtered elements")
ggsave("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Figures/itraxR/Fig4.5_Corr-den_matrix_xrf1_cps.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# cps/inc as Z-scores - filtered elements2 
ggpairs(HER42PB_xrf1_inc_Z, columns = HER42PB_elements2, upper = list(continuous = wrap("cor", size = 3)),
        lower = list(continuous = wrap("points", alpha = 0.5, size=0.2)),
        title="Correlation-density plot: cps/inc as Z-scores for ACF filtered elements")
ggsave("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Figures/itraxR/Fig4.7_Corr-den_matrix_xrf1_Z.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# %cps sum - filtered elements2 
ggpairs(HER42PB_xrf1_pc_cps, columns = HER42PB_elements2, upper = list(continuous = wrap("cor", size = 3)),
        lower = list(continuous = wrap("points", alpha = 0.5, size=0.2)),
        title="Correlation-density plot: %cps sum for ACF filtered elements")
ggsave("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Figures/itraxR/Fig4.6_Corr-den_matrix_xrf1_pc_cps_sum.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# clr - filtered elements2
ggpairs(HER42PB_xrf1_clr, columns = HER42PB_elements2, upper = list(continuous = wrap("cor", size = 3)),
        lower = list(continuous = wrap("points", alpha = 0.5, size=0.2)),
        title="Correlation-density plot: clr of ACF filtered elements")
ggsave("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Figures/itraxR/Fig4.8_Corr-den_matrix_xrf1_clr.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

library(tidypaleo)
theme_set(theme_paleo(8))


# -------------------------------------------------------------------------

# SECTION 5: Summary depth plots

# -------------------------------------------------------------------------

library(tidypaleo)
theme_set(theme_paleo(base_size=8))

# set plot text themes for image/rad plots
HER42PB_data_theme <- theme(plot.title = element_text(size = 8),
                             axis.title.x=element_text(size = 8),
                             axis.text.x=element_text(size = 8))

# cps
HER42PB_xrfStrat_cps <- HER42PB_xrf1_long  %>%
  filter(param %in% HER42PB_elements2) %>%
  mutate(param = fct_relevel(param, HER42PB_elements2)) %>%
  ggplot(aes(x = element, y = depth)) +
  geom_lineh(aes(color = label)) +
  #geom_point(size = 0.01) + #don't use for ITRAX - too many datapoints 
  #geom_lineh(size = 0.5) + #this will make a single black line plot
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(param)) +
  labs(x = "Peak area [cps]", y = "Depth [mm]") +
  ggtitle("HER42PB-itraxR: cps (filtered elements)") +
  HER42PB_data_theme

# cps/inc as Z scores
HER42PB_xrfStrat_inc_Z <- HER42PB_xrf1_inc_Z_long  %>% 
  filter(param %in% HER42PB_elements2) %>%
  mutate(param = fct_relevel(param, HER42PB_elements2)) %>%
  ggplot(aes(x = element, y = depth)) +
  geom_lineh(aes(color = label)) +
  #geom_point(size = 0.01) +
  #geom_lineh(size = 0.5) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(param)) +
  labs(x = "cps/inc [Z-score]", y = "Depth [mm]") +
  ggtitle("HER42PB-itraxR: cps/inc as Z-scores (filtered elements)") +
  HER42PB_data_theme

# %cps_sum
HER42PB_xrfStrat_pc_cps <- HER42PB_xrf1_pc_cps_long  %>% 
  filter(param %in% HER42PB_elements2) %>%
  mutate(param = fct_relevel(param, HER42PB_elements2)) %>%
  ggplot(aes(x = element, y = depth)) +
  geom_lineh(aes(color = label)) +
  #geom_point(size = 0.01) +
  #geom_lineh(size = 0.5) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(param)) +
  labs(x = "%cps sum [%]", y = "Depth [mm]") +
  ggtitle("HER42PB-itraxR: %cps sum (filtered elements)") +
  HER42PB_data_theme

# clr
HER42PB_xrfStrat_clr <- HER42PB_xrf1_clr_long  %>% 
  filter(param %in% HER42PB_elements2) %>%
  mutate(param = fct_relevel(param, HER42PB_elements2)) %>%
  ggplot(aes(x = element, y = depth)) +
  geom_lineh(aes(color = label)) +
  #geom_point(size = 0.01) +
  #geom_lineh(size = 0.5) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(param)) +
  labs(x = "clr", y = "Depth [mm]") +
  ggtitle("HER42PB-itraxR: clr - centred log ratio (filtered elements)") +
  HER42PB_data_theme

# Summary plots - data only 
ggarrange(HER42PB_xrfStrat_cps, HER42PB_xrfStrat_inc_Z, nrow = 2)
ggsave("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Figures/itraxR/Fig5.1_Summary_cps_incZ_depth_plots.pdf",
       height = c(30), width = c(30), dpi = 600, units = "cm")

ggarrange(HER42PB_xrfStrat_pc_cps, HER42PB_xrfStrat_clr, nrow = 2)
ggsave("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Figures/itraxR/Fig5.2_Summary_cpssum_clr_depth_plots.pdf",
       height = c(30), width = c(30), dpi = 600, units = "cm")


# CONISS - this takes a long time -----------------------------------------------------

# cps_sum + CONISS
coniss1_HER42PB_pc_cps <- HER42PB_xrf1_pc_cps_long %>%
  nested_data(qualifiers = c(depth, label), key = param, value = element, trans = scale) %>%
  nested_chclust_coniss()

HER42PB_xrfStrat_pc_cps_CONISS <- HER42PB_xrfStrat_pc_cps +
  layer_dendrogram(coniss1_HER42PB_pc_cps, aes(y = depth), param = "CONISS") +
  layer_zone_boundaries(coniss1_HER42PB_pc_cps, aes(y = depth, colour = "grey", lty = 2, alpha = 0.4)) +
  ggtitle("HER42PB-itraxR: %cps sum (filtered elements) with CONISS")

# cps_sum + CONISS
coniss2_HER42PB_clr <- HER42PB_xrf1_clr_long %>%
  nested_data(qualifiers = c(depth, label), key = param, value = element, trans = scale) %>%
  nested_chclust_coniss()

HER42PB_xrfStrat_clr_CONISS <- HER42PB_xrfStrat_clr +
  layer_dendrogram(coniss2_HER42PB_clr, aes(y = depth), param = "CONISS") +
  layer_zone_boundaries(coniss2_HER42PB_clr, aes(y = depth, colour = "grey", lty = 2, alpha = 0.4)) +
  ggtitle("HER42PB-itraxR: clr (filtered elements) with CONISS") +
  HER42PB_data_theme

ggarrange(HER42PB_xrfStrat_pc_cps_CONISS, HER42PB_xrfStrat_clr_CONISS, nrow = 2)
ggsave("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Figures/itraxR/Fig5.3_Summary_cpssum_clr_CONISS_depth_plots.pdf",
       height = c(30), width = c(30), dpi = 600, units = "cm")

# -------------------------------------------------------------------------

# SECTION 6:  Optical and radiograph images

# -------------------------------------------------------------------------

# First run ggItraxImage and ggItraxRadio functions 

# to display images alternating side by side with radiography overlaid on top of the optical image
# create ggItraxImage and ggItraxRadio functions to reduce size of code with annotation_custom etc ##

# setting up ggItraxImage function
ggItraxImage <- function(section, xposition){
  
  if(xposition == "l"){
    annotation_custom(rasterGrob(section$image$image,
                                 width = unit(1, "npc"),
                                 height = unit(1, "npc")),
                      ymax = max(section$xrf$depth),
                      ymin = min(section$xrf$depth),
                      xmin = min(as.numeric(colnames(section$image$image))),
                      xmax = max(as.numeric(colnames(section$image$image)))
    )} else if(xposition == "r"){
      annotation_custom(rasterGrob(section$image$image,
                                   width = unit(1, "npc"),
                                   height = unit(1, "npc")),
                        ymax = max(section$xrf$depth),
                        ymin = min(section$xrf$depth),
                        xmin = max(as.numeric(colnames(section$image$image))),
                        xmax = max(as.numeric(colnames(section$image$image)))*2)
    } else{stop(`l`)}
}

# setting up ggItraxRadio function
ggItraxRadio <- function(section, xposition){
  if(xposition == "l"){
    annotation_custom(rasterGrob(section$radiograph$image,
                                 width = unit(1, "npc"),
                                 height = unit(1, "npc")),
                      ymin = min(section$xrf$depth),
                      ymax = max(section$xrf$depth),
                      xmin = mean(as.numeric(colnames(section$image$image))) -
                        mean(as.numeric(colnames(section$radiograph$image))),
                      xmax = mean(as.numeric(colnames(section$image$image))) +
                        mean(as.numeric(colnames(section$radiograph$image)))
    )
  } else if(xposition == "r"){
    annotation_custom(rasterGrob(section$radiograph$image,
                                 width = unit(1, "npc"),
                                 height = unit(1, "npc")),
                      ymin = min(section$xrf$depth),
                      ymax = max(section$xrf$depth),
                      xmin = max(as.numeric(colnames(section$image$image))) +
                        mean(as.numeric(colnames(section$image$image))) -
                        mean(as.numeric(colnames(section$radiograph$image))),
                      xmax = max(as.numeric(colnames(section$image$image))) +
                        mean(as.numeric(colnames(section$image$image))) +
                        mean(as.numeric(colnames(section$radiograph$image)))
    )
  } else{stop(`l`)}
}

# set plot text themes for image/rad plots
HER42PB_image_theme <- theme(plot.title = element_text(size = 8),
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank())

# Plot optical image  -------------------------------------------------

# ggplot using ggItraxImage of overlapping cores (side-by-side)
theme_set(theme_paleo(8))
HER42PB_Opt <- ggplot() +
  scale_y_continuous(limits = rev(range(HER42PB_xrf_all$depth))) +
  scale_x_continuous(breaks = round(c(0, max(as.numeric(colnames(HER42PB1A$image$image)))*2)),
                     limits = c(0, max(as.numeric(colnames(HER42PB1A$image$image)))*2)) +
  labs(y = "Depth [mm]", x = "[mm]") +
  ggtitle("Optical") +
  ggItraxImage(HER42PB1A, "l") +
  ggItraxImage(HER42PB1B, "r") +
  ggItraxImage(HER42PB1C, "l") +
  ggItraxImage(HER42PB1D, "r") +
  ggItraxImage(HER42PB1E, "l") +
  ggItraxImage(HER42PB1F, "r") +
  ggItraxImage(HER42PB1G, "l") +
  ggItraxImage(HER42PB1H, "r") +
  ggItraxImage(HER42PB1I, "l") +
  ggItraxImage(HER42PB1J, "r") +
  HER42PB_image_theme

print(HER42PB_Opt + coord_fixed(ratio = 1))
ggsave("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Figures/itraxR/Fig6.1_HER42PB_OptPlot.pdf",
       height = c(30), width = c(5), dpi = 600, units = "cm")

# Plot radiograph image only  -------------------------------------------

# ggplot using ggItraxImage of overlapping cores (side-by-side)
theme_set(theme_classic(8))
HER42PB_XRAD <- ggplot() +
  scale_y_continuous(limits = rev(range(HER42PB_xrf_all$depth))) +
  scale_x_continuous(breaks = round(c(0, max(as.numeric(colnames(HER42PB1A$image$image)))*2)),
                     limits = c(0, max(as.numeric(colnames(HER42PB1A$image$image)))*2)) +
  labs(y = "Depth [mm]", x = "[mm]") +
  ggtitle("XRAD-") +
  ggItraxRadio(HER42PB1A, "l") +
  ggItraxRadio(HER42PB1B, "r") +
  ggItraxRadio(HER42PB1C, "l") +
  ggItraxRadio(HER42PB1D, "r") +
  ggItraxRadio(HER42PB1E, "l") +
  ggItraxRadio(HER42PB1F, "r") +
  ggItraxRadio(HER42PB1G, "l") +
  ggItraxRadio(HER42PB1H, "r") +
  ggItraxRadio(HER42PB1I, "l") +
  ggItraxRadio(HER42PB1J, "r") +
  HER42PB_image_theme

print(HER42PB_XRAD + coord_fixed(ratio = 0.1))
ggsave("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Figures/itraxR/Fig6.2_HER42PB_XRADPlot.pdf",
       height = c(30), width = c(5), dpi = 600, units = "cm")

# Plot radiograph on top of optical image -------------------------------

# plot overlappping optical images side-by-side using imagePlot function
theme_set(theme_classic(8))
HER42PB_Opt1 <- ggplot() +
  scale_y_continuous(limits = rev(range(HER42PB_xrf_all$depth))) +
  scale_x_continuous(breaks = round(c(0, max(as.numeric(colnames(HER42PB1A$image$image)))*2)),
                     limits = c(0, max(as.numeric(colnames(HER42PB1A$image$image)))*2)) +
  labs(y = "Depth [mm]", x = "[mm]") +
  ggtitle("Opt/XRAD-") +
  ggItraxImage(HER42PB1A, "l") +
  ggItraxImage(HER42PB1B, "r") +
  ggItraxImage(HER42PB1C, "l") +
  ggItraxImage(HER42PB1D, "r") +
  ggItraxImage(HER42PB1E, "l") +
  ggItraxImage(HER42PB1F, "r") +
  ggItraxImage(HER42PB1G, "l") +
  ggItraxImage(HER42PB1H, "r") +
  ggItraxImage(HER42PB1I, "l") +
  ggItraxImage(HER42PB1J, "r")

# plot radiographs in top using imagePlot function
HER42PB_OptXRAD <- HER42PB_Opt1 +
  ggItraxRadio(HER42PB1A, "l") +
  ggItraxRadio(HER42PB1B, "r") +
  ggItraxRadio(HER42PB1C, "l") +
  ggItraxRadio(HER42PB1D, "r") +
  ggItraxRadio(HER42PB1E, "l") +
  ggItraxRadio(HER42PB1F, "r") +
  ggItraxRadio(HER42PB1G, "l") +
  ggItraxRadio(HER42PB1H, "r") +
  ggItraxRadio(HER42PB1I, "l") +
  ggItraxRadio(HER42PB1J, "r") +
  HER42PB_image_theme

print(HER42PB_OptXRAD + coord_fixed(ratio = 1))
ggsave("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Figures/itraxR/Fig6.3_HER42PB_OptXRAD.pdf",
       height = c(30), width = c(5), dpi = 600, units = "cm")

# -------------------------------------------------------------------------

#  SECTION 7: Combine strat figures & optical/radiograph images

# -------------------------------------------------------------------------

# choose core images to plot
# HER42PB_Opt
# HER42PB_XRAD
# HER42PB_OptXRAD

# choose strat data to plot
# HER42PB_xrfStrat_cps
# HER42PB_xrfStrat_inc_Z
# HER42PB_xrfStrat_pc_cps
# HER42PB_xrfStrat_clr
# HER42PB_xrfStrat_pc_cps_CONISS
# HER42PB_xrfStrat_clr_CONISS

# Optical + clr data as an exmaple below

# plot data and images with patchwork to align y-axis using wrap_plots
library(tidypaleo)
library(patchwork)

HER42_theme <- theme_set(theme_paleo(base_size=8)) + 
  theme(plot.title = element_text(size = 10, face = "bold"))

remove_x <- theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.x = element_blank()
)

remove_y <- theme(
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.title.y = element_blank()
)

# Choose a final plot 
HER42PB_p1 <- list(HER42PB_Opt + remove_x, HER42PB_xrfStrat_cps + remove_y)
HER42PB_p2 <- list(HER42PB_Opt + remove_x, HER42PB_xrfStrat_inc_Z + remove_y)
HER42PB_p3 <- list(HER42PB_Opt + remove_x, HER42PB_xrfStrat_pc_cps + remove_y)
HER42PB_p4 <- list(HER42PB_Opt + remove_x, HER42PB_xrfStrat_clr + remove_y)
HER42PB_p5 <- list(HER42PB_Opt + remove_x, HER42PB_xrfStrat_pc_cps_CONISS + remove_y)
HER42PB_p6 <- list(HER42PB_Opt + remove_x, HER42PB_xrfStrat_clr_CONISS + remove_y)

HER42PB_p7 <- list(HER42PB_XRAD + remove_x, HER42PB_xrfStrat_cps + remove_y)
HER42PB_p8 <- list(HER42PB_XRAD + remove_x, HER42PB_xrfStrat_inc_Z + remove_y)
HER42PB_p9 <- list(HER42PB_XRAD + remove_x, HER42PB_xrfStrat_pc_cps + remove_y)
HER42PB_p10 <- list(HER42PB_XRAD + remove_x, HER42PB_xrfStrat_clr + remove_y)
HER42PB_p11 <- list(HER42PB_XRAD + remove_x, HER42PB_xrfStrat_pc_cps_CONISS + remove_y)
HER42PB_p12 <- list(HER42PB_XRAD + remove_x, HER42PB_xrfStrat_clr_CONISS + remove_y)

HER42PB_p13 <- list(HER42PB_OptXRAD + remove_x, HER42PB_xrfStrat_cps + remove_y)
HER42PB_p14 <- list(HER42PB_OptXRAD + remove_x, HER42PB_xrfStrat_inc_Z + remove_y)
HER42PB_p15 <- list(HER42PB_OptXRAD + remove_x, HER42PB_xrfStrat_pc_cps + remove_y)
HER42PB_p16 <- list(HER42PB_OptXRAD + remove_x, HER42PB_xrfStrat_clr + remove_y)
HER42PB_p17 <- list(HER42PB_OptXRAD + remove_x, HER42PB_xrfStrat_pc_cps_CONISS + remove_y)
HER42PB_p18 <- list(HER42PB_OptXRAD + remove_x, HER42PB_xrfStrat_clr_CONISS + remove_y)

# FINAL itraxR plot = Optical&XRAD + clr data
HER42PB_itraxR_FINAL <- wrap_plots(HER42PB_p16) +
  plot_layout(
    nrow = 1, ncol = 2,
    widths = c(0.5, 8)) + # widths controls the relative width of image vs strat
  plot_annotation(
    title = "HER42PB",
    caption = NULL,
    theme = HER42_theme
  )
HER42PB_itraxR_FINAL
ggsave("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Figures/itraxR/Fig7.1_HER42PB_Opt3&XRAD_clr.pdf",
       height = c(20), width = c(30), dpi = 600, units = "cm")

# Add GEOTEK_MSCL data ----------------------------------------------------

HER42PB_MSCL_comp <- read_csv("MSCL/HER42PB_COMPOSITE_RCT2.4.csv") %>% 
  filter(!if_any(everything(), is.na)) %>%
  mutate(SITE = 'HER42PB', .before = SECTION_ID) %>% 
  select(SITE:MSCL_ID, Strat_depth:RES_SAT, WMAR_SAT)
HER42PB_MSCL_comp

# convert to long
HER42PB_comp_long2 <- select(HER42PB_comp, SITE:Accum_rate, Den1_SAT, MS1_SAT, DCMS1_SAT, WMAR_SAT) %>%
  pivot_longer(c(Den1_SAT, MS1_SAT, DCMS1_SAT, WMAR_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HER42PB_comp_long2

HER42PB_MSCL <- HER42PB_comp_long2 %>% 
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
    units = c("GRD" = "g cm-3", "MS" = "SI x10^-5", "DCMS" = "x10^8 m3 kg-1", "WMAR" = "g cm-2  a-1"))
HER42PB_MSCL

# Add age on secondary age y axis
HER42PB_adm <- age_depth_model(
  HER42PB_comp, 
  depth = Strat_depth, age = SH20_mean_age)

HER42PB_MSCL_age <- HER42PB_MSCL + 
  scale_y_depth_age(HER42PB_adm, age_name = "Age (cal a BP)")
HER42PB_MSCL_age

# MSCl with CONISS

# Uses Den1_SAT, DCMS_SAT only in _com_long3, but plotted with MS1_SAT and WMAR
# CONISS using nested_chclust_coniss() and the rioja package for constrained cluster analysis 
# By default, this calculates the number of significant zones based on a broken stick simulation
# TO DO <- add number of groups and their depth/ages
HER42PB_comp_long3 <- select(HER42PB_comp, SITE:Accum_rate,  Den1_SAT, DCMS1_SAT) %>%
  pivot_longer(c(Den1_SAT, DCMS1_SAT), 
               names_to = "param", values_to = "value") %>% 
  relocate(param, .before = Strat_depth) %>% 
  arrange(param)
HER42PB_comp_long3

HER42PB_MSCL_coniss <- HER42PB_comp_long3 %>%
  nested_data(qualifiers = c(SH20_mean_age, Strat_depth), 
              key = param, value = value, trans = scale) %>%
  nested_chclust_coniss()

HER42PB_MSCL_coniss <- HER42PB_MSCL +
  layer_dendrogram(HER42PB_MSCL_coniss, aes(y = Strat_depth), 
                   param = "CONISS", size = 0.2, lwd = 0.2) +
  layer_zone_boundaries(HER42PB_coniss, aes(y = Strat_depth), 
                        size = 0.3, lwd = 0.3, colour = "red")
HER42PB_MSCL_coniss

# FINAL MSCL plot = OptXRAD + clr data

HER42PB_MSCL_p1 <- list(HER42PB_OptXRAD + remove_x, HER42PB_MSCL + remove_y)

HER42PB_MSCL_FINAL <- wrap_plots(HER42PB_p16) +
  plot_layout(
    nrow = 1, ncol = 2,
    widths = c(0.5, 8)) + # widths controls the relative width of image vs strat
  plot_annotation(
    title = "HER42PB",
    caption = NULL,
    theme = HER42_theme
  )
HER42PB_itraxR_FINAL
ggsave("Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Figures/itraxR/Fig7.1_HER42PB_Opt3&XRAD_clr.pdf",
       height = c(20), width = c(30), dpi = 600, units = "cm")


# Add age-depth model to itraxR datasets ---------------------------------------------------
library(tidyverse)

setwd("/Users/Steve/Dropbox/BAS/Data/R/")

HER42PB_itraxR_depth <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/itraxR/HER42PB_xrf1.csv") %>% 
  select(depth) %>% 
  mutate(depth, ./10)
HER42PB_itraxR_depth

# switch to BACON

write.table(HER42PB_itraxR_depth, "RBacon/Bacon_runs/Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB_depths.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE)


setwd("/Users/Steve/Dropbox/BAS/Data/R/RBacon")
library(rbacon)
Bacon("HER42PB",depths.file=TRUE,d.max=500,thick=10, cc=3, 
      postbomb=4,rotate.axes=TRUE,mem.mean=0.4,rounded=2,
      mem.strength=20,acc.mean=10, yr.max=15000,
      C14.border=rgb(0, 0, 1, 1), C14.col=rgb(0, 0, 1, 0.5), 
      cal.col=rgb(0, 0.5, 0.5, 0.5), cal.border=rgb(0, 0.5, 0.5, 1),
      dark=0.5, height=50, plot.pdf = TRUE, 
      title.location='topright', mgp=c(1.5, 0.7, 0))

setwd("/Users/Steve/Dropbox/BAS/Data/R/")
HER42PB_xrf1 <- read_csv("Papers/Hodgson_Hermite_2022/Data/ITRAX/Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/Output/itraxR/HER42PB_xrf1.csv") 
HER42PB_xrf1_ages <- read_tsv("R/RBacon/Bacon_runs/Papers/Hodgson_Hermite_2022/Data/ITRAX/HER42PB/HER42PB_ages.txt") %>% 
  mutate(depth = depth*10) %>% 
  left_join(HER42PB_xrf1, HER42PB_ages, by = c("depth" = "depth")) %>% 
  rename(min_age = min, max_age = max, median_age = median, mean_age = mean) %>% 
  relocate(min_age:mean_age, .after = depth) %>%
  relocate(label, .after = last_col())
HER42PB_xrf1_ages

HER42PB_xrf1_ages %>%
  group_by_(depth) %>%
  summarise(across(.cols=everything()), mean, na.rm = TRUE)
  
HER42PB_xrf1_ages1






