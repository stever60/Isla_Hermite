
# HER24 itrax.R

# install latest github version of ITRAX R
remotes::install_github("tombishop1/itraxR")

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
library(PeriodicTable)
library(errors)
library(chemometrics)

# Raw data are located here: 
# https://www.dropbox.com/sh/xd5wdn18m0qfnjf/AACTGoGqIQofNG3K6y_gMu_4a?dl=0

#set working directory
setwd("/Users/Steve/Dropbox/BAS/Data/R/Papers/Hodgson_Hermite_2022/Data/ITRAX/HER24")
#check working directory
getwd() 

# SECTION 1: DATA WRANGLING ----------------------------------------------------------

# Renaming columns in  old ITRAX output files to match itrax.R formats before using itraximport  ----------------------------------

# result.txt = original txt output file 
# Results = kcps renamed as cps (code not used here - see original)
#df<- read_tsv("HER24_3C_XRF/result.txt", col_names = TRUE, skip = 2) %>% 
#  rename(cps = kcps, `Fe a*2` = D1)
#select(df1,-X55) %>% 
#  write_tsv("HER24_3A/XRF/Results.tsv")
# Some reanalysis data contains extra column callled revalidity; validity is the original run validity

# Results1 = kcps & cps retained; cps calculated as element and scatter sum by code below 
# import, calculate cps_sum and rename as cps, retain kcps column - save as Results1.tsv

df1_rowsums <- c(11:48, 53, 54) # sum all elements + coh and inc scatter

dfA <- read_tsv("HER24_3A_XRF/result.txt", col_names = TRUE, skip = 2) %>%
  mutate(cps_sum = rowSums(.[df1_rowsums])) %>% 
  rename(cps = cps_sum, `Fe a*2` = D1) %>% 
  relocate(cps, .before = MSE)
select(dfA,-X55) %>% 
  write_tsv("HER24_3A_XRF/Results1.tsv")

dfB <- read_tsv("HER24_3B_XRF/result.txt", col_names = TRUE, skip = 2) %>%
  mutate(cps_sum = rowSums(.[df1_rowsums])) %>% 
  rename(cps = cps_sum, `Fe a*2` = D1) %>% 
  relocate(cps, .before = MSE)
select(dfB,-X55) %>% 
  write_tsv("HER24_3B_XRF/Results1.tsv")

dfC <- read_tsv("HER24_3C_XRF/result.txt", col_names = TRUE, skip = 2) %>%
mutate(cps_sum = rowSums(.[df1_rowsums])) %>% 
  rename(cps = cps_sum, `Fe a*2` = D1) %>% 
  relocate(cps, .before = MSE)
select(dfC,-X55) %>% 
  write_tsv("HER24_3C_XRF/Results1.tsv")

# open Results1.tsv in excel 
# add first two rows back into all new txt files using copy/paste
# rename .tsv files to .txt in Finder
# now ready to import into itrax.R and runs as normal
# choose whether to use Results.txt or Results1.txt
# below uses Results1.txt to match Manchester output from new ITRAX software

# remove from global environment
rm(dfA, dfA1, dfB, dfB1, dfC, df1_rowsums)



# SECTION 2 - ITRAX DATA IMPORT ------------------------------------------------

# Import ITRAX XRF data, image, radiographs for each core using itrax_import ##
# ITRAX folders need to be named as e.g., HER24_1A_XRF, HER24_1A_RAD etc
# Use start depth from e.g., HER24_ITRAX_OVERLAPS.csv file
# Depth is different from position on xrf scanner
# itrax.R retains all data in overlaps and merges data - different to e.g., HER24_ITRAX_COMPOSITE.csv file

# Core 3A

HER24_S1 <- list(metadata   = itrax_meta("HER24_3A_XRF/document.txt"),
                    xrf        = itrax_import("HER24_3A/XRF/Results1.txt", 
                                              depth = 50, 
                                              parameters = "all"),
                    image      = itrax_image(file = "HER24_3A_XRF/optical1.tif",
                                             meta = "HER24_3A_XRF/document.txt"),
                    radiograph = itrax_radiograph(file = "HER24_3A_XRAD//radiograph1.tif",
                                                  meta = "HER24_3A_XRAD/document.txt",
                                                  trim = as.numeric(itrax_meta("HER24_3A/XRAD/document.txt")[6:7,2]))
                 )

# Core 3B

HER24_S2 <- list(metadata   = itrax_meta("HER24_3B_XRF/document.txt"),
                 xrf        = itrax_import("HER24_3B_XRF/Results1.txt", 
                                           depth = 310, # modified depth
                                           parameters = "all"),
                 image      = itrax_image(file = "HER24_3B_XRF/optical1.tif",
                                          meta = "HER24_3B_XRF/document.txt"),
                 radiograph = itrax_radiograph(file = "HER24_3B_XRAD//radiograph1.tif",
                                               meta = "HER24_3B_XRAD/document.txt",
                                               trim = as.numeric(itrax_meta("HER24_3B/XRAD/document.txt")[6:7,2]))
                 )

# Core 3C

HER24_S3 <- list(metadata   = itrax_meta("HER24_3C_XRF/document.txt"),
              xrf        = itrax_import("HER24_3C_XRF/Results1.txt", 
                                        depth = 405, #modified depth
                                        parameters = "all"),
              image      = itrax_image(file = "HER24_3C_XRF/optical1.tif",
                                       meta = "HER24_3C_XRF/document.txt"),
              radiograph = itrax_radiograph(file = "HER24_3C_XRAD//radiograph1.tif",
                                            meta = "HER24_3C_XRAD/document.txt",
                                            trim = as.numeric(itrax_meta("HER24_3C/XRAD/document.txt")[6:7,2]))
              )


# join the xrf data for the sections together ---
HER24_xrf <- itrax_join(list(S1 = HER24_S1$xrf,
                                S2 = HER24_S2$xrf, 
                                S3 = HER24_S3$xrf)
)

HER24_xrf
write.csv(HER24_xrf,"HER24_xrf.csv", row.names = FALSE)

# Convert to long format for plotting in ggplot and tidypaleo
HER24_xrf_long <- select(HER24_xrf,  depth, MSE, cps, all_of(cps_elementsList)) %>%
  pivot_longer(all_of(cps_elementsList), names_to = "param", values_to = "value")
HER24_xrf_long


# Figure 1 - Summary overlaps plot - all sections -------------------------

Fig1.1 <- ggplot(data = na.omit(HER24_xrf), mapping = aes(x = depth, y = Ti)) +
  geom_line(aes(color = label)) + 
  coord_flip() +
  scale_x_reverse() +
  labs(x = "Depth [mm]", color = "Core") +
  theme_classic() +
  theme(legend.position = "none")
Fig1.2 <- ggplot(data = na.omit(HER24_xrf), mapping = aes(x = depth, y = Br)) +
  geom_line(aes(color = label)) + 
  coord_flip() +
  scale_x_reverse() +
  labs(x = "Depth [mm]", color = "Core") +
  theme_classic() +
  theme(legend.position = "none")
Fig1.3 <- ggplot(data = na.omit(HER24_xrf), mapping = aes(x = depth, y = `Mo coh`/`Mo inc`)) +
  geom_line(aes(color = label)) + 
  coord_flip() +
  scale_x_reverse() +
  labs(x = "Depth [mm]", color = "Core") +
  theme_classic() +
  theme(legend.position = "none")
ggarrange(Fig1.1, Fig1.2, Fig1.3, ncol = 3)
ggsave("Figures/Fig1_overlaps.pdf", 
       height = c(15), width = c(15), dpi = 600, units = "cm")


# SECTION 3 - QC and FILTERING DATA ----------------------------------------------------------

# cps filtering - Fe a*2 & adjust cps to between 10,000-20,000
Fig2 <- ggplot(data = HER24_xrf, mapping = aes(x = cps, y = `Fe a*2`)) + 
  geom_point(alpha = 0.1) + 
  theme_bw()
Fig2
ggsave("Figures/Fig2_tolerance_Fe_count.pdf", 
       height = c(10), width = c(10), dpi = 600, units = "cm")

# cps - 2sd is too strict for HER24
cps.mean <- mean(HER24_xrf$cps)
cps.sd <- 3*sd(HER24_xrf$cps)
cps.min.thres <- cps.mean - cps.sd 
cps.max.thres <- cps.mean + cps.sd 

#  OR 

# cps tolerance filter - could use >mean+/-2s cps (based on kcps or cps_sum)
cps.min.thres <- 40000
cps.max.thres <- 90000

HER24_xrf  %>%
  mutate(in_cps_tolerance = ifelse(cps <=cps.min.thres | cps >=cps.max.thres | is.na(cps) == TRUE, FALSE, TRUE)) %>%
  ggplot(mapping = aes(x = depth, y = cps, col = in_cps_tolerance)) + 
  geom_line(aes(group = 1)) +
  scale_x_reverse() +
  geom_hline(yintercept = c(cps.min.thres, cps.max.thres)) +
  geom_rug(sides = "b", data = . %>% filter(in_cps_tolerance == FALSE)) + 
  theme_bw()
ggsave("Figures/Fig3_tolerance_cps.pdf", 
       height = c(10), width = c(10), dpi = 600, units = "cm")

# MSE tolerance filter - 2 used here but could use >mean+2s - which is 1.764766 for ARD record 

# MSE tolerance - 2sd is too strict
MSE.mean <- mean(HER24_xrf$MSE)
MSE.sd <- 3*sd(HER24_xrf$MSE)
MSE.thres <- MSE.mean + MSE.sd 
MSE.thres

#  OR

MSE.thres <- 2 # use this for ARD

HER24_xrf %>%
  mutate(in_mse_tolerance = ifelse(MSE >=MSE.thres, FALSE, TRUE)) %>% 
  ggplot(mapping = aes(x = depth,  y = MSE, col = in_mse_tolerance)) +
  geom_line(aes(group = 1)) +
  scale_x_reverse() +
  geom_hline(yintercept = MSE.thres) +
  geom_rug(sides = "b", data = . %>% filter(in_mse_tolerance == FALSE)) + 
  theme_bw()
ggsave("Figures/Fig4_tolerance_MSE.pdf", 
       height = c(10), width = c(10), dpi = 600, units = "cm")

# Surface slope tolerance filter
slope.min.thres = -0.5
slope.max.thres = 0.5

#  OR - used for ARD
slope1 <-  HER24_xrf$`sample surface` - lag(HER24_xrf$`sample surface`)
s1 <- as_tibble(slope1) %>% 
  filter(!if_any(everything(), is.na))
slope.mean <- mean(s1$value)
slope.sd <- 3*sd(s1$value)
slope.min.thres <- slope.mean - slope.sd 
slope.max.thres <- slope.mean + slope.sd 

HER24_xrf %>%
  mutate(slope = `sample surface` - dplyr::lag(`sample surface`)) %>%
  mutate(in_tolerance = ifelse(slope <=slope.min.thres | slope >=slope.max.thres | is.na(slope) == TRUE, FALSE, TRUE)) %>% 
  ggplot(mapping = aes(x = depth, y = slope, col = in_tolerance)) +
  scale_y_continuous(limits = c(-0.55, 0.55), oob = scales::squish) +
  geom_line(aes(group = 1)) +
  geom_hline(yintercept = c(slope.min.thres, slope.max.thres)) +
  geom_rug(data = . %>% filter(validity == FALSE)) +
  scale_x_reverse() +
  theme_bw()
ggsave("Figures/Fig5_tolerance_sur_slope_.pdf", 
         height = c(10), width = c(10), dpi = 600, units = "cm")
  
# Combining 'validity' flags   
HER24_xrf <- HER24_xrf %>%
  mutate(slope = `sample surface` - dplyr::lag(`sample surface`)) %>%
  mutate(in_slope_tolerance = ifelse(slope <=slope.min.thres | slope >=slope.max.thres | is.na(slope) == TRUE, FALSE, TRUE)) %>%
  select(-slope) %>%
  mutate(in_cps_tolerance = ifelse(cps <=cps.min.thres | cps >=cps.max.thres | is.na(cps) == TRUE, FALSE, TRUE)) %>%
  mutate(in_mse_tolerance = ifelse(MSE <=MSE.thres, TRUE, FALSE)) %>%
  rowwise() %>%
  mutate(qc = !any(c(validity, in_slope_tolerance, in_cps_tolerance, in_mse_tolerance) == FALSE)) %>%
  ungroup() %>%
  select(-c(in_slope_tolerance, in_cps_tolerance, in_mse_tolerance)) # %>% filter(qc == TRUE) #to remove from HER24_xrf rows that dont pass QC
# plot summary
theme_set(theme_bw(8))
Fig6.1 <- ggplot(data = HER24_xrf, aes(y = depth, x = `Ti`, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  geom_point(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  scale_color_discrete(name = "pass QC") +
  scale_y_reverse(name = "Depth [mm]")
Fig6.2 <- ggplot(data = HER24_xrf, aes(y = depth, x = `Fe`, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  geom_point(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  scale_color_discrete(name = "pass QC") +
  scale_y_reverse(name = "Depth [mm]")
Fig6.3 <- ggplot(data = HER24_xrf, aes(y = depth, x = `Br`, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  geom_point(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  scale_color_discrete(name = "pass QC") +
  scale_y_reverse(name = "Depth [mm]")
Fig6.4 <- ggplot(data = HER24_xrf, aes(y = depth, x = `Mo coh`/`Mo inc`, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  geom_point(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  scale_color_discrete(name = "pass QC") +
  scale_y_reverse(name = "Depth [mm]")
ggarrange(Fig6.1, Fig6.2, Fig6.3, Fig6.4, ncol = 4, common.legend = TRUE)
ggsave("Figures/Fig6_tolerance_combined.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

HER24_xrf 

# Combining 'validity' flags - removing data that doesn't pass - replotting as line only
HER24_xrf1 <- HER24_xrf %>%
  mutate(slope = `sample surface` - dplyr::lag(`sample surface`)) %>%
  mutate(in_slope_tolerance = ifelse(slope <=-0.3 | slope >=0.3 | is.na(slope) == TRUE, FALSE, TRUE)) %>%
  select(-slope) %>%
  mutate(in_cps_tolerance = ifelse(cps <=20000 | cps >=60000 | is.na(cps) == TRUE, FALSE, TRUE)) %>%
  mutate(in_mse_tolerance = ifelse(MSE <=2, TRUE, FALSE)) %>%
  rowwise() %>%
  mutate(qc = !any(c(validity, in_slope_tolerance, in_cps_tolerance, in_mse_tolerance) == FALSE)) %>%
  ungroup() %>%
  select(-c(in_slope_tolerance, in_cps_tolerance, in_mse_tolerance)) %>%  
  filter(qc == TRUE) #to remove from HER24_xrf rows that dont pass QC

# plot summary
theme_set(theme_bw(8))
Fig7.1 <- ggplot(data = HER24_xrf, aes(y = depth, x = `Ti`, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  scale_color_discrete(name = "pass QC") +
  scale_y_reverse(name = "Depth [mm]")
Fig7.2 <- ggplot(data = HER24_xrf, aes(y = depth, x = `Fe`, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  scale_color_discrete(name = "pass QC") +
  scale_y_reverse(name = "Depth [mm]")
Fig7.3 <- ggplot(data = HER24_xrf, aes(y = depth, x = `Br`, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  scale_color_discrete(name = "pass QC") +
  scale_y_reverse(name = "Depth [mm]")
Fig7.4 <- ggplot(data = HER24_xrf, aes(y = depth, x = `Mo coh`/`Mo inc`, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  scale_color_discrete(name = "pass QC") +
  scale_y_reverse(name = "Depth [mm]")
ggarrange(Fig7.1, Fig7.2, Fig7.3, Fig7.4, ncol = 4, common.legend = TRUE)
ggsave("Figures/Fig7_tolerance_filtered.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Correcting for dead time & Uncertainties ---------------------------------------------- 

# This doesn't work on cores scanned at Aber - no Dt (Dwell time) column in results.txt file

# ggplot(data = HER24_xrf, aes(x = depth, y = Dt)) +
#  scale_x_reverse() + 
#  scale_y_continuous(sec.axis = sec_axis( trans=~(.+(1-mean(HER24_xrf$Dt, na.rm = TRUE))), name="Correction Factor")) +
#  geom_line() +
#  geom_hline(yintercept = mean(HER24_xrf$Dt, na.rm = TRUE), linetype = "dotted")
 

# Noisy data ------------------------------------------------------------

# see also section: # Calculate as % of normalising factors TS (Total Scatter), cps_sum - in Bertrand et al.R

# Using autocorrelation to detect noisy signals - non-noisy data should be AC, higher, outside 95% limits, showing some order/pattern
library(forecast)
library(ggpubr)

# Individual elements
Fig8 <- ggarrange(
  ggAcf(HER24_xrf$Ca) + ylim(c(NA,1)), ggAcf(HER24_xrf$Ti) + ylim(c(NA,1)), 
  ggAcf(HER24_xrf$Fe) + ylim(c(NA,1)),  ggAcf(HER24_xrf$Sr) + ylim(c(NA,1)),
  ggAcf(HER24_xrf$P) + ylim(c(NA,1)),ggAcf(HER24_xrf$Cu) + ylim(c(NA,1)), 
  ggAcf(HER24_xrf$Zn) + ylim(c(NA,1)), ggAcf(HER24_xrf$S) + ylim(c(NA,1)),
  ggAcf(HER24_xrf$Ni) + ylim(c(NA,1)), ggAcf(HER24_xrf$Cs) + ylim(c(NA,1)),
  ggAcf(HER24_xrf$Ba) + ylim(c(NA,1)), ggAcf(HER24_xrf$Pb) + ylim(c(NA,1)),
  nrow = 4, ncol = 3)
Fig8
ggsave("Figures/Fig8_ACF.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")


# all elements in summary plot
elementsList <- select(HER24_xrf, c(Mg:`Mo coh`)) %>% 
  names()
elementsList

apply(HER24_xrf %>% select(any_of(elementsList)), 2, FUN = function(x){round(Acf(x, plot = F)$acf, 3)}) %>%
  as_tibble(rownames = "lag") %>%
  pivot_longer(!c("lag"), names_to = "elements", values_to = "value") %>%
  mutate(lag = as.numeric(lag),
         elements = factor(elements, levels = filter(., lag == 5) %>% arrange(desc(value)) %>% pull(elements))) %>%
  group_by(elements) %>%
  ggplot(aes(x = lag, y = value, col = elements)) +
  geom_line()
ggsave("Figures/Fig9_ACF_all.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# identify acceptable variables
# use coh and inc as stop points for well measured: 0.5 takes down to Mo inc, 0.23 goes down to Mo coh 
apply(HER24_xrf %>% select(any_of(elementsList)), 2, FUN = function(x){round(Acf(x, plot = F)$acf, 3)}) %>%
  as_tibble(rownames = "lag") %>%
  pivot_longer(!c("lag"), names_to = "elements", values_to = "value") %>%
  mutate(lag = as.numeric(lag),
         elements = factor(elements, levels = filter(., lag == 5) %>% arrange(desc(value)) %>% pull(elements))) %>%
  group_by(elements) %>%
  filter(lag == 5) %>%
  filter(value >= 0.23) %>%
  pull(elements) %>% 
  ordered() -> myElements
myElements

# get acceptable rows and variables and make into long format and then plot
# add P for ARD regression analysis as not in selected
# remove Ar, Ta, W (if selected above) - these are detector generated elements
HER24_xrf %>% 
  filter(qc == TRUE) %>% # pivot long
  select(P, any_of(myElements), depth, label) %>% 
  select(-c(Ar)) %>% 
  tidyr::pivot_longer(!c("depth", "label"), names_to = "elements", values_to = "peakarea") %>% 
  mutate(elements = factor(elements, levels = c(elementsList, "coh/inc"))) %>%
  # plot
  ggplot(aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(color = label)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 2) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "peak area", y = "Depth [mm]") +
  tidypaleo::theme_paleo() +
  theme(legend.position = "none")
ggsave("Figures/Fig10_filtered_elements.pdf", 
       height = c(15), width = c(30), dpi = 600, units = "cm")


## Visualising Raw Data



# SECTION 4 - PLOTTING ----------------------------------------------------

allelements <- 5:40

theme_set(theme_paleo(8))
xrfStrat <- HER24_xrf %>% 
  mutate(`coh/inc` = `Mo coh`/`Mo inc`) %>%
  mutate(`inc/coh` = `Mo inc`/`Mo coh`) %>%
  mutate(TS_sum = `Mo inc` + `Mo coh`) %>% # added by sjro to match Aber normalisation method
  mutate(cps_sum = rowSums(.[allelements])) %>% # added by sjro to match Aber normalisation method
  # select(Fe, Ti, Mn,`coh/inc`, `inc/coh`,TS_sum, cps_sum, depth, label) %>% # a smaller set of elements defined manually to test.
  select(P, S, any_of(myElements), `coh/inc`, `cps_sum`, depth, label) %>%
  select(-c(Ar)) %>% 
  tidyr::pivot_longer(!c("depth", "label"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c(elementsList, "coh/inc", "inc/coh", "TS_sum", "cps_sum")))  
  # note that the levels controls the order

ggplot(xrfStrat, aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(color = label)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 2) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "peak area", y = "Depth [mm]") +
  #tidypaleo::theme_paleo() +
  theme(legend.position = "none")
ggsave("Figures/Fig11_filtered_elements_plot.pdf", 
       height = c(15), width = c(30), dpi = 600, units = "cm")

# Plots with tolerance filtered row removed but no colour
ggplot(xrfStrat, aes(x = peakarea, y = depth)) +
  geom_lineh() + #aes(color = label)
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 2) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "peak area", y = "Depth [mm]") +
  #tidypaleo::theme_paleo() +
  theme(legend.position = "none")
ggsave("Figures/Fig12_multi_tolerance_filtered_elements.pdf", 
       height = c(15), width = c(30), dpi = 600, units = "cm")

xrfStrat


# SECTION 5 - TRANSFORMING DATA -------------------------------------------

# validity filtered
HER24_xrfNorm <- HER24_xrf %>% # n = 2005
  mutate(`coh_inc` = `Mo coh`/`Mo inc`) %>%
  mutate(`inc_coh` = `Mo inc`/`Mo coh`) %>%
  mutate(TS_sum = `Mo inc` + `Mo coh`) %>% # added by sjro to match Aber normalisation method
  mutate(cps_sum = rowSums(.[allelements])) %>% # added by sjro to match Aber normalisation method

  # transform
  mutate(across(any_of(elementsList)) /`Mo inc`) %>%
  mutate_if(is.numeric, list(~na_if(., Inf))) %>% # convert all Inf to NA
  
  # identify acceptable observations - validity and/or qc
  # comment in/out to choose either/or, or both
  filter(validity == TRUE) %>% # n = 1996
  #filter(qc == TRUE) %>% # n = 1707
  
  # identify acceptable variables
  select(P, S, any_of(myElements), `coh_inc`, `cps_sum`, depth, label) %>%
  select(-c(Ar))
  
# pivot
HER24_xrfNorm_long <-  tidyr::pivot_longer(HER24_xrfNorm, !c("depth", "label"), names_to = "elements", values_to = "peakarea") %>% 
  mutate(elements = factor(elements, levels = c(elementsList, "coh/inc", "inc/coh", "TS_sum", "cps_sum")))
  
# plot
ggplot(HER24_xrfNorm_long, aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(color = label)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 2) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "peak area / Mo. inc.", y = "Depth [mm]") +
  tidypaleo::theme_paleo() +
  theme(legend.position = "none")
ggsave("Figures/Fig13_Mo inc_normalised.pdf", 
         height = c(15), width = c(30), dpi = 600, units = "cm")

# Validity and qc filtered
HER24_xrfNorm_qc <- HER24_xrf %>% # n = 2005
  mutate(`coh_inc` = `Mo coh`/`Mo inc`) %>%
  mutate(`inc_coh` = `Mo inc`/`Mo coh`) %>%
  mutate(TS_sum = `Mo inc` + `Mo coh`) %>% # added by sjro to match Aber normalisation method
  mutate(cps_sum = rowSums(.[allelements])) %>% # added by sjro to match Aber normalisation method
  
  # transform
  mutate(across(any_of(elementsList)) /`Mo inc`) %>%
  mutate_if(is.numeric, list(~na_if(., Inf))) %>% # convert all Inf to NA
  
  # identify acceptable observations - validity and/or qc
  # comment in/out to choose either/or, or both
  filter(validity == TRUE) %>% # n = 1996
  filter(qc == TRUE) %>% # n = 1707
  
  # identify acceptable variables
  select(P, S, any_of(myElements), `coh_inc`, `cps_sum`, depth, label) %>%
  select(-c(Ar))

# pivot
HER24_xrfNorm_qc_long <-  tidyr::pivot_longer(HER24_xrfNorm_qc, !c("depth", "label"), names_to = "elements", values_to = "peakarea") %>% 
  mutate(elements = factor(elements, levels = c(elementsList, "coh/inc", "inc/coh", "TS_sum", "cps_sum")))

# plot
ggplot(HER24_xrfNorm_long_qc, aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(color = label)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 2) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "peak area / Mo. inc.", y = "Depth [mm]") +
  tidypaleo::theme_paleo() +
  theme(legend.position = "none")
ggsave("Figures/Fig13A_Mo inc_normalised_qc.pdf", 
       height = c(15), width = c(30), dpi = 600, units = "cm")


# SMOOTHING  ---------------------------------------------------------------

HER24_xrfSmooth <- HER24_xrf %>%
  # uses a 10 point running mean (2 cm for this data); 5 before, 5 after - 1 cm i.e., 5 point RM 2.5 before/after doesnt work
  mutate(across(any_of(elementsList), 
                function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
  )
  ) 

ggplot(HER24_xrfSmooth, mapping = aes(x = depth, y = Ca)) + 
  geom_line(data = HER24_xrf, col = "grey80") + 
  geom_line() + 
  scale_x_reverse() +
  theme_paleo()
ggsave("Figures/Fig14_Ca_smoothed.pdf", 
       height = c(15), width = c(30), dpi = 600, units = "cm")

# Smoothed stratigraphic diagram -------------------------------------------------
# smoothed data has to be labelled and combined with the original data so it can be faceted.
# make the xrf plot with running means

# make new dataset HER24_xrf1 with coh/inc and cps_sum included
HER24_xrf1 <- HER24_xrf %>%
  mutate(coh_inc = `Mo coh`/`Mo inc`) %>%
  mutate(inc_coh = `Mo inc`/`Mo coh`) %>%
  mutate(TS_sum = `Mo inc` + `Mo coh`) %>% # added by sjro to match Aber normalisation method
  mutate(cps_sum = rowSums(.[allelements]))# added by sjro to match Aber normalisation method
HER24_xrf1

# make new element list
elementsList1 <- select(HER24_xrf1, c(Mg:`Mo coh`, coh_inc, inc_coh, cps_sum)) %>% names()
elementsList1

# Smoothed cps plot - final join, smooth and plot with elements of most interest
full_join(y = HER24_xrf1 %>%
            as_tibble() %>%
            # uses a 10 point running mean (2 cm for this data); 5 before, 5 after
            mutate(across(any_of(c(elementsList1)), 
                          function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
            )
            ) %>%
            mutate(type = "mean"), 
          x = HER24_xrf1 %>% 
            as_tibble() %>% 
            mutate(type = "raw")
) %>% 
  filter(validity == TRUE) %>%
  #filter(qc == TRUE) %>%
  select(P, Ca, Ti, Cu, Zn, Sr, coh_inc, MSE, depth, label, type) %>%
  tidyr::pivot_longer(!c("depth", "label", "type"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c("MSE", elementsList1))) %>%
  mutate(label = as_factor(label),
         type = as_factor(type)
  ) %>%

  glimpse() %>%
  
  ggplot(aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(group = type, colour = label, alpha = type)) +
  scale_alpha_manual(values = c(0.1, 1)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 2) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "peak area", y = "Depth [mm]") +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Smoothed (10 pt, 2 cm RM) cps; validity filtered")
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank()) 
ggsave("Figures/Fig15_smoothed_cps.pdf", 
       height = c(15), width = c(30), dpi = 600, units = "cm")


# Smoothed Mo inc normalised plot - validity filtered - elements of most interest
full_join(y = HER24_xrfNorm %>%
            as_tibble() %>%
            # uses a 10 point running mean (2 cm for this data); 5 before, 5 after
            mutate(across(any_of(c(elementsList1)), 
                          function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
            )
            ) %>%
            mutate(type = "mean"), 
          x = HER24_xrfNorm %>% 
            as_tibble() %>% 
            mutate(type = "raw")
) %>% 
  #filter(validity == TRUE) %>% not needed because HER24_xrfNorm has already been filtered
  #filter(qc == TRUE) %>%
  select(P, Ca, Ti, Cu, Zn, Sr, coh_inc, depth, label, type) %>%
  tidyr::pivot_longer(!c("depth", "label", "type"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c("MSE", elementsList1))) %>%
  mutate(label = as_factor(label),
         type = as_factor(type)
  ) %>%
  
  glimpse() %>%
  
  ggplot(aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(group = type, colour = label, alpha = type)) +
  scale_alpha_manual(values = c(0.1, 1)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 2) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "peak area / Mo inc.", y = "Depth [mm]") +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Smoothed (10 pt, 2 cm RM) inc. normalised; validity filtered")
#axis.text.x = element_blank(),
#axis.ticks.x = element_blank())
ggsave("Figures/Fig16_smoothed_incNorm.pdf", 
       height = c(15), width = c(30), dpi = 600, units = "cm")

# Smoothed Mo inc normalised plot - validity and qc filtered - elements of most interest
full_join(y = HER24_xrfNorm_qc %>%
            as_tibble() %>%
            # uses a 10 point running mean (2 cm for this data); 5 before, 5 after
            mutate(across(any_of(c(elementsList1)), 
                          function(x){unlist(slider::slide(x, mean, .before = 5, .after = 5))}
            )
            ) %>%
            mutate(type = "mean"), 
          x = HER24_xrfNorm_qc %>% 
            as_tibble() %>% 
            mutate(type = "raw")
) %>% 
  #filter(validity == TRUE) %>% not needed because HER24_xrfNorm has already been filtered
  #filter(qc == TRUE) %>%
  select(P, Ca, Ti, Cu, Zn, Sr, coh_inc, depth, label, type) %>%
  tidyr::pivot_longer(!c("depth", "label", "type"), names_to = "elements", values_to = "peakarea") %>% 
  tidyr::drop_na() %>%
  mutate(elements = factor(elements, levels = c("MSE", elementsList1))) %>%
  mutate(label = as_factor(label),
         type = as_factor(type)
  ) %>%
  
  glimpse() %>%
  
  ggplot(aes(x = peakarea, y = depth)) +
  tidypaleo::geom_lineh(aes(group = type, colour = label, alpha = type)) +
  scale_alpha_manual(values = c(0.1, 1)) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 2) +
  facet_geochem_gridh(vars(elements)) +
  labs(x = "peak area / Mo inc.", y = "Depth [mm]") +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Smoothed (10 pt, 2 cm RM) inc. normalised; validity and qc filtered")
#axis.text.x = element_blank(),
#axis.ticks.x = element_blank())
ggsave("Figures/Fig17_smoothed_incNorm_qc.pdf", 
       height = c(15), width = c(30), dpi = 600, units = "cm")

# Smoothed %cps_sum and Ti-log normalised plots - validity filtered - elements of most interest ***** TO DO **** 
# get code from Bertrand et al.R section 
# log normalised

# Write to file  --------------------------------------------------------
write.csv(HER24_xrf,"Output/HER24_xrf.csv", row.names = FALSE)
write.csv(HER24_xrf1,"Output/HER24_xrf1.csv", row.names = FALSE)
write.csv(HER24_xrfNorm,"Output/HER24_xrfNorm.csv", row.names = FALSE)
write.csv(HER24_xrfNorm_qc,"Output/HER24_xrfNorm_qc.csv", row.names = FALSE)


# SECTION 6 - MULTIVARIATE METHODS -------------------------------------------



# SECTION 7 - CALIBRATING DATA -------------------------------------------

