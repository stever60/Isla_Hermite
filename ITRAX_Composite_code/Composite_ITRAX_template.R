# HER49L ITRAX using a previously made composite file = MSCL composite

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

# SECTION 2: IMPORT DATA

# -------------------------------------------------------------------------

#set working directory
setwd("/Users/Steve/Dropbox/BAS/Data/R/Papers/Hodgson_Hermite_2022/Data/ITRAX/")
#check working directory
getwd() 

# Assign elements to/from PeriodicTable package to 'elementsList' (then unload package)
data(periodicTable)
elementsList <- periodicTable$symb

# list of all possible elements only
data(periodicTable)
allelements <- c(symb(1:117))
rm(periodicTable)

# Import an existing composite record and change column names to match itrax.R names/formatting
HER49L_xrf_comp <- read_csv("HER49L/HER49L_ITRAX_COMPOSITE_cps.csv") %>%
  relocate(Mo_inc, .after = U) %>%
  relocate(Mo_coh, .after = Mo_inc) %>% 
  mutate(cps = rowSums(across(Mg:Mo_coh))) %>% #mutate(cps = rowSums(.[df1_rowsums])) %>%
  select(-c(filename, position_mm, 
                          Section_depth_cm, Field_depth_cm, SH20_min_age_95CI: SH20_median_age, 
                          E_gain:F_offset, S1:S3)) %>% 
  rename(`Fe a*2` = D1, SH20_age = SH20_mean_age, 
         depth = Strat_depth_cm, label = Section) %>% 
  relocate(cps, .before = MSE) %>% 
  filter(validity=='1')

HER49L_xrf_comp

write.csv(HER49L_xrf_comp,"HER49L/Output/Composite/HER49L_xrf_comp.csv", row.names = FALSE)

# SECTION 3: QC and FILTERING DATA --------------------------------------------------------------

# Set up data filtering criteria using quality control measures - following itrax.R principles
# can skip this section if happy to accept validity = TRUE from ITRAX measurements  

#  MSE and cps filtering --------------------------------------------------

# cps filtering - Fe a*2
Fig3.1 <- ggplot(data = HER49L_xrf_comp, mapping = aes(x = cps, y = `Fe a*2`)) + 
  geom_point(alpha = 0.1) + 
  theme_bw()
print(Fig3.1)

# cps - 2 std dev is too strict for HER49L and other peat cores with a combination of low and high count matrices
cps.mean <- mean(HER49L_xrf_comp$cps)
cps.sd <- 4*sd(HER49L_xrf_comp$cps)
cps.min.thres <- cps.mean - cps.sd 
cps.max.thres <- cps.mean + cps.sd 

#  OR use this 

# cps tolerance filter 
# cps.min.thres <- 40000
# cps.max.thres <- 90000

Fig3.2 <-HER49L_xrf_comp  %>%
  mutate(in_cps_tolerance = ifelse(cps <=cps.min.thres | cps >=cps.max.thres | is.na(cps) == TRUE, FALSE, TRUE)) %>%
  ggplot(mapping = aes(x = depth, y = cps, col = in_cps_tolerance)) + 
  geom_line(aes(group = 1)) +
  scale_x_reverse() +
  geom_hline(yintercept = c(cps.min.thres, cps.max.thres)) +
  geom_rug(sides = "b", data = . %>% filter(in_cps_tolerance == FALSE)) + 
  theme_bw()
print(Fig3.2)

# MSE tolerance filter set at 6xSD - need to be careful as MSE can indicate different lithologies 

# MSE tolerance - 2 std dev is too strict when MSE values are very similar but all below 2
MSE.mean <- mean(HER49L_xrf_comp$MSE)
MSE.sd <- 4*sd(HER49L_xrf_comp$MSE)
MSE.thres <- MSE.mean + MSE.sd 
MSE.thres

#  OR

#MSE.thres <- 2 # use this as an established general threshold

Fig3.3 <- HER49L_xrf_comp %>%
  mutate(in_mse_tolerance = ifelse(MSE >=MSE.thres, FALSE, TRUE)) %>% 
  ggplot(mapping = aes(x = depth,  y = MSE, col = in_mse_tolerance)) +
  geom_line(aes(group = 1)) +
  scale_x_reverse() +
  geom_hline(yintercept = MSE.thres) +
  geom_rug(sides = "b", data = . %>% filter(in_mse_tolerance == FALSE)) + 
  theme_bw()
print(Fig3.3)

# Surface slope tolerance filter ------------------------------------------

# Either use this with wide margins for peat cores eg 6 std dev equivalent to +/-0.5 
#slope.min.thres = -0.5
#slope.max.thres = 0.5

#  OR based on mean and SD thresholds 

slope1 <-  HER49L_xrf_comp$surface - lag(HER49L_xrf_comp$surface)
s1 <- as_tibble(slope1) %>% 
  filter(!if_any(everything(), is.na))
slope.mean <- mean(s1$value)
slope.sd <- 4*sd(s1$value)
slope.min.thres <- slope.mean - slope.sd 
slope.max.thres <- slope.mean + slope.sd 

Fig3.4 <- HER49L_xrf_comp %>%
  mutate(slope = surface - dplyr::lag(surface)) %>%
  mutate(in_slope_tolerance = ifelse(slope <=slope.min.thres | slope >=slope.max.thres | is.na(slope) == TRUE, FALSE, TRUE)) %>% 
  ggplot(mapping = aes(x = depth, y = slope, col = in_slope_tolerance)) +
  scale_y_continuous(limits = c(-0.55, 0.55), oob = scales::squish) +
  geom_line(aes(group = 1)) +
  geom_hline(yintercept = c(slope.min.thres, slope.max.thres)) +
  geom_rug(data = . %>% filter(validity == FALSE)) +
  scale_x_reverse() +
  theme_bw()
print(Fig3.4)

# show Figure 3.1 - 3.4 together:
ggarrange(Fig3.1, Fig3.2, Fig3.3, Fig3.4, nrow = 2, ncol = 2, labels = c('a', 'b', 'c', 'd'))
ggsave("HER49L/Figures/Composite/Fig3.1_comp_tolerance_combined.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Combining all 'validity' flags  ---------------------------------------------

HER49L_xrf_comp_qc <- HER49L_xrf_comp %>%
  mutate(slope = surface - dplyr::lag(surface)) %>%
  mutate(in_slope_tolerance = ifelse(slope <=slope.min.thres | slope >=slope.max.thres | is.na(slope) == TRUE, FALSE, TRUE)) %>%
  select(-slope) %>%
  mutate(in_cps_tolerance = ifelse(cps <=cps.min.thres | cps >=cps.max.thres | is.na(cps) == TRUE, FALSE, TRUE)) %>%
  mutate(in_mse_tolerance = ifelse(MSE <=MSE.thres, TRUE, FALSE)) %>%
  rowwise() %>%
  mutate(qc = !any(c(validity, in_slope_tolerance, in_cps_tolerance, in_mse_tolerance) == FALSE)) %>%
  ungroup() %>%
  select(-c(in_slope_tolerance, in_cps_tolerance, in_mse_tolerance)) %>% 
  filter(qc == TRUE) #to remove from HER49L_xrf rows that dont pass QC

# plot summary
theme_set(theme_bw(12))
Fig3.5 <- ggplot(data = HER49L_xrf_comp_qc, aes(y = depth, x = Ti, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.width = unit(0.5, 'cm'),
        legend.box.spacing = unit(0.1, 'cm')) +
  scale_color_discrete(name = "Pass QC") +
  scale_y_reverse(name = "Depth [cm]")

Fig3.6 <- ggplot(data = HER49L_xrf_comp_qc, aes(y = depth, x = Fe, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.width = unit(0.5, 'cm'),
        legend.box.spacing = unit(0.1, 'cm')) +
  scale_color_discrete(name = "Pass QC") +
  scale_y_reverse(name = "Depth [cm]")

Fig3.7 <- ggplot(data = HER49L_xrf_comp_qc, aes(y = depth, x = Br, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.width = unit(0.5, 'cm'),
        legend.box.spacing = unit(0.1, 'cm')) +
  scale_color_discrete(name = "Pass QC") +
  scale_y_reverse(name = "Depth [cm]")

Fig3.8 <- ggplot(data = HER49L_xrf_comp_qc, aes(y = depth, x = Mo_inc/Mo_coh, col = qc)) + 
  geom_lineh(aes(group = 1)) +
  geom_rug(sides = "l", data = . %>% filter(qc == FALSE)) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.width = unit(0.5, 'cm'),
        legend.box.spacing = unit(0.1, 'cm')) +
  scale_color_discrete(name = "Pass QC") +
  scale_y_reverse(name = "Depth [cm]")

ggarrange(Fig3.5, Fig3.6, Fig3.7, Fig3.8, ncol = 2, nrow = 2, common.legend = TRUE)
ggsave("HER49L/Figures/Composite/Fig3.2_tolerance_comp_combined_examples.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Remove machine elements and known matrix effect elements from dataframe ----------------
machine_elements <- select(HER49L_xrf_comp_qc, c(Ar, Ta, W)) %>% 
  names()
machine_elements

# Create list of elements - can remove Zr (matrix effect in peat/organic seds) at this stage 
HER49L_elements <- select(HER49L_xrf_comp_qc, any_of(elementsList), Mo_inc, Mo_coh, -c(all_of(machine_elements), `Fe a*2`)) %>% # , Zr
  names()
HER49L_elements

# REE if in dataframe
#REE <- select(HER49L_xrf_comp_qc, c(La:Ho)) %>% 
#  names()
#REE

# -------------------------------------------------------------------------

# Choose dataset to take forward to element filtering section

# -------------------------------------------------------------------------

# Simple QC filtering based on validity = TRUE only
HER49L_xrf <- HER49L_xrf_comp %>%
  select(Site:MSE, any_of(HER49L_elements)) %>%
  filter(validity == 1)
HER49L_xrf

#  OR use below - use this most of the time

# Complex QC filtering based on multiple parameters above
HER49L_xrf <- HER49L_xrf_comp_qc 
HER49L_xrf

# Summary stats for xrf
HER49L_xrf_stats <- HER49L_xrf %>%
  select(kcps:MSE, any_of(HER49L_elements)) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()
tail(HER49L_xrf_stats)

# Save QC xrf dataset and stats to file
write.csv(HER49L_xrf,"HER49L/Output/Composite/HER49L_xrf.csv", row.names = FALSE)
write.csv(HER49L_xrf_stats,"HER49L/Output/Composite/HER49L_xrf_stats.csv", row.names = FALSE)

# -------------------------------------------------------------------------

# SECTION 3: Measurable element filtering

# -------------------------------------------------------------------------

#  Autocorrelation based filtering of elements ----------------------------------

# Use autocorrelation function (acf) and plots to explore noise in a time-series
library(forecast)
library(ggrepel)
library(directlabels)

# Adjust for any element of interest by changing $
# split into two groups below to visualise the most common elements measured by ITRAX
theme_set(theme_bw(8))
Fig3.10 <- ggarrange(
  ggAcf(HER49L_xrf$Al) + ylim(c(NA,1)), ggAcf(HER49L_xrf$Si) + ylim(c(NA,1)), 
  ggAcf(HER49L_xrf$P) + ylim(c(NA,1)), ggAcf(HER49L_xrf$S) + ylim(c(NA,1)),
  ggAcf(HER49L_xrf$Cl) + ylim(c(NA,1)), ggAcf(HER49L_xrf$K) + ylim(c(NA,1)), 
  ggAcf(HER49L_xrf$Ca) + ylim(c(NA,1)), ggAcf(HER49L_xrf$Ti) + ylim(c(NA,1)),
  ggAcf(HER49L_xrf$V) + ylim(c(NA,1)), ggAcf(HER49L_xrf$Cr) + ylim(c(NA,1)),
  ggAcf(HER49L_xrf$Mn) + ylim(c(NA,1)), ggAcf(HER49L_xrf$Fe) + ylim(c(NA,1)),
  nrow = 4, ncol = 3)
print(Fig3.10)
ggsave("HER49L/Figures/Composite/Fig3.3_ACF_comp_pt1.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

Fig3.11 <- ggarrange(
  ggAcf(HER49L_xrf$Co) + ylim(c(NA,1)), ggAcf(HER49L_xrf$Ni) + ylim(c(NA,1)),
  ggAcf(HER49L_xrf$Cu) + ylim(c(NA,1)), ggAcf(HER49L_xrf$Zn) + ylim(c(NA,1)),
  ggAcf(HER49L_xrf$Se) + ylim(c(NA,1)), ggAcf(HER49L_xrf$Br) + ylim(c(NA,1)),
  ggAcf(HER49L_xrf$Rb) + ylim(c(NA,1)), ggAcf(HER49L_xrf$Sr) + ylim(c(NA,1)),
  ggAcf(HER49L_xrf$Zr) + ylim(c(NA,1)), ggAcf(HER49L_xrf$Ba) + ylim(c(NA,1)),
  ggAcf(HER49L_xrf$Mo_inc) + ylim(c(NA,1)), ggAcf(HER49L_xrf$Mo_coh) + ylim(c(NA,1)),
  nrow = 4, ncol = 3)
print(Fig3.11)
ggsave("HER49L/Figures/Composite/Fig3.4_ACF_comp_pt2.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Filter elements based on acf lag thresholds ------------------------------

# set lag threshold to an arbitary value of 5  for whole dataset
# for 1 mm dataset, this checks autocorrelation +/-5 mm around each measurement 
# can be performed on different lithologies later on
# 0.2 or 0.1 is usual cut off point

# define filter and lag thresholds
acf_thres_min <- 0.1
acf_thres_max <- 0.5
lag_thres <- 20

Fig3.12a <- apply(HER49L_xrf %>% select(any_of(elementsList), Mo_inc, Mo_coh), 
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
ggsave("HER49L/Figures/Composite/Fig3.5a_ACF_comp_all_elements.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# remove machine elements and apply acf based filtering to - leaving acf elements
apply(HER49L_xrf %>% select(any_of(elementsList), Mo_inc, Mo_coh), 2, FUN = function(x){round(Acf(x, plot = F)$acf, 3)}) %>%
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

acfElementsList <- select(HER49L_xrf, any_of(acfElements)) %>% 
  names()
acfElementsList

# Replot with acf filtered elements only
Fig3.12b <- apply(HER49L_xrf %>% select(any_of(acfElements)), 
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
ggsave("HER49L/Figures/Composite/Fig3.5b_ACF_comp_all_elements.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Summary acf element plot vs depth --------------------------------------------------
library(tidypaleo)
theme_set(theme_bw(base_size=8))

Fig3.13 <- HER49L_xrf %>% 
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
  ggtitle("HER49L-ITRAX: cps (ACF-filtered elements)")
print(Fig3.13)
ggsave("HER49L/Figures/Composite/Fig3.6_ACF_key_elements.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# nested ggarrange to produce split level summary plot with different number of plots pe row
ggarrange(Fig3.13, # First row
          ggarrange(Fig3.12a, Fig3.12b, ncol = 2, labels = c("B", "C")), # Second row with two plots
          nrow = 2, 
          labels = "A", common.legend = TRUE)
ggsave("HER49L/Figures/Composite/Fig3.7_ACF_elements_depth.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# ACF elements - further filtering based on mean and max cps and %cps sum values

# check acf element list
acfElementsList

# cps filtering based on mean >50 cps, max>100 cps

# mean cps > 50 function
HER49L_mean50 <- function(x){
  is.numeric(x) && (mean(x, na.rm = TRUE) > 50)
}

# mean cps >50 dataframe 
HER49Lmean50 <- select(HER49L_xrf, Site:MSE,
                        any_of(acfElementsList) &
                          where(HER49L_mean50))
# mean cps >50 element list - remove Zr for peat as it's a matrix effect
HER49Lmean50List <- select(HER49Lmean50, -c(Site:MSE, Zr)) %>%
  names() %>% 
  print()

# max cps >100
HER49L_max100 <- function(x){
  is.numeric(x) && (max(x, na.rm = TRUE) > 100)
}
HER49Lmax100 <- select(HER49L_xrf,Site:MSE,
                       any_of(acfElementsList) &
                         where(HER49L_max100))
HER49Lmax100List <- select(HER49Lmax100, -c(Site:MSE, Zr)) %>%
  names() %>% 
  print()

# mean cps >50 and max cps >100
HER49Lmean50_max100 <- select(HER49Lmean50,Site:MSE,
                              any_of(acfElementsList) &
                                where(HER49L_max100))
HER49Lmean50_max100List <- select(HER49Lmean50_max100, -c(Site:MSE, Zr)) %>%
  names() %>% 
  print()

# Select cps filtered dataframe to take forward
HER49L_xrf_filter <- HER49Lmean50_max100
HER49L_xrf_filter

HER49L_filterList <- HER49Lmean50_max100List
HER49L_filterList

# Generate cps summary stats for filtered dataframe
HER49L_xrf_filter_stats <- HER49L_xrf_filter %>%
  select(kcps:cps, any_of(HER49L_filterList)) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()

# Save qc element filtered dataset & stats to file
write.csv(HER49L_xrf_filter,"HER49L/Output/Composite/HER49L_xrf_filter.csv", row.names = FALSE)
write.csv(HER49L_xrf_filter_stats,"HER49L/Output/Composite/HER49L_xrf_filter_stats.csv", row.names = FALSE)

# Create a new dataframe called xrf1 of qc, acf and filtered elements and add  additional useful parameters ----------------

HER49L_xrf[HER49L_xrf == 0] <- NA #set 0 to NA log value calculations

HER49L_xrf1 <- HER49L_xrf %>%
  select(Site:MSE, label, K, Ti, any_of(HER49L_filterList), Mo_inc, Mo_coh) %>% 
  mutate(TS_sum = Mo_inc + Mo_coh) %>%
  mutate(inc_coh = Mo_inc / Mo_coh) %>%
  mutate(coh_inc = Mo_coh / Mo_inc)%>%
  mutate(LnK_Ti = log(K / Ti)) %>% 
  mutate(LnK_Ca = log(K / Ca)) %>% 
  mutate(LnTi_Ca = log(Ti / Ca)) %>% 
  relocate(label, .after = depth) %>%
  relocate(K, .after = S) %>%
  relocate(Ti, .after = Ca) %>%
  relocate(TS_sum, .after = Mo_coh) %>%
  relocate(inc_coh, .after = TS_sum) %>%
  relocate(coh_inc, .after = inc_coh) %>% 
  relocate(LnK_Ti, .after = coh_inc) %>% 
  relocate(LnK_Ca, .after = LnK_Ti) %>% 
  relocate(LnTi_Ca, .after = LnK_Ca) %>% 
  relocate(label, .after = Site)
HER49L_xrf1

#set NA back to 0 for _xrf dataframe
HER49L_xrf <- HER49L_xrf %>% 
  replace(is.na(.), 0) 
HER49L_xrf

# new element and ratio lists
HER49L_elements1 <- select(HER49L_xrf1, any_of(elementsList), Mo_inc, Mo_coh) %>% 
  names()
HER49L_elements1

HER49L_ratioList1 <- select(HER49L_xrf1, inc_coh, coh_inc, LnK_Ti, LnK_Ca, LnTi_Ca) %>% 
  names()
HER49L_ratioList1

# Summary stats for xrf1
HER49L_xrf1_stats <- HER49L_xrf1 %>%
  select(any_of(HER49L_elements1), any_of(HER49L_ratioList1)) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()
tail(HER49L_xrf1_stats)

write.csv(HER49L_xrf1,"HER49L/Output/Composite/HER49L_xrf1.csv", row.names = FALSE)
write.csv(HER49L_xrf1_stats,"HER49L/Output/Composite/HER49L_xrf1_stats.csv", row.names = FALSE)

# -------------------------------------------------------------------------

# SECTION 4.1: Data Transformation using filtered cps dataframe 

# -------------------------------------------------------------------------

# cps as % of cps_sum

# original cps qc dataframe
HER49L_xrf_pc_cps0 <- HER49L_xrf_comp_qc %>%
  select(Site:MSE, Mg:Mo_coh, qc) %>%
  mutate(across(Mg:Mo_coh) /`cps`) %>%
  mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
  replace(is.na(.), 0) %>%
  mutate(across(Mg:Mo_coh) *100) %>% 
  print()
# check sum = 100%
HER49L_xrf_pc_cps1 <- HER49L_xrf_pc_cps0 %>% 
  mutate(cps_sum_all = rowSums(across(Mg:Mo_coh))) 
HER49L_xrf_pc_cps2 <- HER49L_xrf_pc_cps1 %>% 
  mutate(cps_sum_filtered = rowSums(across(all_of(HER49L_elements1)))) 
head(HER49L_xrf_pc_cps2$cps_sum_all)
head(HER49L_xrf_pc_cps2$cps_sum_filtered)

# %cps_sum summary stats table for all ITRAX elements
HER49L_xrf_pc_cps_stats <- HER49L_xrf_pc_cps1 %>%
  select(any_of(HER49L_elements), cps_sum_all) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()
tail(HER49L_xrf_pc_cps_stats)
write.csv(HER49L_xrf_pc_cps_stats,"HER49L/Output/Composite/HER49L_xrf_pc_cps_stats.csv", row.names = FALSE)

# Filtered elements as % of cps_sum
HER49L_xrf1_pc_cps_filtered <- select(HER49L_xrf_pc_cps2, 
                               Site:MSE,all_of(HER49L_elements1), cps_sum_all, cps_sum_filtered)

HER49L_xrf1_pc_cps <- select(HER49L_xrf1, inc_coh:LnTi_Ca) %>% 
  bind_cols(HER49L_xrf1_pc_cps_filtered) %>%
  relocate(inc_coh:LnTi_Ca, .after = Mo_coh)
HER49L_xrf1_pc_cps

# %cps_sum summary stats table for filtered ITRAX elements
HER49L_xrf1_pc_cps_stats <- HER49L_xrf1_pc_cps %>%
       select(any_of(HER49L_elements1), cps_sum_all, cps_sum_filtered) %>% 
       psych::describe(quant=c(.25,.75)) %>%
       as_tibble(rownames="rowname")  %>%
       print()
write.csv(HER49L_xrf1_pc_cps_stats,"HER49L/Output/Composite/HER49L_xrf1_pc_cps_stats.csv", row.names = FALSE)

# -------------------------------------------------------------------------

# Normalised datasets

# Normalised by inc scatter
HER49L_xrf1_inc <- HER49L_xrf1 %>%
  mutate(across(any_of(HER49L_elements1)) /`Mo_inc`) %>%
  mutate_if(is.numeric, list(~na_if(., Inf)))
HER49L_xrf1_inc

# Normalised by Ti
HER49L_xrf1_Ti <- HER49L_xrf1 %>% 
  mutate(across(any_of(HER49L_elements1), .fns = ~./ Ti))

# Natural log of Ti normalised dataframe with -Inf replaced by NA
HER49L_xrf1_Ln_Ti <- HER49L_xrf1_Ti %>% 
  mutate(across(any_of(HER49L_elements1),log)) 
is.na(HER49L_xrf1_Ti)<-sapply(HER49L_xrf1_Ti, is.infinite)
# Replace NA with 0 - if needed
# HER49L_Ln_Ti_norm[is.na(HER49L_Ln_Ti_norm)]<-0
# can replace Ti normalized with other normalisation parameters above

#Check = 1
head(HER49L_xrf1_Ti$Ti)
#Check = 0
head(HER49L_xrf1_Ln_Ti$Ti)

# -------------------------------------------------------------------------

# Z-scores

# Standardised and centre cps dataframe
HER49L_xrf1_Z <- HER49L_xrf1
HER49L_xrf1_Z[, HER49L_elements1] <- scale(HER49L_xrf1[, HER49L_elements1], center = TRUE, scale = TRUE)
HER49L_xrf1_Z

# Standardise and centre inc normalised dataframe
HER49L_xrf1_inc_Z <- HER49L_xrf1_inc
HER49L_xrf1_inc_Z[, HER49L_elements1] <- scale(HER49L_xrf1_inc[, HER49L_elements1], center = TRUE, scale = TRUE)
HER49L_xrf1_inc_Z

# Standardize and center Ln Ti-normalized data
HER49L_xrf1_Ln_Ti_Z <- HER49L_xrf1_Ln_Ti
HER49L_xrf1_Ln_Ti_Z[, HER49L_elements1] <- scale(HER49L_xrf1_Ln_Ti[, HER49L_elements1], center = TRUE, scale = TRUE)
HER49L_xrf1_Ln_Ti_Z

# -------------------------------------------------------------------------

# Save original (xrf), acf filtered (xrf1) and transfromed cps dataframes to csv
write.csv(HER49L_xrf,"HER49L/Output/Composite/HER49L_xrf.csv", row.names = FALSE)
write.csv(HER49L_xrf1,"HER49L/Output/Composite/HER49L_xrf1.csv", row.names = FALSE)
write.csv(HER49L_xrf1_pc_cps,"HER49L/Output/Composite/HER49L_xrf1_pc_cps.csv", row.names = FALSE)
write.csv(HER49L_xrf1_inc,"HER49L/Output/Composite/HER49L_xrf1_inc.csv", row.names = FALSE)
write.csv(HER49L_xrf1_Ti,"HER49L/Output/Composite/HER49L_xrf1_Ti.csv", row.names = FALSE)
write.csv(HER49L_xrf1_Ln_Ti,"HER49L/Output/Composite/HER49L_xrf1_Ln_Ti.csv", row.names = FALSE)
write.csv(HER49L_xrf1_Z,"HER49L/Output/Composite/HER49L_xrf1_Z.csv", row.names = FALSE)
write.csv(HER49L_xrf1_inc_Z,"HER49L/Output/Composite/HER49L_xrf1_inc_Z.csv", row.names = FALSE)
write.csv(HER49L_xrf1_Ln_Ti_Z,"HER49L/Output/Composite/HER49L_xrf1_Ln_Ti_Z.csv", row.names = FALSE)

# %cps sum filtering  ---------------------------------------------------

# mean %cps sum > 0.05%
HER49L_mean1 <- function(x){
  is.numeric(x) && (mean(x, na.rm = TRUE) > 0.05)
}

HER49Lmean1 <- select(HER49L_xrf1_pc_cps, Site:MSE,
                       any_of(acfElementsList) &
                         where(HER49L_mean1))
HER49Lmean1List <- select(HER49Lmean1, -c(Site:MSE)) %>%
  names() %>% 
  print()

# max %cps sum >0.5%
HER49L_max1 <- function(x){
  is.numeric(x) && (max(x, na.rm = TRUE) > 0.5)
}
HER49Lmax1 <- select(HER49L_xrf1_pc_cps,Site:MSE,
                      any_of(acfElementsList) &
                        where(HER49L_max1))
HER49Lmax1List <- select(HER49Lmax1, -c(Site:MSE)) %>%
  names() %>% 
  print()

# mean %cps sum > 0.05% and max %cps sum >0.5%
HER49Lmean1_max1 <- select(HER49Lmean1,Site:MSE,
                              any_of(acfElementsList) &
                                where(HER49L_max1))
HER49Lmean1_max1List <- select(HER49Lmean1_max1, -c(Site:MSE)) %>%
  names() %>% 
  print()

# Select cps filtered dataframe to take forward
HER49L_xrf_filter_pc <- HER49Lmean1_max1
HER49L_xrf_filter_pc

HER49L_filterList_pc <- HER49Lmean1_max1List
HER49L_filterList_pc

# Filtered elements as % of cps_sum
HER49L_xrf_pc_cps3 <- HER49L_xrf_pc_cps2 %>% 
  mutate(cps_sum_filtered_pc = rowSums(across(all_of(HER49L_filterList_pc)))) 
HER49L_xrf2_pc_cps_filtered <- select(HER49L_xrf_pc_cps3, 
                                      Site:MSE,all_of(HER49L_filterList_pc), cps_sum_all, cps_sum_filtered, cps_sum_filtered_pc) %>% 
  relocate(K, .before = Ca)
HER49L_xrf2_pc_cps_filtered

HER49L_xrf2_pc_cps <- select(HER49L_xrf1, inc_coh:LnTi_Ca) %>% 
  bind_cols(HER49L_xrf2_pc_cps_filtered) %>%
  relocate(inc_coh:LnTi_Ca, .after = Mo_coh)
HER49L_xrf2_pc_cps

# %cps_sum summary stats table for filtered ITRAX elements
HER49L_xrf2_pc_cps_stats <- HER49L_xrf2_pc_cps %>%
  select(any_of(HER49L_filterList_pc), cps_sum_all, cps_sum_filtered, cps_sum_filtered_pc) %>% 
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()
write.csv(HER49L_xrf2_pc_cps,"HER49L/Output/Composite/HER49L_xrf2_pc_cps.csv", row.names = FALSE)
write.csv(HER49L_xrf2_pc_cps_stats,"HER49L/Output/Composite/HER49L_xrf2_pc_cps_stats.csv", row.names = FALSE)

# -------------------------------------------------------------------------

# clr (centered log ratios) dataset - filtered cps pc elements

library(compositions)

# Cannot run clr with zeroes or NAs
# removing rows with zeros creates too many gaps in the record for peat cores 
# so +1 to all cps
HER49L_clr0 <- HER49L_xrf1 %>%
  replace(is.na(.), 0) %>%
  select(K, any_of(HER49L_filterList_pc))
HER49L_clr1 <- HER49L_clr0 + 1 #can run without +1 but this reduces the number of elements to 8
as_tibble(HER49L_clr1)

# Or replace all 0 with NA and run clr
#HER49L_clr1[HER49L_clr1 == 0] <- NA
HER49L_clr1 <- as_tibble(HER49L_clr1) %>%
  clr() 
head(HER49L_clr1)
tail(HER49L_clr1)

# Add columns back into the clr dataframe
HER49L_xrf1_clr <- select(HER49L_xrf1, Site:MSE, inc_coh:LnTi_Ca) %>% 
  bind_cols(HER49L_clr1) %>%
  relocate(inc_coh:LnTi_Ca, .after = Mo_coh) %>% 
  print()

write.csv(HER49L_xrf1_clr,"HER49L/Output/Composite/HER49L_xrf1_clr.csv", row.names = FALSE)

# -------------------------------------------------------------------------

# SECTION 4.2: START HERE ONCE DATASETS ABOVE HAVE BEEN ESTABLISHED 

# -------------------------------------------------------------------------

#set working directory
setwd("/Users/Steve/Dropbox/BAS/Data/R/Papers/Hodgson_Hermite_2022/Data/ITRAX/")
#check working directory
getwd() 

# used below
HER49L_xrf <- read_csv("HER49L/Output/Composite/HER49L_xrf.csv")
HER49L_xrf1 <- read_csv("HER49L/Output/Composite/HER49L_xrf1.csv")
HER49L_xrf1_pc_cps <- read_csv("HER49L/Output/Composite/HER49L_xrf1_pc_cps.csv")
HER49L_xrf2_pc_cps <- read_csv("HER49L/Output/Composite/HER49L_xrf2_pc_cps.csv")
HER49L_xrf1_inc_Z <- read_csv("HER49L/Output/Composite/HER49L_xrf1_inc_Z.csv")
HER49L_xrf1_clr <- read_csv("HER49L/Output/Composite/HER49L_xrf1_clr.csv")

# not used below
#HER49L_xrf1_inc <- read_csv("HER49L/Output/Composite/HER49L_xrf1_inc.csv")
#HER49L_xrf1_Ti <- read_csv("HER49L/Output/Composite/HER49L_xrf1_Ti.csv")
#HER49L_xrf1_Ln_Ti <- read_csv("HER49L/Output/Composite/HER49L_xrf1_Ln_Ti.csv")
#HER49L_xrf1_Z <- read_csv("HER49L/Output/Composite/HER49L_xrf1_Z.csv")
#HER49L_xrf1_Ln_Ti_Z <- read_csv("HER49L/Output/Composite/HER49L_xrf1_Ln_Ti_Z.csv")

# Assign elements to/from PeriodicTable package to 'elementsList' (then unload package)
data(periodicTable)
elementsList <- periodicTable$symb

# list of all possible elements only
data(periodicTable)
allelements <- c(symb(1:117))
rm(periodicTable)

# Create lists of elements
HER49L_elements <- select(HER49L_xrf, c(Mg:Mo_coh)) %>% # , Zr
  names()
HER49L_elements

# acf & cps filtered elements
HER49L_elements1 <- select(HER49L_xrf1, any_of(elementsList), Mo_inc, Mo_coh, inc_coh, coh_inc) %>% 
  names()
HER49L_elements1

# acf, cps & %cps_sum filtered elements
HER49L_elements2 <- select(HER49L_xrf2_pc_cps, any_of(elementsList), Mo_inc, Mo_coh, inc_coh, coh_inc) %>% 
  names()
HER49L_elements2

# additional ratio list
HER49L_ratioList1 <- select(HER49L_xrf1, inc_coh, coh_inc, LnK_Ti, LnK_Ca, LnTi_Ca) %>% 
  names()
HER49L_ratioList1

#set NA back to 0 for _xrf1 dataframe
HER49L_xrf1 <- HER49L_xrf1 %>% 
  replace(is.na(.), 0) 
HER49L_xrf1

HER49L_xrf1_inc_Z <- HER49L_xrf1_inc_Z %>% 
  replace(is.na(.), 0) 
HER49L_xrf1_inc_Z

# Pivot transformed dataframes to long format for ggplot using filetered elements 2 (%cps sum filtered) --------------

# Find/replace from here to change to filtered elements 1 (cps filtered)

# used below
HER49L_xrf_long <- select(HER49L_xrf,  Site, depth, SH20_age, label, any_of(HER49L_elements2)) %>%
  pivot_longer(any_of(HER49L_elements2), names_to = "param", values_to = "element")
HER49L_xrf1_long <- select(HER49L_xrf1,  Site, depth, SH20_age, label, any_of(HER49L_elements2)) %>%
  pivot_longer(any_of(HER49L_elements2), names_to = "param", values_to = "element")
HER49L_xrf1_inc_Z_long <- select(HER49L_xrf1_inc_Z,  Site, depth, SH20_age, label, MSE, kcps, cps, any_of(HER49L_elements2)) %>%
  pivot_longer(any_of(HER49L_elements2), names_to = "param", values_to = "element")
HER49L_xrf1_pc_cps_long <- select(HER49L_xrf1_pc_cps,  Site, depth, SH20_age, label, any_of(HER49L_elements2)) %>%
  pivot_longer(any_of(HER49L_elements2), names_to = "param", values_to = "element")
HER49L_xrf2_pc_cps_long <- select(HER49L_xrf2_pc_cps,  Site, depth, SH20_age, label, any_of(HER49L_elements2)) %>%
  pivot_longer(any_of(HER49L_elements2), names_to = "param", values_to = "element")
HER49L_xrf1_clr_long <- select(HER49L_xrf1_clr,  Site, depth, SH20_age, label, any_of(HER49L_elements2)) %>%
  pivot_longer(any_of(HER49L_elements2), names_to = "param", values_to = "element")

# not used below
#HER49L_xrf1_inc_long <- select(HER49L_xrf1_inc,  Site, depth, SH20_age, label, any_of(HER49L_elements2)) %>%
#  pivot_longer(any_of(HER49L_elements2), names_to = "param", values_to = "element")
#HER49L_xrf1_Z_long <- select(HER49L_xrf1_Z,  Site, depth, SH20_age, label, MSE, kcps, cps, any_of(HER49L_elements2)) %>%
#  pivot_longer(any_of(HER49L_elements2), names_to = "param", values_to = "element")
#HER49L_xrf1_Ln_Ti_Z_long <- select(HER49L_xrf1_Ln_Ti_Z,  Site, depth, SH20_age, label, MSE, kcps, cps, any_of(HER49L_elements2)) %>%
#  pivot_longer(any_of(HER49L_elements2), names_to = "param", values_to = "element")


# Correlation plots -------------------------------------------------------

library(GGally)
library(dplyr)

# Correlation matrices

# cps - all elements
theme_set(theme_bw(base_size=2))
ggcorr(HER49L_xrf[,HER49L_elements], method = c("everything", "pearson"), 
       size = 3, label = TRUE, label_alpha = TRUE, label_round=2, label_size= 2)
ggsave("HER49L/Figures/Composite/Fig4.1a_Corr_matrix_xrf1_cps.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# cps - filtered elements2 
theme_set(theme_bw(base_size=2))
ggcorr(HER49L_xrf1[,HER49L_elements2], method = c("everything", "pearson"), 
       size = 3, label = TRUE, label_alpha = TRUE, label_round=2, label_size = 3)
ggsave("HER49L/Figures/Composite/Fig4.1b_Corr_matrix_xrf1_cps.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# cps/inc as Z-scores - filtered elements2 
theme_set(theme_bw(base_size=2))
ggcorr(HER49L_xrf1_inc_Z[,HER49L_elements2], method = c("everything", "pearson"), 
       size = 3, label = TRUE, label_alpha = TRUE, label_round=2, label_size = 3)
ggsave("HER49L/Figures/Composite/Fig4.3_Corr_matrix_incZ.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# %cps sum - filtered elements2 
theme_set(theme_bw(base_size=2))
ggcorr(HER49L_xrf1_pc_cps[,HER49L_elements2], method = c("everything", "pearson"), 
       size = 3, label = TRUE, label_alpha = TRUE, label_round=2, label_size = 3)
ggsave("HER49L/Figures/Composite/Fig4.2_Corr_matrix_xrf1_pc_cps_sum.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# clr - filtered elements2 
ggcorr(HER49L_xrf1_clr[,HER49L_elements2], method = c("everything", "pearson"), 
       size = 3, label = TRUE, label_alpha = TRUE, label_round=2, label_size = 3) 
ggsave("HER49L/Figures/Composite/Fig4.4_Corr_matrix_xrf1_clr.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Correlation density matrix plots
theme_set(theme_bw(base_size=8))
# cps - filtered elements2 
ggpairs(HER49L_xrf1, columns = HER49L_elements2, upper = list(continuous = wrap("cor", size = 3)),
        lower = list(continuous = wrap("points", alpha = 0.5, size=0.2)),
        title="Correlation-density plot: cps for ACF filtered elements")
ggsave("HER49L/Figures/Composite/Fig4.5_Corr-den_matrix_xrf1_cps.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# cps/inc as Z-scores - filtered elements2 
ggpairs(HER49L_xrf1_inc_Z, columns = HER49L_elements2, upper = list(continuous = wrap("cor", size = 3)),
        lower = list(continuous = wrap("points", alpha = 0.5, size=0.2)),
        title="Correlation-density plot: cps/inc as Z-scores for ACF filtered elements")
ggsave("HER49L/Figures/Composite/Fig4.7_Corr-den_matrix_xrf1_Z.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# %cps sum - filtered elements2 
ggpairs(HER49L_xrf1_pc_cps, columns = HER49L_elements2, upper = list(continuous = wrap("cor", size = 3)),
        lower = list(continuous = wrap("points", alpha = 0.5, size=0.2)),
        title="Correlation-density plot: %cps sum for ACF filtered elements")
ggsave("HER49L/Figures/Composite/Fig4.6_Corr-den_matrix_xrf1_pc_cps_sum.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# clr - filtered elements2
ggpairs(HER49L_xrf1_clr, columns = HER49L_elements2, upper = list(continuous = wrap("cor", size = 3)),
        lower = list(continuous = wrap("points", alpha = 0.5, size=0.2)),
        title="Correlation-density plot: clr of ACF filtered elements")
ggsave("HER49L/Figures/Composite/Fig4.8_Corr-den_matrix_xrf1_clr.pdf", 
       height = c(30), width = c(30), dpi = 600, units = "cm")

# -------------------------------------------------------------------------

# Section 5: Summary depth and age plots 

# -------------------------------------------------------------------------

library(tidypaleo)
theme_set(theme_paleo(8))

# set plot text themes for image/rad plots
HER42PB_data_theme <- theme(plot.title = element_text(size = 8),
                            axis.title.x=element_text(size = 8),
                            axis.text.x=element_text(size = 8))

# Summary plots vs depth --------------------------------------------------

# cps
HER49L_xrfStrat_cps <- HER49L_xrf1_long  %>%
  filter(param %in% HER49L_elements2) %>%
  mutate(param = fct_relevel(param, HER49L_elements2)) %>%
  ggplot(aes(x = element, y = depth)) +
  geom_lineh(aes(color = label)) +
  #geom_point(size = 0.01) + #don't use for ITRAX - too many datapoints 
  #geom_lineh(size = 0.5) + #this will make a single black line plot
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(param)) +
  labs(x = "Peak area [cps]", y = "Depth [cm]") +
  ggtitle("HER49L-ITRAX: cps (filtered elements)") +
  HER42PB_data_theme

# cps/inc as Z scores
HER49L_xrfStrat_inc_Z <- HER49L_xrf1_inc_Z_long  %>% 
  filter(param %in% HER49L_elements2) %>%
  mutate(param = fct_relevel(param, HER49L_elements2)) %>%
  ggplot(aes(x = element, y = depth)) +
  geom_lineh(aes(color = label)) +
  #geom_point(size = 0.01) +
  #geom_lineh(size = 0.5) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(param)) +
  labs(x = "cps/inc [Z-score]", y = "Depth [cm]") +
  ggtitle("HER49L-ITRAX: cps/inc as Z-scores (filtered elements)") +
  HER42PB_data_theme

# %cps_sum
HER49L_xrfStrat_pc_cps <- HER49L_xrf1_pc_cps_long  %>% 
  filter(param %in% HER49L_elements2) %>%
  mutate(param = fct_relevel(param, HER49L_elements2)) %>%
  ggplot(aes(x = element, y = depth)) +
  geom_lineh(aes(color = label)) +
  #geom_point(size = 0.01) +
  #geom_lineh(size = 0.5) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(param)) +
  labs(x = "%cps sum [%]", y = "Depth [cm]") +
  ggtitle("HER49L-ITRAX: %cps sum (filtered elements)") +
  HER42PB_data_theme

# clr
HER49L_xrfStrat_clr <- HER49L_xrf1_clr_long  %>% 
  filter(param %in% HER49L_elements2) %>%
  mutate(param = fct_relevel(param, HER49L_elements2)) %>%
  ggplot(aes(x = element, y = depth)) +
  geom_lineh(aes(color = label)) +
  #geom_point(size = 0.01) +
  #geom_lineh(size = 0.5) +
  scale_y_reverse() +
  scale_x_continuous(n.breaks = 4) +
  facet_geochem_gridh(vars(param)) +
  labs(x = "clr", y = "Depth [cm]") +
  ggtitle("HER49L-ITRAX: clr - centred log ratio (filtered elements)") +
  HER42PB_data_theme

# Summary plots
ggarrange(HER49L_xrfStrat_cps, HER49L_xrfStrat_inc_Z, nrow = 2)
ggsave("HER49L/Figures/Composite/Fig5.1_Summary_cps_incZ_depth_plots.pdf",
       height = c(30), width = c(30), dpi = 600, units = "cm")

ggarrange(HER49L_xrfStrat_pc_cps, HER49L_xrfStrat_clr, nrow = 2)
ggsave("HER49L/Figures/Composite/Fig5.2_Summary_cpssum_clr_depth_plots.pdf",
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Add CONISS Zones - this takes a long time -----------------------------------------------------

# cps_sum + CONISS
coniss1_HER49L_pc_cps <- HER49L_xrf1_pc_cps_long %>%
  nested_data(qualifiers = c(SH20_age, depth), key = param, value = element, trans = scale) %>%
  nested_chclust_coniss()

HER49L_xrfStrat_pc_cps_CONISS <- HER49L_xrfStrat_pc_cps +
  layer_dendrogram(coniss1_HER49L_pc_cps, aes(y = depth), param = "CONISS") +
  layer_zone_boundaries(coniss1_HER49L_pc_cps, aes(y = depth, colour = "grey", lty = 2, alpha = 0.4)) +
  ggtitle("HER49L-ITRAX: %cps sum (filtered elements) with CONISS") +
  HER42PB_data_theme

# clr + CONISS
coniss2_HER49L_clr <- HER49L_xrf1_clr_long %>%
  nested_data(qualifiers = c(SH20_age, depth), key = param, value = element, trans = scale) %>%
  nested_chclust_coniss()

HER49L_xrfStrat_clr_CONISS <- HER49L_xrfStrat_clr +
  layer_dendrogram(coniss2_HER49L_clr, aes(y = depth), param = "CONISS") +
  layer_zone_boundaries(coniss2_HER49L_clr, aes(y = depth, colour = "grey", lty = 2, alpha = 0.4)) +
  ggtitle("clr (filtered elements) with CONISS") +
  HER42PB_data_theme

ggarrange(HER49L_xrfStrat_pc_cps_CONISS, HER49L_xrfStrat_clr_CONISS, nrow = 2)
ggsave("HER49L/Figures/Composite/Fig5.3_Summary_cpssum_clr_CONISS_depth_plots.pdf",
       height = c(30), width = c(30), dpi = 600, units = "cm")

# Summary plots vs age --------------------------------------------------------

# cps
HER49L_xrfAge_cps <- HER49L_xrf1_long  %>%
  filter(param %in% HER49L_elements2) %>%
  mutate(param = fct_relevel(param, HER49L_elements2)) %>%
  ggplot(aes(x = SH20_age, y = element)) +
  #geom_point(size = 0.01) +
  geom_line(size = 0.5) +
  scale_y_continuous(n.breaks = 4) +
  facet_geochem_grid(vars(param)) +
  labs(x = "Age (cal a BP)", y = "Peak area [cps]") +
  geom_vline(xintercept = c(0, 1000, 2000, 3000, 4000, 10000, 13800), col = "blue", lty = 3, alpha = 0.7, lwd = 0.3) +
  geom_vline(xintercept = c(4200, 8200, 7500, 11750, 16000), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  ggtitle("HER49L-ITRAX: cps (filtered elements)")

# cps/cin as Z-scores
HER49L_xrfAge_inc_Z <- HER49L_xrf1_inc_Z_long  %>% 
  filter(param %in% HER49L_elements2) %>%
  mutate(param = fct_relevel(param, HER49L_elements2)) %>%
  ggplot(aes(x = SH20_age, y = element)) +
  #geom_point(size = 0.01) +
  geom_line(size = 0.5) +
  scale_y_continuous(n.breaks = 4) +
  facet_geochem_grid(vars(param)) +
  labs(x = "Age (cal a BP)", y = "cps/inc [Z-score]") +
  geom_vline(xintercept = c(0, 1000, 2000, 3000, 4000, 10000, 13800), col = "blue", lty = 3, alpha = 0.7, lwd = 0.3) +
  geom_vline(xintercept = c(4200, 8200, 7500, 11750, 16000), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  ggtitle("HER49L-ITRAX: cps/inc as Z-scores (filtered elements)")

# %cps sum
HER49L_xrfAge_pc_cps <- HER49L_xrf1_pc_cps_long  %>% 
  filter(param %in% HER49L_elements2) %>%
  mutate(param = fct_relevel(param, HER49L_elements2)) %>%
  ggplot(aes(x = SH20_age, y = element)) +
  #geom_point(size = 0.01) +
  geom_line(size = 0.5) +
  scale_y_continuous(n.breaks = 4) +
  facet_geochem_grid(vars(param)) +
  labs(x = "Age (cal a BP)", y = "%cps sum [%]") +
  geom_vline(xintercept = c(0, 1000, 2000, 3000, 4000, 10000, 13800), col = "blue", lty = 3, alpha = 0.7, lwd = 0.3) +
  geom_vline(xintercept = c(4200, 8200, 7500, 11750, 16000), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  ggtitle("HER49L-ITRAX: %cps sum (filtered elements)")

# clr
HER49L_xrfAge_clr <- HER49L_xrf1_clr_long  %>% 
  filter(param %in% HER49L_elements2) %>%
  mutate(param = fct_relevel(param, HER49L_elements2)) %>%
  ggplot(aes(x = SH20_age, y = element)) +
  #geom_point(size = 0.01) +
  geom_line(size = 0.5) +
  scale_y_continuous(n.breaks = 4) +
  facet_geochem_grid(vars(param)) +
  labs(x = "Age (cal a BP)", y = "clr") +
  geom_vline(xintercept = c(0, 1000, 2000, 3000, 4000, 10000, 13800), col = "blue", lty = 3, alpha = 0.7, lwd = 0.3) +
  geom_vline(xintercept = c(4200, 8200, 7500, 11750, 16000), col = "blue", lty = 1, alpha = 0.7, lwd = 0.3) +
  ggtitle("HER49L-ITRAX: clr - centred log ratio (filtered elements)")

# Summary plots
ggarrange(HER49L_xrfAge_cps, HER49L_xrfAge_inc_Z, ncol = 2)
ggsave("HER49L/Figures/Composite/Fig5.4_Summary_cps_incZ_age_plots.pdf",
       height = c(30), width = c(30), dpi = 600, units = "cm")

ggarrange(HER49L_xrfAge_pc_cps, HER49L_xrfAge_clr, ncol = 2)
ggsave("HER49L/Figures/Composite/Fig5.5_Summary_cpssum_clr_age_plots.pdf",
       height = c(30), width = c(30), dpi = 600, units = "cm")

# To do: add CONISS Zones as vlines from age-depth matching

# -------------------------------------------------------------------------

# SECTION 6:  Add Optical and radiograph images

# -------------------------------------------------------------------------



# END ---------------------------------------------------------------------
