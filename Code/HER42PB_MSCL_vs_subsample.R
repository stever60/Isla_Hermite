# -------------------------------------------------------------------------

# HER42 GRD vs wet density - Approximate correlation

# -------------------------------------------------------------------------

# Subsample vs GRD density datasets

# set working directory
setwd("/Users/Steve/Dropbox/BAS/Data/R/")

# TO DO: use itraxR matching code to take mean+/-SD of min-max depth density values

round_r <- function(x,seed=111, tol=1.e-6) { 
  set.seed(seed) 
  round(ifelse(near(x%%1, 0.5), jitter(x, amount = tol), x))
}

ACE09_density <- read_csv("Papers/Hodgson_Hermite_2022/Data/Subsample/ACE09_density.csv") %>% 
  select(-c(Section, MSCL_ID, depth_min, depth_max)) %>% 
  filter(Site == 'HER42PB') %>% 
  mutate_at(vars(depth), round, 0)
ACE09_density

ACE09_GRD <- read_csv("Papers/Hodgson_Hermite_2022/Data/Subsample/ACE09_GRD.csv") %>% 
  select(-c(Section, MSCL_ID, RCT)) %>% 
  filter(Site == 'HER42PB')
ACE09_GRD

ACE09_GRD_density <- left_join(ACE09_density, ACE09_GRD, by = c("depth" = "depth"))
ACE09_GRD_density

library(ggpubr)

#set up colours
jamaHER42PB <- c("#8F7700FF","#B24745FF")
show_col(jamaHER42PB)

# GLM regression
formula1 <- y ~ poly(x, 1, raw = TRUE)
theme_set(theme_bw(base_size=16) + theme(
  plot.title = element_text(color="black", size=16, face="bold.italic"),
  axis.title.x = element_text(color="black", size=16),
  axis.title.y = element_text(color="black", size=16)
))
HER42_Den_LM <- ggplot(ACE09_GRD_density, aes(x=Wet_den, y=Den1_SAT, color=Site.x)) +
  geom_point() +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, aes(group=Site.x)) +
  stat_cor(label.y = c(1.3), label.x = c(1.8), face = "bold") + 
  stat_regline_equation(label.y = c(1.29), label.x = c(1.8), face = "bold") +
  scale_shape_manual(values = c(21)) +
  scale_fill_manual(values = jamaHER42PB) +
  scale_color_manual(values = jamaHER42PB) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 12, face="bold.italic"), 
        legend.justification = c(0, 1), legend.position = "top", #c(1,0 ),
        legend.background = element_rect(fill=NULL,
                                         size=0.5, linetype=NULL,
                                         colour =NULL)) +
  scale_x_continuous(breaks=seq(0, 16000, 2000)) + 
  scale_y_continuous(breaks=seq(1, 2.5, 0.5))

HER42_Den_LM
ggsave("Papers/Hodgson_Hermite_2022/Data/Subsample/HER42PB_DensityvsGRD.pdf", height = c(20), width = c(30), dpi = 600, units = "cm")
