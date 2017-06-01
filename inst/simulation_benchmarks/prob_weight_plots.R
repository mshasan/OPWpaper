library(ggplot2)
library(reshape2)
library(cowplot)
library(gridExtra)

#==================== Start: Example of prob vs. weight=========================

# load the following data generated earlier------------

m = 10000
ranksProb <- read.csv(".../ranksProb_byEffect_m10000.csv", h=T)
ranksWeight <- read.csv(".../Weight_byEffect_cont_m10000.csv", h=T)

#ey_index <- c(6, 16, 26, 36)  # ey = 1
#ey_index <- c(7, 17, 27, 37)  # ey = 2
#null_index <- c(26, 27, 28)   # null = 50%
#null_index <- c(36, 37, 38)   # null = 90%

ey2_null90 <- porb_weight_plots(ey_index = c(7, 17, 27, 37), null_index = c(36, 37, 38),
                                m = 10000, ey = 2, null = 90, prob = ranksProb,
                                weight = ranksWeight)

ey2_null50 <- porb_weight_plots(ey_index = c(7, 17, 27, 37), null_index = c(26, 27, 28),
                                m = 10000, ey = 2, null = 50, prob = ranksProb,
                                weight = ranksWeight)

ey1_null90 <- porb_weight_plots(ey_index = c(6, 16, 26, 36), null_index = c(36, 37, 38),
                                m = 10000, ey = 1, null = 90, prob = ranksProb,
                                weight = ranksWeight)

ey1_null50 <- porb_weight_plots(ey_index = c(6, 16, 26, 36), null_index = c(26, 27, 28),
                                m = 10000, ey = 1, null = 50, prob = ranksProb,
                                weight = ranksWeight)

# save the plots
save.image("prob_vs_weight.RData")

#==============end: Example of prob vs. weight==================================
