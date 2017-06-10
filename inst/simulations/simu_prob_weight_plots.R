library(ggplot2)
library(reshape2)
library(cowplot)
library(gridExtra)
library(OPWeight)
library(OPWpaper)

#==================== Start: Example of prob vs. weight=========================

# load the following data generated earlier------------

m = 10000
ranksProb <- read.csv(".../ranksProb_byEffect_m10000.csv", h=T)
ranksWeight <- read.csv(".../Weight_byEffect_cont_m10000.csv", h=T)

#ey_index <- c(6, 16, 26, 36)  # ey = 1
#ey_index <- c(7, 17, 27, 37)  # ey = 2
#null_index <- c(26, 27, 28)   # null = 50%
#null_index <- c(36, 37, 38)   # null = 90%

prob_weight_plots(ey_index = c(7, 17, 27, 37), null_index = c(14, 17, 19),
                  m = 10000, ey = 2, null = 50, prob = ranksProb,
                  weight = ranksWeight)

prob_weight_plots(ey_index = c(7, 17, 27, 37), null_index = c(34, 37, 39),
                  m = 10000, ey = 2, null = 90, prob = ranksProb,
                  weight = ranksWeight)

prob_weight_plots(ey_index = c(6, 16, 26, 36), null_index = c(14, 17, 19),
                  m = 10000, ey = 1, null = 50, prob = ranksProb,
                  weight = ranksWeight)

prob_weight_plots(ey_index = c(6, 16, 26, 36), null_index = c(34, 37, 39),
                  m = 10000, ey = 1, null = 90, prob = ranksProb,
                  weight = ranksWeight)


#==============end: Example of prob vs. weight==================================
