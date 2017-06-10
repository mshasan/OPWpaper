## ----setup, echo = FALSE-------------------------------------------------
knitr::opts_chunk$set(fig.width = 8, fig.height = 8)
knitr::opts_chunk$set(tidy = FALSE, cache = TRUE, autodep = TRUE)

## ----loadLib, message=FALSE, warning=FALSE-------------------------------
library(OPWeight)       # library for the proposed method
library(OPWpaper)       
library(ggplot2)
library(reshape2)       # library for the melt function
library(cowplot)        # plot_grid function

## ----load_ranks_weight_data----------------------------------------------
ranksProb_dat <- system.file("simulation_benchmarks/results_files/ranksProb_byEffect_m10000.csv", 
                         package = "OPWpaper")
weight_dat <- system.file("simulation_benchmarks/results_files/ranksProb_byEffect_m10000.csv", 
                         package = "OPWpaper")
m = 10000
ranksProb <- read.csv(ranksProb_dat, h = TRUE)
ranksWeight <- read.csv(weight_dat, h = TRUE)

## ----plots---------------------------------------------------------------
prob_weight_plots(ey_index = c(7, 17, 27, 37), null_index = c(26, 27, 28),
                                m = 10000, ey = 2, null = 50, prob = ranksProb,
                                weight = ranksWeight)
                                

prob_weight_plots(ey_index = c(7, 17, 27, 37), null_index = c(34, 37, 39),
                                m = 10000, ey = 2, null = 90, prob = ranksProb,
                                weight = ranksWeight)

prob_weight_plots(ey_index = c(6, 16, 26, 36), null_index = c(26, 27, 28),
                                m = 10000, ey = 1, null = 50, prob = ranksProb,
                                weight = ranksWeight)

prob_weight_plots(ey_index = c(6, 16, 26, 36), null_index = c(34, 37, 39),
                                m = 10000, ey = 1, null = 90, prob = ranksProb,
                                weight = ranksWeight)

