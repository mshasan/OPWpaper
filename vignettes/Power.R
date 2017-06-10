## ----setup, echo = FALSE-------------------------------------------------
knitr::opts_chunk$set(fig.width = 8.5, fig.height = 8)
knitr::opts_chunk$set(tidy = FALSE, cache = TRUE, autodep = TRUE)

## ----loadLib, message=FALSE, warning=FALSE-------------------------------
library(OPWeight)       # library for the proposed method
library(OPWpaper)       
library(ggplot2)
library(reshape2)       # library for the melt function
library(cowplot)        # plot_grid function

## ----load_fwerPowerFdrPower_cont_data------------------------------------
power_byEffect_cont <- system.file("simulations/results", package = "OPWpaper")
setwd(power_byEffect_cont)
load("simu_fwerPowerFdrPower_cont.RDATA")

