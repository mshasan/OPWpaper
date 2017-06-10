## ----setup, echo = FALSE-------------------------------------------------
knitr::opts_chunk$set(fig.width = 8, fig.height = 8)
knitr::opts_chunk$set(tidy = FALSE, cache = TRUE, autodep = TRUE)

## ----continuous_effects--------------------------------------------------
nullSize <- c(20, 50, 75, 90, 99)    # proportion of the true null tests
lapply(nullSize, ranksProb_compare_plots, effectType = "continuous")

## ----binary_effects------------------------------------------------------
nullSize <- c(20, 50, 75, 90, 99)
lapply(nullSize, ranksProb_compare_plots, effectType = "binary")

