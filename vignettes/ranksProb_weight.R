## ----setup, echo = FALSE-------------------------------------------------
knitr::opts_chunk$set(fig.width = 8, fig.height = 8)
knitr::opts_chunk$set(tidy = FALSE, cache = TRUE, autodep = TRUE)

## ----plots---------------------------------------------------------------
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

