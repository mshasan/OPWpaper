library(ggplot2)
library(reshape2)
library(cowplot)
library(gridExtra)

#==================== Start: Example of prob vs. weight=========================

porb_weight_plots <- function(ey_index, null_index, m, ey, null, prob, weight)
{
    ranks = 1:m

    prob_by_null = data.frame(ranks, prob[ , ey_index])
    names(prob_by_null) <- c("ranks", "20%", "50%", "75%", "90%")
    prob_by_null_melt <- melt(prob_by_null, id.var = "ranks", variable.name = "null prop.")
    prob_plot = ggplot(prob_by_null_melt, aes(x = ranks, y = value, group = `null prop.`, colour = `null prop.`)) +
        geom_line(aes(linetype = `null prop.`), size = 1.5) +
        labs(x = "ranks", y = "p(rank | effect)", title = paste0("et = ey = ", ey))


    weight_by_null = data.frame(ranks, weight[ , ey_index])
    names(weight_by_null) <- c("ranks", "20%", "50%", "75%", "90%")
    weight_by_null_melt <- melt(weight_by_null, id.var = "ranks", variable.name = "null prop.")
    weight_plot = ggplot(weight_by_null_melt, aes(x = ranks, y = value, group = `null prop.`, colour = `null prop.`)) +
        geom_line(aes(linetype = `null prop.`), size = 1.5) +
        labs(x = "ranks", y = "weight", title = paste0("et = ey = ", ey))


    prob_by_effect <- data.frame(ranks, prob[ , null_index])
    colnames(prob_by_effect) <- c("ranks", "1.0", "2.0", "3.0")
    prob_by_effect_melt <- melt(prob_by_effect, id.var = "ranks", variable.name = "effect size")
    prob_plot_by_effect = ggplot(prob_by_effect_melt, aes(x = ranks, y = value, group = `effect size`,
                                                          colour = `effect size`)) +
        geom_line(aes(linetype = `effect size`), size = 1.5) +
        labs(x = "ranks", y = "p(rank | effect)", title = paste0("null = ", null, "%"))


    weight_by_effect <- data.frame(ranks, weight[ , null_index])
    colnames(weight_by_effect) <- c("ranks", "1.0", "2.0", "3.0")
    weight_by_effect_melt <- melt(weight_by_effect, id.var = "ranks", variable.name = "effect size")
    weight_plot_by_effect = ggplot(weight_by_effect_melt, aes(x = ranks, y = value, group = `effect size`,
                                                              colour = `effect size`)) +
        geom_line(aes(linetype = `effect size`), size = 1.5) +
        labs(x = "ranks", y = "weight", title = paste0("null = ", null, "%"))


    plots <- grid.arrange(prob_plot, weight_plot, prob_plot_by_effect, weight_plot_by_effect,
                          ncol = 2,top="Continuous: probabilities and weights vs. ranks")

    return(plots)

}

m = 10000
ranksProb <- read.csv(".../ranksProb_byEffect_m10000.csv", h=T)
ranksWeight <- read.csv(".../Weight_byEffect_cont_m10000.csv",h=T)

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
