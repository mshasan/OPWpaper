## ----setup, echo = FALSE-------------------------------------------------
knitr::opts_chunk$set(tidy = FALSE, cache = TRUE, autodep = TRUE)

## ----loadLib-------------------------------------------------------------
library(OPWeight)       # library for the proposed method
library(OPWpaper)
library(ggplot2)
library(reshape2)       # library for the melt function
library(cowplot)        # plot_grid function

## ----load_ranks_data-----------------------------------------------------
ranks_dat <- system.file("simulation_benchmarks/results_files/simu_prob_rank_givenEffect.RDATA", package = "OPWpaper")
load(ranks_dat)

## ----ranks_Compare_Plots-------------------------------------------------
dat00 <- prob_50_0_cont[[1]]
colnames(dat00) <- c("ranks", "SH0","SH1","EH0","EH1","AH0","AH1")
dat00 <- melt(dat00, id.var="ranks")

dat01 <- prob_50_1_cont[[2]]
colnames(dat01) <- c("ranks", "SH0","SH1","EH0","EH1","AH0","AH1")
dat01 <- melt(dat01, id.var="ranks")

dat12 <- prob_50_1_cont[[3]]
colnames(dat12) <- c("ranks", "SH0","SH1","EH0","EH1","AH0","AH1")
dat12 <- melt(dat12, id.var="ranks")

dat02 <- prob_50_2_cont[[2]]
colnames(dat02) <- c("ranks", "SH0","SH1","EH0","EH1","AH0","AH1")
dat02 <- melt(dat02,id.var="ranks")


p00 <- ggplot(dat00, aes(x = ranks, y = value, group = variable,
                       colour = variable)) +
          geom_line(aes(linetype = variable), size = 1.5) +
          labs(x = "Ranks", y = "p(rank | effect)", title = "ey = 0, e.one = 0") +
          theme(legend.position="none") +
          annotate("text", x=50, y=.011, label="P(rank | effect = 0)")


p01 <- ggplot(dat01, aes(x = ranks, y = value, group = variable,
                       colour = variable)) +
          geom_line(aes(linetype = variable), size = 1.5) +
          labs(x = "Ranks", y = "p(rank | effect)", title = "ey ~ U(0, 1), e.one = 1") +
          theme(legend.position="none") +
          annotate("text", x = c(50, 50), y = c(.005, .018),
                   label = c(paste(sprintf('\u2190'),"P(rank | effect = 0)"),
					paste("P(rank | effect = e.one)", sprintf('\u2192'))))

p12 <- ggplot(dat12, aes(x = ranks, y = value, group = variable,
                       colour = variable)) +
          geom_line(aes(linetype = variable), size = 1.5) +
          labs(x = "Ranks", y = "p(rank | effect)", title = "ey ~ U(1, 2), e.one = 1") +
          theme(legend.position="none") +
          annotate("text", x = c(60, 40), y=c(.001, .02),
                   label = c(paste(sprintf('\u2190'),"P(rank | effect = 0)"),
					paste("P(rank | effect = e.one)", sprintf('\u2192'))))

p02 <- ggplot(dat02, aes(x = ranks, y = value, group = variable,
                       colour = variable)) +
          geom_line(aes(linetype = variable), size = 1.5) +
          labs(x = "Ranks", y = "p(rank | effect)", title = "ey ~ U(0, 1), e.one = 2") +
          theme(legend.title = element_blank(), legend.position="bottom") +
          annotate("text", x = c(50, 40), y=c(.04, .15),
                  label = c("P(rank | effect = 0)",paste(sprintf('\u2190'),"P(rank | effect = e.one)")))


# extract the legend from one of the plots
legend_art <- get_legend(p02 + theme(legend.direction="horizontal",
                                 legend.position="bottom"))

p02 = p02 +    theme(legend.position="none")


p = plot_grid(p00, p01, p12, p02, nrow = 2, labels = letters[1:4], align = 'hv')
# now add the title
title <- ggdraw() + draw_label("Continuous: m0 = 50, m1 = 50")
plot_grid(title, p, legend_art, ncol = 1, rel_heights=c(0.1, 1, .1))

## ----continuous_effects--------------------------------------------------
nullSize <- c(20, 50, 75, 90, 99)
lapply(nullSize, ranksProb_compare_plots, effectType = "continuous")

## ----binary_effects------------------------------------------------------
nullSize <- c(20, 50, 75, 90, 99)
lapply(nullSize, ranksProb_compare_plots, effectType = "binary")

