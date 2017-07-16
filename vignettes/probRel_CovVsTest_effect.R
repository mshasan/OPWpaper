## ----setup, echo = FALSE-------------------------------------------------
knitr::opts_chunk$set(fig.width = 8.5, fig.height = 4)
knitr::opts_chunk$set(tidy = FALSE, cache = FALSE, autodep = TRUE)

## ----loadLib, message=FALSE, warning=FALSE-------------------------------
library(OPWpaper)       
library(ggplot2)
library(reshape2)       # library for the melt function
library(cowplot)        # plot_grid function

## ----load_fwerPowerFdrPower_cont_data------------------------------------
load(system.file("simulations/results", "simu_probRel_CovVsTest_effect.RDATA",
                 package = "OPWpaper"), envir = environment())

## ----legend--------------------------------------------------------------
# nice plots---------
ranks = 1:100

# extract the legend from one of the plots--------------------------------------
datRelaion1 <- data.frame(ranks, prob0_ed2, prob1_ed2, prob_test0_cor.2_ed2, prob_test1_cor.2_ed2)
colnames(datRelaion1) <- c("ranks", "CH0","CH1","TH0","TH1")
datRelaion1_melt1 <- melt(datRelaion1, id.var="ranks")

p_.2 <- ggplot(datRelaion1_melt1, aes(x = ranks, y = value, group = variable, colour = variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ranks", y = "p(rank | effect)", title = "cor = .5") +
    theme(legend.title = element_blank(), legend.position="bottom")

legend_rel <- get_legend(p_.2 + theme(legend.direction="horizontal", legend.position="bottom"))

## ----rel90_ed2-----------------------------------------------------------
# for 90% ---------
dat_.2_.9_ed2<- cbind(prob02_ed2, prob12_ed2, prob_test0_cor.22_ed2, prob_test1_cor.22_ed2)
p_.2_.9_ed2 <-  nice_plots(x_vec = ranks, y_matrix = dat_.2_.9_ed2, cor = .2, figure = "ranksProb")

dat_.5_.9_ed2<- cbind(prob02_ed2, prob12_ed2, prob_test0_cor.52_ed2, prob_test1_cor.52_ed2)
p_.5_.9_ed2 <-  nice_plots(x_vec = ranks, y_matrix = dat_.5_.9_ed2, cor = .5, figure = "ranksProb")

dat_.8_.9_ed2<- cbind(prob02_ed2, prob12_ed2, prob_test0_cor.82_ed2, prob_test1_cor.82_ed2)
p_.8_.9_ed2 <-  nice_plots(x_vec = ranks, y_matrix = dat_.8_.9_ed2, cor = .8, figure = "ranksProb")


p_rel2 = plot_grid(p_.2_.9_ed2, p_.5_.9_ed2, p_.8_.9_ed2, ncol = 3, labels = letters[1:3], align = 'hv')
title2 <- ggdraw() + draw_label("et = 2, m0 = 90, m1 = 10")
plot_grid(title2, p_rel2, legend_rel, ncol = 1, rel_heights=c(.1, 1, .1))

## ----rel50_ed2-----------------------------------------------------------
# for 50% -----------------
dat_.2_.5_ed2<- cbind(prob0_ed2, prob1_ed2, prob_test0_cor.2_ed2, prob_test1_cor.2_ed2)
p_.2_.5_ed2 <-  nice_plots(x_vec = ranks, y_matrix = dat_.2_.5_ed2, cor = .2, figure = "ranksProb")

dat_.5_.5_ed2<- cbind(prob0_ed2, prob1_ed2, prob_test0_cor.5_ed2, prob_test1_cor.5_ed2)
p_.5_.5_ed2 <-  nice_plots(x_vec = ranks, y_matrix = dat_.5_.5_ed2, cor = .5, figure = "ranksProb")

dat_.8_.5_ed2<- cbind(prob0_ed2, prob1_ed2, prob_test0_cor.8_ed2, prob_test1_cor.8_ed2)
p_.8_.5_ed2 <-  nice_plots(x_vec = ranks, y_matrix = dat_.8_.5_ed2, cor = .8, figure = "ranksProb")


p_rel1 = plot_grid(p_.2_.5_ed2, p_.5_.5_ed2, p_.8_.5_ed2, ncol = 3, labels = letters[1:3], align = 'hv')
title1 <- ggdraw() + draw_label("et = 2, m0 = 50, m1 = 50")
plot_grid(title1, p_rel1, legend_rel, ncol = 1, rel_heights=c(.1, 1, .1))

