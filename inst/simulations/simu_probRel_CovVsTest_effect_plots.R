library(ggplot2)
library(grid)
library(gridExtra)      # for multiplots in the same page
library(reshape2)       # library for the melt function
library(cowplot)        # plot_grid function


# this code is to load saved workspace from parallel computing
load(".../simu_probRel_filterVstest_effect.RDATA")


# nice plots---------
ranks = 1:100

# extract the legend from one of the plots--------------------------------------
datRelaion1 <- data.frame(ranks, prob0_ed2, prob1_ed2, prob_test0_cor.2_ed2, prob_test1_cor.2_ed2)
colnames(datRelaion1) <- c("ranks", "FH0","FH1","TH0","TH1")
datRelaion1_melt1 <- melt(datRelaion1, id.var="ranks")

p_.2 <- ggplot(datRelaion1_melt1, aes(x = ranks, y = value, group = variable, colour = variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ranks", y = "p(rank | effect)", title = "cor = .5") +
    theme(legend.title = element_blank(), legend.position="bottom")

legend_rel <- get_legend(p_.2 + theme(legend.direction="horizontal", legend.position="bottom"))
#-------------------------------------------------------------------------------


# for 50% and ed = 2-----------------
dat_.2_.5_ed2<- cbind(prob0_ed2, prob1_ed2, prob_test0_cor.2_ed2, prob_test1_cor.2_ed2)
p_.2_.5_ed2 <-  nice_plots(x_vec = ranks, y_matrix = dat_.2_.5_ed2, cor = .2, figure = "ranksProb")

dat_.5_.5_ed2<- cbind(prob0_ed2, prob1_ed2, prob_test0_cor.5_ed2, prob_test1_cor.5_ed2)
p_.5_.5_ed2 <-  nice_plots(x_vec = ranks, y_matrix = dat_.5_.5_ed2, cor = .5, figure = "ranksProb")

dat_.8_.5_ed2<- cbind(prob0_ed2, prob1_ed2, prob_test0_cor.8_ed2, prob_test1_cor.8_ed2)
p_.8_.5_ed2 <-  nice_plots(x_vec = ranks, y_matrix = dat_.8_.5_ed2, cor = .8, figure = "ranksProb")


p_rel1 = plot_grid(p_.2_.5_ed2, p_.5_.5_ed2, p_.8_.5_ed2, ncol = 3, labels = letters[1:3], align = 'hv')
title1 <- ggdraw() + draw_label("et = 2, m0 = 50, m1 = 50")
plot_grid(title1, p_rel1, legend_rel, ncol = 1, rel_heights=c(.1, 1, .1))




# for 90% and ed = 2 ---------
dat_.2_.9_ed2<- cbind(prob02_ed2, prob12_ed2, prob_test0_cor.22_ed2, prob_test1_cor.22_ed2)
p_.2_.9_ed2 <-  nice_plots(x_vec = ranks, y_matrix = dat_.2_.9_ed2, cor = .2, figure = "ranksProb")

dat_.5_.9_ed2<- cbind(prob02_ed2, prob12_ed2, prob_test0_cor.52_ed2, prob_test1_cor.52_ed2)
p_.5_.9_ed2 <-  nice_plots(x_vec = ranks, y_matrix = dat_.5_.9_ed2, cor = .5, figure = "ranksProb")

dat_.8_.9_ed2<- cbind(prob02_ed2, prob12_ed2, prob_test0_cor.82_ed2, prob_test1_cor.82_ed2)
p_.8_.9_ed2 <-  nice_plots(x_vec = ranks, y_matrix = dat_.8_.9_ed2, cor = .8, figure = "ranksProb")


p_rel2 = plot_grid(p_.2_.9_ed2, p_.5_.9_ed2, p_.8_.9_ed2, ncol = 3, labels = letters[1:3], align = 'hv')
title2 <- ggdraw() + draw_label("et = 2, m0 = 90, m1 = 10")
plot_grid(title2, p_rel2, legend_rel, ncol = 1, rel_heights=c(.1, 1, .1))




# for m0=100 and ed = 0-----------------
dat_.2_.5_ed0<- cbind(prob0_ed0, prob1_ed0, prob_test0_cor.2_ed0, prob_test1_cor.2_ed0)
p_.2_.5_ed0 <-  nice_plots(x_vec = ranks, y_matrix = dat_.2_.5_ed0, cor = .2, figure = "ranksProb")

dat_.5_.5_ed0<- cbind(prob0_ed0, prob1_ed0, prob_test0_cor.5_ed0, prob_test1_cor.5_ed0)
p_.5_.5_ed0 <-  nice_plots(x_vec = ranks, y_matrix = dat_.5_.5_ed0, cor = .5, figure = "ranksProb")

dat_.8_.5_ed0<- cbind(prob0_ed0, prob1_ed0, prob_test0_cor.8_ed0, prob_test1_cor.8_ed0)
p_.8_.5_ed0 <-  nice_plots(x_vec = ranks, y_matrix = dat_.8_.5_ed0, cor = .8, figure = "ranksProb")


p_rel3 = plot_grid(p_.2_.5_ed0, p_.5_.5_ed0, p_.8_.5_ed0, ncol = 3, labels = letters[1:3], align = 'hv')
title3 <- ggdraw() + draw_label("et = 0, m0 = 100, m1 = 0")
plot_grid(title3, p_rel3, legend_rel, ncol = 1, rel_heights=c(.1, 1, .1))
