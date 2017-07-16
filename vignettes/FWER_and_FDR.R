## ----setup, echo = FALSE-------------------------------------------------
knitr::opts_chunk$set(fig.width = 8.5, fig.height = 4)
knitr::opts_chunk$set(tidy = FALSE, cache = FALSE, autodep = TRUE)

## ----loadLib, message=FALSE, warning=FALSE-------------------------------
library(OPWeight)       # library for the proposed method
library(OPWpaper)       
library(ggplot2)
library(reshape2)       # library for the melt function
library(cowplot)        # plot_grid function

## ----load_fwer_data------------------------------------------------------
load(system.file("simulations/results", "simu_fwer.RDATA",
                 package = "OPWpaper"), envir = environment())

## ----fwer----------------------------------------------------------------
fwer_by_alpha <- matrix(apply(fwer_mat, 1, mean), nrow = 4, byrow = FALSE)

alphaVal = seq(.01, .1, .02)
datError <- data.frame(alphaVal, t(fwer_by_alpha))
colnames(datError) <- c("alpha","BON","CRW_bin","CRW_cont", "IHW")
datError2 <- melt(datError, id.var="alpha")

ggplot(datError2, aes(x = alpha, y = value, col=variable)) +
    geom_line(size=1.5) +
    geom_abline(linetype="dashed") +
    xlab(expression(bold(paste("Nominal ",alpha)))) +
    ylab("FWER")+
    scale_x_continuous(limits = c(0.01,0.1), breaks=seq(0.01,0.09,length=5)) +
    #ylim(0,0.9) +
    theme(legend.title = element_blank())+
    theme(axis.title = element_text(face="bold"))+
    theme(panel.background = element_rect(fill = 'white', colour = 'black'))

## ----load_fwerPowerFdrPower_cont_data------------------------------------
load(system.file("simulations/results", "simu_fwerPowerFdrPower_cont.RDATA",
                 package = "OPWpaper"), envir = environment())

## ----legend--------------------------------------------------------------
# this part is for legend------------------------------------------------------
ey_vec <- c(seq(0, 1, .2), 2, 3, 5, 8)

dat_99 <- data.frame(ey_vec, t(FwerPowerFdrPower5e1[13:16,]))
colnames(dat_99) <- c("effectSize", "CRW", "BH", "RDW", "IHW")
dat_99_par <- melt(dat_99[1:6,], id.var = "effectSize")

p_99_par <- ggplot(dat_99_par, aes(x = effectSize, y = value,group = variable, col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "power", title = "null = 99%") +
    theme(legend.title = element_blank())

legend <- get_legend(p_99_par + theme(legend.direction = "horizontal", legend.position = "bottom"))

## ----fwer_ey_equals_et---------------------------------------------------
# plots FWER et = ey (i.e cv =0)
#-------------------------------------------------
p_.5_eq_fwer <- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower2e1, fdr = FALSE, power = FALSE, null = 50, figure = "effectVsFPFP")
p_.9_eq_fwer <- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower4e1, fdr = FALSE, power = FALSE, null = 90, figure = "effectVsFPFP")
p_.99_eq_fwer<- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower5e1, fdr = FALSE, power = FALSE, null = 99, figure = "effectVsFPFP")

p_eq_fwer = plot_grid(p_.5_eq_fwer, p_.9_eq_fwer, p_.99_eq_fwer, ncol = 3, labels = letters[1:3], align = 'hv')
title <- ggdraw() + draw_label(expression(paste("FWER: et = ey, ", alpha, " = .05")))
plot_grid(title, p_eq_fwer, legend, ncol = 1, rel_heights=c(.1, .5, .1))

## ----fwer_ey_notequals_et------------------------------------------------
# plots FWER et ~ Normal(ey, ey/2) (i.e cv = 1/2)
#-------------------------------------------------
p_.5_uneq_fwer <- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower2e2, fdr = FALSE, power = FALSE, null = 50, figure = "effectVsFPFP")
p_.9_uneq_fwer <- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower4e2, fdr = FALSE, power = FALSE, null = 90, figure = "effectVsFPFP")
p_.99_uneq_fwer<- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower5e2, fdr = FALSE, power = FALSE, null = 99, figure = "effectVsFPFP")

p_uneq_fwer = plot_grid(p_.5_uneq_fwer, p_.9_uneq_fwer, p_.99_uneq_fwer, ncol = 3, labels = letters[1:3], align = 'hv')
title <- ggdraw() + draw_label(expression(paste("FWER: et ~ Normal(ey, ey/2), ", alpha, " = .05")))
plot_grid(title, p_uneq_fwer, legend, ncol = 1, rel_heights=c(.1, .5, .1))

## ----fdr_ey_equals_et----------------------------------------------------
# plots FDR et = ey (i.e cv =0)
#-------------------------------------------------
p_.5_eq_fdr <- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower2e1, fdr = TRUE, power = FALSE, null = 50, figure = "effectVsFPFP")
p_.9_eq_fdr <- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower4e1, fdr = TRUE, power = FALSE, null = 90, figure = "effectVsFPFP")
p_.99_eq_fdr<- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower5e1, fdr = TRUE, power = FALSE, null = 99, figure = "effectVsFPFP")

p_eq_fdr = plot_grid(p_.5_eq_fdr, p_.9_eq_fdr, p_.99_eq_fdr, ncol = 3, labels = letters[1:3], align = 'hv')
title <- ggdraw() + draw_label(expression(paste("FDR: et = ey, ", alpha, " = .05")))
plot_grid(title, p_eq_fdr, legend, ncol = 1, rel_heights=c(.1, .5, .1))

## ----fdr_ey_notequals_et-------------------------------------------------
# plots FDR et ~ Normal(ey, ey/2) (i.e cv = 1/2)
#-------------------------------------------------
p_.5_uneq_fdr <- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower2e2, fdr = TRUE, power = FALSE, null = 50, figure = "effectVsFPFP")
p_.9_uneq_fdr <- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower4e2, fdr = TRUE, power = FALSE, null = 90, figure = "effectVsFPFP")
p_.99_uneq_fdr<- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower5e2, fdr = TRUE, power = FALSE, null = 99, figure = "effectVsFPFP")

p_uneq_fdr = plot_grid(p_.5_uneq_fdr, p_.9_uneq_fdr, p_.99_uneq_fdr, ncol = 3, labels = letters[1:3], align = 'hv')
title <- ggdraw() + draw_label(expression(paste("FDR: et ~ Normal(ey, ey/2), ", alpha, " = .05")))
plot_grid(title, p_uneq_fdr, legend, ncol = 1, rel_heights=c(.1, .5, .1))

