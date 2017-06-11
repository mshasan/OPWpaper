# load data for continuous power----------------
# this code is to load saved workspace from parallel computing
load(".../simu_fwerPowerFdrPower_cont.RDATA")
#load(".../simu_fwerPowerFdrPower_bin.RDATA")



# this part is for legend------------------------------------------------------
ey_vec <- c(seq(0, 1, .2), 2, 3, 5, 8)

dat_99 <- data.frame(ey_vec, t(FwerPowerFdrPower5e1[13:16,]))
colnames(dat_99) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_99_par <- melt(dat_99[1:6,], id.var = "effectSize")

p_99_par <- ggplot(dat_99_par, aes(x = effectSize, y = value,group = variable,
                                   col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "power", title = "null = 99%") +
    theme(legend.title = element_blank())

legend <- get_legend(p_99_par + theme(legend.direction = "horizontal",
                                      legend.position = "bottom"))
#--------------end: legend------------------------------------------------------






# plots of power for mean filter effect(ey) = mean test effect(et)
#-------------------------------------------------------------------------------
p_.5_eq_power <- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower2e1, null = 50, figure = "effectVsFPFP")
p_.9_eq_power <- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower4e1, null = 90, figure = "effectVsFPFP")
p_.99_eq_power<- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower5e1, null = 99, figure = "effectVsFPFP")

p_.5_low_ef_eq_power <- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower2e1, null = 50, low_eff_plot = TRUE, figure = "effectVsFPFP")
p_.9_low_ef_eq_power <- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower4e1, null = 90, low_eff_plot = TRUE, figure = "effectVsFPFP")
p_.99_low_ef_eq_power<- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower5e1, null = 99, low_eff_plot = TRUE, figure = "effectVsFPFP")

p_eq_power = plot_grid(p_.5_eq_power, p_.9_eq_power, p_.99_eq_power,
                       p_.5_low_ef_eq_power, p_.9_low_ef_eq_power, p_.99_low_ef_eq_power,
                       ncol = 3, labels = letters[1:3], align = 'hv')
title <- ggdraw() + draw_label("Power: et = ey")
plot_grid(title, p_eq_power, legend, ncol = 1, rel_heights=c(.1, 1, .1))





# plots of power for
# mean test effect(et) ~ Normal (mean filter effect, mean filter effect/2) (i.e cv = 1/2)
#-----------------------------------------------------------------------------------------
p_.5_uneq_power <- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower2e2, null = 50, figure = "effectVsFPFP")
p_.9_uneq_power <- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower4e2, null = 90, figure = "effectVsFPFP")
p_.99_uneq_power<- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower5e2, null = 99, figure = "effectVsFPFP")

p_.5_low_ef_uneq_power <- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower2e2, null = 50, low_eff_plot = TRUE, figure = "effectVsFPFP")
p_.9_low_ef_uneq_power <- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower4e2, null = 90, low_eff_plot = TRUE, figure = "effectVsFPFP")
p_.99_low_ef_uneq_power<- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower5e2, null = 99, low_eff_plot = TRUE, figure = "effectVsFPFP")

p_uneq_power = plot_grid(p_.5_uneq_power, p_.9_uneq_power, p_.99_uneq_power,
                         p_.5_low_ef_uneq_power, p_.9_low_ef_uneq_power, p_.99_low_ef_uneq_power,
                         ncol = 3, labels = letters[1:3], align = 'hv')
title <- ggdraw() + draw_label("Power: et ~ Normal(ey, ey/2)")
plot_grid(title, p_uneq_power, legend, ncol = 1, rel_heights=c(.1, 1, .1))




# for supplementry materials FWER----------
# plots FWER et = ey (i.e cv =0)
#-------------------------------------------------
p_.5_eq_fwer <- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower2e1, fdr = FALSE, power = FALSE, null = 50, figure = "effectVsFPFP")
p_.9_eq_fwer <- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower4e1, fdr = FALSE, power = FALSE, null = 90, figure = "effectVsFPFP")
p_.99_eq_fwer<- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower5e1, fdr = FALSE, power = FALSE, null = 99, figure = "effectVsFPFP")

p_eq_fwer = plot_grid(p_.5_eq_fwer, p_.9_eq_fwer, p_.99_eq_fwer, ncol = 3, labels = letters[1:3], align = 'hv')
title <- ggdraw() + draw_label(expression(paste("FWER: et = ey, ", alpha, " = .05")))
plot_grid(title, p_eq_fwer, legend, ncol = 1, rel_heights=c(.1, .5, .1))




# plots FWER et ~ Normal(ey, ey/2) (i.e cv = 1/2)
#-------------------------------------------------
p_.5_uneq_fwer <- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower2e2, fdr = FALSE, power = FALSE, null = 50, figure = "effectVsFPFP")
p_.9_uneq_fwer <- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower4e2, fdr = FALSE, power = FALSE, null = 90, figure = "effectVsFPFP")
p_.99_uneq_fwer<- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower5e2, fdr = FALSE, power = FALSE, null = 99, figure = "effectVsFPFP")

p_uneq_fwer = plot_grid(p_.5_uneq_fwer, p_.9_uneq_fwer, p_.99_uneq_fwer, ncol = 3, labels = letters[1:3], align = 'hv')
title <- ggdraw() + draw_label(expression(paste("FWER: et ~ Normal(ey, ey/2), ", alpha, " = .05")))
plot_grid(title, p_uneq_fwer, legend, ncol = 1, rel_heights=c(.1, .5, .1))




# for supplementry materials FDR----------
# plots FDR et = ey (i.e cv =0)
#-------------------------------------------------
p_.5_eq_fdr <- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower2e1, fdr = TRUE, power = FALSE, null = 50, figure = "effectVsFPFP")
p_.9_eq_fdr <- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower4e1, fdr = TRUE, power = FALSE, null = 90, figure = "effectVsFPFP")
p_.99_eq_fdr<- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower5e1, fdr = TRUE, power = FALSE, null = 99, figure = "effectVsFPFP")

p_eq_fdr = plot_grid(p_.5_eq_fdr, p_.9_eq_fdr, p_.99_eq_fdr, ncol = 3, labels = letters[1:3], align = 'hv')
title <- ggdraw() + draw_label(expression(paste("FDR: et = ey, ", alpha, " = .05")))
plot_grid(title, p_eq_fdr, legend, ncol = 1, rel_heights=c(.1, .5, .1))




# plots FDR et ~ Normal(ey, ey/2) (i.e cv = 1/2)
#-------------------------------------------------
p_.5_uneq_fdr <- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower2e2, fdr = TRUE, power = FALSE, null = 50, figure = "effectVsFPFP")
p_.9_uneq_fdr <- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower4e2, fdr = TRUE, power = FALSE, null = 90, figure = "effectVsFPFP")
p_.99_uneq_fdr<- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower5e2, fdr = TRUE, power = FALSE, null = 99, figure = "effectVsFPFP")

p_uneq_fdr = plot_grid(p_.5_uneq_fdr, p_.9_uneq_fdr, p_.99_uneq_fdr, ncol = 3, labels = letters[1:3], align = 'hv')
title <- ggdraw() + draw_label(expression(paste("FDR: et ~ Normal(ey, ey/2), ", alpha, " = .05")))
plot_grid(title, p_uneq_fdr, legend, ncol = 1, rel_heights=c(.1, .5, .1))





# see the effect fo the miss variance of the effect on Power
#-----------------------------------------
# load data for continuous power----------------
# this code is to load saved workspace from parallel computing
load(".../simu_fwerPowerFdrPower_missVar_cont.RDATA")
#load(".../simu_fwerPowerFdrPower_missVar_bin.RDATA")


# plots of power for miss variance of the test effect size; et ~ normal(ey, CV*ey)
# CV = coefficient of variance (i.e cv = 1, 3, 10)
# 50% null case
#----------------------------------------------------------------------------
p_cv1 <- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower2a2, cv = 1, figure = "CV")
p_cv3 <- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower2a3, cv = 3, figure = "CV")
p_cv10<- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower2a10,cv = 10,figure = "CV")

p_cv1_low_ef <- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower2a2, cv = 1, low_eff_plot = TRUE, figure = "CV")
p_cv3_low_ef <- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower2a3, cv = 3, low_eff_plot = TRUE, figure = "CV")
p_cv10_low_ef<- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower2a10,cv = 10,low_eff_plot = TRUE, figure = "CV")

p_cv_power_50 = plot_grid(p_cv1, p_cv3, p_cv10, p_cv1_low_ef, p_cv3_low_ef, p_cv10_low_ef,
                          ncol = 3, labels = letters[1:3], align = 'hv')
title <- ggdraw() + draw_label("Power: et ~ Normal(ey, CV*ey), null = 50%")
plot_grid(title, p_cv_power_50, legend, ncol = 1, rel_heights=c(.1, 1, .1))




# CV = coefficient of variance (i.e cv = 1, 3, 10)
# 90% null case------------
#--------------------------------------------------------------------
p_cv1 <- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower4a2, cv = 1, figure = "CV")
p_cv3 <- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower4a3, cv = 3, figure = "CV")
p_cv10<- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower4a10,cv = 10,figure = "CV")

p_cv1_low_ef <- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower4a2, cv = 1, low_eff_plot = TRUE, figure = "CV")
p_cv3_low_ef <- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower4a3, cv = 3, low_eff_plot = TRUE, figure = "CV")
p_cv10_low_ef<- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower4a10,cv = 10,low_eff_plot = TRUE, figure = "CV")

p_cv_power_90 = plot_grid(p_cv1, p_cv3, p_cv10, p_cv1_low_ef, p_cv3_low_ef, p_cv10_low_ef,
                          ncol = 3, labels = letters[1:3], align = 'hv')
title <- ggdraw() + draw_label("Power: et ~ Normal(ey, CV*ey), null = 90%")
plot_grid(title, p_cv_power_90, legend, ncol = 1, rel_heights=c(.1, 1, .1))




# see the influence of the null proportion
#------------------------------------------------
# this code is to load saved workspace from parallel computing
load("U:/Documents/My Research (UGA)/Multiple Hypoetheses/Article-1/simu_fwerPowerFdrPower_cont.RDATA")

nullProp <- c(20, 50, 75, 90, 99)

# corr = .3-------------
mat_ef.6<- rbind(FwerPowerFdrPower1f1[13:16, 4], FwerPowerFdrPower2f1[13:16, 4],
                 FwerPowerFdrPower3f1[13:16, 4], FwerPowerFdrPower4f1[13:16, 4],
                 FwerPowerFdrPower5f1[13:16, 4])
p_ef.6 <- nice_plots(x_vec = nullProp, y_matrix = mat_ef.6, ey = 0.6, figure = "nullPropVsPower")


mat_ef1 <- rbind(FwerPowerFdrPower1f1[13:16, 6], FwerPowerFdrPower2f1[13:16, 6],
                 FwerPowerFdrPower3f1[13:16, 6], FwerPowerFdrPower4f1[13:16, 6],
                 FwerPowerFdrPower5f1[13:16, 6])
p_ef1 <- nice_plots(x_vec = nullProp, y_matrix = mat_ef1, ey = 1.0, figure = "nullPropVsPower")


mat_ef3 <- rbind(FwerPowerFdrPower1f1[13:16, 8], FwerPowerFdrPower2f1[13:16, 8],
                 FwerPowerFdrPower3f1[13:16, 8], FwerPowerFdrPower4f1[13:16, 8],
                 FwerPowerFdrPower5f1[13:16, 8])
p_ef3 <- nice_plots(x_vec = nullProp, y_matrix = mat_ef3, ey = 3.0, figure = "nullPropVsPower")


plot_cor.3 <- plot_grid(p_ef.6, p_ef1, p_ef3, align = 'hv', ncol = 3, labels=letters[1:3])
title <- ggdraw() + draw_label("corr = 0.3", fontface='bold')
plot_cor.3_title <- plot_grid(title, plot_cor.3, ncol=1, rel_heights=c(0.1, 1), align = 'hv')



# corr = .7 -------------
mat_ef.6.7 <- rbind(FwerPowerFdrPower1h1[13:16, 4], FwerPowerFdrPower2h1[13:16, 4],
                    FwerPowerFdrPower3h1[13:16, 4], FwerPowerFdrPower4h1[13:16, 4],
                    FwerPowerFdrPower5h1[13:16, 4])
p_ef.6.7 <- nice_plots(x_vec = nullProp, y_matrix = mat_ef.6.7, ey = 0.6, figure = "nullPropVsPower")


mat_ef1.7 <- rbind(FwerPowerFdrPower1h1[13:16, 6], FwerPowerFdrPower2h1[13:16, 6],
                   FwerPowerFdrPower3h1[13:16, 6], FwerPowerFdrPower4h1[13:16, 6],
                   FwerPowerFdrPower5h1[13:16, 6])
p_ef1.7 <- nice_plots(x_vec = nullProp, y_matrix = mat_ef1.7, ey = 1.0, figure = "nullPropVsPower")


mat_ef3.7 <- rbind(FwerPowerFdrPower1h1[13:16, 8], FwerPowerFdrPower2h1[13:16, 8],
                   FwerPowerFdrPower3h1[13:16, 8], FwerPowerFdrPower4h1[13:16, 8],
                   FwerPowerFdrPower5h1[13:16, 8])
p_ef3.7 <- nice_plots(x_vec = nullProp, y_matrix = mat_ef3.7, ey = 3.0, figure = "nullPropVsPower")


plot_cor.7 <- plot_grid(p_ef.6.7, p_ef1.7, p_ef3.7, align = 'hv', ncol = 3, labels=letters[4:6])
title <- ggdraw() + draw_label("corr = 0.7", fontface='bold')
plot_cor.7_title <- plot_grid(title, plot_cor.7, ncol=1, rel_heights=c(0.1, 1), align = 'hv')


# for the main title------------
p_prop = plot_grid(plot_cor.3_title, plot_cor.7_title, ncol = 1, align = 'hv')
title_main <- ggdraw() + draw_label("Binary: Power vs. prop. of null")
# make sure get legend----------
plot_grid(title_main, p_prop, legend, ncol = 1, rel_heights=c(.1, 1, .1))







# see correaltion effect on Power (i.e cv =0)
#------------------------------------------------------
# load data for continuous power----------------
# this code is to load saved workspace from parallel computing
load("U:/Documents/My Research (UGA)/Multiple Hypoetheses/Article-1/simu_FwerPowerFdrPower_cont.RDATA")


filterEffectVec <- c(seq(0,1,.2),2,3,5,8)

# use 2 for 50% and 4 for 90% nulls-----------
E = FwerPowerFdrPower2e1
F = FwerPowerFdrPower2f1
G = FwerPowerFdrPower2g1
H = FwerPowerFdrPower2h1
I = FwerPowerFdrPower2i1

corr = c(0,.3,.5,.7,.9)       # correlations
r = 13                        # row starts for fdr based power

gplots <- list()
for(e in 3:8)				# effect size index
{
    PRO = c(E[r,e],    F[r,e],    G[r,e],    H[r,e],    I[r,e])
    BH = c(E[(r+1),e],F[(r+1),e],G[(r+1),e],H[(r+1),e],I[(r+1),e])
    RDW = c(E[(r+2),e],F[(r+2),e],G[(r+2),e],H[(r+2),e],I[(r+2),e])
    IHW = c(E[(r+3),e],F[(r+3),e],G[(r+3),e],H[(r+3),e],I[(r+3),e])
    dat = data.frame(corr, PRO, BH, RDW, IHW)
    dat2 = melt(dat, id.var = "corr")
    gplots[[e]] <- ggplot(dat2, aes(x = corr, y = value, group = variable,
                                    col = variable)) +
        geom_line(aes(linetype = variable), size = 1.5) +
        labs(x = "corr", y = "power", title = paste("et = ", filterEffectVec[e])) +
        #theme(legend.title = element_blank())
        theme(legend.position="none")
}

gplots[[8]] <- gplots[[8]] + theme(legend.position="bottom", legend.title = element_blank())

legend_corr <- get_legend(gplots[[8]])

gplots[[8]] <- gplots[[8]] + theme(legend.position="none")


p = plot_grid(gplots[[3]],gplots[[4]],gplots[[5]],gplots[[6]],gplots[[7]],gplots[[8]],
              ncol = 3, labels = letters[1:6], align = 'hv')
title <- ggdraw() + draw_label("Binary: null = 50%, et = ey")
plot_grid(title, p, legend_corr, ncol = 1, rel_heights=c(.1, 1, .1))





