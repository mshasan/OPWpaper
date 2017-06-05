# load data for continuous power----------------
# this code is to load saved workspace from parallel computing
load(".../simu_fwerPowerFdrPower_cont.RDATA")
#load(".../simu_fwerPowerFdrPower_bin.RDATA")


# plots of power for mean filter effect(ey) = mean test effect(et) (i.e cv =0)
#-------------------------------------------------------------------------------
dat_50 <- data.frame(filterEffectVec, t(FwerPowerFdrPower2e1[13:16,]))
colnames(dat_50) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_50_all <- melt(dat_50, id.var = "effectSize")
p_50_all <- ggplot(dat_50_all, aes(x = effectSize, y = value,group = variable,
                             col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "power", title = "null = 50%") +
    theme(legend.position="none")

dat_50_par <- melt(dat_50[1:6,], id.var = "effectSize")
p_50_par <- ggplot(dat_50_par, aes(x = effectSize, y = value,group = variable,
                             col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "power", title = "null = 50%") +
    theme(legend.position="none")

prow1 <- plot_grid(p_50_all, p_50_par, align = 'hv', ncol = 1)



dat_90 <- data.frame(filterEffectVec, t(FwerPowerFdrPower4e1[13:16,]))
colnames(dat_90) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_90_all <- melt(dat_90, id.var = "effectSize")
p_90_all <- ggplot(dat_90_all, aes(x = effectSize, y = value,group = variable,
                                   col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "power", title = "null = 90%") +
    theme(legend.position="none")

dat_90_par <- melt(dat_90[1:6,], id.var = "effectSize")
p_90_par <- ggplot(dat_90_par, aes(x = effectSize, y = value,group = variable,
                                   col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "power", title = "null = 90%") +
    theme(legend.position="none")

prow2 <- plot_grid(p_90_all, p_90_par, align = 'hv', ncol = 1)



dat_99 <- data.frame(filterEffectVec, t(FwerPowerFdrPower5e1[13:16,]))
colnames(dat_99) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_99_all <- melt(dat_99, id.var = "effectSize")
p_99_all <- ggplot(dat_99_all, aes(x = effectSize, y = value,group = variable,
                                   col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "power", title = "null = 99%") +
    theme(legend.position="none")

dat_99_par <- melt(dat_99[1:6,], id.var = "effectSize")
p_99_par <- ggplot(dat_99_par, aes(x = effectSize, y = value,group = variable,
                                   col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "power", title = "null = 99%") +
    theme(legend.position="none")

prow3 <- plot_grid(p_99_all, p_99_par, align = 'hv', ncol = 1)


p_99_par <- ggplot(dat_99_par, aes(x = effectSize, y = value,group = variable,
                                   col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "power", title = "null = 99%") +
    theme(legend.title = element_blank())

legend_power <- get_legend(p_99_par + theme(legend.direction="horizontal",
                                        legend.position="bottom"))
p_99_par <- p_99_par + theme(legend.position="none")

grid.arrange(arrangeGrob(prow1, prow2, prow3, ncol=3), legend_power, nrow=2,
             heights=c(7,1), top = "Power: et = ey")



# plots of power for
# mean test effect(et) ~ Normal (mean filter effect, mean filter effect/2) (i.e cv = 1/2)
#-----------------------------------------------------------------------------------------
dat_50 <- data.frame(filterEffectVec, t(FwerPowerFdrPower2e2[13:16,]))
colnames(dat_50) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_50_all <- melt(dat_50, id.var = "effectSize")
p_50_all <- ggplot(dat_50_all, aes(x = effectSize, y = value,group = variable,
                             col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "power", title = "null = 50%") +
    theme(legend.position="none")

dat_50_par <- melt(dat_50[1:6,], id.var = "effectSize")
p_50_par <- ggplot(dat_50_par, aes(x = effectSize, y = value,group = variable,
                             col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "power", title = "null = 50%") +
    theme(legend.position="none")

prow1 <- plot_grid(p_50_all, p_50_par, align = 'hv', ncol = 1)



dat_90 <- data.frame(filterEffectVec, t(FwerPowerFdrPower4e2[13:16,]))
colnames(dat_90) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_90_all <- melt(dat_90, id.var = "effectSize")
p_90_all <- ggplot(dat_90_all, aes(x = effectSize, y = value,group = variable,
                                   col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "power", title = "null = 90%") +
    theme(legend.position="none")

dat_90_par <- melt(dat_90[1:6,], id.var = "effectSize")
p_90_par <- ggplot(dat_90_par, aes(x = effectSize, y = value,group = variable,
                                   col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "power", title = "null = 90%") +
    theme(legend.position="none")

prow2 <- plot_grid(p_90_all, p_90_par, align = 'hv', ncol = 1)



dat_99 <- data.frame(filterEffectVec, t(FwerPowerFdrPower5e2[13:16,]))
colnames(dat_99) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_99_all <- melt(dat_99, id.var = "effectSize")
p_99_all <- ggplot(dat_99_all, aes(x = effectSize, y = value,group = variable,
                                   col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "power", title = "null = 99%") +
    theme(legend.position="none")

dat_99_par <- melt(dat_99[1:6,], id.var = "effectSize")
p_99_par <- ggplot(dat_99_par, aes(x = effectSize, y = value,group = variable,
                                   col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "power", title = "null = 99%") +
    theme(legend.position="none")

prow3 <- plot_grid(p_99_all, p_99_par, align = 'hv', ncol = 1)


p_99_par <- ggplot(dat_99_par, aes(x = effectSize, y = value,group = variable,
                                   col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "power", title = "null = 99%") +
    theme(legend.title = element_blank())

legend_power <- get_legend(p_99_par + theme(legend.direction="horizontal",
                                        legend.position="bottom"))
p_99_par <- p_99_par + theme(legend.position="none")

grid.arrange(arrangeGrob(prow1, prow2, prow3, ncol=3), legend_power, nrow=2,
             heights=c(7,1), top = "Power: et ~ Normal(ey, ey/2)")






# see correaltion effect on Power (when cv =0)
#------------------------------------------------------
# load data for continuous power----------------
# this code is to load saved workspace from parallel computing
load(".../simu_FwerPowerFdrPower_cont.RDATA")
#load(".../simu_fwerPowerFdrPower_bin.RDATA")


filterEffectVec <- c(seq(0,1,.2),2,3,5,8)
E = FwerPowerFdrPower5e1
F = FwerPowerFdrPower5f1
G = FwerPowerFdrPower5g1
H = FwerPowerFdrPower5h1
I = FwerPowerFdrPower5i1
corr = c(0,.3,.5,.7,.9)
r = 13
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

grid.arrange(arrangeGrob(gplots[[3]],gplots[[4]],gplots[[5]],gplots[[6]],gplots[[7]],gplots[[8]], nrow=2),
		legend_corr, nrow=2, heights=c(7,1),
             top = "null = 50%, et = ey")




# see the influence of the null proportion
#------------------------------------------------

# this code is to load saved workspace from parallel computing
load(".../simu_fwerPowerFdrPower_cont.RDATA")
#load(".../simu_fwerPowerFdrPower_bin.RDATA")

effectVec <- c(seq(0,1,.2),2,3,5,8)
nullProp <- c(20, 50, 75, 90, 99)

# corr = .3-------------
mat_ef4 <- rbind(FwerPowerFdrPower1f1[13:16, 4], FwerPowerFdrPower2f1[13:16, 4],
	FwerPowerFdrPower3f1[13:16, 4], FwerPowerFdrPower4f1[13:16, 4],
	FwerPowerFdrPower5f1[13:16, 4])
dat_ef4 <- data.frame(nullProp, mat_ef4)
colnames(dat_ef4) <- c("nullProp", "PRO", "BH", "RDW", "IHW")
dat_ef4_melt <- melt(dat_ef4, id.var = "nullProp")
p_ef4 <- ggplot(dat_ef4_melt, aes(x = nullProp, y = value, group = variable,
                                   col = variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "prop. of null", y = "power", title = "ey = 0.6") +
    theme(legend.position="none")


mat_ef6 <- rbind(FwerPowerFdrPower1f1[13:16, 6], FwerPowerFdrPower2f1[13:16, 6],
	FwerPowerFdrPower3f1[13:16, 6], FwerPowerFdrPower4f1[13:16, 6],
	FwerPowerFdrPower5f1[13:16, 6])
dat_ef6 <- data.frame(nullProp, mat_ef6)
colnames(dat_ef6) <- c("nullProp", "PRO", "BH", "RDW", "IHW")
dat_ef6_melt <- melt(dat_ef6, id.var = "nullProp")
p_ef6 <- ggplot(dat_ef6_melt, aes(x = nullProp, y = value, group = variable,
                                   col = variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "prop. of null", y = "power", title = "ey = 1.0") +
    theme(legend.position="none")


mat_ef8 <- rbind(FwerPowerFdrPower1f1[13:16, 8], FwerPowerFdrPower2f1[13:16, 8],
	FwerPowerFdrPower3f1[13:16, 8], FwerPowerFdrPower4f1[13:16, 8],
	FwerPowerFdrPower5f1[13:16, 8])
dat_ef8 <- data.frame(nullProp, mat_ef8)
colnames(dat_ef8) <- c("nullProp", "PRO", "BH", "RDW", "IHW")
dat_ef8_melt <- melt(dat_ef8, id.var = "nullProp")
p_ef8 <- ggplot(dat_ef8_melt, aes(x = nullProp, y = value, group = variable,
                                   col = variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "prop. of null", y = "power", title = "ey = 3.0") +
    theme(legend.position="none")

plot_cor.3 <- plot_grid(p_ef4, p_ef6, p_ef8, align = 'hv', ncol = 3)
title <- ggdraw() + draw_label("corr = 0.3", fontface='bold')
plot_cor.3_title <- plot_grid(title, plot_cor.3, ncol=1, rel_heights=c(0.1, 1))



# corr = .7 -------------
mat_ef4 <- rbind(FwerPowerFdrPower1h1[13:16, 4], FwerPowerFdrPower2h1[13:16, 4],
	FwerPowerFdrPower3h1[13:16, 4], FwerPowerFdrPower4h1[13:16, 4],
	FwerPowerFdrPower5h1[13:16, 4])
dat_ef4 <- data.frame(nullProp, mat_ef4)
colnames(dat_ef4) <- c("nullProp", "PRO", "BH", "RDW", "IHW")
dat_ef4_melt <- melt(dat_ef4, id.var = "nullProp")
p_ef4 <- ggplot(dat_ef4_melt, aes(x = nullProp, y = value, group = variable,
                                   col = variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "prop. of null", y = "power", title = "ey = 0.6") +
    theme(legend.position="none")


mat_ef6 <- rbind(FwerPowerFdrPower1h1[13:16, 6], FwerPowerFdrPower2h1[13:16, 6],
	FwerPowerFdrPower3h1[13:16, 6], FwerPowerFdrPower4h1[13:16, 6],
	FwerPowerFdrPower5h1[13:16, 6])
dat_ef6 <- data.frame(nullProp, mat_ef6)
colnames(dat_ef6) <- c("nullProp", "PRO", "BH", "RDW", "IHW")
dat_ef6_melt <- melt(dat_ef6, id.var = "nullProp")
p_ef6 <- ggplot(dat_ef6_melt, aes(x = nullProp, y = value, group = variable,
                                   col = variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "prop. of null", y = "power", title = "ey = 1.0") +
    theme(legend.position="none")


mat_ef8 <- rbind(FwerPowerFdrPower1h1[13:16, 8], FwerPowerFdrPower2h1[13:16, 8],
	FwerPowerFdrPower3h1[13:16, 8], FwerPowerFdrPower4h1[13:16, 8],
	FwerPowerFdrPower5h1[13:16, 8])
dat_ef8 <- data.frame(nullProp, mat_ef8)
colnames(dat_ef8) <- c("nullProp", "PRO", "BH", "RDW", "IHW")
dat_ef8_melt <- melt(dat_ef8, id.var = "nullProp")
p_ef8 <- ggplot(dat_ef8_melt, aes(x = nullProp, y = value, group = variable,
                                   col = variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "prop. of null", y = "power", title = "ey = 3.0") +
    theme(legend.position="none")

plot_cor.7 <- plot_grid(p_ef4, p_ef6, p_ef8, align = 'hv', ncol = 3)
title <- ggdraw() + draw_label("corr = 0.7", fontface='bold')
plot_cor.7_title <- plot_grid(title, plot_cor.7, ncol=1, rel_heights=c(0.1, 1))


# for the main title------------
p_ef8 <- ggplot(dat_ef8_melt, aes(x = nullProp, y = value, group = variable,
                                   col = variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "prop. of null", y = "power", title = "ey = 3.0") +
    theme(legend.title = element_blank())

legend_null <- get_legend(p_ef8 + theme(legend.direction="horizontal",
                                        legend.position="bottom"))
p_ef8 <- p_ef8 + theme(legend.position="none")

grid.arrange(arrangeGrob(plot_cor.3_title, plot_cor.7_title), legend_null, nrow=2,
             heights=c(7,1), top = "Continuous: power vs. prop. of null")





# for supplementry materials FWER----------

# plots FWER et = ey (i.e cv =0)
#-------------------------------------------------
m0 = 5000 # I averaged over m twice, so make correction by multiplying by m
dat_50 <- data.frame(filterEffectVec, t(FwerPowerFdrPower2e1[1:4,])*m0)
colnames(dat_50) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_50_all <- melt(dat_50, id.var = "effectSize")
p_50_all <- ggplot(dat_50_all, aes(x = effectSize, y = value,group = variable,
                             col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "FWER", title = "null = 50%") +
    theme(legend.position="none")

m0 = 9000
dat_90 <- data.frame(filterEffectVec, t(FwerPowerFdrPower4e1[1:4,])*m0)
colnames(dat_90) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_90_all <- melt(dat_90, id.var = "effectSize")
p_90_all <- ggplot(dat_90_all, aes(x = effectSize, y = value,group = variable,
                                   col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "FWER", title = "null = 90%") +
    theme(legend.position="none")

m0 = 9900
dat_99 <- data.frame(filterEffectVec, t(FwerPowerFdrPower5e1[1:4,])*m0)
colnames(dat_99) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_99_all <- melt(dat_99, id.var = "effectSize")
p_99_all <- ggplot(dat_99_all, aes(x = effectSize, y = value,group = variable,
                                   col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "FWER", title = "null = 99%") +
    theme(legend.title = element_blank())

legend_FWER <- get_legend(p_99_all + theme(legend.direction="horizontal",
                                        legend.position="bottom"))
p_99_all <- p_99_all + theme(legend.position="none")

grid.arrange(arrangeGrob(p_50_all, p_90_all, p_99_all, nrow=1), legend_FWER, nrow=2,
             heights=c(7,1), top = "FWER: et = ey, alpha = .05")




# plots FWER et ~ Normal(ey, ey/2) (i.e cv = 1/2)
#-------------------------------------------------
m0 = 5000
dat_50 <- data.frame(filterEffectVec, t(FwerPowerFdrPower2e2[1:4,])*m0)
colnames(dat_50) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_50_all <- melt(dat_50, id.var = "effectSize")
p_50_all <- ggplot(dat_50_all, aes(x = effectSize, y = value,group = variable,
                             col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "FWER", title = "null = 50%") +
    theme(legend.position="none")

m0 = 9000
dat_90 <- data.frame(filterEffectVec, t(FwerPowerFdrPower4e2[1:4,])*m0)
colnames(dat_90) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_90_all <- melt(dat_90, id.var = "effectSize")
p_90_all <- ggplot(dat_90_all, aes(x = effectSize, y = value,group = variable,
                                   col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "FWER", title = "null = 90%") +
    theme(legend.position="none")

m0 = 9900
dat_99 <- data.frame(filterEffectVec, t(FwerPowerFdrPower5e2[1:4,])*m0)
colnames(dat_99) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_99_all <- melt(dat_99, id.var = "effectSize")
p_99_all <- ggplot(dat_99_all, aes(x = effectSize, y = value,group = variable,
                                   col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "FWER", title = "null = 99%") +
    theme(legend.title = element_blank())

legend_FWER <- get_legend(p_99_all + theme(legend.direction="horizontal",
                                        legend.position="bottom"))
p_99_all <- p_99_all + theme(legend.position="none")

grid.arrange(arrangeGrob(p_50_all, p_90_all, p_99_all, nrow=1), legend_FWER, nrow=2,
             heights=c(7,1), top = "FWER: et ~ Normal(ey, ey/2), alpha = .05")






# for supplementry materials FDR----------

# plots FDR et = ey (i.e cv =0)
#-------------------------------------------------
dat_50 <- data.frame(filterEffectVec, t(FwerPowerFdrPower2e1[9:12,]))
colnames(dat_50) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_50_all <- melt(dat_50, id.var = "effectSize")
p_50_all <- ggplot(dat_50_all, aes(x = effectSize, y = value,group = variable,
                             col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "FDR", title = "null = 50%") +
    theme(legend.position="none")

dat_90 <- data.frame(filterEffectVec, t(FwerPowerFdrPower4e1[9:12,]))
colnames(dat_90) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_90_all <- melt(dat_90, id.var = "effectSize")
p_90_all <- ggplot(dat_90_all, aes(x = effectSize, y = value,group = variable,
                                   col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "FDR", title = "null = 90%") +
    theme(legend.position="none")

dat_99 <- data.frame(filterEffectVec, t(FwerPowerFdrPower5e1[9:12,]))
colnames(dat_99) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_99_all <- melt(dat_99, id.var = "effectSize")
p_99_all <- ggplot(dat_99_all, aes(x = effectSize, y = value,group = variable,
                                   col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "FDR", title = "null = 99%") +
    theme(legend.title = element_blank())

legend_fdr <- get_legend(p_99_all + theme(legend.direction="horizontal",
                                        legend.position="bottom"))
p_99_all <- p_99_all + theme(legend.position="none")

grid.arrange(arrangeGrob(p_50_all, p_90_all, p_99_all, nrow=1), legend_fdr, nrow=2,
             heights=c(7,1), top = "FDR: et = ey, alpha = .05")




# plots FDR et ~ Normal(ey, ey/2) (i.e cv = 1/2)
#-------------------------------------------------
dat_50 <- data.frame(filterEffectVec, t(FwerPowerFdrPower2e2[9:12,]))
colnames(dat_50) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_50_all <- melt(dat_50, id.var = "effectSize")
p_50_all <- ggplot(dat_50_all, aes(x = effectSize, y = value,group = variable,
                             col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "FDR", title = "null = 50%") +
    theme(legend.position="none")

dat_90 <- data.frame(filterEffectVec, t(FwerPowerFdrPower4e2[9:12,]))
colnames(dat_90) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_90_all <- melt(dat_90, id.var = "effectSize")
p_90_all <- ggplot(dat_90_all, aes(x = effectSize, y = value,group = variable,
                                   col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "FDR", title = "null = 90%") +
    theme(legend.position="none")

dat_99 <- data.frame(filterEffectVec, t(FwerPowerFdrPower5e2[9:12,]))
colnames(dat_99) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_99_all <- melt(dat_99, id.var = "effectSize")
p_99_all <- ggplot(dat_99_all, aes(x = effectSize, y = value,group = variable,
                                   col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "FDR", title = "null = 99%") +
    theme(legend.title = element_blank())

legend_fdr <- get_legend(p_99_all + theme(legend.direction="horizontal",
                                        legend.position="bottom"))
p_99_all <- p_99_all + theme(legend.position="none")

grid.arrange(arrangeGrob(p_50_all, p_90_all, p_99_all, nrow=1), legend_fdr, nrow=2,
             heights=c(7,1), top = "FDR: et ~ Normal(ey, ey/2), alpha = .05")











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
dat_50_cv1 <- data.frame(filterEffectVec, t(FwerPowerFdrPower2a2[13:16,]))
colnames(dat_50_cv1) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_50_cv1_all <- melt(dat_50_cv1, id.var = "effectSize")
p_50_cv1 <- ggplot(dat_50_cv1_all, aes(x = effectSize, y = value,group = variable,
                             col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "power", title = "CV = 1") +
    theme(legend.position="none")

dat_50_cv3 <- data.frame(filterEffectVec, t(FwerPowerFdrPower2a3[13:16,]))
colnames(dat_50_cv3) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_50_cv3_all <- melt(dat_50_cv3, id.var = "effectSize")
p_50_cv3 <- ggplot(dat_50_cv3_all, aes(x = effectSize, y = value,group = variable,
                             col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "power", title = "CV = 3") +
    theme(legend.position="none")

dat_50_cv10 <- data.frame(filterEffectVec, t(FwerPowerFdrPower2a5[13:16,]))
colnames(dat_50_cv10) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_50_cv10_all <- melt(dat_50_cv10, id.var = "effectSize")
p_50_cv10 <- ggplot(dat_50_cv10_all, aes(x = effectSize, y = value,group = variable,
                             col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "power", title = "CV = 10") +
    theme(legend.position="none")


dat_50_cv1_par <- data.frame(filterEffectVec[1:6], t(FwerPowerFdrPower2a2[13:16,1:6]))
colnames(dat_50_cv1_par) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_50_cv1_par <- melt(dat_50_cv1_par, id.var = "effectSize")
p_50_cv1_par <- ggplot(dat_50_cv1_par, aes(x = effectSize, y = value,group = variable,
                             col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "power") +
    theme(legend.position="none")


dat_50_cv3_par <- data.frame(filterEffectVec[1:6], t(FwerPowerFdrPower2a3[13:16,1:6]))
colnames(dat_50_cv3_par) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_50_cv3_par <- melt(dat_50_cv3_par, id.var = "effectSize")
p_50_cv3_par <- ggplot(dat_50_cv3_par, aes(x = effectSize, y = value,group = variable,
                             col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "power") +
    theme(legend.position="none")

dat_50_cv10_par <- data.frame(filterEffectVec[1:6], t(FwerPowerFdrPower2a5[13:16,1:6]))
colnames(dat_50_cv10_par) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_50_cv10_par <- melt(dat_50_cv10_par, id.var = "effectSize")
p_50_cv10_par <- ggplot(dat_50_cv10_par, aes(x = effectSize, y = value,group = variable,
                             col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "power") +
    theme(legend.position="none")


p_50_cv10_par <- ggplot(dat_50_cv10_par, aes(x = effectSize, y = value,group = variable,
                                   col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "power") +
    theme(legend.title = element_blank())

legend_misVar <- get_legend(p_50_cv10_par + theme(legend.direction="horizontal",
                                        legend.position="bottom"))
p_50_cv10_par <- p_50_cv10_par + theme(legend.position="none")

grid.arrange(arrangeGrob(p_50_cv1, p_50_cv3, p_50_cv10, p_50_cv1_par, p_50_cv3_par, p_50_cv10_par, nrow=2),
   		legend_misVar, heights=c(7,1), top = "Power: et ~ Normal(ey, CV*ey), null = 50%")



# CV = coefficient of variance (i.e cv = 1, 3, 10)
# 90% null case------------
dat_90_cv1 <- data.frame(filterEffectVec, t(FwerPowerFdrPower4a2[13:16,]))
colnames(dat_90_cv1) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_90_cv1_all <- melt(dat_90_cv1, id.var = "effectSize")
p_90_cv1 <- ggplot(dat_90_cv1_all, aes(x = effectSize, y = value,group = variable,
                             col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "power", title = "CV = 1") +
    theme(legend.position="none")

dat_90_cv3 <- data.frame(filterEffectVec, t(FwerPowerFdrPower4a3[13:16,]))
colnames(dat_90_cv3) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_90_cv3_all <- melt(dat_90_cv3, id.var = "effectSize")
p_90_cv3 <- ggplot(dat_90_cv3_all, aes(x = effectSize, y = value,group = variable,
                             col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "power", title = "CV = 3") +
    theme(legend.position="none")

dat_90_cv10 <- data.frame(filterEffectVec, t(FwerPowerFdrPower4a5[13:16,]))
colnames(dat_90_cv10) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_90_cv10_all <- melt(dat_90_cv10, id.var = "effectSize")
p_90_cv10 <- ggplot(dat_90_cv10_all, aes(x = effectSize, y = value,group = variable,
                             col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "power", title = "CV = 10") +
    theme(legend.position="none")


dat_90_cv1_par <- data.frame(filterEffectVec[1:6], t(FwerPowerFdrPower4a2[13:16,1:6]))
colnames(dat_90_cv1_par) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_90_cv1_par <- melt(dat_90_cv1_par, id.var = "effectSize")
p_90_cv1_par <- ggplot(dat_90_cv1_par, aes(x = effectSize, y = value,group = variable,
                             col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "power") +
    theme(legend.position="none")


dat_90_cv3_par <- data.frame(filterEffectVec[1:6], t(FwerPowerFdrPower4a3[13:16,1:6]))
colnames(dat_90_cv3_par) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_90_cv3_par <- melt(dat_90_cv3_par, id.var = "effectSize")
p_90_cv3_par <- ggplot(dat_90_cv3_par, aes(x = effectSize, y = value,group = variable,
                             col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "power") +
    theme(legend.position="none")

dat_90_cv10_par <- data.frame(filterEffectVec[1:6], t(FwerPowerFdrPower4a5[13:16,1:6]))
colnames(dat_90_cv10_par) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_90_cv10_par <- melt(dat_90_cv10_par, id.var = "effectSize")
p_90_cv10_par <- ggplot(dat_90_cv10_par, aes(x = effectSize, y = value,group = variable,
                             col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "power") +
    theme(legend.position="none")


p_90_cv10_par <- ggplot(dat_90_cv10_par, aes(x = effectSize, y = value,group = variable,
                                   col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ey", y = "power") +
    theme(legend.title = element_blank())

legend_misVar <- get_legend(p_90_cv10_par + theme(legend.direction="horizontal",
                                        legend.position="bottom"))
p_90_cv10_par <- p_90_cv10_par + theme(legend.position="none")

grid.arrange(arrangeGrob(p_90_cv1, p_90_cv3, p_90_cv10, p_90_cv1_par, p_90_cv3_par, p_90_cv10_par, nrow=2),
   		legend_misVar, heights=c(7,1), top = "Power: et ~ Normal(ey, CV*ey), null = 90%")




