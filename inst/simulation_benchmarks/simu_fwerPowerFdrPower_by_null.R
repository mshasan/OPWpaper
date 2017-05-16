# see the influence of the null proportion
#------------------------------------------------

# this code is to load saved workspace from parallel computing
load("U:/Documents/My Research (UGA)/Multiple Hypoetheses/Article-1/simu_fwerPowerFdrPower_cont_null_influence.RDATA")

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
             heights=c(7,1), top = "power vs. prop. of null")
