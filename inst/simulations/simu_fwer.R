library(snow)
library("IHW")		# independent hypotheis weight
library(tibble)     # tibble data frame
library(OPWeight)
library(OPWpaper)

# for parallel computing---------------

cl <- makeCluster(10, type = "MPI")		# start zcluster

clusterExport(cl,"ihw")
clusterExport(cl,"adj_pvalues")
clusterExport(cl,"rejections")
clusterExport(cl,"tibble")
clusterExport(cl, "simu_fwer")





alphaVec = seq(.01, .1, .02)
simVal = 1:3
fwer_mat = parSapply(simVal, simu_fwer, m = 10000, alphaVec = alphaVec)



stopCluster(cl)


#------------------------------------------------------------
# save the workspace to the file .RData in the cwd
save.image("smu_fwer.RData")



# this code is to load saved workspace from parallel computing and plot
load(".../smu_fwer.RData")


fwer_by_alpha <- matrix(apply(fwer_mat, 1, mean), nrow = 4, byrow = FALSE)


alphaVal = seq(.01, .1, .02)
datError <- data.frame(alphaVal, t(fwer_by_alpha))
colnames(datError) <- c("alpha","BON","PRO_bin","PRO_cont", "IHW")
datError2 <- melt(datError, id.var="alpha")

panel_d = ggplot(datError2, aes(x = alpha, y = value, col=variable)) +
    geom_line(size=1.5) +
    geom_abline(linetype="dashed") +
    xlab(expression(bold(paste("Nominal ",alpha)))) +
    ylab("FWER")+
    scale_x_continuous(limits = c(0.01,0.1), breaks=seq(0.01,0.09,length=5)) +
    #ylim(0,0.9) +
    theme(legend.title = element_blank())+
    theme(axis.title = element_text(face="bold"))+
    theme(panel.background = element_rect(fill = 'white', colour = 'black'))

#df = data.frame(colour=mycolours, last_vals=c(0.08,   0.085,   0.09), label=c("BON","PRO_bin","PRO_cont"))
#panel_d <- pretty_legend(panel_d, df, .092)
panel_d
