# this codes are for parallel computing
library(snow)
library(mvnfast)	# fast generate multi variate normal
library("IHW")		# independent hypotheis weight
library(tibble)     # tibble data frame
library(OPWeight)
library(OPWpaper)


cl <- makeCluster(20, type = "MPI")		# start zcluster

clusterExport(cl,"rmvn")
clusterExport(cl,"ihw")
clusterExport(cl,"adj_pvalues")
clusterExport(cl,"rejections")
clusterExport(cl,"tibble")
clusterExport(cl,"fwerPowerFdrPower_by_effect")

simuVal = 1000
filterEffectVec <- c(seq(0, 1, .2), 2, 3, 5, 8)

clusterExport(cl, "filterEffectVec")
clusterExport(cl, "simuVal")

# here random = 0 means cv = 0, and random = 1 means cv = 1/2






#------------start: continuous section-----------------------
datWeight_cont <- read.csv("weight_byEffect_cont_m10000.csv", h = TRUE)
clusterExport(cl, "datWeight_cont")


FwerPowerFdrPower2e1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.5, corr= 0, cv= 0,
                                  alpha=.05, groupSize=100, effectType = "continuous", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_cont[,11:20])
FwerPowerFdrPower2e2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.5, corr= 0, cv=.5,
                                  alpha=.05, groupSize=100, effectType = "continuous", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_cont[,11:20])
FwerPowerFdrPower2f1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.5, corr=.3, cv= 0,
                                  alpha=.05, groupSize=100, effectType = "continuous", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_cont[,11:20])
FwerPowerFdrPower2f2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.5, corr=.3, cv=.5,
                                  alpha=.05, groupSize=100, effectType = "continuous", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_cont[,11:20])
FwerPowerFdrPower2g1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.5, corr=.5, cv= 0,
                                  alpha=.05, groupSize=100, effectType = "continuous", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_cont[,11:20])
FwerPowerFdrPower2g2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.5, corr=.5, cv=.5,
                                  alpha=.05, groupSize=100, effectType = "continuous", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_cont[,11:20])
FwerPowerFdrPower2h1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.5, corr=.7, cv= 0,
                                  alpha=.05, groupSize=100, effectType = "continuous", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_cont[,11:20])
FwerPowerFdrPower2h2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.5, corr=.7, cv=.5,
                                  alpha=.05, groupSize=100, effectType = "continuous", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_cont[,11:20])
FwerPowerFdrPower2i1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.5, corr=.9, cv= 0,
                                  alpha=.05, groupSize=100, effectType = "continuous", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_cont[,11:20])
FwerPowerFdrPower2i2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.5, corr=.9, cv=.5,
                                  alpha=.05, groupSize=100, effectType = "continuous", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_cont[,11:20])



FwerPowerFdrPower4e1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.9, corr= 0, cv= 0,
                                  alpha=.05, groupSize=100, effectType = "continuous", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_cont[,31:40])
FwerPowerFdrPower4e2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.9, corr= 0, cv=.5,
                                  alpha=.05, groupSize=100, effectType = "continuous", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_cont[,31:40])
FwerPowerFdrPower4f1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.9, corr=.3, cv= 0,
                                  alpha=.05, groupSize=100, effectType = "continuous", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_cont[,31:40])
FwerPowerFdrPower4f2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.9, corr=.3, cv=.5,
                                  alpha=.05, groupSize=100, effectType = "continuous", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_cont[,31:40])
FwerPowerFdrPower4g1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.9, corr=.5, cv= 0,
                                  alpha=.05, groupSize=100, effectType = "continuous", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_cont[,31:40])
FwerPowerFdrPower4g2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.9, corr=.5, cv=.5,
                                  alpha=.05, groupSize=100, effectType = "continuous", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_cont[,31:40])
FwerPowerFdrPower4h1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.9, corr=.7, cv= 0,
                                  alpha=.05, groupSize=100, effectType = "continuous", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_cont[,31:40])
FwerPowerFdrPower4h2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.9, corr=.7, cv=.5,
                                  alpha=.05, groupSize=100, effectType = "continuous", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_cont[,31:40])
FwerPowerFdrPower4i1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.9, corr=.9, cv= 0,
                                  alpha=.05, groupSize=100, effectType = "continuous", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_cont[,31:40])
FwerPowerFdrPower4i2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.9, corr=.9, cv=.5,
                                  alpha=.05, groupSize=100, effectType = "continuous", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_cont[,31:40])



FwerPowerFdrPower5e1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.99, corr= 0, cv= 0,
                                  alpha=.05, groupSize=100, effectType = "continuous", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_cont[,41:50])
FwerPowerFdrPower5e2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.99, corr= 0, cv=.5,
                                  alpha=.05, groupSize=100, effectType = "continuous", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_cont[,41:50])
FwerPowerFdrPower5f1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.99, corr=.3, cv= 0,
                                  alpha=.05, groupSize=100, effectType = "continuous", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_cont[,41:50])
FwerPowerFdrPower5f2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.99, corr=.3, cv=.5,
                                  alpha=.05, groupSize=100, effectType = "continuous", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_cont[,41:50])
FwerPowerFdrPower5g1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.99, corr=.5, cv= 0,
                                  alpha=.05, groupSize=100, effectType = "continuous", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_cont[,41:50])
FwerPowerFdrPower5g2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.99, corr=.5, cv=.5,
                                  alpha=.05, groupSize=100, effectType = "continuous", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_cont[,41:50])
FwerPowerFdrPower5h1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.99, corr=.7, cv= 0,
                                  alpha=.05, groupSize=100, effectType = "continuous", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_cont[,41:50])
FwerPowerFdrPower5h2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.99, corr=.7, cv=.5,
                                  alpha=.05, groupSize=100, effectType = "continuous", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_cont[,41:50])
FwerPowerFdrPower5i1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.99, corr=.9, cv= 0,
                                  alpha=.05, groupSize=100, effectType = "continuous", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_cont[,41:50])
FwerPowerFdrPower5i2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.99, corr=.9, cv=.5,
                                  alpha=.05, groupSize=100, effectType = "continuous", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_cont[,41:50])

#-------------------------------------------------------------
# save the workspace to the file .RData in the cwd
save.image("simu_fwerPowerFdrPower_cont.RData")







FwerPowerFdrPower2a2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.5, corr=  0, cv= 1,
                                  alpha=.05, groupSize=100, effectType = "continuous", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_cont[,11:20])
FwerPowerFdrPower2a3 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.5, corr=  0, cv= 3,
                                  alpha=.05, groupSize=100, effectType = "continuous", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_cont[,11:20])
FwerPowerFdrPower2a5 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.5, corr=  0, cv=10,
                                  alpha=.05, groupSize=100, effectType = "continuous", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_cont[,11:20])


FwerPowerFdrPower4a2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.9, corr=  0, cv= 1,
                                  alpha=.05, groupSize=100, effectType = "continuous", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_cont[,31:40])
FwerPowerFdrPower4a3 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.9, corr=  0, cv= 3,
                                  alpha=.05, groupSize=100, effectType = "continuous", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_cont[,31:40])
FwerPowerFdrPower4a5 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.9, corr=  0, cv=10,
                                  alpha=.05, groupSize=100, effectType = "continuous", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_cont[,31:40])

#-------------------------------------------------------------
# save the workspace to the file .RData in the cwd
save.image("simu_fwerPowerFdrPower_missVar_cont.RData")



#------------end: continuous section-----------------------






#------------start: binary section-----------------------
datWeight_bin <- read.csv("weight_byEffect_bin_m10000.csv", h = TRUE)
clusterExport(cl, "datWeight_bin")


FwerPowerFdrPower2e1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.5, corr= 0, cv= 0,
                                  alpha=.05, groupSize=100, effectType = "binary", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_bin[,11:20])
FwerPowerFdrPower2e2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.5, corr= 0, cv=.5,
                                  alpha=.05, groupSize=100, effectType = "binary", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_bin[,11:20])
FwerPowerFdrPower2f1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.5, corr=.3, cv= 0,
                                  alpha=.05, groupSize=100, effectType = "binary", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_bin[,11:20])
FwerPowerFdrPower2f2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.5, corr=.3, cv=.5,
                                  alpha=.05, groupSize=100, effectType = "binary", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_bin[,11:20])
FwerPowerFdrPower2g1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.5, corr=.5, cv= 0,
                                  alpha=.05, groupSize=100, effectType = "binary", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_bin[,11:20])
FwerPowerFdrPower2g2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.5, corr=.5, cv=.5,
                                  alpha=.05, groupSize=100, effectType = "binary", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_bin[,11:20])
FwerPowerFdrPower2h1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.5, corr=.7, cv= 0,
                                  alpha=.05, groupSize=100, effectType = "binary", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_bin[,11:20])
FwerPowerFdrPower2h2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.5, corr=.7, cv=.5,
                                  alpha=.05, groupSize=100, effectType = "binary", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_bin[,11:20])
FwerPowerFdrPower2i1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.5, corr=.9, cv= 0,
                                  alpha=.05, groupSize=100, effectType = "binary", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_bin[,11:20])
FwerPowerFdrPower2i2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.5, corr=.9, cv=.5,
                                  alpha=.05, groupSize=100, effectType = "binary", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_bin[,11:20])



FwerPowerFdrPower4e1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.9, corr= 0, cv= 0,
                                  alpha=.05, groupSize=100, effectType = "binary", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_bin[,31:40])
FwerPowerFdrPower4e2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.9, corr= 0, cv=.5,
                                  alpha=.05, groupSize=100, effectType = "binary", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_bin[,31:40])
FwerPowerFdrPower4f1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.9, corr=.3, cv= 0,
                                  alpha=.05, groupSize=100, effectType = "binary", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_bin[,31:40])
FwerPowerFdrPower4f2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.9, corr=.3, cv=.5,
                                  alpha=.05, groupSize=100, effectType = "binary", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_bin[,31:40])
FwerPowerFdrPower4g1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.9, corr=.5, cv= 0,
                                  alpha=.05, groupSize=100, effectType = "binary", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_bin[,31:40])
FwerPowerFdrPower4g2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.9, corr=.5, cv=.5,
                                  alpha=.05, groupSize=100, effectType = "binary", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_bin[,31:40])
FwerPowerFdrPower4h1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.9, corr=.7, cv= 0,
                                  alpha=.05, groupSize=100, effectType = "binary", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_bin[,31:40])
FwerPowerFdrPower4h2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.9, corr=.7, cv=.5,
                                  alpha=.05, groupSize=100, effectType = "binary", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_bin[,31:40])
FwerPowerFdrPower4i1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.9, corr=.9, cv= 0,
                                  alpha=.05, groupSize=100, effectType = "binary", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_bin[,31:40])
FwerPowerFdrPower4i2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.9, corr=.9, cv=.5,
                                  alpha=.05, groupSize=100, effectType = "binary", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_bin[,31:40])



FwerPowerFdrPower5e1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.99, corr= 0, cv= 0,
                                  alpha=.05, groupSize=100, effectType = "binary", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_bin[,41:50])
FwerPowerFdrPower5e2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.99, corr= 0, cv=.5,
                                  alpha=.05, groupSize=100, effectType = "binary", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_bin[,41:50])
FwerPowerFdrPower5f1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.99, corr=.3, cv= 0,
                                  alpha=.05, groupSize=100, effectType = "binary", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_bin[,41:50])
FwerPowerFdrPower5f2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.99, corr=.3, cv=.5,
                                  alpha=.05, groupSize=100, effectType = "binary", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_bin[,41:50])
FwerPowerFdrPower5g1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.99, corr=.5, cv= 0,
                                  alpha=.05, groupSize=100, effectType = "binary", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_bin[,41:50])
FwerPowerFdrPower5g2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.99, corr=.5, cv=.5,
                                  alpha=.05, groupSize=100, effectType = "binary", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_bin[,41:50])
FwerPowerFdrPower5h1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.99, corr=.7, cv= 0,
                                  alpha=.05, groupSize=100, effectType = "binary", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_bin[,41:50])
FwerPowerFdrPower5h2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.99, corr=.7, cv=.5,
                                  alpha=.05, groupSize=100, effectType = "binary", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_bin[,41:50])
FwerPowerFdrPower5i1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.99, corr=.9, cv= 0,
                                  alpha=.05, groupSize=100, effectType = "binary", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_bin[,41:50])
FwerPowerFdrPower5i2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.99, corr=.9, cv=.5,
                                  alpha=.05, groupSize=100, effectType = "binary", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_bin[,41:50])

#-------------------------------------------------------------
# save the workspace to the file .RData in the cwd
save.image("simu_fwerPowerFdrPower_bin.RData")







FwerPowerFdrPower2a2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.5, corr=  0, cv= 1,
                                  alpha=.05, groupSize=100, effectType = "binary", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_bin[,11:20])
FwerPowerFdrPower2a3 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.5, corr=  0, cv= 3,
                                  alpha=.05, groupSize=100, effectType = "binary", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_bin[,11:20])
FwerPowerFdrPower2a5 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.5, corr=  0, cv=10,
                                  alpha=.05, groupSize=100, effectType = "binary", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_bin[,11:20])


FwerPowerFdrPower4a2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.9, corr=  0, cv= 1,
                                  alpha=.05, groupSize=100, effectType = "binary", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_bin[,31:40])
FwerPowerFdrPower4a3 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.9, corr=  0, cv= 3,
                                  alpha=.05, groupSize=100, effectType = "binary", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_bin[,31:40])
FwerPowerFdrPower4a5 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_by_effect, simu=simuVal, null=.9, corr=  0, cv=10,
                                  alpha=.05, groupSize=100, effectType = "binary", filterEffectVec=filterEffectVec,
                                  datWeightByNull=datWeight_bin[,31:40])

#-------------------------------------------------------------
# save the workspace to the file .RData in the cwd
save.image("simu_fwerPowerFdrPower_missVar_bin.RData")

#------------end: binary section-----------------------


stopCluster(cl)
