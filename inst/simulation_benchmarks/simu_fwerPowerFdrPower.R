# part - I (continuous: for different effects)
#===============================================================================
library(snow)
library(mvnfast)		# fast generate multi variate normal
library("IHW")		# independent hypotheis weight
library(tibble)     # tibble data frame
library(OPWeight)
library(OPWpaper)


cl <- makeCluster(10, type = "MPI")		# start zcluster

clusterExport(cl,"rmvn")
clusterExport(cl,"ihw")
clusterExport(cl,"adj_pvalues")
clusterExport(cl,"rejections")
clusterExport(cl,"tibble")
clusterExport(cl, "weightSum_by_c")
clusterExport(cl, "roeder_wasserman_weight")
clusterExport(cl, "runif_by_mean")
clusterExport(cl, "test_by_block")
clusterExport(cl, "fwerPowerFdrPower")


filterEffectVec <- c(seq(0, 1, .2), 2, 3, 5, 8)
datWeight_cont <- read.csv("weight_byEffect_cont_m10000.csv", h = TRUE)
simuVal = 1000       # should use at least 1000

clusterExport(cl, "datWeight_cont")
clusterExport(cl, "filterEffectVec")
clusterExport(cl, "simuVal")


# FwerPowerFdrPower1e1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.2, corr= 0, cv= 0, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,1:10])
# FwerPowerFdrPower1e2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.2, corr= 0, cv=.5, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,1:10])
FwerPowerFdrPower1f1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.2, corr=.3, cv= 0, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,1:10])
# FwerPowerFdrPower1f2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.2, corr=.3, cv=.5, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,1:10])
# FwerPowerFdrPower1g1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.2, corr=.5, cv= 0, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,1:10])
# FwerPowerFdrPower1g2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.2, corr=.5, cv=.5, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,1:10])
FwerPowerFdrPower1h1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.2, corr=.7, cv= 0, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,1:10])
# FwerPowerFdrPower1h2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.2, corr=.7, cv=.5, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,1:10])
# FwerPowerFdrPower1i1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.2, corr=.9, cv= 0, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,1:10])
# FwerPowerFdrPower1i2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.2, corr=.9, cv=.5, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,1:10])


FwerPowerFdrPower2e1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.5, corr= 0, cv= 0, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,11:20])
FwerPowerFdrPower2e2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.5, corr= 0, cv=.5, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,11:20])
FwerPowerFdrPower2f1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.5, corr=.3, cv= 0, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,11:20])
# FwerPowerFdrPower2f2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.5, corr=.3, cv=.5, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,11:20])
FwerPowerFdrPower2g1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.5, corr=.5, cv= 0, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,11:20])
# FwerPowerFdrPower2g2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.5, corr=.5, cv=.5, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,11:20])
FwerPowerFdrPower2h1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.5, corr=.7, cv= 0, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,11:20])
# FwerPowerFdrPower2h2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.5, corr=.7, cv=.5, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,11:20])
FwerPowerFdrPower2i1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.5, corr=.9, cv= 0, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,11:20])
# FwerPowerFdrPower2i2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.5, corr=.9, cv=.5, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,11:20])


# FwerPowerFdrPower3e1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.75, corr= 0, cv= 0, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,21:30])
# FwerPowerFdrPower3e2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.75, corr= 0, cv=.5, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,21:30])
FwerPowerFdrPower3f1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.75, corr=.3, cv= 0, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,21:30])
# FwerPowerFdrPower3f2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.75, corr=.3, cv=.5, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,21:30])
# FwerPowerFdrPower3g1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.75, corr=.5, cv= 0, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,21:30])
# FwerPowerFdrPower3g2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.75, corr=.5, cv=.5, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,21:30])
FwerPowerFdrPower3h1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.75, corr=.7, cv= 0, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,21:30])
# FwerPowerFdrPower3h2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.75, corr=.7, cv=.5, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,21:30])
# FwerPowerFdrPower3i1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.75, corr=.9, cv= 0, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,21:30])
# FwerPowerFdrPower3i2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.75, corr=.9, cv=.5, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,21:30])


FwerPowerFdrPower4e1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.9, corr= 0, cv= 0, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,31:40])
FwerPowerFdrPower4e2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.9, corr= 0, cv=.5, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,31:40])
FwerPowerFdrPower4f1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.9, corr=.3, cv= 0, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,31:40])
# FwerPowerFdrPower4f2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.9, corr=.3, cv=.5, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,31:40])
FwerPowerFdrPower4g1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.9, corr=.5, cv= 0, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,31:40])
# FwerPowerFdrPower4g2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.9, corr=.5, cv=.5, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,31:40])
FwerPowerFdrPower4h1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.9, corr=.7, cv= 0, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,31:40])
# FwerPowerFdrPower4h2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.9, corr=.7, cv=.5, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,31:40])
FwerPowerFdrPower4i1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.9, corr=.9, cv= 0, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,31:40])
# FwerPowerFdrPower4i2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.9, corr=.9, cv=.5, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,31:40])


FwerPowerFdrPower5e1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.99, corr= 0, cv= 0, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,41:50])
FwerPowerFdrPower5e2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.99, corr= 0, cv=.5, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,41:50])
FwerPowerFdrPower5f1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.99, corr=.3, cv= 0, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,41:50])
# FwerPowerFdrPower5f2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.99, corr=.3, cv=.5, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,41:50])
FwerPowerFdrPower5g1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.99, corr=.5, cv= 0, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,41:50])
# FwerPowerFdrPower5g2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.99, corr=.5, cv=.5, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,41:50])
FwerPowerFdrPower5h1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.99, corr=.7, cv= 0, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,41:50])
# FwerPowerFdrPower5h2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.99, corr=.7, cv=.5, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,41:50])
FwerPowerFdrPower5i1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.99, corr=.9, cv= 0, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,41:50])
# FwerPowerFdrPower5i2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.99, corr=.9, cv=.5, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,41:50])



stopCluster(cl)



#-------------------------------------------------------------
# save the workspace to the file .RData in the cwd
save.image("simu_fwerPowerFdrPower_cont.RData")
#===============================================================================






#===============================================================================
# part - II (continuous: for different coefficient of variances)
library(snow)
library(mvnfast)		# fast generate multi variate normal
library("IHW")		# independent hypotheis weight
library(tibble)     # tibble data frame
library(OPWeight)
library(OPWpaper)


cl <- makeCluster(10, type = "MPI")		# start zcluster

clusterExport(cl,"rmvn")
clusterExport(cl,"ihw")
clusterExport(cl,"adj_pvalues")
clusterExport(cl,"rejections")
clusterExport(cl,"tibble")
clusterExport(cl, "weightSum_by_c")
clusterExport(cl, "roeder_wasserman_weight")
clusterExport(cl, "runif_by_mean")
clusterExport(cl, "test_by_block")
clusterExport(cl, "fwerPowerFdrPower")



filterEffectVec <- c(seq(0, 1, .2), 2, 3, 5, 8)
datWeight_cont <- read.csv("weight_byEffect_cont_m10000.csv", h = TRUE)
simuVal = 1000       # should use at least 1000

clusterExport(cl, "datWeight")
clusterExport(cl, "filterEffectVec")
clusterExport(cl, "simuVal")




FwerPowerFdrPower2a2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.5, cv= 1, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,11:20])
FwerPowerFdrPower2a3 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.5, cv= 3, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,11:20])
FwerPowerFdrPower2a10<- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.5, cv= 10,effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,11:20])


FwerPowerFdrPower4a2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.9, cv= 1, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,31:40])
FwerPowerFdrPower4a3 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.9, cv= 3, effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,31:40])
FwerPowerFdrPower4a10<- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.9, cv= 10,effectType = "continuous", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_cont[,31:40])




stopCluster(cl)


#-------------------------------------------------------------
# save the workspace to the file .RData in the cwd
save.image("simu_fwerPowerFdrPower_missVar_cont.RData")
#===============================================================================








#------------end: continuous section-----------------------








#===============================================================================
# part - III (Binary: for different effects)
library(snow)
library(mvnfast)		# fast generate multi variate normal
library("IHW")		# independent hypotheis weight
library(tibble)     # tibble data frame
library(OPWeight)
library(OPWpaper)


cl <- makeCluster(10, type = "MPI")		# start zcluster

clusterExport(cl,"rmvn")
clusterExport(cl,"ihw")
clusterExport(cl,"adj_pvalues")
clusterExport(cl,"rejections")
clusterExport(cl,"tibble")
clusterExport(cl, "weightSum_by_c")
clusterExport(cl, "roeder_wasserman_weight")
clusterExport(cl, "runif_by_mean")
clusterExport(cl, "test_by_block")
clusterExport(cl, "fwerPowerFdrPower")


filterEffectVec <- c(seq(0, 1, .2), 2, 3, 5, 8)
datWeight_bin <- read.csv("weight_byEffect_bin_m10000.csv", h = TRUE)
simuVal = 1000       # should use at least 1000

clusterExport(cl, "datWeight_bin")
clusterExport(cl, "filterEffectVec")
clusterExport(cl, "simuVal")


# FwerPowerFdrPower1e1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.2, corr= 0, cv= 0, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,1:10])
# FwerPowerFdrPower1e2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.2, corr= 0, cv=.5, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,1:10])
FwerPowerFdrPower1f1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.2, corr=.3, cv= 0, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,1:10])
# FwerPowerFdrPower1f2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.2, corr=.3, cv=.5, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,1:10])
# FwerPowerFdrPower1g1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.2, corr=.5, cv= 0, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,1:10])
# FwerPowerFdrPower1g2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.2, corr=.5, cv=.5, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,1:10])
FwerPowerFdrPower1h1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.2, corr=.7, cv= 0, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,1:10])
# FwerPowerFdrPower1h2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.2, corr=.7, cv=.5, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,1:10])
# FwerPowerFdrPower1i1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.2, corr=.9, cv= 0, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,1:10])
# FwerPowerFdrPower1i2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.2, corr=.9, cv=.5, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,1:10])





FwerPowerFdrPower2e1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.5, corr= 0, cv= 0, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,11:20])
FwerPowerFdrPower2e2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.5, corr= 0, cv=.5, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,11:20])
FwerPowerFdrPower2f1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.5, corr=.3, cv= 0, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,11:20])
# FwerPowerFdrPower2f2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.5, corr=.3, cv=.5, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,11:20])
FwerPowerFdrPower2g1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.5, corr=.5, cv= 0, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,11:20])
# FwerPowerFdrPower2g2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.5, corr=.5, cv=.5, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,11:20])
FwerPowerFdrPower2h1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.5, corr=.7, cv= 0, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,11:20])
# FwerPowerFdrPower2h2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.5, corr=.7, cv=.5, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,11:20])
FwerPowerFdrPower2i1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.5, corr=.9, cv= 0, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,11:20])
# FwerPowerFdrPower2i2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.5, corr=.9, cv=.5, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,11:20])





# FwerPowerFdrPower3e1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.75, corr= 0, cv= 0, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,21:30])
# FwerPowerFdrPower3e2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.75, corr= 0, cv=.5, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,21:30])
FwerPowerFdrPower3f1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.75, corr=.3, cv= 0, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,21:30])
# FwerPowerFdrPower3f2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.75, corr=.3, cv=.5, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,21:30])
# FwerPowerFdrPower3g1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.75, corr=.5, cv= 0, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,21:30])
# FwerPowerFdrPower3g2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.75, corr=.5, cv=.5, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,21:30])
FwerPowerFdrPower3h1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.75, corr=.7, cv= 0, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,21:30])
# FwerPowerFdrPower3h2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.75, corr=.7, cv=.5, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,21:30])
# FwerPowerFdrPower3i1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.75, corr=.9, cv= 0, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,21:30])
# FwerPowerFdrPower3i2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.75, corr=.9, cv=.5, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,21:30])



FwerPowerFdrPower4e1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.9, corr= 0, cv= 0, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,31:40])
FwerPowerFdrPower4e2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.9, corr= 0, cv=.5, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,31:40])
FwerPowerFdrPower4f1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.9, corr=.3, cv= 0, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,31:40])
# FwerPowerFdrPower4f2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.9, corr=.3, cv=.5, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,31:40])
FwerPowerFdrPower4g1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.9, corr=.5, cv= 0, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,31:40])
# FwerPowerFdrPower4g2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.9, corr=.5, cv=.5, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,31:40])
FwerPowerFdrPower4h1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.9, corr=.7, cv= 0, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,31:40])
# FwerPowerFdrPower4h2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.9, corr=.7, cv=.5, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,31:40])
FwerPowerFdrPower4i1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.9, corr=.9, cv= 0, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,31:40])
# FwerPowerFdrPower4i2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.9, corr=.9, cv=.5, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,31:40])



FwerPowerFdrPower5e1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.99, corr= 0, cv= 0, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,41:50])
FwerPowerFdrPower5e2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.99, corr= 0, cv=.5, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,41:50])
FwerPowerFdrPower5f1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.99, corr=.3, cv= 0, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,41:50])
# FwerPowerFdrPower5f2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.99, corr=.3, cv=.5, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,41:50])
FwerPowerFdrPower5g1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.99, corr=.5, cv= 0, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,41:50])
# FwerPowerFdrPower5g2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.99, corr=.5, cv=.5, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,41:50])
FwerPowerFdrPower5h1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.99, corr=.7, cv= 0, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,41:50])
# FwerPowerFdrPower5h2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.99, corr=.7, cv=.5, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,41:50])
FwerPowerFdrPower5i1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.99, corr=.9, cv= 0, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,41:50])
# FwerPowerFdrPower5i2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.99, corr=.9, cv=.5, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,41:50])



stopCluster(cl)



#-------------------------------------------------------------
# save the workspace to the file .RData in the cwd
save.image("simu_fwerPowerFdrPower_bin.RData")
#===============================================================================






#===============================================================================
# part - IV (binary: for different coefficient of variances)
library(snow)
library(mvnfast)		# fast generate multi variate normal
library("IHW")		# independent hypotheis weight
library(tibble)     # tibble data frame
library(OPWeight)
library(OPWpaper)


cl <- makeCluster(10, type = "MPI")		# start zcluster

clusterExport(cl,"rmvn")
clusterExport(cl,"ihw")
clusterExport(cl,"adj_pvalues")
clusterExport(cl,"rejections")
clusterExport(cl,"tibble")
clusterExport(cl, "weightSum_by_c")
clusterExport(cl, "roeder_wasserman_weight")
clusterExport(cl, "runif_by_mean")
clusterExport(cl, "test_by_block")
clusterExport(cl, "fwerPowerFdrPower")



filterEffectVec <- c(seq(0, 1, .2), 2, 3, 5, 8)
datWeight <- read.csv("weight_byEffect_bin_m10000.csv", h = TRUE)
simuVal = 1000       # should use at least 1000

clusterExport(cl, "datWeight")
clusterExport(cl, "filterEffectVec")
clusterExport(cl, "simuVal")




FwerPowerFdrPower2a2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.5, cv= 1, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,11:20])
FwerPowerFdrPower2a3 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.5, cv= 3, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,11:20])
FwerPowerFdrPower2a10<- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.5, cv= 10,effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,11:20])


FwerPowerFdrPower4a2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.9, cv= 1, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,31:40])
FwerPowerFdrPower4a3 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.9, cv= 3, effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,31:40])
FwerPowerFdrPower4a10<- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower, simu=simuVal, null=.9, cv= 10,effectType = "binary", filterEffectVec=filterEffectVec, datWeightByNull=datWeight_bin[,31:40])




stopCluster(cl)


#-------------------------------------------------------------
# save the workspace to the file .RData in the cwd
save.image("simu_fwerPowerFdrPower_missVar_bin.RData")
#===============================================================================


