# these codes are for parallel computing

library(snow)
library(OPWeight)
library(OPWpaper)


cl <- makeCluster(15, type = "MPI")		# start zcluster

clusterExport(cl, "weight_byEffect_bin")


dat <- read.csv("ranksProb_byEffect_m10000.csv",h=T)		# this prob. are saved in an excel file
filterEffectVec <- c(seq(0, 1, .2), 2, 3, 5, 8)
mVal = 10000
del = .000001
aVal = .05

clusterExport(cl,"dat")
clusterExport(cl,"filterEffectVec")
clusterExport(cl,"mVal")
clusterExport(cl,"del")
clusterExport(cl,"aVal")



weightByEffect_null.2  <- parSapply(cl,1:length(filterEffectVec),weight_byEffect_bin,alpha=aVal,null=.2,m=mVal,delInterval=del,filterEffectVec=filterEffectVec, datByNull=dat[,1:10])
weightByEffect_null.5  <- parSapply(cl,1:length(filterEffectVec),weight_byEffect_bin,alpha=aVal,null=.5,m=mVal,delInterval=del,filterEffectVec=filterEffectVec, datByNull=dat[,11:20])
weightByEffect_null.75 <- parSapply(cl,1:length(filterEffectVec),weight_byEffect_bin,alpha=aVal,null=.75,m=mVal,delInterval=del,filterEffectVec=filterEffectVec, datByNull=dat[,21:30])
weightByEffect_null.9  <- parSapply(cl,1:length(filterEffectVec),weight_byEffect_bin,alpha=aVal,null=.9,m=mVal,delInterval=del,filterEffectVec=filterEffectVec, datByNull=dat[,31:40])
weightByEffect_null.99 <- parSapply(cl,1:length(filterEffectVec),weight_byEffect_bin,alpha=aVal,null=.99,m=mVal,delInterval=del,filterEffectVec=filterEffectVec, datByNull=dat[,41:50])


# # save data
# #---------------
 datWeightByEffect <- data.frame(weightByEffect_null.2,weightByEffect_null.5,weightByEffect_null.75,weightByEffect_null.9,weightByEffect_null.99)
 colNames <- c(paste("null.2_ef",filterEffectVec ,sep=""),paste("null.5_ef",filterEffectVec ,sep=""),paste("null.75_ef",filterEffectVec ,sep=""),
               paste("null.9_ef",filterEffectVec ,sep=""),paste("null.99_ef",filterEffectVec ,sep=""))
 colnames(datWeightByEffect) <- colNames
 write.csv(datWeightByEffect,file="Weight_byEffect_bin_m10000.csv",row.names = F)


#####################################################################################################################################

save.image("weight_byEffect_bin.RData")




