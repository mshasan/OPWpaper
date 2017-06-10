library(snow)
library(OPWeight)
library(OPWpaper)


cl <- makeCluster(10, type = "MPI")		# start zcluster


clusterExport(cl, "prob_rank_givenEffect")
clusterExport(cl, "ranksProb_byEffect")


filterEffectVec <- c(seq(0, 1, .2), 2, 3, 5, 8)
nrep = 100000
mVal = 10000
clusterExport(cl, "filterEffectVec")
clusterExport(cl, "nrep")
clusterExport(cl, "mVal")


ranksProb_byEffect_null.2  <- parSapply(cl, 1:length(filterEffectVec), ranksProb_byEffect, null= .2, m=mVal, nrep=nrep, filterEffectVec=filterEffectVec)
ranksProb_byEffect_null.5  <- parSapply(cl, 1:length(filterEffectVec), ranksProb_byEffect, null= .5, m=mVal, nrep=nrep, filterEffectVec=filterEffectVec)
ranksProb_byEffect_null.75 <- parSapply(cl, 1:length(filterEffectVec), ranksProb_byEffect, null=.75, m=mVal, nrep=nrep, filterEffectVec=filterEffectVec)
ranksProb_byEffect_null.9  <- parSapply(cl, 1:length(filterEffectVec), ranksProb_byEffect, null= .9, m=mVal, nrep=nrep, filterEffectVec=filterEffectVec)
ranksProb_byEffect_null.99 <- parSapply(cl, 1:length(filterEffectVec), ranksProb_byEffect, null=.99, m=mVal, nrep=nrep, filterEffectVec=filterEffectVec)


stopCluster(cl)


# # save data
# #---------------------------
# colNames <- c(paste("null.2_ef", filterEffectVec, sep=""), paste("null.5_ef", filterEffectVec,sep=""),
#               paste("null.75_ef",filterEffectVec, sep=""), paste("null.9_ef",filterEffectVec,sep=""),
#               paste("null.99_ef",filterEffectVec, sep=""))
# data_ProbRanks <- data.frame(ranksProb_byEffect_dat1, ranksProb_byEffect_dat2,
#                              ranksProb_byEffect_dat3, ranksProb_byEffect_dat4,
#                              ranksProb_byEffect_dat5)
# colnames(data_ProbRanks) <- colNames
# write.csv(data_ProbRanks, file = "ranksProb_byEffect_m10000.csv", row.names = F)
#
# #dat = read.csv("ranksProb_byEffect_m10000.csv",h=T)
# #ranksProb_byEffect_dat = dat[,31:40]
#
# windows()
# ranksProb_byEffect_dat <- ranksProb_byEffect_dat2
# apply(ranksProb_byEffect_dat,2,sum)
# par(mfrow=c(3,3))
# for(i in 1:9)
# {
#     plot(1:10000,ranksProb_byEffect_dat[,i],type="l")
# }

#################################################################################################################################


save.image("simu_ranksProb_byEffect.RData")

