# these codes are for parallel computing

library(snow)
library(OPWpaper)



cl <- makeCluster(15, type = "MPI")		# start zcluster

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

# # save data
# #---------------------------
colNames <- c(paste("null.2_ef", filterEffectVec, sep=""), paste("null.5_ef", filterEffectVec,sep=""),
              paste("null.75_ef",filterEffectVec, sep=""), paste("null.9_ef",filterEffectVec,sep=""),
              paste("null.99_ef",filterEffectVec, sep=""))
data_ProbRanks <- data.frame(ranksProb_byEffect_null.2, ranksProb_byEffect_null.5,
                             ranksProb_byEffect_null.75, ranksProb_byEffect_null.9,
                             ranksProb_byEffect_null.99)
colnames(data_ProbRanks) <- colNames
write.csv(data_ProbRanks, file = "ranksProb_byEffect_m10000.csv", row.names = F)


# save R objects
save.image("ranksProb_byEffect.RData")

