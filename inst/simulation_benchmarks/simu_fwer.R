library(snow)
library(qvalue)
library(tibble)
library(IHW)
library(OPWeight)
library(OPWpaper)


cl <- makeCluster(15, type = "MPI")		# start zcluster


clusterExport(cl, "qvalue")
clusterExport(cl, "tibble")
clusterExport(cl, "ihw")
clusterExport(cl,"rejections")



alphaVec = seq(.01, .1, .02)
simVal = 1:3
fwer_mat = parSapply(simVal, simu_fwer, m = 10000, alphaVec = alphaVec)



stopCluster(cl)


#------------------------------------------------------------
# save the workspace to the file .RData in the cwd
save.image("smu_fwer.RData")

