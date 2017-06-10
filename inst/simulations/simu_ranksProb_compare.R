library(snow)
library(OPWeight)
library(OPWpaper)


cl <- makeCluster(10, type = "MPI")		# start zcluster

clusterExport(cl, "ranksProb_compare")

#================================start of parrallelcomputing====================
# this data are generated using parallel cmputing system
# null test size = c(20, 50, 75, 90, 99)
sampleSize = 1000000 #( use atleast 1,000,000)
clusterExport(cl, "sampleSize")
effect <- c(0, 1, 2)

# continuous cases==================================================================
prob_20_0_cont = parLapply(cl, effect, ranksProb_compare, e.one = 0, m0=20, m1=80, sampleSize = sampleSize, effectType = "continuous")
prob_20_1_cont = parLapply(cl, effect, ranksProb_compare, e.one = 1, m0=20, m1=80, sampleSize = sampleSize, effectType = "continuous")
prob_20_2_cont = parLapply(cl, effect, ranksProb_compare, e.one = 2, m0=20, m1=80, sampleSize = sampleSize, effectType = "continuous")

prob_50_0_cont = parLapply(cl, effect, ranksProb_compare, e.one = 0, m0=50, m1=50, sampleSize = sampleSize, effectType = "continuous")
prob_50_1_cont = parLapply(cl, effect, ranksProb_compare, e.one = 1, m0=50, m1=50, sampleSize = sampleSize, effectType = "continuous")
prob_50_2_cont = parLapply(cl, effect, ranksProb_compare, e.one = 2, m0=50, m1=50, sampleSize = sampleSize, effectType = "continuous")

prob_75_0_cont = parLapply(cl, effect, ranksProb_compare, e.one = 0, m0=75, m1=25, sampleSize = sampleSize, effectType = "continuous")
prob_75_1_cont = parLapply(cl, effect, ranksProb_compare, e.one = 1, m0=75, m1=25, sampleSize = sampleSize, effectType = "continuous")
prob_75_2_cont = parLapply(cl, effect, ranksProb_compare, e.one = 2, m0=75, m1=25, sampleSize = sampleSize, effectType = "continuous")

prob_90_0_cont = parLapply(cl, effect, ranksProb_compare, e.one = 0, m0=90, m1=10, sampleSize = sampleSize, effectType = "continuous")
prob_90_1_cont = parLapply(cl, effect, ranksProb_compare, e.one = 1, m0=90, m1=10, sampleSize = sampleSize, effectType = "continuous")
prob_90_2_cont = parLapply(cl, effect, ranksProb_compare, e.one = 2, m0=90, m1=10, sampleSize = sampleSize, effectType = "continuous")

prob_99_0_cont = parLapply(cl, effect, ranksProb_compare, e.one = 0, m0=99, m1=1, sampleSize = sampleSize, effectType = "continuous")
prob_99_1_cont = parLapply(cl, effect, ranksProb_compare, e.one = 1, m0=99, m1=1, sampleSize = sampleSize, effectType = "continuous")
prob_99_2_cont = parLapply(cl, effect, ranksProb_compare, e.one = 2, m0=99, m1=1, sampleSize = sampleSize, effectType = "continuous")


# binary cases==================================================================
prob_20_0_bin = parLapply(cl, effect, ranksProb_compare, e.one = 0, m0=20, m1=80, sampleSize = sampleSize, effectType = "binary")
prob_20_1_bin = parLapply(cl, effect, ranksProb_compare, e.one = 1, m0=20, m1=80, sampleSize = sampleSize, effectType = "binary")
prob_20_2_bin = parLapply(cl, effect, ranksProb_compare, e.one = 2, m0=20, m1=80, sampleSize = sampleSize, effectType = "binary")

prob_50_0_bin = parLapply(cl, effect, ranksProb_compare, e.one = 0, m0=50, m1=50, sampleSize = sampleSize, effectType = "binary")
prob_50_1_bin = parLapply(cl, effect, ranksProb_compare, e.one = 1, m0=50, m1=50, sampleSize = sampleSize, effectType = "binary")
prob_50_2_bin = parLapply(cl, effect, ranksProb_compare, e.one = 2, m0=50, m1=50, sampleSize = sampleSize, effectType = "binary")

prob_75_0_bin = parLapply(cl, effect, ranksProb_compare, e.one = 0, m0=75, m1=25, sampleSize = sampleSize, effectType = "binary")
prob_75_1_bin = parLapply(cl, effect, ranksProb_compare, e.one = 1, m0=75, m1=25, sampleSize = sampleSize, effectType = "binary")
prob_75_2_bin = parLapply(cl, effect, ranksProb_compare, e.one = 2, m0=75, m1=25, sampleSize = sampleSize, effectType = "binary")

prob_90_0_bin = parLapply(cl, effect, ranksProb_compare, e.one = 0, m0=90, m1=10, sampleSize = sampleSize, effectType = "binary")
prob_90_1_bin = parLapply(cl, effect, ranksProb_compare, e.one = 1, m0=90, m1=10, sampleSize = sampleSize, effectType = "binary")
prob_90_2_bin = parLapply(cl, effect, ranksProb_compare, e.one = 2, m0=90, m1=10, sampleSize = sampleSize, effectType = "binary")

prob_99_0_bin = parLapply(cl, effect, ranksProb_compare, e.one = 0, m0=99, m1=1, sampleSize = sampleSize, effectType = "binary")
prob_99_1_bin = parLapply(cl, effect, ranksProb_compare, e.one = 1, m0=99, m1=1, sampleSize = sampleSize, effectType = "binary")
prob_99_2_bin = parLapply(cl, effect, ranksProb_compare, e.one = 2, m0=99, m1=1, sampleSize = sampleSize, effectType = "binary")



stopCluster(cl)


#------------------------------------------------------------
# save the workspace to the file .RData in the cwd
save.image("simu_prob_rank_givenEffect.RData")

#===========================end of parallel computing===========================
