
library(snow)
library(OPWeight)
library(OPWpaper)


cl <- makeCluster(10, type = "MPI")		# start zcluster


clusterExport(cl, "prob_rank_givenEffect")
clusterExport(cl, "probRel_filterVstest_effect")


# null = 50% and ed = 2
prob_test0_cor.2_ed2 <- parSapply(cl, 1:100, probRel_filterVstest_effect, rho=.2, H0=0, ed=2, m0=50, m1=50)
prob_test1_cor.2_ed2 <- parSapply(cl, 1:100, probRel_filterVstest_effect, rho=.2, H0=1, ed=2, m0=50, m1=50)

prob_test0_cor.5_ed2 <- parSapply(cl, 1:100, probRel_filterVstest_effect, rho=.5, H0=0, ed=2, m0=50,m1=50)
prob_test1_cor.5_ed2 <- parSapply(cl, 1:100, probRel_filterVstest_effect, rho=.5, H0=1, ed=2, m0=50,m1=50)

prob_test0_cor.8_ed2 <- parSapply(cl, 1:100, probRel_filterVstest_effect, rho=.8, H0=0, ed=2,m0=50,m1=50)
prob_test1_cor.8_ed2 <- parSapply(cl, 1:100, probRel_filterVstest_effect, rho=.8, H0=1, ed=2,m0=50,m1=50)

prob0_ed2 <- parSapply(cl, 1:100, prob_rank_givenEffect, et=0, ey=2, nrep = 10000, m0=50, m1=50)
prob1_ed2 <- parSapply(cl, 1:100, prob_rank_givenEffect, et=2, ey=2, nrep = 10000, m0=50, m1=50)


# null = 90% and ed = 2
prob_test0_cor.22_ed2 <- parSapply(cl, 1:100, probRel_filterVstest_effect, rho=.2, H0=0, ed=2, m0=90, m1=10)
prob_test1_cor.22_ed2 <- parSapply(cl, 1:100, probRel_filterVstest_effect, rho=.2, H0=1, ed=2, m0=90, m1=10)

prob_test0_cor.52_ed2 <- parSapply(cl, 1:100, probRel_filterVstest_effect, rho=.5, H0=0, ed=2, m0=90,m1=10)
prob_test1_cor.52_ed2 <- parSapply(cl, 1:100, probRel_filterVstest_effect, rho=.5, H0=1, ed=2, m0=90,m1=10)

prob_test0_cor.82_ed2 <- parSapply(cl, 1:100, probRel_filterVstest_effect, rho=.8, H0=0, ed=2,m0=90,m1=10)
prob_test1_cor.82_ed2 <- parSapply(cl, 1:100, probRel_filterVstest_effect, rho=.8, H0=1, ed=2,m0=90,m1=10)

prob02_ed2 <- parSapply(cl, 1:100, prob_rank_givenEffect, et=0, ey=2, nrep = 10000, m0=90, m1=10)
prob12_ed2 <- parSapply(cl, 1:100, prob_rank_givenEffect, et=2, ey=2, nrep = 10000, m0=90, m1=10)



# null = 50% and ed = 0
prob_test0_cor.2_ed0 <- parSapply(cl, 1:100, probRel_filterVstest_effect, rho=.2, H0=0, ed=0, m0=50, m1=50)
prob_test1_cor.2_ed0 <- parSapply(cl, 1:100, probRel_filterVstest_effect, rho=.2, H0=1, ed=0, m0=50, m1=50)

prob_test0_cor.5_ed0 <- parSapply(cl, 1:100, probRel_filterVstest_effect, rho=.5, H0=0, ed=0, m0=50,m1=50)
prob_test1_cor.5_ed0 <- parSapply(cl, 1:100, probRel_filterVstest_effect, rho=.5, H0=1, ed=0, m0=50,m1=50)

prob_test0_cor.8_ed0 <- parSapply(cl, 1:100, probRel_filterVstest_effect, rho=.8, H0=0, ed=0,m0=50,m1=50)
prob_test1_cor.8_ed0 <- parSapply(cl, 1:100, probRel_filterVstest_effect, rho=.8, H0=1, ed=0,m0=50,m1=50)

prob0_ed0 <- parSapply(cl, 1:100, prob_rank_givenEffect, et=0, ey=0, nrep = 10000, m0=50, m1=50)
prob1_ed0 <- parSapply(cl, 1:100, prob_rank_givenEffect, et=0, ey=0, nrep = 10000, m0=50, m1=50)


# null = 90% and ed = 0
prob_test0_cor.22_ed0 <- parSapply(cl, 1:100, probRel_filterVstest_effect, rho=.2, H0=0, ed=0, m0=90, m1=10)
prob_test1_cor.22_ed0 <- parSapply(cl, 1:100, probRel_filterVstest_effect, rho=.2, H0=1, ed=0, m0=90, m1=10)

prob_test0_cor.52_ed0 <- parSapply(cl, 1:100, probRel_filterVstest_effect, rho=.5, H0=0, ed=0, m0=90,m1=10)
prob_test1_cor.52_ed0 <- parSapply(cl, 1:100, probRel_filterVstest_effect, rho=.5, H0=1, ed=0, m0=90,m1=10)

prob_test0_cor.82_ed0 <- parSapply(cl, 1:100, probRel_filterVstest_effect, rho=.8, H0=0, ed=0,m0=90,m1=10)
prob_test1_cor.82_ed0 <- parSapply(cl, 1:100, probRel_filterVstest_effect, rho=.8, H0=1, ed=0,m0=90,m1=10)

prob02_ed0 <- parSapply(cl, 1:100, prob_rank_givenEffect, et=0, ey=0, nrep = 10000, m0=90, m1=10)
prob12_ed0 <- parSapply(cl, 1:100, prob_rank_givenEffect, et=0, ey=0, nrep = 10000, m0=90, m1=10)


stopCluster(cl)


#------------------------------------------------------------
# save the workspace to the file .RData in the cwd
save.image("simu_probRel_filterVstest_effect.RData")







