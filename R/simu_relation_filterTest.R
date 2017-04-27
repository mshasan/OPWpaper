library(snow)


cl <- makeCluster(10, type = "MPI")		# start zcluster



# prob = p(rank=k|effect=ey)
#===============================================================================
prob_rank_givenEffect <- function(k, et, ey, nrep = 10000, m0, m1)
{
  t <- rnorm(nrep, et, 1)
  p0 <- pnorm(-t)
  p1 <- pnorm(ey - t)
  mean0 <- (m0 - 1)*p0 + m1*p1 + 1
  mean1 <- m0*p0 + (m1 - 1)*p1 + 1
  var0 <- (m0 - 1)*p0*(1 - p0) + m1*p1*(1 - p1)
  var1 <- m0*p0*(1 - p0) + (m1 - 1)*p1*(1 - p1)
  prob <- ifelse(et == 0, mean(dnorm(k, mean0, sqrt(var0))),
                 mean(dnorm(k, mean1, sqrt(var1))))
  return(prob)
}

clusterExport(cl, "prob_rank_givenEffect")



prob_relation_filterTestEffect <- function(r, rho, H0, ed, m0, m1)
{
  mean_ey = rho*ed
  sd_ey = sqrt(1 - rho^2)
  ey_val = rnorm(5000, mean_ey, sd_ey)
  prob_condition_ey <- function(ey)
  {
    et <- ifelse(H0==0, 0, ey)
    probs_per_ey = prob_rank_givenEffect(k=r, et=et, ey=ey, nrep = 10000,
                                         m0=m0, m1=m1)
    return(probs_per_ey)
  }
  prob_per_r = mean(sapply(ey_val,  prob_condition_ey))
  return(prob_per_r)
}


prob_test0_cor.2 <- parSapply(cl, 1:100, prob_relation_filterTestEffect, rho=.2, H0=0, 
                           ed=2, m0=50, m1=50)
prob_test1_cor.2 <- parSapply(cl, 1:100, prob_relation_filterTestEffect, rho=.2, H0=1, 
                           ed=2, m0=50, m1=50)

prob_test0_cor.5 <- parSapply(cl, 1:100, prob_relation_filterTestEffect, rho=.5, H0=0, 
                           ed=2, m0=50,m1=50)
prob_test1_cor.5 <- parSapply(cl, 1:100, prob_relation_filterTestEffect, rho=.5, H0=1, 
                           ed=2, m0=50,m1=50)

prob_test0_cor.8 <- parSapply(cl, 1:100, prob_relation_filterTestEffect, rho=.8, H0=0, 
                           ed=2,m0=50,m1=50)
prob_test1_cor.8 <- parSapply(cl, 1:100, prob_relation_filterTestEffect, rho=.8, H0=1, 
                           ed=2,m0=50,m1=50)

prob0 <- parSapply(cl, 1:100, prob_rank_givenEffect, et=0, ey=2, nrep = 10000, m0=50, m1=50)
prob1 <- parSapply(cl, 1:100, prob_rank_givenEffect, et=2, ey=2, nrep = 10000, m0=50, m1=50)


prob_test0_cor.22 <- parSapply(cl, 1:100, prob_relation_filterTestEffect, rho=.2, H0=0, 
                              ed=2, m0=90, m1=10)
prob_test1_cor.22 <- parSapply(cl, 1:100, prob_relation_filterTestEffect, rho=.2, H0=1, 
                              ed=2, m0=90, m1=10)

prob_test0_cor.52 <- parSapply(cl, 1:100, prob_relation_filterTestEffect, rho=.5, H0=0, 
                              ed=2, m0=90,m1=10)
prob_test1_cor.52 <- parSapply(cl, 1:100, prob_relation_filterTestEffect, rho=.5, H0=1, 
                              ed=2, m0=90,m1=10)

prob_test0_cor.82 <- parSapply(cl, 1:100, prob_relation_filterTestEffect, rho=.8, H0=0, 
                              ed=2,m0=90,m1=10)
prob_test1_cor.82 <- parSapply(cl, 1:100, prob_relation_filterTestEffect, rho=.8, H0=1, 
                              ed=2,m0=90,m1=10)

prob02 <- parSapply(cl, 1:100, prob_rank_givenEffect, et=0, ey=2, nrep = 10000, m0=90, m1=10)
prob12 <- parSapply(cl, 1:100, prob_rank_givenEffect, et=2, ey=2, nrep = 10000, m0=90, m1=10)

stopCluster(cl)


#------------------------------------------------------------
# save the workspace to the file .RData in the cwd
save.image("smu_relation_filterTest.RData")




