library(snow)
library(qvalue)
library(mvnfast)		# fast generate multi variate normal
library(IHW)
library("tibble")       # data table


cl <- makeCluster(10, type = "MPI")		# start zcluster


clusterExport(cl, "qvalue")
clusterExport(cl, "tibble")
clusterExport(cl, "ihw")
clusterExport(cl,"rmvn")
clusterExport(cl,"adj_pvalues")
clusterExport(cl,"rejections")




rw_weight <- function(testStat, gamma, alpha, group=5L, tail=2L)
{
    m = length(testStat)
    groupSize <- m/group
    testGroup <- rep(1:group, each = groupSize)
    testMeans <- tapply(testStat, testGroup, mean)
    testSd <- tapply(testStat, testGroup, sd)
    pi_hat <- testMeans^2/(testMeans^2 + testSd^2 - 1)
    effect_hat <- testMeans/pi_hat
    effect_hat[pi_hat <= 1/groupSize] <- 0

    if(sum(effect_hat, na.rm=T) == 0) {normWeight.w <- rep(1,m)
    }
    else{
        theta = as.vector(effect_hat)
        c0 = qnorm(alpha/(tail*m), lower.tail = FALSE)^2/2

        findc = function(thet, alpha, c0, m)
        {
            cc = seq(.1, c0, .005)
            tot <- sapply(cc, function(c) return(sum((m/alpha)*pnorm((thet/2 + c/thet),
                                                                     lower.tail = FALSE))))
            cout = cc[min(abs(tot - m)) == abs(tot - m)]
            coutAdj = ifelse(length(cout) > 1, sample(cout, 1), cout)
            return(coutAdj)
        }

        c <- vapply(theta, findc, 1, alpha, c0, m)
        weight_k <- (m/alpha)*pnorm((effect_hat/2 + c/effect_hat), lower.tail = FALSE)
        weight_k_smooth <- (1-gamma)*weight_k + gamma*sum(weight_k)/group
        weight_tests <- rep(weight_k_smooth, each=groupSize)
        normWeight.w <- weight_tests/sum(weight_tests)*m
    }
    return(normWeight.w)
}

clusterExport(cl, "rw_weight")



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


# function to compute weight continuous case
#--------------------------------------
weight_continuous <- function(alpha, et, m, tail=2L, delInterval=.0001, prob)
{
    prob <- prob/sum(prob, na.rm = T)
    delta <- seq(0, 1, delInterval)
    findDelta <- function(delta)
    {
        weight <- tail*(m/alpha)*pnorm(et/2 + 1/et*log(delta/(alpha*prob)),
                                       lower.tail=FALSE)
        return(sum(weight, na.rm = TRUE))
    }
    weightSumVec <- vapply(delta, findDelta, 1)
    deltaOut <- delta[min(abs(weightSumVec - m)) == abs(weightSumVec - m)]
    deltaOut <- ifelse(length(deltaOut) > 1, .0001, deltaOut)
    weight.out <- tail*(m/alpha)*pnorm(et/2 + 1/et*log(deltaOut/(alpha*prob)),
                                       lower.tail=FALSE)
    sumWeight <- sum(weight.out, na.rm = TRUE)
    normWeight <- if(sumWeight == 0) {rep(1, m)} else {weight.out/sumWeight*m}
    return(normWeight)
}

clusterExport(cl, "weight_continuous")


runif_by_mean <- function(n, mean)
{
    sd = mean/2
    uni_rv <- mean + sd*scale(runif(n, 0, 1))
    return(as.vector(uni_rv))
}

clusterExport(cl, "runif_by_mean")



test_by_block <- function(r, eVec, groupSize, Sigma)
{
    eSub <- eVec[(groupSize*r + 1 - groupSize):(groupSize*r)]
    test <- as.vector(rmvn(1, eSub, Sigma, ncores = 50))
    return(test)
}



clusterExport(cl, "test_by_block")


fwerPowerFdrPower_by_effect <- function(filterEffect, m , null, testCorr=0,
                                        random=0, groupSize=100, alpha=.05, Sigma)
{
    testEffect = if(random == 0){filterEffect
    } else {(c(2, 3, 5, 1/2, 1/3, 1/5)*filterEffect)[random]}

    xf = runif_by_mean(n = m, mean = filterEffect)         # covariate
    xt = runif_by_mean(n = m, mean = testEffect)

    H = rbinom(m, 1 , 1-null)          	# alternative hypothesis true or false
    ef <- H*xf  				# filter effect vector (mixture of null and alt)
    et <- H*xt					# test effect vector (mixture of null and alt)

    mGrp = m/groupSize				# subgroup of tests.

    filter <- if(testCorr == 0) {rnorm(m, ef, 1)
    } else {as.vector(vapply(1:mGrp, test_by_block, 1, eVec=ef, groupSize=groupSize, Sigma = Sigma))}	# filter test stat

    test <- if(testCorr == 0) {rnorm(m, et, 1)
    } else {as.vector(vapply(1:mGrp, test_by_block, 1, eVec=et, groupSize=groupSize, Sigma=Sigma))}	# actual test stat

    pval = pnorm(test, lower.tail = FALSE)

    dat = cbind(test, pval, et, filter)
    OD = dat[order(dat[,4], decreasing=T), ]			# odered by covariate for full data set

    null_est = qvalue(pval, pi0.method="bootstrap")$pi0
    m0 = ceiling(m*null_est)
    m1 = m-m0

    model = lm(filter ~ test)

    test_effect = if(m1 == 0) {0
    } else {sort(test, decreasing = T)[1:m1]}		# two-tailed test

    et_cont = mean(test_effect, na.rm = T)
    ey_cont = model$coef[[1]] + model$coef[[2]]*et_cont

    ranksProb <-vapply(1:m, prob_rank_givenEffect, 1, et = ey_cont,
                       ey = ey_cont, nrep = 10000, m0 = m0, m1 = m1)

    ranksProb <- ranksProb/sum(ranksProb, na.ram = T)

    weight_pro = weight_continuous(alpha = alpha, et = et_cont, m = m, tail = 2,
                                   delInterval=.0001 , prob = ranksProb)

    # pro=proposed,bon=bonferroni,rdw=roeder and wasserman,IHW=independent Hyp Weight
    #----------------------------------------------------------------------------------------
    weight_rdw <- as.vector(rw_weight(testStat = dat[,1], gamma=.05, alpha=alpha, group=5, tail=1))		# roeder wasserman weight
    ihw_fwer <- ihw(dat[,2], dat[,4], alpha=alpha, adjustment_type = "bonferroni")		# IHW method for FWER
    ihw_fdr <-  ihw(dat[,2], dat[,4], alpha=alpha, adjustment_type = "BH")			# IHW method for FDR

    rej_pro <- OD[,2] <= alpha*weight_pro/m			# total rejections of all methods
    rej_bon <- dat[,2] <= alpha/m
    rej_rdw <- dat[,2] <= alpha*weight_rdw/m
    rej_ihwFwer <- adj_pvalues(ihw_fwer) <= alpha

    n_null <- max(1, sum(OD[,3] == 0, na.rm = T))
    n_alt <-  max(1, sum(OD[,3] != 0, na.rm = T))

    FWER_pro <- sum(rej_pro[OD[,3] == 0])/n_null			# FWER of proposed method
    FWER_bon <- sum(rej_bon[dat[,3] == 0])/n_null			# FWER of bonferroni method
    FWER_rdw <- sum(rej_rdw[dat[,3] == 0])/n_null			# FWER of Roeder Wasserman method
    FWER_ihw <- sum(rej_ihwFwer[dat[,3] == 0])/n_null			# FWER of IHW method

    POWER_pro <- sum(rej_pro[OD[,3] != 0])/n_alt		# power of proposed
    POWER_bon <- sum(rej_bon[dat[,3] != 0])/n_alt		# power of bonferroni
    POWER_rdw <- sum(rej_rdw[dat[,3] != 0])/n_alt		# power of Roeder Wasserman method
    POWER_ihw <- sum(rej_ihwFwer[dat[,3] !=0 ])/n_alt	# power of IHW method

    adjPval_pro <- p.adjust(OD[,2]/weight_pro, method="BH", n=m)	# adjusted pvalue to compute FDR
    adjPval_bon <- p.adjust(dat[,2], method="BH", n=m)
    adjPval_rdw <- p.adjust(dat[,2]/weight_rdw, method="BH", n=m)
    adjPval_ihw <- adj_pvalues(ihw_fdr)

    FDR_pro <- sum(adjPval_pro[OD[,3] == 0] <= alpha)/max(1, sum(adjPval_pro <= alpha))	# FDR of proposed
    FDR_bh  <- sum(adjPval_bon[dat[,3] == 0] <= alpha)/max(1, sum(adjPval_bon <= alpha))	# FDR of benjaminin and hochberg
    FDR_rdw <- sum(adjPval_rdw[dat[,3] == 0] <= alpha)/max(1, sum(adjPval_rdw <= alpha))	# FDR of wasserman
    FDR_ihw <- sum(adjPval_ihw[dat[,3] == 0] <= alpha)/max(1, rejections(ihw_fdr))		# FDR of IHW method

    FDR_POWER_pro <- sum(adjPval_pro[OD[,3] != 0] <= alpha)/n_alt	# FDR of proposed
    FDR_POWER_bh  <- sum(adjPval_bon[dat[,3] != 0] <= alpha)/n_alt	# FDR of benjaminin and hochberg
    FDR_POWER_rdw <- sum(adjPval_rdw[dat[,3] != 0] <= alpha)/n_alt	# FDR of wasserman
    FDR_POWER_ihw <- sum(adjPval_ihw[dat[,3] != 0] <= alpha)/n_alt		# FDR of IHW method

    return(c(FWER_pro, FWER_bon, FWER_rdw, FWER_ihw,
             POWER_pro, POWER_bon, POWER_rdw, POWER_ihw,
             FDR_pro, FDR_bh, FDR_rdw, FDR_ihw,
             FDR_POWER_pro, FDR_POWER_bh, FDR_POWER_rdw, FDR_POWER_ihw))
}


clusterExport(cl, "fwerPowerFdrPower_by_effect")

SigmaVal <- diag(100)
clusterExport(cl, "SigmaVal")

# function for simulations--------------------------
simu_fwerPowerFdrPower_by_effect <- function(s, filterEffect, m, null, testCorr, random, groupSize, alpha)
{
    fwerPowerFdrPower_mat = sapply(effecSize, fwerPowerFdrPower_by_effect, m = m,
                                   null = null, testCorr = testCorr, random = random,
                                   groupSize = 100, alpha = .05, Sigma = SigmaVal)
    return(fwerPowerFdrPower_mat)
}



#================parallel computing start=======================================

mVal = 10000
effecSize = c(.3, .5, .7, 1, 2, 3)

clusterExport(cl, "mVal")
clusterExport(cl, "effecSize")

simuVal = 1:3

FwerPowerFdrPower_simu_b1 <- parSapply(cl, simuVal, simu_fwerPowerFdrPower_by_effect, filterEffect=effecSize, m = mVal, null=.9, testCorr=0, random=1, groupSize = 100, alpha = .05, Sigma = SigmaVal)
FwerPowerFdrPower_simu_b2 <- parSapply(cl, simuVal, simu_fwerPowerFdrPower_by_effect, filterEffect=effecSize, m = mVal, null=.9, testCorr=0, random=2, groupSize = 100, alpha = .05, Sigma = SigmaVal)
FwerPowerFdrPower_simu_b3 <- parSapply(cl, simuVal, simu_fwerPowerFdrPower_by_effect, filterEffect=effecSize, m = mVal, null=.9, testCorr=0, random=3, groupSize = 100, alpha = .05, Sigma = SigmaVal)
FwerPowerFdrPower_simu_b4 <- parSapply(cl, simuVal, simu_fwerPowerFdrPower_by_effect, filterEffect=effecSize, m = mVal, null=.9, testCorr=0, random=4, groupSize = 100, alpha = .05, Sigma = SigmaVal)
FwerPowerFdrPower_simu_b5 <- parSapply(cl, simuVal, simu_fwerPowerFdrPower_by_effect, filterEffect=effecSize, m = mVal, null=.9, testCorr=0, random=5, groupSize = 100, alpha = .05, Sigma = SigmaVal)
FwerPowerFdrPower_simu_b6 <- parSapply(cl, simuVal, simu_fwerPowerFdrPower_by_effect, filterEffect=effecSize, m = mVal, null=.9, testCorr=0, random=6, groupSize = 100, alpha = .05, Sigma = SigmaVal)



stopCluster(cl)


#------------------------------------------------------------
# save the workspace to the file .RData in the cwd
save.image("simu_fwerPowerFdrPower_random.RData")



