library(snow)
library(mvnfast)		# fast generate multi variate normal
library("IHW")		# independent hypotheis weight
library(qvalue)


cl <- makeCluster(10, type = "MPI")		# start zcluster


clusterExport(cl, "rmvn")
clusterExport(cl, "ihw")
clusterExport(cl, "adj_pvalues")
clusterExport(cl, "rejections")
clusterExport(cl, "qvalue")



# define various functions-------------
# normWeight.w = normalized weight
#---------------------------------------------------------------------------------------------
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
        c0 = (qnorm(1-alpha/(tail*m)))^2/2

        findc = function(thet, alpha, c0, m)
        {
            cc = seq(.1, c0, .005)
            tot <- sapply(cc, function(c) return(sum((m/alpha)*pnorm((thet/2 + c/thet),
                                                                     lower.tail = FALSE))))
            cout = cc[min(abs(tot - m)) == abs(tot - m)]
            coutAdj = ifelse(length(cout) > 1, sample(cout, 1), cout)
            return(coutAdj)
        }

        c <- sapply(theta, findc, alpha, c0, m)
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


# function to compute weight binary case
#--------------------------------------
weight_binary <- function(alpha, et, m, m1, tail = 2L, delInterval = .0001, prob)
{
    prob <- prob/sum(prob)
    delta <- seq(0, 1, delInterval)
    findDelta <- function(delta)
    {
        weight <- tail*(m/alpha)*pnorm(et/2 + 1/et*log(delta*m/(alpha*m1*prob)),
                                       lower.tail = FALSE)
        return(sum(weight, nam.rm = TRUE))
    }
    weightSumVec <- sapply(delta, findDelta)
    deltaOut <- delta[min(abs(weightSumVec - m)) == abs(weightSumVec - m)]
    deltaOut <- ifelse(length(deltaOut) > 1, .0001, deltaOut)
    weight.out <- tail*(m/alpha)*pnorm(et/2 + 1/et*log(deltaOut*m/(alpha*m1*prob)),
                                       lower.tail = FALSE)
    sumWeight <- sum(weight.out, na.rm = TRUE)
    normWeight <- if(sumWeight == 0) {rep(1, m)} else {weight.out/sumWeight*m}
    return(normWeight)
}

clusterExport(cl, "weight_binary")



# function to compute weight continuous case
#--------------------------------------
weight_continuous <- function(alpha, et, m, tail=2L, delInterval=.0001, prob)
{
    prob <- prob/sum(prob)
    delta <- seq(0, 1, delInterval)
    findDelta <- function(delta)
    {
        weight <- tail*(m/alpha)*pnorm(et/2 + 1/et*log(delta/(alpha*prob)),
                                       lower.tail=FALSE)
        return(sum(weight, nam.rm = TRUE))
    }
    weightSumVec <- sapply(delta, findDelta)
    deltaOut <- delta[min(abs(weightSumVec - m)) == abs(weightSumVec - m)]
    deltaOut <- ifelse(length(deltaOut) > 1, .0001, deltaOut)
    weight.out <- tail*(m/alpha)*pnorm(et/2 + 1/et*log(deltaOut/(alpha*prob)),
                                       lower.tail=FALSE)
    sumWeight <- sum(weight.out, na.rm = TRUE)
    normWeight <- if(sumWeight == 0) {rep(1, m)} else {weight.out/sumWeight*m}
    return(normWeight)
}

clusterExport(cl, "weight_continuous")




# function to conduct FWER simulations
#----------------------------------------
fun_Fwer_Simu <- function(simu, corr)
{
    fun.Fwer <- function(corr=0, alpha=.05, m=10000)
    {
        effect <- rep(0, 2)							# effect size
        Sigma <- matrix(corr, 2, 2) + diag(2)*(1-corr)		# test correlation matrix
        bivar_test <-rmvn(m, effect, Sigma)
        pval.et <- 2*pnorm(abs(bivar_test[ , 1]), lower.tail = FALSE)		# actual test pvalues
        Data <- data.frame(test = bivar_test[ , 1], pval = pval.et, filter = bivar_test[ , 2] )	# data of effect,test, and pvlaues
        OD <- Data[order(Data$filter, decreasing=T),]			# ordered data ordered by filter tests

        null = qvalue(p = Data$pval, pi0.method="bootstrap")$pi0
        m0 <- ceiling(m*null)							# no. of null hypothesis
        m1 <- m-m0									# no. of alt hyp.

        model <- lm(Data$filter ~ Data$test)
        test_effect = sort(abs(Data$test), decreasing = T)[1:m1]		# two-tailed test
        et_bin = median(test_effect, na.rm = T)		# median from the summary for the binary case test effect
        et_cont = mean(test_effect, na.rm = T)
        ey_bin = model$coef[[1]] + model$coef[[2]]*et_bin
        ey_cont = model$coef[[1]] + model$coef[[2]]*et_cont

        ranks = 1:m
        ranksProb_bin = sapply(ranks, prob_rank_givenEffect, et=et_bin, ey=ey_bin,
                               nrep=10000, m0=m0, m1=m1)
        ranksProb_cont = sapply(ranks, prob_rank_givenEffect, et=et_cont, ey=ey_cont,
                                nrep=10000, m0=m0, m1=m1)

        W_bin = weight_binary(alpha=alpha, et=et_bin, m=m, m1=m1, tail=2,
                              delInterval=.0001, prob=ranksProb_bin)
        w_cont = weight_continuous(alpha=alpha, et=et_cont, m=m, tail=2,
                                   delInterval=.0001, prob=ranksProb_cont)

        # pro=proposed,bon=bonferroni,rdw=roeder and wasserman,IHW=independent Hyp Weight
        #----------------------------------------------------------------------------------------
        weight_rdw <- as.vector(rw_weight(testStat=Data$test, gamma=.05, alpha=alpha,
                                          group=5L, tail=2))		# roeder wasserman weight
        ihw_fwer <- ihw(Data$pval, Data$filter, alpha=alpha, adjustment_type = "bonferroni")		# IHW method for FWER
        ihw_fdr <-  ihw(Data$pval, Data$filter, alpha=alpha, adjustment_type = "BH")			# IHW method for FDR

        rej_pro_bin <- OD$pval <= alpha*W_bin/m			# total rejections of all methods
        rej_pro_cont <- OD$pval <= alpha*w_cont/m
        rej_bon <- OD$pval <= alpha/m
        rej_rdw <- OD$pval <= alpha*weight_rdw/m
        rej_ihwFwer <- adj_pvalues(ihw_fwer) <= alpha

        FWER_pro_bin <- sum(rej_pro_bin)/m			# FWER of proposed method
        FWER_pro_cont <- sum(rej_pro_cont)/m			# FWER of proposed method
        FWER_bon <- sum(rej_bon)/m			# FWER of bonferroni method
        FWER_rdw <- sum(rej_rdw)/m			# FWER of Roeder Wasserman method
        FWER_ihw <- sum(rej_ihwFwer)/m			# FWER of IHW method
        return(c(FWER_pro_bin,  FWER_pro_cont, FWER_bon, FWER_rdw, FWER_ihw))
    }
    rejections <- fun.Fwer(corr=corr, alpha=.05, m=10000)
    return(rejections)
}


sim = 1:1000

Fwer1 <- parSapply(cl, sim, fun_Fwer_Simu, corr = 0)
Fwer2 <- parSapply(cl, sim, fun_Fwer_Simu, corr = .3)
Fwer3 <- parSapply(cl, sim, fun_Fwer_Simu, corr = .5)
Fwer4 <- parSapply(cl, sim, fun_Fwer_Simu, corr = .7)
Fwer5 <- parSapply(cl, sim, fun_Fwer_Simu, corr = .9)


stopCluster(cl)


#------------------------------------------------------------
# save the workspace to the file .RData in the cwd
save.image("FWER.RData")



