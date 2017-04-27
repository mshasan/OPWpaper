library(snow)
library(mvnfast)		# fast generate multi variate normal
library("IHW")		# independent hypotheis weight
library(tibble)     # tibble data frame


cl <- makeCluster(15, type = "MPI")		# start zcluster

clusterExport(cl,"rmvn")
clusterExport(cl,"ihw")
clusterExport(cl,"adj_pvalues")
clusterExport(cl,"rejections")
clusterExport(cl,"tibble")


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



# function to generate multivariate tests statistics by block
# input: r=no. of test groups,eVec=effect vector,Sigma=corr matrix
# output: test = multivariate test statistics
#-------------------------------------------------------------------
test_by_block <- function(r, eVec, groupSize, Sigma)
{
    eSub <- eVec[(groupSize*r + 1 - groupSize):(groupSize*r)]
    test <- as.vector(rmvn(1, eSub, Sigma, ncores = 50))
    return(test)
}


clusterExport(cl, "test_by_block")




################################################################################################################
#----------------------2b-3:fwerPowerFdrPower_bin----------------------------
#function to compute Simulated FWER, POWER, and FDR by effect size
#-------------------------------------------------------------------------------------------------------
# Input:
#=============
# i = i-th effect
# simu = number of simulation
# null =  proportion of null hypothesis
# corr = test correlation
# random = 0, test effect = filter effect; 1 if test effect~N(filter effect, filter effect/2)
# alpha = FWER level
# filterEffectVec =  different effect size (10 values)
# datWeightByNull = weight matrix (10000 by 40) for various null = c(20,50,90,99)% computed before;
# 			m=10,000 and effect = filterEffectVec*4 = (10 effects)*(4 null scenerios) = 40

#output:
#=============
# FwerPowerFdr = A 16 by 10 matrix consits of simulated FWER(first 4 rows), POWER (2nd 4 rows),
# FDR (3rd 4 rows) and FDRPower(last 4 rows) for 10 different effect sizes
#-----------------------------------------------------------------------------------------------------

fwerPowerFdrPower_bin <- function(i, simu, null, corr=0, random, alpha, groupSize,
                                  filterEffectVec, datWeightByNull)
    {
    	W = datWeightByNull[ , i]							# weight vector for a specific effect size
    	m = length(W)								# test size
    	weight_pro <- if(sum(W)==0){rep(1, m)} else {W/sum(W)*m}			# normalizing proposed weight
    	ey <- filterEffectVec[i]							# effect size
    	m0 <- ceiling(m*null)							# no. of null hypothesis
    	m1 <- m - m0									# no. of alt hyp.
    	xf <- rep(ey, m)								# only alt. filter effect vector
    	xt <- if(random == 0){rep(ey, m)} else {rnorm(m, ey, ey/2)} 			# alt. test effect vector
    	Sigma <- matrix(corr, 100, 100) + diag(100)*(1 - corr)			# test correlation matrix


	# function to do replications
	#----------------------------------
	fwerPowerFdrPower_simu <- function(s) # input: s=replication; output: vector of FWER,POWER,FDR(4*4=16 obs.)
		{
    		H <- rbinom(m, 1, 1 - null)          	# alternative hypothesis true or false
    		ef <- H*xf  				# filter effect vector (mixture of null and alt)
    		et <- H*xt					# test effect vector (mixture of null and alt)
    		mGrp = m/groupSize				# subgroup of tests.

    		test <- if(corr == 0) {rnorm(m, et, 1)
    		    } else {as.vector(sapply(1:mGrp, test_by_block, eVec = et,
    		                             groupSize = groupSize, Sigma = Sigma))}	# filter test stat

    		filter <- if(corr == 0) {rnorm(m, ef, 1)
    		    } else {as.vector(sapply(1:mGrp, test_by_block, eVec = ef,
    		                             groupSize = groupSize, Sigma = Sigma))}	# actual test stat

    		pval <- pnorm(test, lower.tail = F)		# filter test pvalues

    		dat = tibble(test, pval, et, filter)
    		OD = dat[order(dat$filter, decreasing = T), ]

    		# pro=proposed,bon=bonferroni,rdw=roeder and wasserman,IHW=independent Hyp Weight
    		#----------------------------------------------------------------------------------------
    		weight_rdw <- as.vector(rw_weight(testStat = dat$test, gamma=.05, alpha=alpha, group=5, tail=1))		# roeder wasserman weight
    		ihw_fwer <- ihw(dat$pval, dat$filter, alpha = alpha, adjustment_type = "bonferroni")		# IHW method for FWER
    		ihw_fdr <-  ihw(dat$pval, dat$filter, alpha = alpha, adjustment_type = "BH")			# IHW method for FDR

    		rej_pro <- OD$pval <= alpha*weight_pro/m			# total rejections of all methods
    		rej_bon <- dat$pval <= alpha/m
    		rej_rdw <- dat$pval <= alpha*weight_rdw/m
    		rej_ihwFwer <- adj_pvalues(ihw_fwer) <= alpha

    		n_null <- max(1, sum(dat$et == 0, na.rm = T))
    		n_alt <-  max(1, sum(dat$et != 0, na.rm = T))

    		FWER_pro <- sum(rej_pro[OD$et == 0])			# FWER of proposed method
    		FWER_bon <- sum(rej_bon[dat$et == 0])			# FWER of bonferroni method
    		FWER_rdw <- sum(rej_rdw[dat$et == 0])			# FWER of Roeder Wasserman method
    		FWER_ihw <- sum(rej_ihwFwer[dat$et == 0])			# FWER of IHW method

    		POWER_pro <- sum(rej_pro[OD$et != 0])/n_alt		# power of proposed
    		POWER_bon <- sum(rej_bon[dat$et != 0])/n_alt		# power of bonferroni
    		POWER_rdw <- sum(rej_rdw[dat$et != 0])/n_alt		# power of Roeder Wasserman method
    		POWER_ihw <- sum(rej_ihwFwer[dat$et != 0])/n_alt	# power of IHW method

    		adjPval_pro <- p.adjust(OD$pval/weight_pro, method="BH")	# adjusted pvalue to compute FDR
    		adjPval_bon <- p.adjust(dat$pval, method="BH")
    		adjPval_rdw <- p.adjust(dat$pval/weight_rdw, method="BH")
    		adjPval_ihw <- adj_pvalues(ihw_fdr)

    		FDR_pro <- sum(adjPval_pro[OD$et == 0] <= alpha)/max(1,sum(adjPval_pro <= alpha))	# FDR of proposed
    		FDR_bh  <- sum(adjPval_bon[dat$et == 0] <= alpha)/max(1,sum(adjPval_bon <= alpha))	# FDR of benjaminin and hochberg
    		FDR_rdw <- sum(adjPval_rdw[dat$et == 0] <= alpha)/max(1,sum(adjPval_rdw <= alpha))	# FDR of wasserman
    		FDR_ihw <- sum(adjPval_ihw[dat$et == 0] <= alpha)/max(1,rejections(ihw_fdr))		# FDR of IHW method

    		FDR_POWER_pro <- sum(adjPval_pro[OD$et!=0] <= alpha)/n_alt	# FDR of proposed
    		FDR_POWER_bh  <- sum(adjPval_bon[OD$et!=0] <= alpha)/n_alt	# FDR of benjaminin and hochberg
    		FDR_POWER_rdw <- sum(adjPval_rdw[OD$et!=0] <= alpha)/n_alt	# FDR of wasserman
    		FDR_POWER_ihw <- sum(adjPval_ihw[OD$et!=0] <= alpha)/n_alt		# FDR of IHW method

    		return(c(FWER_pro, FWER_bon, FWER_rdw, FWER_ihw,
    		         POWER_pro, POWER_bon, POWER_rdw, POWER_ihw,
    			     FDR_pro, FDR_bh, FDR_rdw, FDR_ihw, FDR_POWER_pro,
    			     FDR_POWER_bh, FDR_POWER_rdw, FDR_POWER_ihw))
    	}

	fwerPowerFdrPower_bysimu <- sapply(1:simu, fwerPowerFdrPower_simu)		# result of each replication
	fwerPowerFdrPower <- apply(fwerPowerFdrPower_bysimu, 1, mean, na.rm=T)			# mean over replications

	return(fwerPowerFdrPower)

	}


filterEffectVec <- c(seq(0, 1, .2), 2, 3, 5, 8)
datWeight <- read.csv("Binary_WeightByEffect_m10000.csv", h = T)
simuVal = 1

clusterExport(cl, "datWeight")
clusterExport(cl, "filterEffectVec")
clusterExport(cl, "simuVal")


FwerPowerFdrPower2e1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_bin, simu=simuVal, null=.5, corr= 0, random=0, alpha=.05, groupSize=100, filterEffectVec=filterEffectVec, datWeightByNull=datWeight[,11:20])
FwerPowerFdrPower2e2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_bin, simu=simuVal, null=.5, corr= 0, random=1, alpha=.05, groupSize=100, filterEffectVec=filterEffectVec, datWeightByNull=datWeight[,11:20])
FwerPowerFdrPower2f1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_bin, simu=simuVal, null=.5, corr=.3, random=0, alpha=.05, groupSize=100, filterEffectVec=filterEffectVec, datWeightByNull=datWeight[,11:20])
FwerPowerFdrPower2f2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_bin, simu=simuVal, null=.5, corr=.3, random=1, alpha=.05, groupSize=100, filterEffectVec=filterEffectVec, datWeightByNull=datWeight[,11:20])
FwerPowerFdrPower2g1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_bin, simu=simuVal, null=.5, corr=.5, random=0, alpha=.05, groupSize=100, filterEffectVec=filterEffectVec, datWeightByNull=datWeight[,11:20])
FwerPowerFdrPower2g2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_bin, simu=simuVal, null=.5, corr=.5, random=1, alpha=.05, groupSize=100, filterEffectVec=filterEffectVec, datWeightByNull=datWeight[,11:20])
FwerPowerFdrPower2h1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_bin, simu=simuVal, null=.5, corr=.7, random=0, alpha=.05, groupSize=100, filterEffectVec=filterEffectVec, datWeightByNull=datWeight[,11:20])
FwerPowerFdrPower2h2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_bin, simu=simuVal, null=.5, corr=.7, random=1, alpha=.05, groupSize=100, filterEffectVec=filterEffectVec, datWeightByNull=datWeight[,11:20])
FwerPowerFdrPower2i1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_bin, simu=simuVal, null=.5, corr=.9, random=0, alpha=.05, groupSize=100, filterEffectVec=filterEffectVec, datWeightByNull=datWeight[,11:20])
FwerPowerFdrPower2i2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_bin, simu=simuVal, null=.5, corr=.9, random=1, alpha=.05, groupSize=100, filterEffectVec=filterEffectVec, datWeightByNull=datWeight[,11:20])



FwerPowerFdrPower4e1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_bin, simu=simuVal, null=.9, corr= 0, random=0, alpha=.05, groupSize=100, filterEffectVec=filterEffectVec, datWeightByNull=datWeight[,31:40])
FwerPowerFdrPower4e2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_bin, simu=simuVal, null=.9, corr= 0, random=1, alpha=.05, groupSize=100, filterEffectVec=filterEffectVec, datWeightByNull=datWeight[,31:40])
FwerPowerFdrPower4f1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_bin, simu=simuVal, null=.9, corr=.3, random=0, alpha=.05, groupSize=100, filterEffectVec=filterEffectVec, datWeightByNull=datWeight[,31:40])
FwerPowerFdrPower4f2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_bin, simu=simuVal, null=.9, corr=.3, random=1, alpha=.05, groupSize=100, filterEffectVec=filterEffectVec, datWeightByNull=datWeight[,31:40])
FwerPowerFdrPower4g1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_bin, simu=simuVal, null=.9, corr=.5, random=0, alpha=.05, groupSize=100, filterEffectVec=filterEffectVec, datWeightByNull=datWeight[,31:40])
FwerPowerFdrPower4g2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_bin, simu=simuVal, null=.9, corr=.5, random=1, alpha=.05, groupSize=100, filterEffectVec=filterEffectVec, datWeightByNull=datWeight[,31:40])
FwerPowerFdrPower4h1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_bin, simu=simuVal, null=.9, corr=.7, random=0, alpha=.05, groupSize=100, filterEffectVec=filterEffectVec, datWeightByNull=datWeight[,31:40])
FwerPowerFdrPower4h2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_bin, simu=simuVal, null=.9, corr=.7, random=1, alpha=.05, groupSize=100, filterEffectVec=filterEffectVec, datWeightByNull=datWeight[,31:40])
FwerPowerFdrPower4i1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_bin, simu=simuVal, null=.9, corr=.9, random=0, alpha=.05, groupSize=100, filterEffectVec=filterEffectVec, datWeightByNull=datWeight[,31:40])
FwerPowerFdrPower4i2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_bin, simu=simuVal, null=.9, corr=.9, random=1, alpha=.05, groupSize=100, filterEffectVec=filterEffectVec, datWeightByNull=datWeight[,31:40])



FwerPowerFdrPower5e1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_bin, simu=simuVal, null=.99, corr= 0, random=0, alpha=.05, groupSize=100, filterEffectVec=filterEffectVec, datWeightByNull=datWeight[,41:50])
FwerPowerFdrPower5e2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_bin, simu=simuVal, null=.99, corr= 0, random=1, alpha=.05, groupSize=100, filterEffectVec=filterEffectVec, datWeightByNull=datWeight[,41:50])
FwerPowerFdrPower5f1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_bin, simu=simuVal, null=.99, corr=.3, random=0, alpha=.05, groupSize=100, filterEffectVec=filterEffectVec, datWeightByNull=datWeight[,41:50])
FwerPowerFdrPower5f2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_bin, simu=simuVal, null=.99, corr=.3, random=1, alpha=.05, groupSize=100, filterEffectVec=filterEffectVec, datWeightByNull=datWeight[,41:50])
FwerPowerFdrPower5g1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_bin, simu=simuVal, null=.99, corr=.5, random=0, alpha=.05, groupSize=100, filterEffectVec=filterEffectVec, datWeightByNull=datWeight[,41:50])
FwerPowerFdrPower5g2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_bin, simu=simuVal, null=.99, corr=.5, random=1, alpha=.05, groupSize=100, filterEffectVec=filterEffectVec, datWeightByNull=datWeight[,41:50])
FwerPowerFdrPower5h1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_bin, simu=simuVal, null=.99, corr=.7, random=0, alpha=.05, groupSize=100, filterEffectVec=filterEffectVec, datWeightByNull=datWeight[,41:50])
FwerPowerFdrPower5h2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_bin, simu=simuVal, null=.99, corr=.7, random=1, alpha=.05, groupSize=100, filterEffectVec=filterEffectVec, datWeightByNull=datWeight[,41:50])
FwerPowerFdrPower5i1 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_bin, simu=simuVal, null=.99, corr=.9, random=0, alpha=.05, groupSize=100, filterEffectVec=filterEffectVec, datWeightByNull=datWeight[,41:50])
FwerPowerFdrPower5i2 <- parSapply(cl, 1:length(filterEffectVec), fwerPowerFdrPower_bin, simu=simuVal, null=.99, corr=.9, random=1, alpha=.05, groupSize=100, filterEffectVec=filterEffectVec, datWeightByNull=datWeight[,41:50])



stopCluster(cl)



#-------------------------------------------------------------
# save the workspace to the file .RData in the cwd
save.image("simu_fwerPowerFdrPower_bin.RData")








