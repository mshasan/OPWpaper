getwd()
setwd("U:/Documents/My Research (UGA)/Multiple Hypoetheses/Article-1")
#setwd("C:/Users/Apu-Jerrica/Google Drive/UGA/My R Packages/OPWeight")

library(OPWeight)       # library for the proposed method
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(gridExtra)      # for multiplots in the same page
library(xtable)
library(reshape2)       # library for the melt function
library(cowplot)        # plot_grid function
library(mvnfast)		# fast generate multi variate normal
library(qvalue)         # qvalue
library(IHW)
library("tibble")       # data table




################################################################################
# function to compare rank probbailities from three approahes: 1) sikmulation,
# 2) exact fromual, and 3) normal approximation

# inpout:----------------
# ey = effect size
# e.one = vary one test effect across all tests
# m0 = number of null tests
# m1 = number of alternative tests
# effectType = type of effect size c("binary","continuous")
# sampleSize = total number of sample generated (use sample size at least 100,000)

# internal parameters:-----
# m = total number of tests
# ranks = sequesnce of ranks or index numbers
# rank = rank of the tests per sample
# probbailities from various approaches (see corresponding help files)

# output:---------------
# Data = proabbilities
# or plots of the probabilities of ranks vs. ranks of the tests

#===============================================================================
ranksProb_compare <- function(ey, e.one, m0, m1, sampleSize,
                              effectType = c("binary", "continuous"))
{
    m = m0 +m1
    ranks <- 1:m
    et = e.one

    # simulation approach=======================================================
    # compute rank of the tests
    rank <- sapply(1:sampleSize, prob_rank_givenEffect_simu, ey = ey, e.one = e.one,
                                    m0=m0, m1=m1, effectType = effectType)

    # rank may generate missing valuse because of the large effcet size,
    # therefore, to make a matplot equal vectors size are needed. This procedure
    # will replace the missing value to make equal sized vector
    # probability of rank of a null test
    prob0 <- rep(NA, m)
    prob0_x <- tapply(rank[1,], rank[1,], length)/sampleSize
    prob0[as.numeric(names(prob0_x))] <- as.vector(prob0_x)

    # probability of rank of an alternative test
    prob1 <- rep(NA, m)
    prob1_x <- tapply(rank[2,], rank[2,], length)/sampleSize
    prob1[as.numeric(names(prob1_x))] <- as.vector(prob1_x)

    # exact approach============================================================
    # do not compute probability for a large number of tests
    prob0_exact <- sapply(ranks, prob_rank_givenEffect_exact, et=0, ey=ey,
                          nrep=10000, m0=m0, m1=m1, effectType = effectType)
    prob1_exact <- sapply(ranks, prob_rank_givenEffect_exact, et=et, ey=ey,
                          nrep=10000,m0=m0, m1=m1, effectType = effectType)

    # normal approximation approcah=============================================
    prob0_approx <- sapply(ranks, prob_rank_givenEffect_approx, et=0, ey=ey,
                           nrep=10000, m0=m0, m1=m1, effectType = effectType)
    prob1_approx <- sapply(ranks, prob_rank_givenEffect_approx, et=et, ey=ey,
                           nrep=10000, m0=m0, m1=m1, effectType = effectType)

    # nice plots================================================================
    Data <- data.frame(ranks, prob0, prob1, prob0_exact, prob1_exact,
                       prob0_approx, prob1_approx)

    return(Data)
}

# applying function------------
probData <- ranksProb_compare(ey = 1, e.one = 1, m0=90, m1=10,
                  sampleSize = sampleSize, effectType = "binary")


# plots------------
colnames(probData) <- c("ranks", "SH0","SH1","EH0","EH1","AH0","AH1")
dat <- melt(probData, id.var = "ranks")
ey_plot <- if(effectType == "binary"){ey
} else {if(ey == 0){0} else {paste0(ey-1, "-", ey)}}
ggplot(dat, aes(x = ranks, y = value, group = variable, colour = variable)) +
    geom_line(aes(color = variable), size = 1) +
    labs(x = "ranks", y = "p(rank | effect)") +
    guides(fill = guide_legend(title = NULL)) +
    theme(legend.title = element_blank()) +
    annotate("text", x = 50, y = .011,
             label = paste0("P(rank | effect = ", ey_plot,")"))





#================================start of parrallelcomputing====================
# null test size = c(20, 50, 75, 90, 99)
sampleSize = 1000 #( use atleast 1,000,000)
effect <- c(0, 1, 2)

# continuous cases==============================================================
prob_20_0_cont = parLapply(cl, effect, ranksProb_compare, e.one = 0, m0=20, m1=80, sampleSize = sampleSize, effectType = "binary")
prob_20_1_cont = parLapply(cl, effect, ranksProb_compare, e.one = 1, m0=20, m1=80, sampleSize = sampleSize, effectType = "binary")
prob_20_2_cont = parLapply(cl, effect, ranksProb_compare, e.one = 2, m0=20, m1=80, sampleSize = sampleSize, effectType = "binary")

prob_50_0_cont = parLapply(cl, effect, ranksProb_compare, e.one = 0, m0=50, m1=50, sampleSize = sampleSize, effectType = "binary")
prob_50_1_cont = parLapply(cl, effect, ranksProb_compare, e.one = 1, m0=50, m1=50, sampleSize = sampleSize, effectType = "binary")
prob_50_2_cont = parLapply(cl, effect, ranksProb_compare, e.one = 2, m0=50, m1=50, sampleSize = sampleSize, effectType = "binary")

prob_75_0_cont = parLapply(cl, effect, ranksProb_compare, e.one = 0, m0=75, m1=25, sampleSize = sampleSize, effectType = "binary")
prob_75_1_cont = parLapply(cl, effect, ranksProb_compare, e.one = 1, m0=75, m1=25, sampleSize = sampleSize, effectType = "binary")
prob_75_2_cont = parLapply(cl, effect, ranksProb_compare, e.one = 2, m0=75, m1=25, sampleSize = sampleSize, effectType = "binary")

prob_90_0_cont = parLapply(cl, effect, ranksProb_compare, e.one = 0, m0=90, m1=10, sampleSize = sampleSize, effectType = "binary")
prob_90_1_cont = parLapply(cl, effect, ranksProb_compare, e.one = 1, m0=90, m1=10, sampleSize = sampleSize, effectType = "binary")
prob_90_2_cont = parLapply(cl, effect, ranksProb_compare, e.one = 2, m0=90, m1=10, sampleSize = sampleSize, effectType = "binary")

prob_99_0_cont = parLapply(cl, effect, ranksProb_compare, e.one = 0, m0=99, m1=1, sampleSize = sampleSize, effectType = "binary")
prob_99_1_cont = parLapply(cl, effect, ranksProb_compare, e.one = 1, m0=99, m1=1, sampleSize = sampleSize, effectType = "binary")
prob_99_2_cont = parLapply(cl, effect, ranksProb_compare, e.one = 2, m0=99, m1=1, sampleSize = sampleSize, effectType = "binary")



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





#===========================end of parallel computing===========================

# this code is to load saved workspace from parallel computing
load("U:/Documents/My Research (UGA)/Multiple Hypoetheses/Article-1/simu_prob_rank_givenEffect.RDATA")



#=======================probability plot for the main article===================
# probability plols for the main article----------------------
#prob_50_0_cont           # to get (ey=0, and e.one = 0)
#prob_50_1_cont           # to get (ey~U(0, 1) and ey~U(1, 2) for e.one = 1)
#prob_50_2_cont           # ey~U(0, 1) and e.one = 2

dat00 <- prob_50_0_cont[[1]]
colnames(dat00) <- c("ranks", "SH0","SH1","EH0","EH1","AH0","AH1")
dat00 <- melt(dat00, id.var="ranks")

dat01 <- prob_50_1_cont[[2]]
colnames(dat01) <- c("ranks", "SH0","SH1","EH0","EH1","AH0","AH1")
dat01 <- melt(dat01, id.var="ranks")

dat12 <- prob_50_1_cont[[3]]
colnames(dat12) <- c("ranks", "SH0","SH1","EH0","EH1","AH0","AH1")
dat12 <- melt(dat12, id.var="ranks")

dat02 <- prob_50_2_cont[[2]]
colnames(dat02) <- c("ranks", "SH0","SH1","EH0","EH1","AH0","AH1")
dat02 <- melt(dat02,id.var="ranks")


p00 <- ggplot(dat00, aes(x = ranks, y = value, group = variable,
                       colour = variable)) +
          geom_line(aes(linetype = variable), size = 1) +
          labs(x = "ranks", y = "p(rank | effect)", title = "(a) ey = 0, e.one = 0") +
          theme(legend.position="none") +
          annotate("text", x=50, y=.011, label="P(rank | effect = 0)")

p01 <- ggplot(dat01, aes(x = ranks, y = value, group = variable,
                       colour = variable)) +
          geom_line(aes(linetype = variable), size = 1) +
          labs(x = "ranks", y = "p(rank | effect)", title = "(b) ey ~ U(0, 1), e.one = 1") +
          theme(legend.position="none") +
          annotate("text", x = c(50, 50), y = c(.005, .018),
                   label = c("<---P(rank | effect = 0)","P(rank | effect = e.one)--->"))

p12 <- ggplot(dat12, aes(x = ranks, y = value, group = variable,
                       colour = variable)) +
          geom_line(aes(linetype = variable), size = 1) +
          labs(x = "ranks", y = "p(rank | effect)", title = "(c) ey ~ U(1, 2), e.one = 1") +
          theme(legend.position="none") +
          annotate("text", x = c(60, 40), y=c(.001, .02),
                  label = c("<---P(rank | effect = 0)","P(rank | effect = e.one)--->"))

p02 <- ggplot(dat02, aes(x = ranks, y = value, group = variable,
                       colour = variable)) +
          geom_line(aes(linetype = variable), size = 1) +
          labs(x = "ranks", y = "p(rank | effect)", title = "(d) ey ~ U(0, 1), e.one = 2") +
          theme(legend.title = element_blank()) +
          annotate("text", x = c(50, 50), y=c(.04, .18),
                  label = c("P(rank | effect = 0)","<---P(rank | effect = e.one)"))


# extract the legend from one of the plots
legend <- get_legend(p02 + theme(legend.direction="horizontal",
                                 legend.position="bottom"))
# arrange the plots
grid.arrange(p00, p01, p12, p02, nrow=2,
               top= "Continuous: m0 = 50, m1 = 50")

#==========================end==================================================





#=======================probability plot for the supplementry materials=========

# Function to plots from the parallel computing outputs=========================
# probability plols for the supplementry materials
prob_plots <- function(m0, effectType = c("binary", "continuous"))
{
    effType <- ifelse(effectType == "binary", "_bin", "_cont")
    prow <- list()
    for(j in 1:3)
    {
        e.one <- j-1
        # make text into varriable name
        probData_com <- eval(parse(text=paste0("prob_",m0,"_",e.one, effType)))
        g <- list()
        for(i in 1:3)
            {
            probData <- probData_com[[i]]
            colnames(probData) <- c("ranks", "SH0","SH1","EH0","EH1","AH0","AH1")
            dat <- melt(probData, id.var = "ranks")

            ey_val <- ifelse(effectType == "binary", i-1,
                             ifelse(i==1, 0, paste0("U(",i-2,", ", i-1,")")))
            eySm <- ifelse(effectType == "binary", "ey = ", "ey ~ ")

            g[[i]] = ggplot(dat, aes(x = ranks, y = value, group = variable,
                                     colour = variable)) +
                geom_line(aes(linetype = variable), size = 1) +
                labs(x = "ranks", y = "p(rank | effect)",
                     title = paste0(eySm, ey_val, ", e.one = ", e.one)) +
                theme(legend.title = element_blank())
           }

        # extract the legend from one of the plots
        legend <- get_legend(g[[1]] + theme(legend.direction="horizontal",
                                            legend.position="bottom"))

        # arrange the three plots in a single row
        prow[[j]] <- plot_grid(g[[1]] + theme(legend.position="none"),
                          g[[2]] + theme(legend.position="none"),
                          g[[3]] + theme(legend.position="none"),
                          align = 'hv', nrow = 1)
    }

    effType2 <- ifelse(effectType == "binary", "Binary: m0 = ", "Continuous: m0 = ")
    plots = grid.arrange(prow[[1]], prow[[2]], prow[[3]], legend, ncol=1,
            heights=c(7,7,7,1), top= paste(effType2, m0, ", m1 = ", 100-m0))
    return(list(plots))
}

# applying plot function========================================================
nullSize <- c(20, 50, 75, 90, 99)
lapply(nullSize, prob_plots, effectType = "continuous")
lapply(nullSize, prob_plots, effectType = "binary")

#================================end============================================





#==================== Example prob vs. weight===================================
ranks = 1:10000
nrep = 10000  # orginally I used nrep = 1,000,000
prob <- sapply(ranks, prob_rank_givenEffect, et=1, ey=1, nrep = nrep , m0=5000, m1=5000)
#prob <- read.csv("dataBin_ProbRanksByEffect_m10000.csv",h=T)[,6]

weight <- weight_continuous(alpha=.05, et=1, m=10000, tail=1, delInterval=.0001,
                            prob=prob)

dat2 = data.frame(ranks, prob, weight)
p_prob = ggplot(dat2, aes(x = ranks, y = prob)) +  geom_line() +
                labs(y = "p(rank | effect)")
p_weight = ggplot(dat2, aes(x = ranks, y = weight)) + geom_line()
grid.arrange(p_prob, p_weight, ncol = 2,
    top=paste0("Continuous: E(", expression(epsilon),") = 1, m = 10000, m0 = 5000, m1 = 5000"))




# ##############################################################################
#----------------------fun.Fwer----------------------------
#function to compute Simulated FWER, POWER, and FDR by effect size
#-------------------------------------------------------------------------------
# Input:
#=============
# i = i-th effect
# simu = number of simulation
# null =  proportion of null hypothesis
# corr = test correlation
# random = 0, test effect = filter effect; 1 if test effect~N(filter effect, filter effect/2)
# alpha = FWER level
# effectVec =  different effect size (10 values)
# datWeightByNull = weight matrix (10000 by 40) for various null=c(20,50,90,99)% computed before;
# 			m=10,000 and effect = effectVec*4 = (10 effects)*(4 null scenerios)=40

#output:
#=============
# FwerPowerFdr = A 16 by 10 matrix consits of simulated FWER(first 4 rows), POWER (2nd 4 rows),
# FDR (3rd 4 rows) and FDRPower(last 4 rows) for 10 different effect sizes
#-----------------------------------------------------------------------------------------------------



# function to conduct FWER simulations
#----------------------------------------
simu_TypeI_error <- function(s)
{
  TypeI_error <- function(alpha, m, nrep)
  {
    # create data set--------------
    filter = rnorm(m)
    test = rnorm(m)
    pval = 2*pnorm(abs(test), lower.tail = FALSE)
    dat = tibble(test, pval, filter)
    
    OD = dat[order(dat$filter, decreasing=T), ]			# odered by covariate for full data set
    odered.pvalue = OD$pval
    
    null = qvalue(pval)$pi0
    m0 = ceiling(m*null) 
    m1 = m-m0
    
    model = lm(filter ~ test)
    test_effect = sort(abs(test), decreasing = T)[1:m1]		# two-tailed test
    et_bin = median(test_effect, na.rm = T)		# median from the summary for the binary case test effect
    et_cont = mean(test_effect, na.rm = T)
    ey_bin = model$coef[[1]] + model$coef[[2]]*et_bin
    ey_cont = model$coef[[1]] + model$coef[[2]]*et_cont
    
    ranks = 1:m
    prob_bin <-sapply(ranks, prob_rank_givenEffect, et = abs(ey_bin), 
                             ey =abs(ey_bin), nrep = nrep, m0 = m0, m1 = m1)
    prob_cont <-sapply(ranks, prob_rank_givenEffect, et = abs(ey_cont),
                              ey = abs(ey_cont), nrep = nrep, m0 = m0, m1 = m1)
    
    w_bin <- weight_binary(alpha = alpha, et = et_bin, m = m, m1 = m1, tail = 2, 
                           delInterval = .0001, prob = prob_bin) 
    w_cont = weight_continuous(alpha= alpha, et = et_cont, m = m, tail = 2, 
                               delInterval=.0001 , prob = prob_cont)
    
    # complete pvalues
    #-----------------------
    bon = sum(odered.pvalue <= alpha/m, na.rm = T)
    pro_bin = sum(odered.pvalue <= alpha*w_bin/m, na.rm = T)
    pro_cont = sum(odered.pvalue <= alpha*w_cont/m, na.rm = T)

    return(c(bon, pro_bin, pro_cont))
  }

  alphaVal = seq(.01, .1, .02)
  TypeI_error_mat = sapply(alphaVal, TypeI_error, m=100, nrep = 10000)

  return(TypeI_error_mat)
}

simuVal = 1:2
simu_TypeI_error_mat = sapply(simuVal, simu_TypeI_error)


typeIerror_mat = matrix(apply(simu_TypeI_error_mat, 1, mean), nrow=5, byrow = T)

alphaVal = seq(.01, .1, .02)
datError <- data.frame(alphaVal, typeIerror_mat)
colnames(datError) <- c("alpha","BON","PRO_bin","PRO_cont")
datError2 <- melt(datError, id.var="alpha")

ggplot(datError2, aes(x = alpha, y = value, col=variable)) +
    geom_line(size=1.2) +
    geom_abline(linetype="dashed") + 
    xlab(expression(bold(paste("Nominal ",alpha)))) + 
    ylab("FWER")+
    scale_x_continuous(limits= c(0.01,0.1), breaks=seq(0.01,0.09,length=5)) +
    #ylim(0,0.9) +
    theme(legend.title = element_blank())+
    theme(axis.title = element_text(face="bold") )


#================end of FWER====================================================






################################################################################
#----------------------fun.FwerPowerFdrPower----------------------------
#function to compute Simulated FWER, POWER, and FDR by effect size
#-------------------------------------------------------------------------------
# Input:
#=============
# i = i-th effect
# simu = number of simulation
# null =  proportion of null hypothesis
# corr = test correlation
# random = 0, test effect = filter effect; 1 if test effect~N(filter effect, filter effect/2)
# alpha = FWER level
# effectVec =  different effect size (10 values)
# datWeightByNull = weight matrix (10000 by 40) for various null=c(20,50,90,99)% computed before;
# 			m=10,000 and effect = effectVec*4 = (10 effects)*(4 null scenerios)=40

#output:
#=============
# FwerPowerFdr = A 16 by 10 matrix consits of simulated FWER(first 4 rows), POWER (2nd 4 rows),
# FDR (3rd 4 rows) and FDRPower(last 4 rows) for 10 different effect sizes
#-----------------------------------------------------------------------------------------------------

fun.PowerFdrPower <- function(i,simu,null,corr=0,random,alpha,effectVec,datWeightByNull)
{
    W = datWeightByNull[,i]							# weight vector for a specific effect size
    m = length(W)								# test size
    weight_pro <- if(sum(W)==0){rep(1,m)} else {W/sum(W)*m}			# normalizing proposed weight
    ey <- effectVec[i]							# effect size
    m0 <- ceiling(m*null)							# no. of null hypothesis
    m1 <- m-m0									# no. of alt hyp.
    xf <- rep(ey,m)								# only alt. filter effect vector
    xt <- if(random==0){rep(ey,m)} else {rnorm(m,ey,ey/2)} 			# alt. test effect vector
    #xt <- ifelse(random==0,1,4)						# to see the effect of false weight
    Sigma <- matrix(corr, 100, 100) + diag(100)*(1-corr)			# test correlation matrix


    # function to do replications
    #----------------------------------
    fun.FwerPowerFdrPower.simu <- function(s) # input: s=replication; output: vector of FWER,POWER,FDR(4*4=16 obs.)
    {
        H <- rbinom(m,1,1-null)          	# alternative hypothesis true or false
        ef <- H*xf  				# filter effect vector (mixture of null and alt)
        et <- H*xt					# test effect vector (mixture of null and alt)
        mGrp = m/100				# subgroup of tests.

        # function to generate test
        # input: r=no. of test groups,eVec=effect vector,Sigma=corr matrix
        # output: test = multivariate test statistics
        #-------------------------------------------------------------------
        fun.test <- function(r,eVec,Sigma)
        {
            eSub <- eVec[(100*r+1-100):(100*r)]
            test <- as.vector(rmvn(1,eSub,Sigma))
            return(test)
        }
        test.ef <- if(corr==0) {rnorm(m,ef,1)} else {as.vector(sapply(1:mGrp,fun.test,eVec=ef,Sigma=Sigma))}	# filter test stat
        test.et <- if(corr==0) {rnorm(m,et,1)} else {as.vector(sapply(1:mGrp,fun.test,eVec=et,Sigma=Sigma))}	# actual test stat
        pval.ef <- 1-pnorm(test.ef)		# filter test pvalues
        pval.et <- 1-pnorm(test.et)		# actual test pvalues
        Data <- data.frame(ef,tf=test.ef,pf=pval.ef,et,tt=test.et,pt=pval.et)	# data of effect,test, and pvlaues
        OD <- Data[order(Data$tf,decreasing=T),]			# ordered data ordered by filter tests

        # pro=proposed,bon=bonferroni,rdw=roeder and wasserman,IHW=independent Hyp Weight
        #----------------------------------------------------------------------------------------
        weight_rdw <- as.vector(RoederWasermanWeight(OD$tt,m=m,gamma=.05,alpha=alpha,rk=1000))		# roeder wasserman weight
        ihw_fwer <- ihw(OD$pt,OD$tf,alpha=alpha,adjustment_type = "bonferroni")		# IHW method for FWER
        ihw_fdr <-  ihw(OD$pt,OD$tf,alpha=alpha,adjustment_type = "BH")			# IHW method for FDR

        rej_pro <- OD$pt <= alpha*weight_pro/m			# total rejections of all methods
        rej_bon <- OD$pt <= alpha/m
        rej_rdw <- OD$pt <= alpha*weight_rdw/m
        rej_ihwFwer <- adj_pvalues(ihw_fwer) <= alpha

        FWER_pro <- sum(rej_pro[OD$et==0])/sum(OD$et==0)			# FWER of proposed method
        FWER_bon <- sum(rej_bon[OD$et==0])/sum(OD$et==0)			# FWER of bonferroni method
        FWER_rdw <- sum(rej_rdw[OD$et==0])/sum(OD$et==0)			# FWER of Roeder Wasserman method
        FWER_ihw <- sum(rej_ihwFwer[OD$et==0])/sum(OD$et==0)			# FWER of IHW method

        POWER_pro <- sum(rej_pro[OD$et!=0])/max(1,sum(OD$et!=0))		# power of proposed
        POWER_bon <- sum(rej_bon[OD$et!=0])/max(1,sum(OD$et!=0))		# power of bonferroni
        POWER_rdw <- sum(rej_rdw[OD$et!=0])/max(1,sum(OD$et!=0))		# power of Roeder Wasserman method
        POWER_ihw <- sum(rej_ihwFwer[OD$et!=0])/max(1,sum(OD$et!=0))	# power of IHW method

        adjPval_pro <- p.adjust(OD$pt/weight_pro, method="BH")	# adjusted pvalue to compute FDR
        adjPval_bon <- p.adjust(OD$pt, method="BH")
        adjPval_rdw <- p.adjust(OD$pt/weight_rdw, method="BH")
        adjPval_ihw <- adj_pvalues(ihw_fdr)

        FDR_pro <- sum(adjPval_pro[OD$et==0] <= alpha)/max(1,sum(adjPval_pro <= alpha))	# FDR of proposed
        FDR_bh  <- sum(adjPval_bon[OD$et==0] <= alpha)/max(1,sum(adjPval_bon <= alpha))	# FDR of benjaminin and hochberg
        FDR_rdw <- sum(adjPval_rdw[OD$et==0] <= alpha)/max(1,sum(adjPval_rdw <= alpha))	# FDR of wasserman
        FDR_ihw <- sum(adjPval_ihw[OD$et==0] <= alpha)/max(1,rejections(ihw_fdr))		# FDR of IHW method

        FDR_POWER_pro <- sum(adjPval_pro[OD$et!=0] <= alpha)/max(1,sum(OD$et!=0))	# FDR of proposed
        FDR_POWER_bh  <- sum(adjPval_bon[OD$et!=0] <= alpha)/max(1,sum(OD$et!=0))	# FDR of benjaminin and hochberg
        FDR_POWER_rdw <- sum(adjPval_rdw[OD$et!=0] <= alpha)/max(1,sum(OD$et!=0))	# FDR of wasserman
        FDR_POWER_ihw <- sum(adjPval_ihw[OD$et!=0] <= alpha)/max(1,sum(OD$et!=0))		# FDR of IHW method

        return(c(FWER_pro,FWER_bon,FWER_rdw,FWER_ihw,POWER_pro,POWER_bon,POWER_rdw,POWER_ihw,
                 FDR_pro,FDR_bh,FDR_rdw,FDR_ihw,FDR_POWER_pro,FDR_POWER_bh,FDR_POWER_rdw,FDR_POWER_ihw))
    }
    FwerPowerFdrPower.simu <- sapply(1:simu,fun.FwerPowerFdrPower.simu)		# result of each replication
    FwerPowerFdrPower <- apply(FwerPowerFdrPower.simu,1,mean,na.rm=T)			# mean over replications
    return(FwerPowerFdrPower)
}



effectVec <- c(seq(0,1,.2),2,3)
datWeight <- read.csv("Binary_WeightByEffect_m10000.csv",h=T)

clusterExport(cl,"datWeight")
clusterExport(cl,"effectVec")
sim=1000
clusterExport(cl,"sim")


FwerPowerFdrPower1e1 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.2,corr=0,random=0,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,1:10])
FwerPowerFdrPower1e2 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.2,corr=0,random=1,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,1:10])
FwerPowerFdrPower1f1 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.2,corr=.3,random=0,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,1:10])
FwerPowerFdrPower1f2 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.2,corr=.3,random=1,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,1:10])
FwerPowerFdrPower1g1 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.2,corr=.5,random=0,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,1:10])
FwerPowerFdrPower1g2 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.2,corr=.5,random=1,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,1:10])
FwerPowerFdrPower1h1 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.2,corr=.7,random=0,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,1:10])
FwerPowerFdrPower1h2 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.2,corr=.7,random=1,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,1:10])
FwerPowerFdrPower1i1 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.2,corr=.9,random=0,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,1:10])
FwerPowerFdrPower1i2 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.2,corr=.9,random=1,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,1:10])



FwerPowerFdrPower2e1 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.5,corr=0,random=0,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,11:20])
FwerPowerFdrPower2e2 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.5,corr=0,random=1,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,11:20])
FwerPowerFdrPower2f1 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.5,corr=.3,random=0,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,11:20])
FwerPowerFdrPower2f2 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.5,corr=.3,random=1,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,11:20])
FwerPowerFdrPower2g1 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.5,corr=.5,random=0,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,11:20])
FwerPowerFdrPower2g2 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.5,corr=.5,random=1,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,11:20])
FwerPowerFdrPower2h1 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.5,corr=.7,random=0,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,11:20])
FwerPowerFdrPower2h2 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.5,corr=.7,random=1,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,11:20])
FwerPowerFdrPower2i1 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.5,corr=.9,random=0,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,11:20])
FwerPowerFdrPower2i2 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.5,corr=.9,random=1,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,11:20])


FwerPowerFdrPower3e1 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.75,corr=0,random=0,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,21:30])
FwerPowerFdrPower3e2 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.75,corr=0,random=1,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,21:30])
FwerPowerFdrPower3f1 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.75,corr=.3,random=0,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,21:30])
FwerPowerFdrPower3f2 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.75,corr=.3,random=1,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,21:30])
FwerPowerFdrPower3g1 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.75,corr=.5,random=0,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,21:30])
FwerPowerFdrPower3g2 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.75,corr=.5,random=1,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,21:30])
FwerPowerFdrPower3h1 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.75,corr=.7,random=0,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,21:30])
FwerPowerFdrPower3h2 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.75,corr=.7,random=1,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,21:30])
FwerPowerFdrPower3i1 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.75,corr=.9,random=0,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,21:30])
FwerPowerFdrPower3i2 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.75,corr=.9,random=1,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,21:30])


FwerPowerFdrPower4e1 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.9,corr=0,random=0,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,31:40])
FwerPowerFdrPower4e2 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.9,corr=0,random=1,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,31:40])
FwerPowerFdrPower4f1 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.9,corr=.3,random=0,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,31:40])
FwerPowerFdrPower4f2 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.9,corr=.3,random=1,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,31:40])
FwerPowerFdrPower4g1 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.9,corr=.5,random=0,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,31:40])
FwerPowerFdrPower4g2 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.9,corr=.5,random=1,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,31:40])
FwerPowerFdrPower4h1 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.9,corr=.7,random=0,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,31:40])
FwerPowerFdrPower4h2 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.9,corr=.7,random=1,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,31:40])
FwerPowerFdrPower4i1 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.9,corr=.9,random=0,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,31:40])
FwerPowerFdrPower4i2 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.9,corr=.9,random=1,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,31:40])


FwerPowerFdrPower5e1 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.99,corr=0,random=0,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,41:50])
FwerPowerFdrPower5e2 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.99,corr=0,random=1,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,41:50])
FwerPowerFdrPower5f1 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.99,corr=.3,random=0,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,41:50])
FwerPowerFdrPower5f2 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.99,corr=.3,random=1,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,41:50])
FwerPowerFdrPower5g1 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.99,corr=.5,random=0,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,41:50])
FwerPowerFdrPower5g2 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.99,corr=.5,random=1,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,41:50])
FwerPowerFdrPower5h1 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.99,corr=.7,random=0,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,41:50])
FwerPowerFdrPower5h2 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.99,corr=.7,random=1,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,41:50])
FwerPowerFdrPower5i1 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.99,corr=.9,random=0,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,41:50])
FwerPowerFdrPower5i2 <- parSapply(cl,1:length(effectVec),fun.FwerPowerFdrPower,simu=sim,null=.99,corr=.9,random=1,alpha=.05,effectVec=effectVec,datWeightByNull=datWeight[,41:50])



stopCluster(cl)













#======================== relationship between filter and test effect ==========


prob_relation_filterTestEffect <- function(r, rho,H0, ed, m0, m1)
{
    mean_ey = rho*ed
    sd_ey = sqrt(1 - rho^2)
    ey_val = rnorm(100, mean_ey, sd_ey)
    prob_condition_ey <- function(ey)
    {
	et <- ifelse(H0==0,0,ey)
      probs_per_ey = prob_rank_givenEffect(k=r, et=et, ey=ey, nrep = 10000,
                                             m0=m0, m1=m1)
      return(probs_per_ey)
    }
    prob_per_r = mean(sapply(ey_val,  prob_condition_ey))
    return(prob_per_r)
}
prob_test0 <- sapply(1:100, prob_relation_filterTestEffect, rho=0,H0=0, ed=2,m0=50,m1=50)
prob_test1 <- sapply(1:100, prob_relation_filterTestEffect, rho=0,H0=1, ed=2,m0=50,m1=50)

prob_test00 <- sapply(1:100, prob_relation_filterTestEffect, rho=.9,H0=0, ed=2,m0=50,m1=50)
prob_test11 <- sapply(1:100, prob_relation_filterTestEffect, rho=.9,H0=1, ed=2,m0=50,m1=50)

prob0 <- sapply(1:100, prob_rank_givenEffect, et=0, ey=2, nrep = 10000, m0=50, m1=50)
prob1 <- sapply(1:100, prob_rank_givenEffect, et=2, ey=2, nrep = 10000, m0=50, m1=50)

par(mfrow=c(2,2))
matplot(1:100, cbind(prob_test0,prob_test1, prob0, prob1))
matplot(1:100, cbind(prob_test00, prob_test11, prob0, prob1))


dev.off()
par(oma = c(1, 0, 0,0),mfrow=c(1,3))
matplot(1:100, cbind(prob0,prob1,prob_test0,prob_test1),
		type="l",lwd=2,col=6:1,xlab="ranks",ylab="P(rank|effect)")
matplot(1:100, cbind(prob0,prob1,prob_test00,prob_test00),
		type="l",lwd=2,col=6:1,xlab="ranks",ylab="P(rank|effect)")
matplot(1:100, cbind(prob0,prob1,rankProbH0_exact[,3],rankProbH1_exact[,3],rankProbH0_aprox[,3],rankProbH1_aprox[,3]),
		type="l",lwd=2,col=6:1,xlab="ranks",ylab="P(rank|effect)")
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 2, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n",main="et=2, m0 = 50, m1 = 50")
legend("bottom", c("DH0","DH1","FH0","FH1"), xpd = TRUE, horiz = TRUE, inset = c(0, 
    0), bty = "n", lty=1:4, ,col=4:1, cex = 1,lwd=2)











