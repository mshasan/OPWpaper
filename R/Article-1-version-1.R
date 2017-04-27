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
#source("https://bioconductor.org/biocLite.R")
#biocLite("qvalue")
library(qvalue)         # qvalue
#biocLite("IHW")
library(IHW)
#biocLite("IHWpaper")
library(IHWpaper)
library(tibble)       # data table
library(MASS)           # boc-cox transforamtion
#devtools::install_github("hadley/lineprof")
library(lineprof)       # check code performance
library(Rcpp)           # C++ library
library(wesanderson)  # for plot colors
#colors <- wes_palette("Cavalcanti")[1:4]
#biocLite("DESeq2")
library("DESeq2")

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
          geom_line(aes(linetype = variable), size = 1.5) +
          labs(x = "ranks", y = "p(rank | effect)", title = "(a) ey = 0, e.one = 0") +
          theme(legend.position="none") +
          annotate("text", x=50, y=.011, label="P(rank | effect = 0)")


p01 <- ggplot(dat01, aes(x = ranks, y = value, group = variable,
                       colour = variable)) +
          geom_line(aes(linetype = variable), size = 1.5) +
          labs(x = "ranks", y = "p(rank | effect)", title = "(b) ey ~ U(0, 1), e.one = 1") +
          theme(legend.position="none") +
          annotate("text", x = c(50, 50), y = c(.005, .018),
                   label = c(paste(sprintf('\u2190'),"P(rank | effect = 0)"),
					paste("P(rank | effect = e.one)", sprintf('\u2192'))))

p12 <- ggplot(dat12, aes(x = ranks, y = value, group = variable,
                       colour = variable)) +
          geom_line(aes(linetype = variable), size = 1.5) +
          labs(x = "ranks", y = "p(rank | effect)", title = "(c) ey ~ U(1, 2), e.one = 1") +
          theme(legend.position="none") +
          annotate("text", x = c(60, 40), y=c(.001, .02),
                   label = c(paste(sprintf('\u2190'),"P(rank | effect = 0)"),
					paste("P(rank | effect = e.one)", sprintf('\u2192'))))

p02 <- ggplot(dat02, aes(x = ranks, y = value, group = variable,
                       colour = variable)) +
          geom_line(aes(linetype = variable), size = 1.5) +
          labs(x = "ranks", y = "p(rank | effect)", title = "(d) ey ~ U(0, 1), e.one = 2") +
          theme(legend.title = element_blank(), legend.position="bottom") +
          annotate("text", x = c(50, 40), y=c(.04, .15),
                  label = c("P(rank | effect = 0)",paste(sprintf('\u2190'),"P(rank | effect = e.one)")))


# extract the legend from one of the plots
legend_art <- get_legend(p02 + theme(legend.direction="horizontal",
                                 legend.position="bottom"))

p02 = p02 +    theme(legend.position="none")


# arrange the plots
grid.arrange(arrangeGrob(p00, p01, p12, p02, nrow = 2),legend_art, nrow=2, heights=c(7,1),
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
prob <- sapply(ranks, prob_rank_givenEffect, et=1, ey=1, nrep = nrep , m0=9000, m1=1000)
#prob <- read.csv("dataBin_ProbRanksByEffect_m10000.csv",h=T)[,16]

weight <- weight_continuous(alpha=.05, et=1, m=10000, tail=1, delInterval=.0001,
                            prob=prob)

dat2 = data.frame(ranks, prob, weight)
p_prob = ggplot(dat2, aes(x = ranks, y = prob)) +  geom_line(size=1.5, col="firebrick4") +
                   labs(y = "p(rank | effect)")
p_weight = ggplot(dat2, aes(x = ranks, y = weight)) + geom_line(size=1.5, col="firebrick4")
 
grid.arrange(p_prob, p_weight, ncol = 2,
                      top=paste0("Continuous: E(", expression(epsilon),") = 1, m = 10000, m0 = 9000, m1 = 1000"))





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

simu_TypeI_error <- function(s)
{
    TypeI_error <- function(alpha, m, nrep)
    {
        # create data set--------------
        pval <- runif(m)
        pval_filter <- runif(m)
        test = qnorm(pval, lower.tail=FALSE)
        filter = qnorm(pval_filter, lower.tail=FALSE)

        #dat = tibble(test, pval, filter)

        #OD = dat[order(dat$filter, decreasing=T), ]
        #odered.pvalue = OD$pval
        odered.pvalue = sort(pval)

        null = qvalue(pval)$pi0
        m0 = ceiling(m*null)
        m1 = m-m0

        #model = lm(filter ~ test)
        test_effect <- if(m1 == 0) {0
        } else {sort(test, decreasing = T)[1:m1]}

        et_bin = median(test_effect, na.rm = T)		        # median from the summary for the binary case test effect
        et_cont = mean(test_effect, na.rm = T)
        #ey_bin = model$coef[[1]] + model$coef[[2]]*et_bin
        #ey_cont = model$coef[[1]] + model$coef[[2]]*et_cont

        #prob_bin <-vapply(1:m, prob_rank_givenEffect, 1, et = ey_bin,
        #                  ey = ey_bin, nrep = nrep, m0 = m0, m1 = m1)
        #prob_cont <-vapply(1:m, prob_rank_givenEffect, 1, et = abs(ey_cont),
        #                   ey = abs(ey_cont), nrep = nrep, m0 = m0, m1 = m1)
        prob_bin <- runif(m)
        prob_cont <- runif(m)

        w_bin <- weight_binary(alpha = alpha, et = et_bin, m = m, m1 = m1, tail = 1,
                               delInterval = .0001, prob = prob_bin)
        w_cont = weight_continuous(alpha= alpha, et = et_cont, m = m, tail = 1,
                                   delInterval=.0001 , prob = prob_cont)

        #ihw_fwer <- ihw(dat$pval, dat$filter, alpha=alpha, adjustment_type = "bonferroni")
        ihw_fwer <- ihw(pval, filter, alpha=alpha, adjustment_type = "bonferroni")# IHW method for FWER

        # complete pvalues
        #-----------------------
        bon = sum(pval <= alpha/m, na.rm = T)/m
        pro_bin = sum(odered.pvalue <= alpha*w_bin/m, na.rm = T)/m
        pro_cont = sum(odered.pvalue <= alpha*w_cont/m, na.rm = T)/m
        IHW <- rejections(ihw_fwer)/m

        return(c(bon, pro_bin, pro_cont, IHW))
    }

    # function to conduct FWER simulations
    #----------------------------------------
    alphaVal = seq(.01, .1, .02)
    TypeI_error_mat = sapply(alphaVal, TypeI_error, m=1000, nrep = 10000)
    return(TypeI_error_mat)
}

simuVal = 1:1000
simu_TypeI_error_mat = sapply(simuVal, simu_TypeI_error)

# load saved data---------------------
# this code is to load saved workspace from parallel computing
load("U:/Documents/My Research (UGA)/Multiple Hypoetheses/Article-1/smu_Type_I_error5.RDATA")

# IHW paper------------
null_folder <- system.file("simulation_benchmarks/result_files/ihw_bonf_null",
                           package = "IHWpaper")
null_sim <- bind_rows(lapply(file.path(null_folder,list.files(null_folder)), function(x) readRDS(x))) %>%
    mutate(fdr_method = ifelse(fdr_method == "IHW-Bonferroni E3", "IHW-Bonferroni", fdr_method))
null_sim_val = null_sim$FWER[c(1,7,11,15,19)]
#-----------end------------IHWpaper------------------

typeIerror_mat = matrix(apply(simu_TypeI_error_mat, 1, mean), nrow=4, byrow = FALSE)

alphaVal = seq(.01, .1, .02)
datError <- data.frame(alphaVal, t(typeIerror_mat[-4,])*10000,null_sim_val)
colnames(datError) <- c("alpha","BON","PRO_bin","PRO_cont", "IHW")
datError2 <- melt(datError, id.var="alpha")

panel_d = ggplot(datError2, aes(x = alpha, y = value, col=variable)) +
      geom_line(size=1.2) +
      geom_abline(linetype="dashed") +
      xlab(expression(bold(paste("Nominal ",alpha)))) +
      ylab("FWER")+
      scale_x_continuous(limits= c(0.01,0.1), breaks=seq(0.01,0.09,length=5)) +
      #ylim(0,0.9) +
      theme(legend.title = element_blank())+
      theme(axis.title = element_text(face="bold") )
      #theme(panel.background = element_rect(fill = 'white', colour = 'black'))

#df = data.frame(colour=mycolours, last_vals=c(0.08,   0.085,   0.09), label=c("BON","PRO_bin","PRO_cont"))
#panel_d <- pretty_legend(panel_d, df, .092)
panel_d




#================end of FWER====================================================











# function to generate uniform random number  for fixed mean
runif_by_mean <- function(n, mean)
  {
  sd = mean/2
  uni_rv <- mean + sd*scale(runif(n, 0, 1))
  return(as.vector(uni_rv))
  }

x = runif_by_mean(n = 100, mean = 3)
summary(x)


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

eVec = rep(2,50)
S <- matrix(.9, 50, 50) + diag(50)*(1-.1)		# test correlation matrix
testsStat = test_by_block(r=1, eVec, groupSize=50,  Sigma)

  ################################################################################################################
  #----------------------2b-3:fun.FwerPowerFdrPower----------------------------
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
# effectVec =  different effect size (10 values)
# datWeightByNull = weight matrix (10000 by 40) for various null=c(20,50,90,99)% computed before;
# 			m=10,000 and effect = effectVec*4 = (10 effects)*(4 null scenerios)=40

#output:
#=============
# FwerPowerFdr = A 16 by 10 matrix consits of simulated FWER(first 4 rows), POWER (2nd 4 rows),
# FDR (3rd 4 rows) and FDRPower(last 4 rows) for 10 different effect sizes
#-----------------------------------------------------------------------------------------------------


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

# an example-----------------

effecSize = c(.3, .5, .7, 1, 2, 3)
SigmaVal <- matrix(testCorr, groupSize, groupSize) + diag(groupSize)*(1-testCorr)		# test correlation matrix

fwerPowerFdrPower_mat = sapply(effecSize, fwerPowerFdrPower_by_effect, m = 1000,
                               null = .9, testCorr = 0, random=0, groupSize = 100,
                               alpha = .05, Sigma = SigmaVal)




# function for simulations--------------------------
simu_fwerPowerFdrPower_by_effect <- function(s, filterEffect, m, null, testCorr, random, groupSize, alpha)
  {
    fwerPowerFdrPower_mat = sapply(effecSize, fwerPowerFdrPower_by_effect, m = m,
                                   null = null, testCorr = testCorr, random = random,
                                   groupSize = 100, alpha = .05, Sigma = SigmaVal)
    return(fwerPowerFdrPower_mat)
  }


simuVal = 1:3
FwerPowerFdrPower_simu <- sapply(simuVal, simu_fwerPowerFdrPower_by_effect, filterEffect=effecSize,
                                 m = 1000, null=.9, testCorr=0, random=1, groupSize = 100, alpha = .05)

# load data for continuous power----------------
# this code is to load saved workspace from parallel computing
load("U:/Documents/My Research (UGA)/Multiple Hypoetheses/Article-1/simu_Continuous_FwerPowerFdrPower.RDATA")



# plots
#-----------
dat_50 <- data.frame(effectVec, t(FwerPowerFdrPower2i1[13:16,]))
colnames(dat_50) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_50_all <- melt(dat_50, id.var = "effectSize")
p_50_all <- ggplot(dat_50_all, aes(x = effectSize, y = value,group = variable,
                             col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "effect size", y = "power", title = "null = 50%") +
    theme(legend.position="none")

dat_50_par <- melt(dat_50[1:6,], id.var = "effectSize")
p_50_par <- ggplot(dat_50_par, aes(x = effectSize, y = value,group = variable,
                             col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "effect size", y = "power", title = "null = 50%") +
    theme(legend.position="none")

prow1 <- plot_grid(p_50_all, p_50_par, align = 'hv', ncol = 1)



dat_90 <- data.frame(effectVec, t(FwerPowerFdrPower4i1[13:16,]))
colnames(dat_90) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_90_all <- melt(dat_90, id.var = "effectSize")
p_90_all <- ggplot(dat_90_all, aes(x = effectSize, y = value,group = variable,
                                   col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "effect size", y = "power", title = "null = 90%") +
    theme(legend.position="none")

dat_90_par <- melt(dat_90[1:6,], id.var = "effectSize")
p_90_par <- ggplot(dat_90_par, aes(x = effectSize, y = value,group = variable,
                                   col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "effect size", y = "power", title = "null = 90%") +
    theme(legend.position="none")

prow2 <- plot_grid(p_90_all, p_90_par, align = 'hv', ncol = 1)



dat_99 <- data.frame(effectVec, t(FwerPowerFdrPower5i1[13:16,]))
colnames(dat_99) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_99_all <- melt(dat_99, id.var = "effectSize")
p_99_all <- ggplot(dat_99_all, aes(x = effectSize, y = value,group = variable,
                                   col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "effect size", y = "power", title = "null = 99%") +
    theme(legend.position="none")

dat_99_par <- melt(dat_99[1:6,], id.var = "effectSize")
p_99_par <- ggplot(dat_99_par, aes(x = effectSize, y = value,group = variable,
                                   col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "effect size", y = "power", title = "null = 99%") +
    theme(legend.title = element_blank())

prow3 <- plot_grid(p_99_all, p_99_par, align = 'hv', ncol = 1)



legend_power <- get_legend(p_99_par + theme(legend.direction="horizontal",
                                        legend.position="bottom"))
p_99_par <- p_99_par + theme(legend.position="none")

grid.arrange(arrangeGrob(prow1, prow2, prow3, ncol=3), legend_power, nrow=2, 
             heights=c(7,1), top = "Power: mean test effect = mean filter effect")


# if test and filter effect size are not same----------------------------------
# this code is to load saved workspace from parallel computing
load("U:/Documents/My Research (UGA)/Multiple Hypoetheses/Article-1/simu_fwerPowerFdrPower_random.RDATA")


FwerPowerFdrPower_simu <- FwerPowerFdrPower_simu_b1
FwerPowerFdrPower = apply(FwerPowerFdrPower_simu, 1, mean)
FwerPowerFdrPower_by_effectSize <- matrix(FwerPowerFdrPower, nrow=16, byrow = FALSE)
effecSize = c(.3, .5, .7, 1, 2, 3)
rownames(FwerPowerFdrPower_by_effectSize) <- c('FWER_pro', 'FWER_bon', 'FWER_rdw', 'FWER_ihw',
                                               'POWER_pro', 'POWER_bon', 'POWER_rdw', 'POWER_ihw',
                                               'FDR_pro', 'FDR_bh', 'FDR_rdw', 'FDR_ihw',
                                               'FDR_POWER_pro', 'FDR_POWER_bh', 'FDR_POWER_rdw', 'FDR_POWER_ihw')
FwerPowerFdrPower_by_effectSize2 <- data.frame(effecSize, t(FwerPowerFdrPower_by_effectSize))

datFwer = FwerPowerFdrPower_by_effectSize2[,c(1, 14:17)]
datFwer2 <- melt(datFwer, id.var="effecSize")

panel_d = ggplot(datFwer2, aes(x = effecSize, y = value, col=variable)) +
    geom_line(size=1.2) +
    xlab("effect size") +
    ylab("FWER")+
    scale_x_continuous(breaks=effecSize) +
    #scale_x_continuous(breaks=effecSize) +
    #ylim(0,1) +
    theme(legend.title = element_blank())+
    theme(axis.title = element_text(face="bold") )
#theme(panel.background = element_rect(fill = 'white', colour = 'black'))

#df = data.frame(colour=mycolours, last_vals=c(0.08,   0.085,   0.09), label=c("BON","PRO_bin","PRO_cont"))
#panel_d <- pretty_legend(panel_d, df, .092)
panel_d


# plots
# 90% null for all cases, cor =0, test effect != filter effect
#----------------------------------------------------------------
effecSize = c(.3, .5, .7, 1, 2, 3)

rand_power_b1 <- matrix(apply(FwerPowerFdrPower_simu_b1, 1, mean), nrow=16, byrow = FALSE)
dat_rand_b1 <- data.frame(effecSize, t(rand_power_b1[13:16,]))
colnames(dat_rand_b1) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_rand_b1_melt <- melt(dat_rand_b1, id.var = "effectSize")

p_rand_b1 <- ggplot(dat_rand_b1_melt, aes(x = effectSize, y = value,group = variable,
                                   col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "effect size", y = "power", title = "et = 2*ef") +
    theme(legend.position="none")


rand_power_b2 <- matrix(apply(FwerPowerFdrPower_simu_b2, 1, mean), nrow=16, byrow = FALSE)
dat_rand_b2 <- data.frame(effecSize, t(rand_power_b2[13:16,]))
colnames(dat_rand_b2) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_rand_b2_melt <- melt(dat_rand_b2, id.var = "effectSize")

p_rand_b2 <- ggplot(dat_rand_b2_melt, aes(x = effectSize, y = value,group = variable,
                                          col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "effect size", y = "power", title = "et = 3*ef") +
    theme(legend.position="none")


rand_power_b3 <- matrix(apply(FwerPowerFdrPower_simu_b3, 1, mean), nrow=16, byrow = FALSE)
dat_rand_b3 <- data.frame(effecSize, t(rand_power_b3[13:16,]))
colnames(dat_rand_b3) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_rand_b3_melt <- melt(dat_rand_b3, id.var = "effectSize")

p_rand_b3 <- ggplot(dat_rand_b3_melt, aes(x = effectSize, y = value,group = variable,
                                          col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "effect size", y = "power", title = "et = 5*ef") +
    theme(legend.position="none")


rand_power_b4 <- matrix(apply(FwerPowerFdrPower_simu_b4, 1, mean), nrow=16, byrow = FALSE)
dat_rand_b4 <- data.frame(effecSize, t(rand_power_b4[13:16,]))
colnames(dat_rand_b4) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_rand_b4_melt <- melt(dat_rand_b4, id.var = "effectSize")

p_rand_b4 <- ggplot(dat_rand_b4_melt, aes(x = effectSize, y = value,group = variable,
                                          col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "effect size", y = "power", title = "et = 1/2*ef") +
    theme(legend.position="none")


rand_power_b5 <- matrix(apply(FwerPowerFdrPower_simu_b5, 1, mean), nrow=16, byrow = FALSE)
dat_rand_b5 <- data.frame(effecSize, t(rand_power_b5[13:16,]))
colnames(dat_rand_b5) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_rand_b5_melt <- melt(dat_rand_b5, id.var = "effectSize")

p_rand_b5 <- ggplot(dat_rand_b5_melt, aes(x = effectSize, y = value,group = variable,
                                          col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "effect size", y = "power", title = "et = 1/3*ef") +
    theme(legend.position="none")


rand_power_b6 <- matrix(apply(FwerPowerFdrPower_simu_b6, 1, mean), nrow=16, byrow = FALSE)
dat_rand_b6 <- data.frame(effecSize, t(rand_power_b6[13:16,]))
colnames(dat_rand_b6) <- c("effectSize", "PRO", "BH", "RDW", "IHW")
dat_rand_b6_melt <- melt(dat_rand_b6, id.var = "effectSize")

p_rand_b6 <- ggplot(dat_rand_b6_melt, aes(x = effectSize, y = value,group = variable,
                                          col=variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "effect size", y = "power", title = "et = 1/5*ef") +
    theme(legend.title = element_blank())




legend_power_rand <- get_legend(p_rand_b6 + theme(legend.direction="horizontal",
                                            legend.position="bottom"))
p_rand_b6 <- p_rand_b6 + theme(legend.position="none")

grid.arrange(arrangeGrob(p_rand_b1, p_rand_b2, p_rand_b3, p_rand_b4, p_rand_b5,
                         p_rand_b6, nrow=2), legend_power_rand, nrow=2, heights=c(7,1), 
             top = "Power: mean test effect != mean filter effect")






# see correaltion effect FWER/POWER/FDR
#-----------------------------------------
# load data for continuous power----------------
# this code is to load saved workspace from parallel computing
load("U:/Documents/My Research (UGA)/Multiple Hypoetheses/Article-1/simu_Continuous_FwerPowerFdrPower.RDATA")


effectVec <- c(seq(0,1,.2),2,3,5,8)
E = FwerPowerFdrPower4e1
F = FwerPowerFdrPower4f1
G = FwerPowerFdrPower4g1
H = FwerPowerFdrPower4h1
I = FwerPowerFdrPower4i1
corr = c(0,.3,.5,.7,.9)
r = 13
gplots <- list()
for(e in 3:8)				# effect size index
{
    PRO = c(E[r,e],    F[r,e],    G[r,e],    H[r,e],    I[r,e])
    BH = c(E[(r+1),e],F[(r+1),e],G[(r+1),e],H[(r+1),e],I[(r+1),e])
    RDW = c(E[(r+2),e],F[(r+2),e],G[(r+2),e],H[(r+2),e],I[(r+2),e])
    IHW = c(E[(r+3),e],F[(r+3),e],G[(r+3),e],H[(r+3),e],I[(r+3),e])
    dat = data.frame(corr, PRO, BH, RDW, IHW)
    dat2 = melt(dat, id.var = "corr")
    gplots[[e]] <- ggplot(dat2, aes(x = corr, y = value, group = variable,
                                          col = variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "corr", y = "power", title = paste("et = ", effectVec[e])) +
    #theme(legend.title = element_blank())
    theme(legend.position="none")
}

gplots[[8]] <- gplots[[8]] + theme(legend.position="bottom", legend.title = element_blank())

legend_corr <- get_legend(gplots[[8]])

gplots[[8]] <- gplots[[8]] + theme(legend.position="none")

grid.arrange(arrangeGrob(gplots[[3]],gplots[[4]],gplots[[5]],gplots[[6]],gplots[[7]],gplots[[8]], nrow=2), 
		legend_corr, nrow=2, heights=c(7,1), 
             top = "null = 90%, mean test effect = mean filter effect")




                                                                                 







#======================== relationship between filter and test effect ==========


prob_relation_filterTestEffect <- function(r, rho, H0, ed, m0, m1)
{
    mean_ey = rho*ed
    sd_ey = sqrt(1 - rho^2)
    ey_val = rnorm(100, mean_ey, sd_ey)
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

m = 100
m0 = 50
m1 = 50
ed=0

prob_test0_cor.2 <- sapply(1:100, prob_relation_filterTestEffect, rho=.2, H0=0,
                           ed=ed, m0=m0, m1=m1)
prob_test1_cor.2 <- sapply(1:100, prob_relation_filterTestEffect, rho=.2, H0=1,
                           ed=ed, m0=m0, m1=m1)

prob_test0_cor.5 <- sapply(1:100, prob_relation_filterTestEffect, rho=.5, H0=0,
                           ed=ed, m0=m0,m1=m1)
prob_test1_cor.5 <- sapply(1:100, prob_relation_filterTestEffect, rho=.5, H0=1,
                           ed=ed, m0=m0,m1=m1)

prob_test0_cor.8 <- sapply(1:100, prob_relation_filterTestEffect, rho=.8, H0=0,
                           ed=ed,m0=m0,m1=m1)
prob_test1_cor.8 <- sapply(1:100, prob_relation_filterTestEffect, rho=.8, H0=1,
                           ed=ed,m0=m0,m1=m1)

prob0 <- sapply(1:100, prob_rank_givenEffect, et=0, ey=ed, nrep = 10000, m0=m0, m1=m1)
prob1 <- sapply(1:100, prob_rank_givenEffect, et=ed, ey=ed, nrep = 10000, m0=m0, m1=m1)

par(mfrow=c(1,3))
matplot(1:100, cbind(prob_test0_cor.2, prob_test1_cor.2, prob0, prob1))
matplot(1:100, cbind(prob_test0_cor.5, prob_test1_cor.5, prob0, prob1))
matplot(1:100, cbind(prob_test0_cor.8, prob_test1_cor.8, prob0, prob1))



# this code is to load saved workspace from parallel computing
load("U:/Documents/My Research (UGA)/Multiple Hypoetheses/Article-1/smu_relation_filterTest.RDATA")


# nice plots---------
ranks=1:100
datRelaion1 <- data.frame(ranks, prob0, prob1, prob_test0_cor.2, prob_test1_cor.2)
colnames(datRelaion1) <- c("ranks", "FH0","FH1","TH0","TH1")
datRelaion1_melt1 <- melt(datRelaion1, id.var="ranks")

datRelaion2 <- data.frame(ranks, prob0, prob1, prob_test0_cor.5, prob_test1_cor.5)
colnames(datRelaion2) <- c("ranks", "FH0","FH1","TH0","TH1")
datRelaion1_melt2 <- melt(datRelaion2, id.var="ranks")

datRelaion3 <- data.frame(ranks, prob0, prob1, prob_test0_cor.8, prob_test1_cor.8)
colnames(datRelaion3) <- c("ranks", "FH0","FH1","TH0","TH1")
datRelaion1_melt3 <- melt(datRelaion3, id.var="ranks")


p_.2 <- ggplot(datRelaion1_melt1, aes(x = ranks, y = value, group = variable,
                         colour = variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ranks", y = "p(rank | effect)", title = "cor = .2") +
    theme(legend.position="none")

p_.5 <- ggplot(datRelaion1_melt2, aes(x = ranks, y = value, group = variable,
                                    colour = variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ranks", y = "p(rank | effect)", title = "cor = .5") +
    theme(legend.title = element_blank(), legend.position="bottom")

p_.8 <- ggplot(datRelaion1_melt3, aes(x = ranks, y = value, group = variable,
                         colour = variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ranks", y = "p(rank | effect)", title = "cor = .8") +
    theme(legend.position="none")

# extract the legend from one of the plots
legend_rel <- get_legend(p_.5 + theme(legend.direction="horizontal",
                                 legend.position="bottom"))
p_.5 = p_.5 + theme(legend.position="none")

# arrange the plots
grid.arrange(p_.2, p_.5, p_.8,ggplot(NULL),legend_rel,ggplot(NULL),ncol=3, heights=c(7,1),
             top= "et = 2, m0 = 50, m1 = 50")


# nice plots---------
ranks=1:100
datRelaion1 <- data.frame(ranks, prob02, prob12, prob_test0_cor.22, prob_test1_cor.22)
colnames(datRelaion1) <- c("ranks", "FH0","FH1","TH0","TH1")
datRelaion1_melt1 <- melt(datRelaion1, id.var="ranks")

datRelaion2 <- data.frame(ranks, prob02, prob12, prob_test0_cor.52, prob_test1_cor.52)
colnames(datRelaion2) <- c("ranks", "FH0","FH1","TH0","TH1")
datRelaion1_melt2 <- melt(datRelaion2, id.var="ranks")

datRelaion3 <- data.frame(ranks, prob02, prob12, prob_test0_cor.82, prob_test1_cor.82)
colnames(datRelaion3) <- c("ranks", "FH0","FH1","TH0","TH1")
datRelaion1_melt3 <- melt(datRelaion3, id.var="ranks")


p_.2 <- ggplot(datRelaion1_melt1, aes(x = ranks, y = value, group = variable,
                                      colour = variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ranks", y = "p(rank | effect)", title = "cor = .2") +
    theme(legend.position="none")

p_.5 <- ggplot(datRelaion1_melt2, aes(x = ranks, y = value, group = variable,
                                      colour = variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ranks", y = "p(rank | effect)", title = "cor = .5") +
    theme(legend.title = element_blank(), legend.position="bottom")

p_.8 <- ggplot(datRelaion1_melt3, aes(x = ranks, y = value, group = variable,
                                      colour = variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ranks", y = "p(rank | effect)", title = "cor = .8") +
    theme(legend.position="none")

# extract the legend from one of the plots
legend_rel <- get_legend(p_.5 + theme(legend.direction="horizontal",
                                      legend.position="bottom"))
p_.5 = p_.5 + theme(legend.position="none")

# arrange the plots
#dev.new(width=8, height=4)
grid.arrange(p_.2, p_.5, p_.8,ggplot(NULL),legend_rel,ggplot(NULL),ncol=3, heights=c(7,1),
             top= "et = 2, m0 = 90, m1 = 10")
#=============================end===============================================



# nice plots---------
ranks=1:100
datRelaion1 <- data.frame(ranks, prob0, prob1, prob_test0_cor.2, prob_test1_cor.2)
colnames(datRelaion1) <- c("ranks", "FH0","FH1","TH0","TH1")
datRelaion1_melt1 <- melt(datRelaion1, id.var="ranks")

datRelaion2 <- data.frame(ranks, prob0, prob1, prob_test0_cor.5, prob_test1_cor.5)
colnames(datRelaion2) <- c("ranks", "FH0","FH1","TH0","TH1")
datRelaion1_melt2 <- melt(datRelaion2, id.var="ranks")

datRelaion3 <- data.frame(ranks, prob0, prob1, prob_test0_cor.8, prob_test1_cor.8)
colnames(datRelaion3) <- c("ranks", "FH0","FH1","TH0","TH1")
datRelaion1_melt3 <- melt(datRelaion3, id.var="ranks")


p_.2 <- ggplot(datRelaion1_melt1, aes(x = ranks, y = value, group = variable,
                                      colour = variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ranks", y = "p(rank | effect)", title = "cor = .2") +
    theme(legend.position="none")

p_.5 <- ggplot(datRelaion1_melt2, aes(x = ranks, y = value, group = variable,
                                      colour = variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ranks", y = "p(rank | effect)", title = "cor = .5") +
    theme(legend.title = element_blank(), legend.position="bottom")

p_.8 <- ggplot(datRelaion1_melt3, aes(x = ranks, y = value, group = variable,
                                      colour = variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = "ranks", y = "p(rank | effect)", title = "cor = .8") +
    theme(legend.position="none")

# extract the legend from one of the plots
legend_rel <- get_legend(p_.5 + theme(legend.direction="horizontal",
                                      legend.position="bottom"))
p_.5 = p_.5 + theme(legend.position="none")

# arrange the plots
#dev.new(width=8, height=4)
grid.arrange(p_.2, p_.5, p_.8,ggplot(NULL),legend_rel,ggplot(NULL),ncol=3, heights=c(7,1),
             top= "et = 2, m0 = 99, m1 = 1")
#=============================end===============================================



#========================data application=======================================

##########################################################################################################

#---------------Example-1: bottomly data (RNA-seq)-----------------------------


########################################################################################################

# data processing
#-------------------------
bottomly_count_table <- read.table("bottomly_count_table.txt",h=T)
bottomly_phenodata <- read.table("bottomly_phenodata.txt",h=T)
countData <- as.matrix(bottomly_count_table[,-1])		# counts
condition <- factor(bottomly_phenodata[,3])				# strain as factor
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), ~ condition)
dds <- DESeq(dds)
bottomly <- as.data.frame(results(dds))
colnames(bottomly)

hist(pval)    # two-tailed test 

test = bottomly$stat
pval = bottomly$pvalue
filter = bottomly$baseMean

Data <- tibble(test, pval, filter)	# data of filter covariate and pvlaues

# fit regression to obtain filter and test effect sizes
#-------------------------------------------------------------
nrep = 10000
m = length(pval)
null = qvalue(pval, pi0.method="bootstrap")$pi0
m0 = ceiling(null*m)
m1 = m-m0

bc <- boxcox(filter ~ test)
trans <- bc$x[which.max(bc$y)]
model_bot <- lm(filter^trans ~ test)


test_effect = if(m1 == 0) {0
} else {sort(abs(test), decreasing = T)[1:m1]}		# two-tailed test

et_cont = mean(test_effect, na.rm = T)
ey_cont = model_bot$coef[[1]] + model_bot$coef[[2]]*et_cont
et_bin = median(test_effect, na.rm = T)
ey_bin = model_bot$coef[[1]] + model_bot$coef[[2]]*et_bin

prob_bin <-vapply(1:m, prob_rank_givenEffect, 1, et = ey_bin,
                  ey = ey_bin, nrep = nrep, m0 = m0, m1 = m1)
prob_cont <-vapply(1:m, prob_rank_givenEffect, 1, et = ey_cont,
                   ey = ey_cont, nrep = nrep, m0 = m0, m1 = m1)

alphaVec = seq(.05,.1,length.out = 5)
w_bin <- sapply(alphaVec, weight_binary, et = et_bin, m = m, m1 = m1, tail = 1,
                delInterval = .0001, prob = prob_bin)
w_cont = sapply(alphaVec, weight_continuous, et = et_cont, m = m, tail = 1,
                delInterval=.0001 , prob = prob_cont)


# function to compute number of rejections
# input:
#=======
# i=index number
# alpha = significance level
# Data = data to be analyzed composed of P=pvalues and Xf=filter covariate
# W_bin_mat = Binary weight matrix for diffecrent alpha
# W_cont_mat = Continuous weight matrix for different alpha
# output:
#========
# rej_mat = a rejection matrix composed of rejections from variaous methods
#---------------------------------------------------------------
fun.rejections <- function(i,alphaVec,Data,W_bin_mat,W_cont_mat)
{
    alpha=alphaVec[i]
    W_bin=as.vector(W_bin_mat[,i])
    W_cont=as.vector(W_cont_mat[,i])
    
    m = length(W_bin)
    
    OD <- Data[order(Data$filter,decreasing=T),]		# odered by covariate
    Ordered.pvalue <- OD$pval					# odered pvalues for all tests
    
    # preprocesing before counting rejections
    #----------------------------------------
    ihw_res_bon <- ihw(Data$pval,Data$filter, alpha=alpha, nbins=4,nsplits_internal=5,
                       lambdas=seq(0,3,length=20),adjustment_type = "bonferroni")
    padj_Pro_bin <-p.adjust(Ordered.pvalue/W_bin, method = "BH")		# proposed method based on right-tailed pvalue
    padj_Pro_cont <-p.adjust(Ordered.pvalue/W_cont, method = "BH")
    padj_BH <- p.adjust(Data$pval, method = "BH")
    ihw_res_fdr <- ihw(Data$pval,Data$filter, alpha=alpha, nbins=13,nsplits_internal=5L, 
                       nfolds_internal=4L)
    
    # rejections by FWER
    #-------------------
    Pro_bon_bin = sum(Ordered.pvalue <= alpha*W_bin/m, na.rm = TRUE)
    Pro_bon_cont = sum(Ordered.pvalue <= alpha*W_cont/m, na.rm = TRUE)
    bon = sum(Data$pval <= alpha/m, na.rm = TRUE)
    ihw_bon = rejections(ihw_res_bon)
    
    # rejections by FDR
    #------------------
    Pro_bh_bin = sum(padj_Pro_bin <= alpha, na.rm = TRUE)
    Pro_bh_cont = sum(padj_Pro_cont <= alpha, na.rm = TRUE)
    bh = sum(padj_BH <= alpha, na.rm = TRUE)
    ihw_bh = rejections(ihw_res_fdr)
    
    return(c(Pro_bon_bin,Pro_bon_cont,bon,ihw_bon,Pro_bh_bin,Pro_bh_cont,bh,ihw_bh))
}
rej_mat_bot = sapply(1:length(alphaVec),fun.rejections,alphaVec,Data=Data,W_bin_mat=w_bin,W_cont_mat=w_cont)


# from IHW paper-------------
rnaseq_file <- system.file("real_data_examples/result_files", "RNAseq_benchmark.Rds", package = "IHWpaper")
rnaseq_data <- readRDS(file=rnaseq_file)
panel_a_data <- group_by(rnaseq_data$alpha_df, alpha) %>% summarize(BH = max(bh_rejections), IHW=max(rejections)) %>% 
    gather(method, rejections, BH, IHW)


# nice plots---------------
# rej_mat_bot_FWER <- data.frame(alphaVec, t(rej_mat_bot[1:4,]))
# colnames(rej_mat_bot_FWER) <- c("alpha", "PRO_bin","PRO_cont","BON","IHW" )
# rej_mat_bot_FWER2 <- melt(rej_mat_bot_FWER, id.var = "alpha")

IHW = c(1312, 1434, 1535, 1665, 1783)
BH = c(1211, 1339, 1414, 1507, 1618)
rej_mat_bot_FDR <- data.frame(alphaVec, t(rej_mat_bot[5:6,]),BH, IHW)
colnames(rej_mat_bot_FDR) <- c("alpha", "PRO_bin","PRO_cont","BH","IHW" )
rej_mat_bot_FDR2 <- melt(rej_mat_bot_FDR, id.var = "alpha")


# p_fwer <- ggplot(rej_mat_bot_FWER2, aes(x = alpha, y = value, group = variable,
#                                       colour = variable)) +
#     geom_line(aes(linetype = variable), size = 1.2) +
#     labs(x = "alpha", y = "discoveries", title = "FWER based") +
#     theme(legend.title = element_blank())
#     #theme(panel.background = element_rect(fill = 'white', colour = 'black'))

#df = data.frame(colour=colors, last_vals=c(370, 382, 348, 391),
                label=c("PRO_bin","PRO_cont","BH","IHW"))
#panel_d <- pretty_legend(p_fwer, df, .1)
#panel_d

p_fdr_bot <- ggplot(rej_mat_bot_FDR2, aes(x = alpha, y = value, group = variable,
                                       colour = variable)) +
    geom_line(aes(linetype = variable), size = 1.2) +
    labs(x = "alpha", y = "discoveries", title = "FDR based") +
    theme(legend.title = element_blank())


# # arrange the plots
# grid.arrange(p_fwer, p_fdr, ncol=2,
#     top= paste0("Bottomly: et_bin = ", round(et_bin,1), ", et_cont = ", round(et_cont,1),
#                 ", ey_bin = ", round(ey_bin,1), ", ey_cont = ", round(ey_cont,1),"\n",
#                          "m = ",m,", null = ", round(null*100),"%"))





#########################################################################################################

#--------------Example-2: proteomics------------------------------

##########################################################################################################

# data processing
#-------------------------
proteomics_file <- system.file("extdata/real_data","science_signaling.csv", package = "IHWpaper")
proteomics_df <- read.csv(proteomics_file, stringsAsFactors = F)
# pvalues were adjusted by BH method so rewrite to obtain orginal pvlaues
proteomics_df$pvalue <- rank(proteomics_df$p1, ties.method="first")*proteomics_df$p1/nrow(proteomics_df)
proteomics_df$test = qnorm(proteomics_df$pvalue, lower.tail = F)
names(proteomics_df)


test = proteomics_df$test
test[test == -Inf] <- NA
test[test == Inf] <- NA
pval = proteomics_df$pvalue
filter = proteomics_df$X..peptides
hist(pval)  # one-tailed pvalue

Data <- tibble(test, pval, filter)	# data of filter covariate and pvlaues

#fit simple regression to obtain filter and test effect sizes
#-------------------------------------------------------------
nrep = 10000
m = length(pval)
#null=propTrueNull(proteomics_df$pvalue)
null = qvalue(pval, pi0.method = "bootstrap")$pi0
m0 = ceiling(null*m)
m1 = m-m0

bc2 <- boxcox(filter ~ test)
trans2 <- bc2$x[which.max(bc2$y)]
model_prot <- lm(filter^trans2 ~ test)

test_effect = if(m1 == 0) {0
} else {sort(test, decreasing = T)[1:m1]}		# two-tailed test

et_cont = mean(test_effect, na.rm = T)
ey_cont = model_prot$coef[[1]] + model_prot$coef[[2]]*et_cont
et_bin = median(test_effect, na.rm = T)
ey_bin = model_prot$coef[[1]] + model_prot$coef[[2]]*et_bin

prob_bin <-vapply(1:m, prob_rank_givenEffect, 1, et = ey_bin,
                  ey = ey_bin, nrep = nrep, m0 = m0, m1 = m1)
prob_cont <-vapply(1:m, prob_rank_givenEffect, 1, et = ey_cont,
                   ey = ey_cont, nrep = nrep, m0 = m0, m1 = m1)

alphaVec = seq(.05,.1,length.out = 5)
w_bin <- sapply(alphaVec, weight_binary, et = et_bin, m = m, m1 = m1, tail = 1,
                       delInterval = .0001, prob = prob_bin)
w_cont = sapply(alphaVec, weight_continuous, et = et_cont, m = m, tail = 1,
                           delInterval=.0001 , prob = prob_cont)


# function to compute number of rejections
# input:
#=======
# i=index number
# alpha = significance level
# Data = data to be analyzed composed of P=pvalues and Xf=filter covariate
# W_bin_mat = Binary weight matrix for diffecrent alpha
# W_cont_mat = Continuous weight matrix for different alpha
# output:
#========
# rej_mat = a rejection matrix composed of rejections from variaous methods
#---------------------------------------------------------------
fun.rejections <- function(i,alphaVec,Data,W_bin_mat,W_cont_mat)
{
    alpha=alphaVec[i]
    W_bin=as.vector(W_bin_mat[,i])
    W_cont=as.vector(W_cont_mat[,i])

    m = length(W_bin)

    OD <- Data[order(Data$filter,decreasing=T),]		# odered by covariate
    Ordered.pvalue <- OD$pval					# odered pvalues for all tests

    # preprocesing before counting rejections
    #----------------------------------------
    ihw_res_bon <- ihw(Data$pval,Data$filter, alpha=alpha, nbins=4,nsplits_internal=5,
                       lambdas=seq(0,3,length=20),adjustment_type = "bonferroni")
    padj_Pro_bin <-p.adjust(Ordered.pvalue/W_bin, method = "BH")		# proposed method based on right-tailed pvalue
    padj_Pro_cont <-p.adjust(Ordered.pvalue/W_cont, method = "BH")
    padj_BH <- p.adjust(Data$pval, method = "BH")
    ihw_res_fdr <- ihw(Data$pval,Data$filter, alpha=alpha, nbins=4,nsplits_internal=5, lambdas=seq(0,3,length=20))

    # rejections by FWER
    #-------------------
    Pro_bon_bin = sum(Ordered.pvalue <= alpha*W_bin/m, na.rm = TRUE)
    Pro_bon_cont = sum(Ordered.pvalue <= alpha*W_cont/m, na.rm = TRUE)
    bon = sum(Data$pval <= alpha/m, na.rm = TRUE)
    ihw_bon = rejections(ihw_res_bon)

    # rejections by FDR
    #------------------
    Pro_bh_bin = sum(padj_Pro_bin <= alpha, na.rm = TRUE)
    Pro_bh_cont = sum(padj_Pro_cont <= alpha, na.rm = TRUE)
    bh = sum(padj_BH <= alpha, na.rm = TRUE)
    ihw_bh = rejections(ihw_res_fdr)

    return(c(Pro_bon_bin,Pro_bon_cont,bon,ihw_bon,Pro_bh_bin,Pro_bh_cont,bh,ihw_bh))
}
rej_mat_prot = sapply(1:length(alphaVec),fun.rejections,alphaVec,Data=Data,W_bin_mat=w_bin,W_cont_mat=w_cont)


# from IHW paper-------------
proteomics_file <- system.file("real_data_examples/result_files", "proteomics_benchmark.Rds", package = "IHWpaper")
proteomics_data <- readRDS(file=proteomics_file)
panel_c_data <- group_by(proteomics_data$alpha_df, alpha) %>% summarize(BH = max(bh_rejections), IHW=max(rejections)) %>% 
  gather(method, rejections, BH, IHW)


# # nice plots-----------
# rej_mat_prot_FWER <- data.frame(alphaVec, t(rej_mat_prot[1:4,]))
# colnames(rej_mat_prot_FWER) <- c("alpha", "PRO_bin","PRO_cont","BON","IHW" )
# rej_mat_prot_FWER2 <- melt(rej_mat_prot_FWER, id.var = "alpha")

rej_mat_prot_FDR <- data.frame(alphaVec, t(rej_mat_prot[5:7,]),IHW=c(192, 216, 238, 262, 271))
colnames(rej_mat_prot_FDR) <- c("alpha", "PRO_bin","PRO_cont","BH","IHW" )
rej_mat_prot_FDR2 <- melt(rej_mat_prot_FDR, id.var = "alpha")


# p_fwer_prot <- ggplot(rej_mat_prot_FWER2, aes(x = alpha, y = value, group = variable,
#                                         colour = variable)) +
#     geom_line(aes(linetype = variable), size = 1.2) +
#     labs(x = "alpha", y = "discoveries", title = "FWER based") +
#     theme(legend.title = element_blank())
# #theme(panel.background = element_rect(fill = 'white', colour = 'black'))
# 
# #df = data.frame(colour=colors, last_vals=c(370, 382, 348, 391),
# label=c("PRO_bin","PRO_cont","BH","IHW"))
# #panel_d <- pretty_legend(p_fwer, df, .1)
# #panel_d



# final plots----------------------

p_fdr_bot <- ggplot(rej_mat_bot_FDR2, aes(x = alpha, y = value, group = variable,
                                          colour = variable)) +
  geom_line(aes(linetype = variable), size = 1.5) +
  labs(x = "alpha", y = "discoveries", title = "(a) Bottomly") +
  theme(legend.title = element_blank(), legend.position="bottom")

p_fdr_prot <- ggplot(rej_mat_prot_FDR2, aes(x = alpha, y = value, group = variable,
                                            colour = variable)) +
  geom_line(aes(linetype = variable), size = 1.5) +
  labs(x = "alpha", y = "discoveries", title = "(b) Proteomics") +
  theme(legend.position="none")


# extract the legend from one of the plots
legend_example <- get_legend(p_fdr_bot + theme(legend.direction="horizontal",
                                      legend.position="bottom"))
p_fdr_bot = p_fdr_bot + theme(legend.position="none")

# arrange the plots
#dev.new(width=8, height=4)
grid.arrange(arrangeGrob(p_fdr_bot, p_fdr_prot, nrow=1),legend_example, nrow=2, heights=c(7,1))
           





