#' @title Simulate FWER, POWER, FDR, and POWER
#'
#' @description This function simulate Family Wise Error Rate (FWER) and corresponding
#' Power, and False Discovery Rate (FDR) and the corresponding Power for different
#' effect sizes
#'
#' @param i i-th filter effect
#' @param simu number of replications
#' @param null proportion of the true null hypothesis
#' @param corr correlation between test statistics
#' @param cv determine whether the test mean effect and the filter mean effects are the same
#' @param alpha significance threshold
#' @param groupSize number of test statistics per group
#' @param effectType type of effect sizes, c("continuous", "binary")
#' @param filterEffectVec a vector of different effect size
#' @param datWeightByNull a matrix of weights, each column corresponds to a effect size
#'
#' @details
#' This function simulate Family Wise Error Rate (FWER) and corresponding
#' Power, and False Discovery Rate (FDR) and the corresponding Power for different
#' effect sizes
#'
#' @author Mohamad S. Hasan, \email{mshasan@uga.edu}
#' @export
#'
#' @seealso \code{\link{weight_byEffect_cont}}
#' \code{\link{ranksProb_byEffect}}
#'
#' @return a matrix of 16 rows containing information about FWER, POWER, FDR, and
#' POWER (4 rows for each)
#'
#' @references Hasan and Schliekelman (2017)
#'
#' @examples
#' # vector of effect sizes
#' filterEffectVec <- c(1, 1.5, 2)
#'
#' # compute probability matrix
#' ranksProb_byEffect <- sapply(1:length(filterEffectVec), ranksProb_byEffect,
#'              null = .9, m = 100, filterEffectVec = filterEffectVec)
#'
#' # compute weights
#' weightByEffect <- sapply(1:length(filterEffectVec), weight_byEffect_cont,
#'                    alpha = .05, null = .9, m = 100, delInterval = .0001,
#'                    filterEffectVec = filterEffectVec,
#'                    datByNull = ranksProb_byEffect)
#'
#' simuVal = 3  # in actual case use at least simVal = 1000
#' result <- sapply(1:length(filterEffectVec), fwerPowerFdrPower, simu = simuVal,
#'              null = .5, corr = 0, cv = 0, alpha = .05, groupSize = 100,
#'              effectType = "continuous", filterEffectVec = filterEffectVec,
#'              datWeightByNull = weightByEffect)
#'
#===============================================================================
#----------------------fwerPowerFdrPower----------------------------
# function to compute Simulated FWER, POWER, and FDR by effect size
#
# inpout:----------------
# i = i-th filter effect
# simu = number of replications
# null = proportion of the true null hypothesis
# corr = correlation between test statistics
# cv = determine whether the test mean effect and the filter mean effects are the same
# random =1 mean cv=1/2
# alpha = significance threshold
# groupSize = number of test statistics per group
# filterEffectVec = a vector of different effect size
# datWeightByNull = a matrix of weights, each column corresponds to a effect size
#
# internal parameters:-----
# W = weight vector for a specific effect size
# m = test size
# weight_pro = normalizing proposed weight
# ey = filter effect size
# m0 = no. of null hypothesis
# m1 = no. of alt hyp.
# xf = only alt. filter effect vector
# xt = alt. test effect vector
# Sigma = test correlation matrix
# H = alternative hypothesis true or false
# ef = filter effect vector (mixture of null and alt)
# et = test effect vector (mixture of null and alt)
# mGrp = subgroup of tests.
# test = filter test stat
# filter = actual test stat
# pval = filter test pvalues
# pro = proposed, bon = bonferroni, rdw = roeder and wasserman, IHW=  independent Hyp. Weight
#
# output:---------------
#  a matrix of 16 rows containing information about FWER, POWER, FDR, and
#  POWER (4 rows for each in respective order)
#
#===============================================================================

fwerPowerFdrPower <- function(i, simu, null, corr = 0, cv = 0, alpha = .05,
                                        groupSize = 100, effectType = c("continuous", "binary"),
                                        filterEffectVec, datWeightByNull)
    {
        W = datWeightByNull[ , i]
        m = length(W)
        weight_pro <- if(sum(W)==0){rep(1, m)} else {W/sum(W)*m}
        ey <- filterEffectVec[i]
        m0 <- ceiling(m*null)
        m1 <- m - m0

        if(effectType == "continuous"){
            xf <- as.vector(runif_by_mean(n = m, mean = ey))
        } else {
            xf <- rep(ey, m)
        }

        xt <- if(cv == 0){xf} else {rnorm(m, ey, cv*ey)}
        Sigma <- matrix(corr, groupSize, groupSize) + diag(groupSize)*(1 - corr)


        fwerPowerFdrPower_simu <- function(s)
            {
                H <- rbinom(m, 1, 1 - null)
                ef <- H*xf
                et <- H*xt
                mGrp = m/groupSize

                test <- if(corr == 0) {rnorm(m, et, 1)
                } else {as.vector(sapply(1:mGrp, test_by_block, eVec = et,
                                         groupSize = groupSize, Sigma = Sigma))}

                filter <- if(corr == 0) {rnorm(m, ef, 1)
                } else {as.vector(sapply(1:mGrp, test_by_block, eVec = ef,
                                         groupSize = groupSize, Sigma = Sigma))}

                pval <- pnorm(test, lower.tail = FALSE)

                dat = tibble(test, pval, et, filter)
                OD = dat[order(dat$filter, decreasing = TRUE), ]

                weight_rdw <- roeder_wasserman_weight(pvalue = pval, filter = filter,
                                        gamma = .05, alpha = alpha, group = 5L, tail = 1L)
                ihw_fwer <- ihw(OD$pval, OD$filter, alpha = alpha, adjustment_type = "bonferroni")
                ihw_fdr <-  ihw(OD$pval, OD$filter, alpha = alpha, adjustment_type = "BH")

                rej_pro <- OD$pval <= alpha*weight_pro/m
                rej_bon <- OD$pval <= alpha/m
                rej_rdw <- OD$pval <= alpha*weight_rdw/m
                rej_ihwFwer <- adj_pvalues(ihw_fwer) <= alpha

                n_null <- max(1, sum(OD$et == 0, na.rm = TRUE))
                n_alt <-  max(1, sum(OD$et != 0, na.rm = TRUE))

                FWER_pro <- sum(rej_pro[OD$et == 0])
                FWER_bon <- sum(rej_bon[OD$et == 0])
                FWER_rdw <- sum(rej_rdw[OD$et == 0])
                FWER_ihw <- sum(rej_ihwFwer[OD$et == 0])

                POWER_pro <- sum(rej_pro[OD$et != 0])/n_alt
                POWER_bon <- sum(rej_bon[OD$et != 0])/n_alt
                POWER_rdw <- sum(rej_rdw[OD$et != 0])/n_alt
                POWER_ihw <- sum(rej_ihwFwer[OD$et != 0])/n_alt

                adjPval_pro <- p.adjust(OD$pval/weight_pro, method="BH")
                adjPval_bon <- p.adjust(OD$pval, method="BH")
                adjPval_rdw <- p.adjust(OD$pval/weight_rdw, method="BH")
                adjPval_ihw <- adj_pvalues(ihw_fdr)

                FDR_pro <- sum(adjPval_pro[OD$et == 0] <= alpha)/max(1, sum(adjPval_pro <= alpha))
                FDR_bh  <- sum(adjPval_bon[OD$et == 0] <= alpha)/max(1, sum(adjPval_bon <= alpha))
                FDR_rdw <- sum(adjPval_rdw[OD$et == 0] <= alpha)/max(1, sum(adjPval_rdw <= alpha))
                FDR_ihw <- sum(adjPval_ihw[OD$et == 0] <= alpha)/max(1, rejections(ihw_fdr))

                FDR_POWER_pro <- sum(adjPval_pro[OD$et != 0] <= alpha)/n_alt
                FDR_POWER_bh  <- sum(adjPval_bon[OD$et != 0] <= alpha)/n_alt
                FDR_POWER_rdw <- sum(adjPval_rdw[OD$et != 0] <= alpha)/n_alt
                FDR_POWER_ihw <- sum(adjPval_ihw[OD$et != 0] <= alpha)/n_alt

                return(c(FWER_pro, FWER_bon, FWER_rdw, FWER_ihw,
                         POWER_pro, POWER_bon, POWER_rdw, POWER_ihw,
                         FDR_pro, FDR_bh, FDR_rdw, FDR_ihw, FDR_POWER_pro,
                         FDR_POWER_bh, FDR_POWER_rdw, FDR_POWER_ihw))
            }

        fwerPowerFdrPower_bysimu <- sapply(1:simu, fwerPowerFdrPower_simu)
        fwerPowerFdrPower <- apply(fwerPowerFdrPower_bysimu, 1, mean, na.rm=TRUE)

        return(fwerPowerFdrPower)

    }





