#' @title Simulate FWER, POWER, FDR, and POWER
#'
#' @description This function simulate the Family Wise Error Rate (FWER) and the
#' corresponding Power, and the False Discovery Rate (FDR) and the corresponding
#' Power for the different effect sizes
#'
#' @param i Integer, i-th effect size of a vector of effects
#' @param simu Integer, number of replications
#' @param null Numeric, proportion of the true null hypothesis
#' @param corr Numeric, correlation between the test statistics
#' @param cv Numeric, coefficient of variation of the test statistics
#' @param alpha Numeric value of the significance threshold
#' @param groupSize Integer, number of test statistics per group
#' @param effectType Character ("continuous" or "binary"), type of effect sizes
#' @param covariateEffectVec A numeric vector of different covariate-effect size
#' @param datWeightByNull A numeric matrix of weights, each column corresponds
#' to an covariate-effect size
#'
#' @details
#' This function simulate Family Wise Error Rate (FWER) and corresponding
#' Power, and False Discovery Rate (FDR) and the corresponding Power for the
#' different effect sizes
#'
#' @author Mohamad S. Hasan, \email{shakilmohamad7@gmail.com}
#' @export
#'
#' @seealso \code{\link{weight_byEffect_cont}}
#' \code{\link{ranksProb_byEffect}}
#'
#' @return A matrix of 16 rows containing information about FWER, POWER, FDR,
#' and POWER (4 rows for each item)
#'
#' @references Hasan and Schliekelman (2017)
#'
#' @examples
#' # vector of covariate-effect sizes
#' covariateEffectVec <- c(1, 1.5, 2)
#'
#' # compute probability matrix
#' ranksProb_byEffect <- sapply(1:length(covariateEffectVec), ranksProb_byEffect,
#'              null = .9, m = 100, covariateEffectVec = covariateEffectVec)
#'
#' # compute weights
#' weightByEffect <- sapply(1:length(covariateEffectVec), weight_byEffect_cont,
#'                    alpha = .05, null = .9, m = 100, delInterval = .01,
#'                    covariateEffectVec = covariateEffectVec,
#'                    datByNull = ranksProb_byEffect)
#'
#' simuVal = 3  # in the actual case use at least simVal = 1000
#' result <- sapply(1:length(covariateEffectVec), fwerPowerFdrPower, simu = simuVal,
#'              null = .5, corr = 0, cv = 0, alpha = .05, groupSize = 100,
#'              effectType = "continuous", covariateEffectVec = covariateEffectVec,
#'              datWeightByNull = weightByEffect)
#'
#===============================================================================
#----------------------fwerPowerFdrPower----------------------------
# internal parameters:-----
# W = weight vector for a specific effect size
# m = test size
# weight_CRW = normalizing proposed weight
# ey = covariate effect size
# m0 = no. of null hypothesis
# m1 = no. of alt hyp.
# xf = only alt. covariate effect vector
# xt = alt. test effect vector
# Sigma = test correlation matrix
# H = alternative hypothesis true or false
# ef = covariate effect vector (mixture of null and alt)
# et = test effect vector (mixture of null and alt)
# mGrp = subgroup of tests.
# test = covariate test stat
# covariate = actual test stat
# pval = covariate test pvalues
# CRW = proposed, bon = bonferroni, rdw = roeder and wasserman,
# IHW=  independent Hyp. Weight
#===============================================================================

fwerPowerFdrPower <- function(i, simu, null, corr = 0, cv = 0, alpha = .05,
                        groupSize = 100, effectType = c("continuous", "binary"),
                        covariateEffectVec, datWeightByNull)
{
    W = datWeightByNull[ , i]
    m = length(W)
    weight_CRW <- if(sum(W)==0){rep(1, m)} else {W/sum(W)*m}
    ey <- covariateEffectVec[i]
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

        covariate <- if(corr == 0) {rnorm(m, ef, 1)
        } else {as.vector(sapply(1:mGrp, test_by_block, eVec = ef,
                                 groupSize = groupSize, Sigma = Sigma))}

        pval <- pnorm(test, lower.tail = FALSE)

        dat = tibble(test, pval, et, covariate)
        OD = dat[order(dat$covariate, decreasing = TRUE), ]

        weight_rdw <- roeder_wasserman_weight(pvalue = OD$pval, alpha = alpha)
        ihw_fwer <- ihw(OD$pval, OD$covariate, alpha = alpha, adjustment_type = "bonferroni")
        ihw_fdr <-  ihw(OD$pval, OD$covariate, alpha = alpha, adjustment_type = "BH")

        rej_CRW <- OD$pval <= alpha*weight_CRW/m
        rej_bon <- OD$pval <= alpha/m
        rej_rdw <- OD$pval <= alpha*weight_rdw/m
        rej_ihwFwer <- adj_pvalues(ihw_fwer) <= alpha

        n_null <- max(1, sum(OD$et == 0, na.rm = TRUE))
        n_alt <-  max(1, sum(OD$et != 0, na.rm = TRUE))

        FWER_CRW <- sum(rej_CRW[OD$et == 0])
        FWER_bon <- sum(rej_bon[OD$et == 0])
        FWER_rdw <- sum(rej_rdw[OD$et == 0])
        FWER_ihw <- sum(rej_ihwFwer[OD$et == 0])

        POWER_CRW <- sum(rej_CRW[OD$et != 0])/n_alt
        POWER_bon <- sum(rej_bon[OD$et != 0])/n_alt
        POWER_rdw <- sum(rej_rdw[OD$et != 0])/n_alt
        POWER_ihw <- sum(rej_ihwFwer[OD$et != 0])/n_alt

        adjPval_CRW <- p.adjust(OD$pval/weight_CRW, method="BH")
        adjPval_bon <- p.adjust(OD$pval, method="BH")
        adjPval_rdw <- p.adjust(OD$pval/weight_rdw, method="BH")
        adjPval_ihw <- adj_pvalues(ihw_fdr)

        FDR_CRW <- sum(adjPval_CRW[OD$et == 0] <= alpha)/max(1, sum(adjPval_CRW <= alpha))
        FDR_bh  <- sum(adjPval_bon[OD$et == 0] <= alpha)/max(1, sum(adjPval_bon <= alpha))
        FDR_rdw <- sum(adjPval_rdw[OD$et == 0] <= alpha)/max(1, sum(adjPval_rdw <= alpha))
        FDR_ihw <- sum(adjPval_ihw[OD$et == 0] <= alpha)/max(1, rejections(ihw_fdr))

        FDR_POWER_CRW <- sum(adjPval_CRW[OD$et != 0] <= alpha)/n_alt
        FDR_POWER_bh  <- sum(adjPval_bon[OD$et != 0] <= alpha)/n_alt
        FDR_POWER_rdw <- sum(adjPval_rdw[OD$et != 0] <= alpha)/n_alt
        FDR_POWER_ihw <- sum(adjPval_ihw[OD$et != 0] <= alpha)/n_alt

        return(c(FWER_CRW, FWER_bon, FWER_rdw, FWER_ihw,
                 POWER_CRW, POWER_bon, POWER_rdw, POWER_ihw,
                 FDR_CRW, FDR_bh, FDR_rdw, FDR_ihw, FDR_POWER_CRW,
                 FDR_POWER_bh, FDR_POWER_rdw, FDR_POWER_ihw))
    }

    fwerPowerFdrPower_bysimu <- sapply(1:simu, fwerPowerFdrPower_simu)
    fwerPowerFdrPower <- rowMeans(fwerPowerFdrPower_bysimu, na.rm=TRUE)

    return(fwerPowerFdrPower)

}





