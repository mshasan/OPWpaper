#' @title Simulate Family Wise Error Rate (FWER)
#'
#' @description This function simulate family wise error rate or test type I
#' error rate
#'
#' @param s Integer, number of replications in a simulation
#' @param m Integer, total number of hypothesis test
#' @param alphaVec A numeric vector of the significance levels
#'
#' @details
#' This function generate pvalues from the \code{uniform(0, 1)} and then
#' apply the function \code{opw} from the R package \code{OPWeight} to obtain
#' the Familly Wise Error Rate (FWER)
#'
#' @author Mohamad S. Hasan, \email{shakilmohamad7@gmail.com}
#' @export
#'
#' @seealso \code{\link{qvalue}}
#' \code{\link{prob_rank_givenEffect}}
#' \code{\link{weight_binary}}
#' \code{\link{weight_continuous}}

#'
#' @return A numeric matrix of FWER for the different p-value weighting methods
#'
#' @references Hasan and Schliekelman (2017)
#'
#' @examples
#' alphaVec = .05
#' simVal = 1:3  # in actual case use at least simVal = 1000
#' typeIerror_mat = sapply(simVal, simu_fwer, m = 100, alphaVec = alphaVec)
#'
#===============================================================================
# internal parameters:-----
# pval = pvalues from null tests
# pval_covariate = covariate pvalues from null tests
# test = test statistics
# covariate = covariate test statistics
# dat = a data frame
# OD = ordered data by covariate
# odered.pvalue = ordered pvalue by covariate
# nullprop = prportion of null
# m0 = true null test size
# m1 = true alternative test size
# test_effect =  estimated true alternative test effects
# prob_bin = binary ranks probablity
# prob_cont = continuous ranks probability
# w_bin = binary weight
# w_cont = continuous weight
#
#===============================================================================

simu_fwer <- function(s, m, alphaVec)
    {
    fwer_per_rep <- function(alpha)
        {
            pval <- runif(m)
            pval_covariate <- runif(m)
            test = qnorm(pval, lower.tail = FALSE)
            covariate = qnorm(pval_covariate, lower.tail = FALSE)

            dat = tibble(test, pval, covariate)

            OD = dat[order(dat$covariate, decreasing=TRUE), ]
            odered.pvalue = OD$pval

            nullprop = qvalue(pval)$pi0
            m0 = ceiling(m*nullprop)
            m1 = m - m0

            model = lm(covariate ~ test)

            test_effect <- if(m1 == 0) {0
                           } else {sort(test, decreasing = TRUE)[1:m1]}

            et_bin = median(test_effect, na.rm = TRUE)
            et_cont = mean(test_effect, na.rm = TRUE)

            ey_bin = model$coef[[1]] + model$coef[[2]]*et_bin
            ey_cont = model$coef[[1]] + model$coef[[2]]*et_cont

            prob_bin <-sapply(1:m, prob_rank_givenEffect, et = ey_bin,
                              ey = ey_bin,m0 = m0, m1 = m1)
            prob_cont <-sapply(1:m, prob_rank_givenEffect, et = ey_cont,
                               ey = ey_cont, m0 = m0, m1 = m1)

            w_bin <- weight_binary(alpha = alpha, et = et_bin, m = m, m1 = m1,
                            tail = 1, delInterval = .0001, ranksProb = prob_bin)
            w_cont = weight_continuous(alpha = alpha, et = et_cont, m = m,
                            tail = 1, delInterval = .0001 , ranksProb = prob_cont)

            ihw_fwer <- ihw(dat$pval, dat$covariate, alpha = alpha,
                                            adjustment_type = "bonferroni")

            bon = sum(pval <= alpha/m, na.rm = TRUE)
            CRW_bin = sum(odered.pvalue <= alpha*w_bin/m, na.rm = TRUE)
            CRW_cont = sum(odered.pvalue <= alpha*w_cont/m, na.rm = TRUE)
            IHW <- rejections(ihw_fwer)

            return(c(bon, CRW_bin, CRW_cont, IHW))
        }

        fwer_per_rep_mat = sapply(alphaVec, fwer_per_rep)
        return(fwer_per_rep_mat)
    }







