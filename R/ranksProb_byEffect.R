#' @title Ranks probability of the covariate given the test-effect size
#'
#' @description Comnpute the ranks probability for different test-effect sizes
#'
#' @param i i Integer, i-th effect size of a vector of effects
#' @param null Numeric, proportion of the true null hypothesis
#' @param m Integer, total number of hypothesis test
#' @param nrep Integer, number of replications for the importance sampling
#' @param covariateEffectVec A numeric vector of the covariate-effect sizes
#'
#' @details This function compute ranks probabilities for the different effect
#' sizes. It apply the function \code{prob_ranks_givenEffect} from the
#' \code{OPWeight} package and compute the probabilities.
#'
#' @author Mohamad S. Hasan, \email{shakilmohamad7@gmail.com}
#' @export
#' @import stats
#' @seealso \code{\link{prob_rank_givenEffect}}
#'
#' @return A numeric matrix of the ranks pobabilities in which each column
#' corresponds to an effect size
#'
#' @examples
#' # vector of effect sizes
#' covariateEffectVec <- c(1, 1.5, 2)
#'
#' # compute ranks probability matrix
#' ranksProb_byEffect <- sapply(1:length(covariateEffectVec), ranksProb_byEffect,
#'              null = .9, m = 100, covariateEffectVec = covariateEffectVec)
#'
#===============================================================================
# internal parameters:----------
# et = effect of the targeted test for importance sampling
# ey = mean covariate efffect from external information
# m0 = number of true null hypothesis
# m1 = number of true alternative hypothesis
#
#===============================================================================

ranksProb_byEffect <- function(i, null, m, nrep = 10000, covariateEffectVec)
{
    ey <- covariateEffectVec[i]
    et <- ey
    m0 <- ceiling(m*null)
    m1 <- m - m0
    rankProbsH1 <- sapply(1:m, prob_rank_givenEffect, et = ey, ey = ey,
                          nrep = nrep, m0 = m0, m1 = m1)
    return(rankProbsH1)
}


