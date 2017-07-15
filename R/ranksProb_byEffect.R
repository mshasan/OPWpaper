#' @title Ranks probability of the tests given the effect size
#'
#' @description Comnpute the probability of the ranks for different effect sizes
#'
#' @param i Integer, i-th effect size
#' @param null Numeric,proportion of the true null hypothesis
#' @param m Integer, total number of hypothesis test
#' @param nrep Integer, number of replications for the importance sampling
#' @param filterEffectVec A numeric vector of filter effect sizes
#'
#' @details This function compute ranks probabilities for the different effect
#' sizes. It apply the function \code{prob_ranks_givenEffect} from the
#' \code{OPWeight} package and comute the probabilities.
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
#' filterEffectVec <- c(1, 1.5, 2)
#'
#' # compute ranks probability matrix
#' ranksProb_byEffect <- sapply(1:length(filterEffectVec), ranksProb_byEffect,
#'              null = .9, m = 100, filterEffectVec = filterEffectVec)
#'
#===============================================================================
# internal parameters:----------
# et = effect of the targeted test for importance sampling
# ey = mean filter efffect from external information
# m0 = number of true null hypothesis
# m1 = number of true alternative hypothesis
#
#===============================================================================

ranksProb_byEffect <- function(i, null, m, nrep = 10000, filterEffectVec)
{
    ey <- filterEffectVec[i]
    et <- ey
    m0 <- ceiling(m*null)
    m1 <- m - m0
    rankProbsH1 <- sapply(1:m, prob_rank_givenEffect, et = ey, ey = ey,
                          nrep = nrep, m0 = m0, m1 = m1)
    return(rankProbsH1)
}


