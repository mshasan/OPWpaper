#' @title Probability of rank of test given effect size
#'
#' @description Comnpute the probability of the ranks for different effect sizes
#'
#' @param i i-th effect size
#' @param null proportion of the true null hypothesis
#' @param m total number of hypothesis test
#' @param nrep number of replication for importance sampling
#' @param filterEffectVec a vector of filter effect sizes
#'
#' @details This function compute ranks probabilities for different effect sizes.
#' It apply the function \code{prob_ranks_givenEffect} from the \code{OPWeight}
#' package and comute the probabilities.
#'
#' @author Mohamad S. Hasan and Paul Schliekelman
#' @export
#' @import stats
#' @seealso \code{\link{prob_rank_givenEffect}}
#'
#' @return a matrix of ranks pobabilities each column corresponds to an effect size
#' @examples
#' # vector of effect sizes
#' filterEffectVec <- c(1, 1.5, 2)
#'
#' # compute probability matrix
#' ranksProb_byEffect <- sapply(1:length(filterEffectVec), ranksProb_byEffect,
#'              null = .9, m = 100, filterEffectVec = filterEffectVec)
#'
#===============================================================================
# function to compute "P(rank|effect) by effect size for m=10,000"
#
# Input:----------
# i = i-th effect
# null =  proportion of null hypothesis
# m = test size
# nrep = number of replication for importance sampling
# filterEffectVec = a vector of filter effect
#
# internal parameters:----------
# et = effect of the targeted test for importance sampling
# ey = mean filter efffect from external information
# m0 = number of true null hypothesis
# m1 = number of true alternative hypothesis
#
#output:-----------
# P(rank | effect = ey) for the alternaive case,
# because we need only alterantive case to compute weight
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


