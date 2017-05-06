#' @title Weight for different effect sizes in binary case
#'
#' @description Compute weight from the probability of the rank given the effect
#' size for different effect sizes in a binary effect size situation
#'
#' @param i i-th effect size
#' @param alpha significance level
#' @param null proportion of true true null tests
#' @param m total number of hypotheis test
#' @param tail one-tailed or two-tailed hypothesis test
#' @param delInterval interval between the \code{delta} values of a sequence.
#' @param filterEffectVec a vector of filter effect sizes
#' Note that, \code{delta} is a LaGrange multiplier, necessary to normalize the weight
#' @param datByNull a matrix of ranks pobabilities each column corresponds to
#'                    an effect size
#'
#' @details This function compute the weights in a binary settings by applying
#' the ranks probabilities for the different effect sizes. It applies the function
#' function \code{weight_binary} from the package 'OPWeight' to
#' comute the weights from a probability matirx.
#'
#' @author Mohamad S. Hasan, mshasan@uga.edu
#' @export
#'
#' @import stats
#'
#' @seealso \code{\link{ranksProb_byEffect}}
#' \code{\link{weight_binary}}
#'
#' @return a matrix of weights each column corresponds to an effect size
#' @examples
#' # vector of effect sizes
#' filterEffectVec <- c(1, 1.5, 2)
#'
#' # compute probability matrix
#' ranksProb_byEffect <- sapply(1:length(filterEffectVec), ranksProb_byEffect,
#'              null = .9, m = 100, filterEffectVec = filterEffectVec)
#'
#' # compute weights
#'weightByEffect <- sapply(1:length(filterEffectVec), weight_byEffect_bin,
#'                    alpha = .05, null = .9, m = 100, delInterval = .0001,
#'                    filterEffectVec = filterEffectVec,
#'                    datByNull = ranksProb_byEffect)
#'
#===============================================================================
# function to compute  weight for the binary case
#-----------------------------------------------------------------------
#
# Input:-----------
# i = i-th effect
# alpha = significance level
# null = proportion of true null test
# m = test size
# tail = one-tailed or two-tailed hypothesis test
# delInterval = increasing rate of the delta (lagrange multiplier)
# datByNull = P(rank|effect) data generated before will call here
#
# internal parameters:-----
# et = effect of the targeted test for importance sampling
# m0 = number of true null hypothesis
# m1 = number of true alternative hypothesis
# prob = probability of rank given the only effect size
# delta = sequene of delta (lagrange multiplier) values
# findDelta = function to compute sum of weight for each dleta
# deltaOut = optimal delta value
# sumWeight = sum of the weights
#
# Output:-----------
# normWeight = a matrix of weights each column corresponds to an effect size
#===============================================================================

weight_byEffect_bin <- function(i, alpha, null, m, tail = 1L, delInterval,
                                filterEffectVec, datByNull)
    {
        et <- filterEffectVec [i]
        m0 <- ceiling(m*null)
        m1 <- m-m0
        prob <- datByNull[,i]
        prob <- prob/sum(prob, na.rm = TRUE)
        delta <- seq(0, 1, delInterval)
        findDelta <- function(delta)
        {
            weight <- tail*(m/alpha)*pnorm(et/2 + 1/et*log(delta*m/(alpha*m1*prob)),
                                           lower.tail = FALSE)
            return(sum(weight, na.rm = TRUE))
        }
        weightSumVec <- vapply(delta, findDelta, 1)
        deltaOut <- delta[min(abs(weightSumVec - m)) == abs(weightSumVec - m)]
        deltaOut <- ifelse(length(deltaOut) > 1, .0001, deltaOut)
        weight.out <- tail*(m/alpha)*pnorm(et/2 + 1/et*log(deltaOut*m/(alpha*m1*prob)),
                                           lower.tail = FALSE)
        sumWeight <- sum(weight.out, na.rm = TRUE)
        normWeight <- if(sumWeight == 0) {rep(1, m)} else {weight.out/sumWeight*m}
        return(normWeight)
    }


