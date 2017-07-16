#' @title Weight for the different continuous effect sizes
#'
#' @description Compute weight from the the ranks probability given the effect
#' sizes when the effect sizes are continuous.
#'
#' @param i Integer, i-th effect size of a vector of effects
#' @param alpha Numeric, a significance level of the hypothesis tests
#' @param null Numeric, proportion of true true null tests
#' @param m Integer, total number of hypotheis tests
#' @param tail Integer (1 or 2), right-tailed or two-tailed hypothesis test.
#' @param delInterval Numeric, interval between the \code{delta} values of a
#' sequence. Note that, \code{delta} is a Langrangian multiplier, necessary to
#' normalize the weight
#' @param covariateEffectVec A numeric vector of covariate effect sizes
#' @param datByNull A numerics matrix of ranks pobabilities in which each column
#' corresponds to an effect size
#'
#' @details This function compute the weights when the effect sizes are
#' continuous by applying the ranks probabilities of the different effect sizes.
#' It applies the function function \code{weight_continuous} from the R package
#' \code{OPWeight} to compute the weights from a ranks probability matirx.
#'
#' @author Mohamad S. Hasan, \email{shakilmohamad7@gmail.com}
#' @export
#'
#' @import stats
#'
#' @seealso \code{\link{ranksProb_byEffect}}
#' \code{\link{weight_continuous}}
#'
#' @return A numeics matrix of weights in which each column corresponds to an
#' effect size
#'
#' @examples
#' # vector of effect sizes
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
#===============================================================================
# function to compute  weight for the continuous case
#-----------------------------------------------------------------------
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
#===============================================================================

weight_byEffect_cont <- function(i, alpha, null, m, tail = 1L, delInterval,
                                 covariateEffectVec, datByNull)
{					#input:effect index; output: weight
    et <- covariateEffectVec[i]
    m0 <- ceiling(m*null)
    m1 <- m - m0
    ranksProb <- datByNull[ , i]
    prob <- ranksProb/sum(ranksProb, na.rm = TRUE)
    delta <- seq(0, 1, delInterval)

    weightSumVec <- sapply(delta, weight_by_delta, alpha = alpha, et = et, m = m,
                           m1 = NULL, tail = tail, ranksProb = prob,
                           effectType = "continuous")

    deltaOut <- delta[min(abs(weightSumVec - m)) == abs(weightSumVec - m)]
    deltaOut <- ifelse(length(deltaOut) > 1, .0001, deltaOut)
    weight.out <- tail*(m/alpha)*pnorm(et/2 + 1/et*log(deltaOut/(alpha*prob)),
                                       lower.tail=FALSE)
    sumWeight <- sum(weight.out, na.rm = TRUE)
    normWeight <- if(sumWeight == 0) {rep(1, m)} else {weight.out/sumWeight*m}
    return(normWeight)
}














