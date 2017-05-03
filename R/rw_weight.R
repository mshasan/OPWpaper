#' @title Weight from the Roeder and Waserman (2009) paper
#'
#' @description Compute weight by splitting test statistics. This method is taken
#' Roeder and Waserman (2009) paper.
#' @param testStat test statistics
#' @param gamma smooting parameter
#' @param alpha Significance level
#' @param group number of group
#' @param tail one-tailed or two-tailed hypothesis test
#'
#' @details
#' NOne
#' @author Mohamad S. Hasan, \email{mshasan@uga.edu}
#' @export
#' @seealso \code{\link{weight_binary}} \code{\link{weight_continuous}}
#' @return \code{weight} normalized weights of the tests
#' @references Roeder, Kathryn, and Larry Wasserman. "Genome-wide significance
#' levels and weighted hypothesis testing." Statistical science: a review journal
#' of the Institute of Mathematical Statistics 24.4 (2009): 398.
#' \url{www.stat.cmu.edu/~roeder/publications/statsci.pdf}
#' @examples
#' # generate test statistics
#' testStat <- rnorm(100000,2,1)
#'
#' # Compute wiehgt
#' weight = rw_weight(testStat=testStat, gamma=.05, alpha=.05, group=10, tail=2)
#'
#' # plot the weight
#' plot(testStat, weight)
#'
#===============================================================================
# Function to estiamte weight from data by using Roeder-Waserman algorithm

# Input:-----
# testStat = test statistics
# gamma = smooting parameter
# alpha = Significance level
# group = number of group
# tail = one-tailed or two-tailed hypothesis test

# internal parameters:-----
# testGroup = patitioning 'm' tests into groups
# testMeans = calculate mean of each group
# testSd = calculate standard deviation of each group
# pi_hat = estiamte pi parameter
# effect_hat = estiamte effect size
# theta = vector of the effect size
# findc = find log-constant of the formula so that average weight equals to 1
# weight_k = calculate weight
# weight_k_smooth = smoothing weight by smotthing parameter gamma
# weight_tests = distribute weight over all tests
# normWeight.w = little adjustment to obtain sum of weight equals to m

# output:-----
# normWeight.w = normalized weight
#---------------------------------------------------------------------------------------------
rw_weight <- function(testStat, gamma, alpha, group=5L, tail=2L)
{
    m = length(testStat)
    groupSize <- m/group
    testGroup <- rep(1:group, each = groupSize)
    testMeans <- tapply(testStat, testGroup, mean)
    testSd <- tapply(testStat, testGroup, sd)
    pi_hat <- testMeans^2/(testMeans^2 + testSd^2 - 1)
    effect_hat <- testMeans/pi_hat
    effect_hat[pi_hat <= 1/groupSize] <- 0

    if(sum(effect_hat, na.rm=T) == 0) {normWeight.w <- rep(1,m)
    }
    else{
        theta = as.vector(effect_hat)
        c0 = qnorm(alpha/(tail*m), lower.tail = FALSE)^2/2

        findc = function(thet, alpha, c0, m)
        {
            cc = seq(.1, c0, .005)
            tot <- sapply(cc, function(c) return(sum((m/alpha)*pnorm((thet/2 + c/thet),
                                                                     lower.tail = FALSE))))
            cout = cc[min(abs(tot - m)) == abs(tot - m)]
            coutAdj = ifelse(length(cout) > 1, sample(cout, 1), cout)
            return(coutAdj)
        }

        c <- vapply(theta, findc, 1, alpha, c0, m)
        weight_k <- (m/alpha)*pnorm((effect_hat/2 + c/effect_hat), lower.tail = FALSE)
        weight_k_smooth <- (1-gamma)*weight_k + gamma*sum(weight_k)/group
        weight_tests <- rep(weight_k_smooth, each=groupSize)
        normWeight.w <- weight_tests/sum(weight_tests)*m
    }
    return(normWeight.w)
}

