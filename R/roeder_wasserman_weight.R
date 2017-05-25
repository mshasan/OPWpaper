#' @title Compute sum of normal cummulative values for different constant
#'
#' @description Compute sum of normal cummulative values for the different constant, c.
#'  This constant c is defined in the paper Roeder and Waserman (2009).
#'
#' @param c a constant value
#' @param alpha Significance level
#' @param thet effect size
#' @param m number of hypothesis tests
#'
#' @details
#' NOne
#' @author Mohamad S. Hasan, \email{mshasan@uga.edu}
#'
#' @export
#'
#' @references Roeder, Kathryn, and Larry Wasserman. "Genome-wide significance
#' levels and weighted hypothesis testing." Statistical science: a review journal
#' of the Institute of Mathematical Statistics 24.4 (2009): 398.
#' \url{www.stat.cmu.edu/~roeder/publications/statsci.pdf}
#'
#' @examples
#' cc = seq(.1, 1, .005)
#' tot <- sapply(cc, fun_c, alpha = .05, thet = .3, m = 100)
#===============================================================================

fun_c <- function(c, alpha = .05, thet, m)
{
    cdf_sum <- sum((m/alpha)*pnorm((thet/2 + c/thet), lower.tail = FALSE))
    return(cdf_sum)
}

#===============================================================================

#' @title Compute constant, c, which makes the sum of weights equal to
#' the number of tests
#'
#' @description Compute constant, c, which makes the sum of weights equal to
#' the number of tests
#'
#' @param thet effect size
#' @param alpha significance level
#' @param c0 square of the tests' critical value
#' @param m number of hypothesis tests
#'
#' @details
#' NOne
#' @author Mohamad S. Hasan, \email{mshasan@uga.edu}
#'
#' @export
#'
#' @references Roeder, Kathryn, and Larry Wasserman. "Genome-wide significance
#' levels and weighted hypothesis testing." Statistical science: a review journal
#' of the Institute of Mathematical Statistics 24.4 (2009): 398.
#' \url{www.stat.cmu.edu/~roeder/publications/statsci.pdf}
#'
#' @examples
#' theta = seq(0, 3, .5)
#' c <- sapply(theta, findc, alpha = .05, c0 = 1, m = 100)
#'
#===============================================================================

findc = function(thet, alpha = .05, c0, m)
{
    cc = seq(.1, c0, .005)
    tot <- sapply(cc, fun_c, alpha = alpha, thet = thet, m = m)

    cout = cc[min(abs(tot - m)) == abs(tot - m)]
    coutAdj = ifelse(length(cout) > 1, sample(cout, 1), cout)

    return(coutAdj)
}

#===============================================================================

#' @title Weight from the Roeder and Wasserman (2009) paper
#'
#' @description Compute weights by splitting test statistics raked by the filter statistics.
#'  This method is taken Roeder and Waserman (2009) paper.
#' @param pvalue a vector of pvalues of the test statistics
#' @param filter a vector of filter statistics
#' @param gamma smooting parameter
#' @param alpha Significance level
#' @param group number of groups
#' @param tail right-tailed or two-tailed hypothesis test. default is right-tailed test
#'
#' @details
#' NOne
#' @author Mohamad S. Hasan, \email{mshasan@uga.edu}
#' @export
#'
#' @import tibble tibble
#'
#' @seealso \code{\link{weight_binary}} \code{\link{weight_continuous}}
#' @return \code{weight} normalized weights of the tests
#' @references Roeder, Kathryn, and Larry Wasserman. "Genome-wide significance
#' levels and weighted hypothesis testing." Statistical science: a review journal
#' of the Institute of Mathematical Statistics 24.4 (2009): 398.
#' \url{www.stat.cmu.edu/~roeder/publications/statsci.pdf}
#' @examples
#'
#' #generate pvalues and filter statistics
#' m = 10000
#' set.seed(123)
#' filters = runif(m, min = 0, max = 2.5)          # filter statistics
#' H = rbinom(m, size = 1, prob = 0.1)             # hypothesis true or false
#' tests = rnorm(m, mean = H * filters)            # Z-score
#' pvals = 1 - pnorm(tests)                        # pvalue
#'
#' # Compute wiehgt
#' weight = roeder_wasserman_weight(pvalue = pvals, filter = filters)
#'
#' # plot the weight
#' plot(weight)
#'
#' # compute number of rejections
#' Data <- tibble(tests, pvals, filters)
#' OD <- Data[order(Data$filters, decreasing = TRUE), ]
#' alpha = .05
#' rwd <- sum(OD$pvals <= alpha*weight/m)
#' bon <- sum(pvals <= alpha/m)
#'
#===============================================================================
# Function to estiamte weight from data by using Roeder-Waserman algorithm

# Input:-----
# tests = test statistics
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
# cr_val = critical value
# findc = find log-constant of the formula so that average weight equals to 1
# weight_k = calculate weight
# weight_k_smooth = smoothing weight by smotthing parameter gamma
# weight_tests = distribute weight over all tests
# normWeight.w = little adjustment to obtain sum of weight equals to m

# output:-----
# normWeight.w = normalized weight
#-------------------------------------------------------------------------------

roeder_wasserman_weight <- function(pvalue, filter, gamma = .05, alpha = .05,
                                    group = 5L, tail = 1L)
    {
        m = length(pvalue)
        tests <- qnorm(pvalue/tail, lower.tail = FALSE)

        Data <- tibble(tests, pvalue, filter)
        OD <- Data[order(Data$filter, decreasing = TRUE), ]
        rankedtests <- OD$tests

        groupSize <- m/group
        testGroup <- rep(1:group, each = groupSize)

        testMeans <- as.vector(tapply(rankedtests, testGroup, mean))
        testSd <- as.vector(tapply(rankedtests, testGroup, sd))

        pi_hat <- testMeans*testMeans/(testMeans*testMeans + testSd*testSd - 1)
        effect_hat <- testMeans/pi_hat
        effect_hat[pi_hat <= 1/groupSize] <- 0

        if(sum(effect_hat, na.rm = TRUE) == 0){
            normWeight.w <- rep(1, m)
        } else {
            theta = effect_hat
            cr_val = qnorm(alpha/(tail*m), lower.tail = FALSE)
            c0 = cr_val*cr_val/2

            c <- sapply(theta, findc, alpha = alpha, c0 = c0, m = m)
            wgt_per_grp <- (m/alpha)*pnorm((effect_hat/2 + c/effect_hat),
                                           lower.tail = FALSE)

            wgt_smooth <- (1 - gamma)*wgt_per_grp + gamma*sum(wgt_per_grp)/group
            wgt_per_test <- rep(wgt_smooth, each = groupSize)

            norm_wgt <- wgt_per_test/sum(wgt_per_test)*m
        }

        return(norm_wgt)
    }











