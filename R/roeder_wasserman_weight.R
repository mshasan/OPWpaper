#' @title Find sum of weights for c, function of Langragian multiplier
#'
#' @description Compute sum of weights for a given value of c, the function of
#' LaGrange multiplier
#'
#' @param c Numeric, a constant value (function of Langrangian multiplier)
#' @param m Integer, number of the hypothesis tests
#' @param gamma Numeric value of smooting parameter
#' @param alpha Numeric value of the Significance level
#' @param group Integer, number of groups
#' @param tail Integer (1 or 2), right-tailed or two-tailed hypothesis test.
#' @param effect_hat Numeric, estimated effect size
#'
#' @details
#' NOne
#' @author Mohamad S. Hasan, \email{shakilmohamad7@gmail.com}
#'
#' @export
#'
#' @return \code{sum weight} Numeric, sum of the weights for each C
#'
#' @references Roeder, Kathryn, and Larry Wasserman. "Genome-wide significance
#' levels and weighted hypothesis testing." Statistical science: a review journal
#' of the Institute of Mathematical Statistics 24.4 (2009): 398.
#' \url{www.stat.cmu.edu/~roeder/publications/statsci.pdf}
#'
#' @examples
#' cc = seq(-10, 10, .05)
#' et <- seq(0, 3, .5)
#' wgtSum_by_c <- sapply(cc, weightSum_by_c, m = 10000, effect_hat = et)
#'
#===============================================================================

weightSum_by_c <- function(c, m, gamma = .05,  alpha = .05, group = 5L,
                                  tail = 1L, effect_hat)
    {
        groupSize <- m/group

        weight_per_c <- tail*(m/alpha)*pnorm((effect_hat/2 + c/effect_hat),
                                                 lower.tail = FALSE)
        wgt_smooth_c <- (1 - gamma)*weight_per_c + gamma*sum(weight_per_c)/group
        wgt_per_test_c <- rep(wgt_smooth_c, each = groupSize)
        sumWeight_per_c <- sum(wgt_per_test_c, na.rm = TRUE)

        return(sumWeight_per_c)
    }

#===============================================================================

#' @title Weight from the Roeder and Wasserman (2009) paper
#'
#' @description Compute weights by splitting test statistics raked by the
#' filter statistics. This method is taken from Roeder and Waserman (2009).
#'
#' @param pvalue A numeric vector of ordered p-values sorted  by the covariate
#' @param gamma Numeric value of smooting parameter
#' @param alpha Numeric value of the Significance level
#' @param group Integer, number of groups
#' @param tail Integer (1 or 2), right-tailed or two-tailed hypothesis test.
#' @param c_interval A nuumeric vector of a sequence of intervals between the
#' \code{c}. Note that, \code{c} is a function of Langrangian multiplier,
#' necessary to normalize the weight
#'
#' @details
#' NOne
#' @author Mohamad S. Hasan, \email{shakilmohamad7@gmail.com}
#' @export
#'
#' @import tibble tibble
#'
#' @seealso \code{\link{weight_binary}} \code{\link{weight_continuous}}
#'
#' @return \code{weight} A numeric vector of the normalized weights of the tests
#'
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
#' library(tibble)
#' Data <- tibble(tests, pvals, filters)
#' OD <- Data[order(Data$filters, decreasing = TRUE), ]
#'
#' weight = roeder_wasserman_weight(pvalue = OD$pvals)
#'
#' # plot the weight
#' plot(weight)
#'
#' # compute number of rejections
#' alpha = .05
#' rwd <- sum(OD$pvals <= alpha*weight/m)
#' bon <- sum(pvals <= alpha/m)
#'
#===============================================================================
# Function to estiamte weight from data by using Roeder-Waserman algorithm
# internal parameters:-----
# testGroup = patitioning 'm' tests into groups
# testMeans = calculate mean of each group
# testSd = calculate standard deviation of each group
# pi_hat = estiamte pi parameter
# effect_hat = estiamte effect size
# wgtSum_by_c = sum of the weights per c
# c_Out = optimal c value
# wgt_per_grp = calculate weight per group
# weight.out = final weight
# sumWeight = sum of the weights
#===============================================================================

roeder_wasserman_weight <- function(pvalue, gamma = .05, alpha = .05,
                                    group = 5L, tail = 1L, c_interval = .01)
    {
        m = length(pvalue)
        # ordered pvalues, thus ordered tests------------
        rankedtests <- qnorm(pvalue/tail, lower.tail = FALSE)

        groupSize <- m/group
        testGroup <- rep(1:group, each = groupSize)

        testMeans <- as.vector(tapply(rankedtests, testGroup, mean))
        testSd <- as.vector(tapply(rankedtests, testGroup, sd))

        pi_hat <- testMeans*testMeans/(testMeans*testMeans + testSd*testSd - 1)
        effect_hat <- testMeans/pi_hat
        effect_hat[pi_hat <= 1/groupSize] <- 0

        if(sum(effect_hat, na.rm = TRUE) == 0){
            norm_wgt <- rep(1, m)
        } else {
            c <- seq(-10, 10, c_interval)
            wgtSum_by_c <- sapply(c, weightSum_by_c, m, gamma = .05,  alpha = .05,
                                     group = 5L, tail = 1L, effect_hat)

            c_out <- c[min(abs(wgtSum_by_c - m)) == abs(wgtSum_by_c - m)]
            c_out <- ifelse(length(c_out) > 1, -1, c_out)
            weight.out <- tail*(m/alpha)*pnorm((effect_hat/2 + c_out/effect_hat),
                                                               lower.tail = FALSE)

            wgt_smooth_cOut <- (1 - gamma)*weight.out + gamma*sum(weight.out)/group
            wgt_per_test_cOut <- rep(wgt_smooth_cOut, each = groupSize)
            sumWeight <- sum(wgt_per_test_cOut, na.rm = TRUE)

            norm_wgt <- if(sumWeight == 0) {
                            rep(1, m)
                        } else {
                            wgt_per_test_cOut/sumWeight*m
                        }
        }

        return(norm_wgt)
    }

