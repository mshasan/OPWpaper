#' @title Relationship between filter and test effect sizes
#'
#' @description Compute the relationship between filter and test effect sizes interms of
#' the probability of the rank of the test given the test effect sizes
#'
#' @param r rank of the test statistics
#' @param rho correlation between the filter and test efect sizes
#' @param H0 determine null or alternative hypotheisis; H0 = 0 for null and H0 = 1 for alternative
#' @param ed mean effect size of the test statistics
#' @param m0 number of the true null hypothesis
#' @param m1 number of the true alternative hypothesis
#'
#' @details Compute the relationship between filter and test effect sizes interms of
#' the probability of the rank of the tes given the test effect sizes. The weight
#' identity is based on the test effect size; however, the ranks probability needs
#' to compute from the filter effect. Therefore, it is expected that there is a
#' relationship between the filter and test effect sizes.
#'
#' @author Mohamad S. Hasan and Paul Schliekelman
#' @export
#' @import stats
#' @seealso \code{\link{prob_rank_givenEffect}}
#'
#' @return \code{prob} probability of the rank given the mean test effect
#' @examples
#' ranks = 1:100
#' prob_test0 <- sapply(ranks, probRel_filterVstest_effect, rho = .8,
#'                        H0 = 0, ed = 2, m0 = 90, m1 = 10)
#'
#' prob_test1 <- sapply(ranks, probRel_filterVstest_effect, rho=.8,
#'                        H0 = 1, ed = 2, m0 = 90, m1 = 10)
#'
#' prob0 <- sapply(ranks, prob_rank_givenEffect, et = 0, ey = 2,
#'                                            m0 = 90, m1 = 10)
#'
#' prob1 <- sapply(ranks, prob_rank_givenEffect, et = 2, ey = 2,
#'                                            m0 = 90, m1 = 10)
#'
#' matplot(1:100, cbind(prob_test0, prob_test1, prob0, prob1),
#'                      type = "l", xlab = "ranks", ylab = "p(rank | effect)")
#' legend("topright", legend = c("prob_test0", "prob_test1", "prob0", "prob1"),
#'                        col = 1:4, lty = 1:4, lwd = 2)
#'
#===============================================================================
# function to compute "P(rank|effect) by effect size for m=10,000"
#
# Input:----------
# r = rank of the test statistics
# rho = correlation between the filter and test efect sizes
# H0 = determine null or alternative hypotheisis; H0 = 0 for null and H0 = 1 for
# alternative
# ed = mean effect size of the test statistics
# m0 = number of the true null hypothesis
# m1 = number of the true alternative hypothesis
#
# internal parameters:----------
# mean_ey = conditional mean of the filter effect
# sd_ey = conditional standard deviation of the filter effect
# ey_val = a vector of the filter effects
# et = mean test effect
# probs_per_ey = probability of the rank of the test for each mean filter effect
#
# output:-----------
# prob = probability of the rank given the mean test effect
#===============================================================================

probRel_filterVstest_effect <- function(r, rho, H0, ed, m0, m1)
    {
        mean_ey = rho*ed
        sd_ey = sqrt(1 - rho*rho)
        ey_val = rnorm(100, mean_ey, sd_ey)

        prob_condition_ey <- function(ey)
            {
                et <- ifelse(H0 == 0, 0, ey)
                probs_per_ey = prob_rank_givenEffect(k = r, et = et, ey = ey,
                                                     m0 = m0, m1 = m1)
                return(probs_per_ey)
            }

        prob_per_r = mean(sapply(ey_val,  prob_condition_ey))
        return(prob_per_r)
    }


