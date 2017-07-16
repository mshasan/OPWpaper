#' @title Relationship between covariate and test effect sizes
#'
#' @description Compute the relationship between the covariate and test effect
#' sizes in terms of the ranks probability of the covariate given the test effect
#' sizes
#'
#' @param r Integer, rank of the covariate
#' @param rho Numeric, correlation between the covariate and the test efect sizes
#' @param H0 Binary 0 or 1, determine the null or the alternative hypotheisis;
#' H0 = 0 if null and H0 = 1 if alternative
#' @param ed Numeric, mean effect size of the test statistics
#' @param m0 Integer, number of true null hypothesis
#' @param m1 Integer, number of the true alternative hypothesis
#' @param n_ey Integer, number of covariate-effects to be generated
#'
#' @details Compute the relationship between the covariate and the test effect
#' sizes in terms of the ranks probability of the covariate given the test effect sizes.
#' The weight identity is based on the test effect size; however, the ranks
#' probability needs to compute from the covariate effects. Therefore, it is
#' expected that there is a relationship between the covariate and test effect sizes.
#'
#' @author Mohamad S. Hasan, \email{shakilmohamad7@gmail.com}
#' @export
#' @import stats
#' @seealso \code{\link{prob_rank_givenEffect}}
#'
#' @return \code{prob} A numeric value of the ranks probability of the test
#' given the mean test effect
#' @examples
#' ranks = 1:10
#' prob_test0 <- sapply(ranks, probRel_CovVsTest_effect, rho = .8,
#'                        H0 = 0, ed = 2, m0 = 9, m1 = 1)
#'
#' # prob_test1 <- sapply(ranks, probRel_CovVsTest_effect, rho=.8,
#' #                       H0 = 1, ed = 2, m0 = 9, m1 = 1)
#'
#' # prob0 <- sapply(ranks, prob_rank_givenEffect, et = 0, ey = 2,
#' #                                           m0 = 9, m1 = 1)
#'
#' # prob1 <- sapply(ranks, prob_rank_givenEffect, et = 2, ey = 2,
#' #                                           m0 = 9, m1 = 1)
#'
#' # matplot(1:10, cbind(prob_test0, prob_test1, prob0, prob1),
#' #                     type = "l", xlab = "ranks", ylab = "p(rank | effect)")
#' # legend("topright", legend = c("prob_test0", "prob_test1", "prob0", "prob1"),
#' #                       col = 1:4, lty = 1:4, lwd = 2)
#'
#===============================================================================
# internal parameters:----------
# mean_ey = conditional mean of the covariate effect
# sd_ey = conditional standard deviation of the covariate effect
# ey_val = a vector of the covariate effects
# et = mean test effect
# probs_per_ey = probability of the rank of the test for each mean covariate effect
#
#===============================================================================

probRel_CovVsTest_effect <- function(r, rho, H0, ed, m0, m1, n_ey = 100)
    {
        mean_ey = rho*ed
        sd_ey = sqrt(1 - rho*rho)
        ey_val = rnorm(n_ey, mean_ey, sd_ey)

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


