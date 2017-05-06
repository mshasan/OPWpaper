#' @title Compare rank probabilities
#'
#' @description \code{OPWeight} package proposed methods to compute the probabilities
#' of the rank of test given the effect size. This funciton uses the methods to
#' compare the rank probabilities from three approahes: 1) simulation,
#' 2) exact formula, and 3) normal approximation
#'
#' @param ey mean filter effect sizevary
#' @param e.one one test effect that will vary across all tests
#' @param m0 number of true null tests
#' @param m1 number of true alternative tests
#' @param sampleSize total number of sample generated (use sample size at least 100,000)
#' @param effectType type of effect sizes; c("continuous", "binary")
#'
#' @details
#' \code{OPWeight} package proposed methods to compute the probabilities
#' of the rank of test given the effect size. This funciton uses the methods to
#' compare the rank probabilities from three approahes: 1) simulation,
#' 2) exact formula, and 3) normal approximation\cr
#'
#' rank may generate missing valuse because of the large effcet size,
#' therefore, to make a matplot equal vectors size are needed. This procedure
#' will replace the missing value to make equal sized vector
#' probability of rank of a null test
#'
#' @author Mohamad S. Hasan, \email{mshasan@uga.edu}
#' @export
#'
#' @seealso \code{\link{prob_rank_givenEffect_simu}}
#' \code{\link{prob_rank_givenEffect_exact}}
#' \code{\link{prob_rank_givenEffect_approx}}
#'
#' @return \code{Data} A data frame containing seven columns; ranks and null and
#' alternative probabilities of the test from the three approaches
#'
#' @references Hasan and Schliekelman (2017)
#'
#' @examples
#' # use sample size at least 100,000 for better result
#' sampleSize = 1000
#' probData <- ranksProb_compare(ey = 1, e.one = 2, m0 = 90, m1 = 10,
#'                          sampleSize = sampleSize, effectType = "binary")
#'
#' # plots------------
#' colnames(probData) <- c("ranks", "SH0","SH1","EH0","EH1","AH0","AH1")
#' matplot(probData[, 1], probData[, 2:5], type = "l", lty = 1:6, col =1:6,
#' lwd = 2, xlab = "ranks", ylab = "P(rank | effect)")
#' legend("topright", legend = c("SH0","SH1","EH0","EH1","AH0","AH1"),
#'                 lty = 1:6, col =1:6, lwd = 2)
#'
#===============================================================================
# inpout:----------------
# ey = effect size
# e.one = vary one test effect across all tests
# m0 = number of null tests
# m1 = number of alternative tests
# effectType = type of effect size c("binary","continuous")
# sampleSize = total number of sample generated (use sample size at least 100,000)
# effectType = type of effect sizes; c("continuous", "binary")
#
# internal parameters:-----
# m = total number of tests
# ranks = sequesnce of ranks or index numbers
# rank = rank of the tests per sample
#
# prob0 = probability of rank of a null test by simulation
# prob1 = probability of rank of an alternative test by simulation
# prob0, prob1 =  rank may generate missing valuse because of the large effcet size,
# therefore, to make a matplot equal vectors size are needed. This procedure
# will replace the missing value to make equal sized vector
# probability of rank of a null test
#
# prob0_exact = probability of rank of a null test by exact mehtod
# prob1_exact = probability of rank of an alternative test by exact mehtod
# do not compute probability for a large number of tests
#
# prob0_approx = probability of rank of a null test by normal approximaiton
# prob1_exact = probability of rank of an alternative test by normal approximaiton
#
# output:---------------
# Data = data frame of the probabilities of the tests
#
#===============================================================================
ranksProb_compare <- function(ey, e.one, m0, m1, sampleSize,
                              effectType = c("continuous", "binary"))
    {
        m = m0 +m1
        ranks <- 1:m
        et = e.one

        rank <- sapply(1:sampleSize, prob_rank_givenEffect_simu, ey = ey,
                       e.one = e.one, m0 = m0, m1 = m1, effectType = effectType)

        prob0 <- rep(NA, m)
        prob0_x <- tapply(rank[1,], rank[1,], length)/sampleSize
        prob0[as.numeric(names(prob0_x))] <- as.vector(prob0_x)

        prob1 <- rep(NA, m)
        prob1_x <- tapply(rank[2,], rank[2,], length)/sampleSize
        prob1[as.numeric(names(prob1_x))] <- as.vector(prob1_x)

        prob0_exact <- sapply(ranks, prob_rank_givenEffect_exact, et = 0, ey = ey,
                        nrep = 10000, m0 = m0, m1 = m1, effectType = effectType)
        prob1_exact <- sapply(ranks, prob_rank_givenEffect_exact, et = et, ey = ey,
                        nrep = 10000, m0 = m0, m1 = m1, effectType = effectType)

        prob0_approx <- sapply(ranks, prob_rank_givenEffect_approx, et = 0, ey = ey,
                         nrep = 10000, m0 = m0, m1 = m1, effectType = effectType)
        prob1_approx <- sapply(ranks, prob_rank_givenEffect_approx, et = et, ey = ey,
                         nrep = 10000, m0 = m0, m1 = m1, effectType = effectType)

        Data <- tibble(ranks, prob0, prob1, prob0_exact, prob1_exact,
                           prob0_approx, prob1_approx)

        return(Data)
    }
