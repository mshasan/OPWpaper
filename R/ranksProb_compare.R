#' @title Compare rank probabilities
#'
#' @description \code{OPWeight} package proposed a method to compute the
#' ranks probabilities of the covariate given the test-effect sizes from three
#' approaches: simualation, exact formula, and normal approximation. This
#' funciton uses the methods to compare the ranks probabilities from the three
#' approahes
#'
#' @param ey Numerics, mean covariate-effect size
#' @param e.one Numeric, one test effect which will vary across all tests
#' @param m0 Integer, number of true null tests
#' @param m1 Integer, number of true alternative tests
#' @param sampleSize Integer, total number of sample generated (use sample size
#' at least 100,000)
#' @param effectType Character ("continuous" or "binary"), type of effect sizes
#'
#' @details
#' The \code{OPWeight} package proposed methods to compute the ranks
#' probabilitiesof the covariate given the test effect size. This funciton uses
#' the methods to compare the rank probabilities from three approahes:
#' 1) simulation, 2) exact formula, and 3) normal approximation\cr
#'
#' The lower rank may generate missing values because of the large effcet sizes.
#' This is particularly true for the simulaiton approach. however,
#' \code{matplot} function requires equal sized vectors. This procedure will
#' replace the missing values by NA so that the vectors size become equal.
#'
#' @author Mohamad S. Hasan, \email{shakilmohamad7@gmail.com}
#' @export
#'
#' @seealso \code{\link{prob_rank_givenEffect_simu}}
#' \code{\link{prob_rank_givenEffect_exact}}
#' \code{\link{prob_rank_givenEffect_approx}}
#'
#' @return \code{Data} A data frame containing the seven columns; the ranks and
#' the corresponding ranks probability of the true null and the true alternative
#' hypothesis of the three approaches.
#'
#' @references Hasan and Schliekelman (2017)
#'
#' @examples
#' # use sample size at least 100,000 for better result
#' # This is just an example
#' sampleSize = 1000
#' probData <- ranksProb_compare(ey = 1, e.one = 2, m0 = 5, m1 = 5,
#'                          sampleSize = sampleSize, effectType = "binary")
#'
#' # plots------------
#' # colnames(probData) <- c("ranks", "SH0","SH1","EH0","EH1","AH0","AH1")
#' # matplot(probData[, 1], probData[, 2:5], type = "l", lty = 1:6, col =1:6,
#' # lwd = 2, xlab = "ranks", ylab = "P(rank | effect)")
#' # legend("topright", legend = c("SH0","SH1","EH0","EH1","AH0","AH1"),
#' #               lty = 1:6, col =1:6, lwd = 2)
#'
#===============================================================================
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
