#' @title Generate multivariate test statistics by block wise
#'
#' @description This function generate multivariate normal test statistics by block wise.
#'
#' @param r number of test groups
#' @param groupSize number of test statistics per group
#' @param eVec a vector of means
#' @param Sigma correlation matrix
#'
#' @details Generating large number of tests from multivariate normal distribution
#' is computionally very slow. Therefore, it is someitme convenient to genrate
#' parts of tests and then concatenate all tests to obtian the desired number of
#' tests statistics. This is called block wise generation of the test statistics.
#'
#' @author Mohamad S. Hasan and Paul Schliekelman
#' @export
#'
#' @import stats
#' @seealso \code{\link{rmvn}}
#'
#' @return \code{test} a vector of multivariate test statistics
#' @examples
#' # mean vector
#' eVec = runif(1000, 0, 2.5)
#'
#' # correaltion matrix
#' Sigma <- matrix(.9, 50, 50) + diag(50)*(1 - .1)
#'
#' # test statistics
#' testsStat = test_by_block(r=1, eVec, groupSize = 50,  Sigma)
#'
#===============================================================================
# function to generate uniform random numbers for a fixed mean
#
# Input:----------
# r = number of test groups
# groupSize = number of test statistics per group
# eVec = a vector of means
# Sigma = correlation matrix
#
# internal parameters:----------
# eSub = block of test statistics
#
#output:-----------
# a vector of multivariate test statistics
#===============================================================================

test_by_block <- function(r, eVec, groupSize, Sigma)
    {
        eSub <- eVec[(groupSize*r + 1 - groupSize):(groupSize*r)]
        test <- as.vector(rmvn(1, eSub, Sigma, ncores = 15))
        return(test)
    }


