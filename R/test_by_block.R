#' @title Generate multivariate test statistics
#'
#' @description This function generate multivariate normal test statistics by
#' block wise.
#'
#' @param r Integer, number of the test groups
#' @param groupSize Integer, number of test statistics per group
#' @param eVec A numeric vector of the means
#' @param Sigma A numeric matrix correlations
#'
#' @details Generating a large number of tests statistics from the multivariate normal
#' distribution is computionally very expensive. Therefore, it is someitme convenient
#' to genrate parts of tests and then concatenate all the tests to obtian the
#' desired number of the tests statistics. This is called block wise generation
#' of the test statistics.
#'
#' @author Mohamad S. Hasan, \email{shakilmohamad7@gmail.com}
#' @export
#'
#' @import stats
#'
#' @seealso \code{\link{rmvn}}
#'
#' @return \code{test} A numeric vector of multivariate test statistics
#'
#' @examples
#' # mean vector
#' eVec = runif(1000, 0, 2.5)
#'
#' # correaltion matrix
#' Sigma <- matrix(.9, 50, 50) + diag(50)*(1 - .1)
#'
#' # test statistics
#' testsStat = test_by_block(r = 1, eVec, groupSize = 50,  Sigma)
#'
#===============================================================================
# internal parameters:----------
# eSub = block of test statistics
#
#===============================================================================

test_by_block <- function(r, eVec, groupSize, Sigma)
    {
        eSub <- eVec[(groupSize*r + 1 - groupSize):(groupSize*r)]
        test <- as.vector(rmvn(1, eSub, Sigma, ncores = 15))
        return(test)
    }


