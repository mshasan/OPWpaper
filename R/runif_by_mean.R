#' @title Generate uniform random numbers for a fixed mean
#'
#' @description A function to generate uniform random numbers from a fixed mean
#'
#' @param n number of observations to be generated
#' @param mean mean of the random numbers
#'
#' @details Sometime it is necessary to generate random numbers with a prespecified
#' mean from the uniform distribution. This function will fullfill the goal
#'
#' @author Mohamad S. Hasan and Paul Schliekelman
#' @export
#'
#' @import stats
#' @seealso \code{\link{scale}} \code{\link{runif}}
#'
#' @return a vector of random numbers
#' @examples
#' x = runif_by_mean(n = 100, mean = 3)
#' summary(x)
#'
#===============================================================================
# function to generate uniform random numbers for a fixed mean
#
# Input:----------
# n = number of observations to be generated
# mean = mean of the random numbers
#
# internal parameters:----------
# sd = standard deviation
#
#output:-----------
# a vector of random numbers
#===============================================================================

runif_by_mean <- function(n, mean)
    {
        sd = mean/2
        uni_rv <- mean + sd*scale(runif(n, 0, 1))
        return(as.vector(uni_rv))
    }

