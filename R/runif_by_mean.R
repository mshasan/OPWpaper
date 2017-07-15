#' @title Generate uniform random numbers for a fixed mean
#'
#' @description A function to generate uniform random numbers from a fixed mean
#'
#' @param n Integer, number of observations to be generated
#' @param mean Numeric, mean of the random numbers
#'
#' @details Sometime it is necessary to generate random numbers with a prespecified
#' mean from the uniform distribution. This function will fullfill the goal
#'
#' @author Mohamad S. Hasan, \email{shakilmohamad7@gmail.com}
#' @export
#'
#' @import stats
#' @seealso \code{\link{scale}} \code{\link{runif}}
#'
#' @return A numeric vector of uniform random numbers
#'
#' @examples
#' x = runif_by_mean(mean = 3, n = 100)
#' summary(x)
#'
#===============================================================================
# internal parameters:----------
# sd = standard deviation
#===============================================================================

runif_by_mean <- function(mean, n)
    {
        sd = mean/2
        uni_rv <- mean + sd*scale(runif(n, 0, 1))
        return(as.vector(uni_rv))
    }

