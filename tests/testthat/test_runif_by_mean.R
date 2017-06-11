#===============================================================================

runif_by_mean <- function(mean, n)
    {
        sd = mean/2
        uni_rv <- mean + sd*scale(runif(n, 0, 1))
        return(as.vector(uni_rv))
    }

