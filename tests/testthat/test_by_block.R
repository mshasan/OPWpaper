test_by_block <- function(r, eVec, groupSize, Sigma)
    {
        eSub <- eVec[(groupSize*r + 1 - groupSize):(groupSize*r)]
        test <- as.vector(rmvn(1, eSub, Sigma, ncores = 15))
        return(test)
    }


