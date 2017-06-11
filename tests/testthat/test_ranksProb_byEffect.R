
#===============================================================================

ranksProb_byEffect <- function(i, null, m, nrep = 10000, filterEffectVec)
{
    ey <- filterEffectVec[i]
    et <- ey
    m0 <- ceiling(m*null)
    m1 <- m - m0
    rankProbsH1 <- sapply(1:m, prob_rank_givenEffect, et = ey, ey = ey,
                          nrep = nrep, m0 = m0, m1 = m1)
    return(rankProbsH1)
}


