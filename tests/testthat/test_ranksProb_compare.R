
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
