library(snow)


cl <- makeCluster(10, type = "MPI")		# start zcluster

# define various functions-------------

#===============================================================================
# function to compute p(rank=k|filterEffect=ey) by simulation
# we used only uniform effects for continuous case.
prob_rank_givenEffect_simu <- function(s, ey, e.one, m0, m1,
                                       effectType = c("binary", "continuous"))
{
    m = m0 + m1
    ey0 <- rep(0, m0)

    if(effectType == "binary"){ey1 <- rep(ey, m1)
    } else { if(ey == 0){ey1 <- rep(0, m1)
    } else {ey1 <- runif(m1, ey-1, ey)
    }
    }
    Ey <- c(ey0, ey1)
    Ey[m0+1] <- e.one
    t01 <- rnorm(m, Ey, 1)
    r0 <- rank(-t01)[1]
    r1 <- rank(-t01)[m0+1]
    cbind(r0, r1)
}

# call the function=============================================================
clusterExport(cl, "prob_rank_givenEffect_simu")


#===============================================================================
# function to compute p(rank=k|filterEffect=ey) by exact method
# we used only uniform effects for continuous case.
prob_rank_givenEffect_exact <- function(k, et, ey, nrep = 100000, m0, m1,
                                        effectType = c("binary", "continuous"))
{
    k0 <- 1:k
    fun.k0 <- function(k0)
    {
        t <- rnorm(nrep, et, 1)
        p0 <- pnorm(-t)

        if(effectType == "binary"){p1 <- pnorm(ey - t)
        } else { if(ey == 0){p1 <- pnorm(ey - t)
        } else {
            m = m0 + m1
            a = ey - 1
            b = ey
            xb = b - t
            xa = a - t
            p1 = (xb*pnorm(xb) - xa*pnorm(xa) + dnorm(xb) - dnorm(xa))/(b-a)
        }
        }

        E.T <- ifelse(et == 0, mean(dbinom(k0-1, m0-1, p0)*dbinom(k-k0, m1, p1)),
                      mean(dbinom(k0-1, m0, p0)*dbinom(k-k0, m1-1, p1)))
        return(E.T)
    }
    prob <- sum(sapply(k0,fun.k0))
    return(prob)
}

# call the function=============================================================
clusterExport(cl, "prob_rank_givenEffect_exact")


#===============================================================================
# function to compute p(rank=k|filterEffect=ey) by normal approximation
# we used only uniform effects for continuous case.
prob_rank_givenEffect_approx <- function(k, et, ey, nrep = 100000, m0, m1,
                                         effectType = c("binary", "continuous"))
{
    t <- rnorm(nrep, et, 1)
    p0 <- pnorm(-t)

    if(effectType == "binary"){p1 <- pnorm(ey - t)
    } else { if(ey == 0){p1 <- pnorm(ey - t)
    } else {
        m = m0 + m1
        a = ey - 1
        b = ey
        xb = b - t
        xa = a - t
        p1 = (xb*pnorm(xb) - xa*pnorm(xa) + dnorm(xb) - dnorm(xa))/(b-a)
    }
    }

    mean0 <- (m0 - 1)*p0 + m1*p1 + 1
    mean1 <- m0*p0 + (m1 - 1)*p1 + 1
    var0 <- (m0 - 1)*p0*(1 - p0) + m1*p1*(1 - p1)
    var1 <- m0*p0*(1 - p0) + (m1 - 1)*p1*(1 - p1)
    prob <- ifelse(et == 0, mean(dnorm(k, mean0, sqrt(var0))),
                   mean(dnorm(k, mean1, sqrt(var1))))
    return(prob)
}

# call the function=============================================================
clusterExport(cl, "prob_rank_givenEffect_approx")


#===============================================================================
# function to compare rank probbailities from three approahes: 1) sikmulation,
# 2) exact fromual, and 3) normal approximation
ranksProb_compare <- function(ey, e.one, m0, m1, sampleSize,
                              effectType = c("binary", "continuous"))
{
    m = m0 +m1
    ranks <- 1:m
    et = e.one

    # simulation approach=======================================================
    # compute rank of the tests
    rank <- sapply(1:sampleSize, prob_rank_givenEffect_simu, ey = ey, e.one = e.one,
                   m0=m0, m1=m1, effectType = effectType)

    # rank may generate missing valuse because of the large effcet size,
    # therefore, to make a matplot equal vectors size are needed. This procedure
    # will replace the missing value to make equal sized vector
    # probability of rank of a null test
    prob0 <- rep(NA, m)
    prob0_x <- tapply(rank[1,], rank[1,], length)/sampleSize
    prob0[as.numeric(names(prob0_x))] <- as.vector(prob0_x)

    # probability of rank of an alternative test
    prob1 <- rep(NA, m)
    prob1_x <- tapply(rank[2,], rank[2,], length)/sampleSize
    prob1[as.numeric(names(prob1_x))] <- as.vector(prob1_x)

    # exact approach============================================================
    # do not compute probability for a large number of tests
    prob0_exact <- sapply(ranks, prob_rank_givenEffect_exact, et=0, ey=ey,
                          nrep=10000, m0=m0, m1=m1, effectType = effectType)
    prob1_exact <- sapply(ranks, prob_rank_givenEffect_exact, et=et, ey=ey,
                          nrep=10000,m0=m0, m1=m1, effectType = effectType)

    # normal approximation approcah=============================================
    prob0_approx <- sapply(ranks, prob_rank_givenEffect_approx, et=0, ey=ey,
                           nrep=10000, m0=m0, m1=m1, effectType = effectType)
    prob1_approx <- sapply(ranks, prob_rank_givenEffect_approx, et=et, ey=ey,
                           nrep=10000, m0=m0, m1=m1, effectType = effectType)

    # nice plots================================================================
    Data <- data.frame(ranks, prob0, prob1, prob0_exact, prob1_exact,
                       prob0_approx, prob1_approx)

    return(Data)
}


# For parallel computing========================================================
# null test size = c(20, 50, 75, 90, 99)
sampleSize = 10000000 #( use atleast 1,000,000)
clusterExport(cl, "sampleSize")

effect <- c(0, 1, 2)



# binary cases==================================================================
prob_20_0_bin = parLapply(cl, effect, ranksProb_compare, e.one = 0, m0=20, m1=80, sampleSize = sampleSize, effectType = "binary")
prob_20_1_bin = parLapply(cl, effect, ranksProb_compare, e.one = 1, m0=20, m1=80, sampleSize = sampleSize, effectType = "binary")
prob_20_2_bin = parLapply(cl, effect, ranksProb_compare, e.one = 2, m0=20, m1=80, sampleSize = sampleSize, effectType = "binary")

prob_50_0_bin = parLapply(cl, effect, ranksProb_compare, e.one = 0, m0=50, m1=50, sampleSize = sampleSize, effectType = "binary")
prob_50_1_bin = parLapply(cl, effect, ranksProb_compare, e.one = 1, m0=50, m1=50, sampleSize = sampleSize, effectType = "binary")
prob_50_2_bin = parLapply(cl, effect, ranksProb_compare, e.one = 2, m0=50, m1=50, sampleSize = sampleSize, effectType = "binary")

prob_75_0_bin = parLapply(cl, effect, ranksProb_compare, e.one = 0, m0=75, m1=25, sampleSize = sampleSize, effectType = "binary")
prob_75_1_bin = parLapply(cl, effect, ranksProb_compare, e.one = 1, m0=75, m1=25, sampleSize = sampleSize, effectType = "binary")
prob_75_2_bin = parLapply(cl, effect, ranksProb_compare, e.one = 2, m0=75, m1=25, sampleSize = sampleSize, effectType = "binary")

prob_90_0_bin = parLapply(cl, effect, ranksProb_compare, e.one = 0, m0=90, m1=10, sampleSize = sampleSize, effectType = "binary")
prob_90_1_bin = parLapply(cl, effect, ranksProb_compare, e.one = 1, m0=90, m1=10, sampleSize = sampleSize, effectType = "binary")
prob_90_2_bin = parLapply(cl, effect, ranksProb_compare, e.one = 2, m0=90, m1=10, sampleSize = sampleSize, effectType = "binary")

prob_99_0_bin = parLapply(cl, effect, ranksProb_compare, e.one = 0, m0=99, m1=1, sampleSize = sampleSize, effectType = "binary")
prob_99_1_bin = parLapply(cl, effect, ranksProb_compare, e.one = 1, m0=99, m1=1, sampleSize = sampleSize, effectType = "binary")
prob_99_2_bin = parLapply(cl, effect, ranksProb_compare, e.one = 2, m0=99, m1=1, sampleSize = sampleSize, effectType = "binary")


# continuous cases==================================================================
prob_20_0_cont = parLapply(cl, effect, ranksProb_compare, e.one = 0, m0=20, m1=80, sampleSize = sampleSize, effectType = "binary")
prob_20_1_cont = parLapply(cl, effect, ranksProb_compare, e.one = 1, m0=20, m1=80, sampleSize = sampleSize, effectType = "binary")
prob_20_2_cont = parLapply(cl, effect, ranksProb_compare, e.one = 2, m0=20, m1=80, sampleSize = sampleSize, effectType = "binary")

prob_50_0_cont = parLapply(cl, effect, ranksProb_compare, e.one = 0, m0=50, m1=50, sampleSize = sampleSize, effectType = "binary")
prob_50_1_cont = parLapply(cl, effect, ranksProb_compare, e.one = 1, m0=50, m1=50, sampleSize = sampleSize, effectType = "binary")
prob_50_2_cont = parLapply(cl, effect, ranksProb_compare, e.one = 2, m0=50, m1=50, sampleSize = sampleSize, effectType = "binary")

prob_75_0_cont = parLapply(cl, effect, ranksProb_compare, e.one = 0, m0=75, m1=25, sampleSize = sampleSize, effectType = "binary")
prob_75_1_cont = parLapply(cl, effect, ranksProb_compare, e.one = 1, m0=75, m1=25, sampleSize = sampleSize, effectType = "binary")
prob_75_2_cont = parLapply(cl, effect, ranksProb_compare, e.one = 2, m0=75, m1=25, sampleSize = sampleSize, effectType = "binary")

prob_90_0_cont = parLapply(cl, effect, ranksProb_compare, e.one = 0, m0=90, m1=10, sampleSize = sampleSize, effectType = "binary")
prob_90_1_cont = parLapply(cl, effect, ranksProb_compare, e.one = 1, m0=90, m1=10, sampleSize = sampleSize, effectType = "binary")
prob_90_2_cont = parLapply(cl, effect, ranksProb_compare, e.one = 2, m0=90, m1=10, sampleSize = sampleSize, effectType = "binary")

prob_99_0_cont = parLapply(cl, effect, ranksProb_compare, e.one = 0, m0=99, m1=1, sampleSize = sampleSize, effectType = "binary")
prob_99_1_cont = parLapply(cl, effect, ranksProb_compare, e.one = 1, m0=99, m1=1, sampleSize = sampleSize, effectType = "binary")
prob_99_2_cont = parLapply(cl, effect, ranksProb_compare, e.one = 2, m0=99, m1=1, sampleSize = sampleSize, effectType = "binary")



stopCluster(cl)


#------------------------------------------------------------
# save the workspace to the file .RData in the cwd
save.image("par_prob_rank_givenEffect.RData")


