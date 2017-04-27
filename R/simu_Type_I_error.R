library(snow)
library(qvalue)
library(tibble)
library(IHW)


cl <- makeCluster(75, type = "MPI")		# start zcluster


clusterExport(cl, "qvalue")
clusterExport(cl, "tibble")
clusterExport(cl, "ihw")
clusterExport(cl,"adj_pvalues")
clusterExport(cl,"rejections")



# prob = p(rank=k|effect=ey)
#===============================================================================
prob_rank_givenEffect <- function(k, et, ey, nrep = 10000, m0, m1)
{
    t <- rnorm(nrep, et, 1)
    p0 <- pnorm(-t)
    p1 <- pnorm(ey - t)

    mean0 <- (m0 - 1)*p0 + m1*p1 + 1
    mean1 <- m0*p0 + (m1 - 1)*p1 + 1

    var0 <- (m0 - 1)*p0*(1 - p0) + m1*p1*(1 - p1)
    var1 <- m0*p0*(1 - p0) + (m1 - 1)*p1*(1 - p1)

    prob <- ifelse(et == 0, mean(dnorm(k, mean0, sqrt(var0))),
                   mean(dnorm(k, mean1, sqrt(var1))))
    return(prob)
}


clusterExport(cl, "prob_rank_givenEffect")


# function to compute weight binary case
#--------------------------------------
weight_binary <- function(alpha, et, m, m1, tail = 2L, delInterval = .0001, prob)
{
    prob <- prob/sum(prob, na.rm = T)
    delta <- seq(0, 1, delInterval)
    findDelta <- function(delta)
    {
        weight <- tail*(m/alpha)*pnorm(et/2 + 1/et*log(delta*m/(alpha*m1*prob)),
                                       lower.tail = FALSE)
        return(sum(weight, na.rm = TRUE))
    }
    weightSumVec <- vapply(delta, findDelta, 1)
    deltaOut <- delta[min(abs(weightSumVec - m)) == abs(weightSumVec - m)]
    deltaOut <- ifelse(length(deltaOut) > 1, .0001, deltaOut)
    weight.out <- tail*(m/alpha)*pnorm(et/2 + 1/et*log(deltaOut*m/(alpha*m1*prob)),
                                       lower.tail = FALSE)
    sumWeight <- sum(weight.out, na.rm = TRUE)
    normWeight <- if(sumWeight == 0) {rep(1, m)} else {weight.out/sumWeight*m}
    return(normWeight)
}

clusterExport(cl, "weight_binary")



# function to compute weight continuous case
#--------------------------------------
weight_continuous <- function(alpha, et, m, tail=2L, delInterval=.0001, prob)
{
    prob <- prob/sum(prob, na.rm = T)
    delta <- seq(0, 1, delInterval)
    findDelta <- function(delta)
    {
        weight <- tail*(m/alpha)*pnorm(et/2 + 1/et*log(delta/(alpha*prob)),
                                       lower.tail=FALSE)
        return(sum(weight, na.rm = TRUE))
    }
    weightSumVec <- vapply(delta, findDelta, 1)
    deltaOut <- delta[min(abs(weightSumVec - m)) == abs(weightSumVec - m)]
    deltaOut <- ifelse(length(deltaOut) > 1, .0001, deltaOut)
    weight.out <- tail*(m/alpha)*pnorm(et/2 + 1/et*log(deltaOut/(alpha*prob)),
                                       lower.tail=FALSE)
    sumWeight <- sum(weight.out, na.rm = TRUE)
    normWeight <- if(sumWeight == 0) {rep(1, m)} else {weight.out/sumWeight*m}
    return(normWeight)
}

clusterExport(cl, "weight_continuous")


simu_Type_I_error <- function(alpha, m, simu, nrep = 1000, seed = NULL)
  {
    if (!is.null(seed)) set.seed(seed)
  
    pvals = matrix(runif(m*simu), nrow = m, ncol = simu, byrow = FALSE)
    pvals_filter = matrix(runif(m*simu), nrow = m, ncol = simu, byrow = FALSE)
    tests = qnorm(pvals, lower.tail = FALSE)
    filters = qnorm(pvals_filter, lower.tail = FALSE)
    
    Type_I_error <- function(simu)
      {
        test = tests[ , simu]
        pval = pvals[ , simu]
        filter = filters[ , simu]
      
        dat = tibble(test, pval, filter)
        OD = dat[order(dat$filter, decreasing=T), ]
        odered.pvalue = OD$pval
        
        model = lm(filter ~ test)
        
        null = qvalue(pval)$pi0
        m0 = ceiling(m*null)
        m1 = m-m0
        
        test_effect <- if(m1 == 0) {0
        } else {sort(test, decreasing = T)[1:m1]
          }
        
        et_bin = median(test_effect, na.rm = T)		        # median from the summary for the binary case test effect
        et_cont = mean(test_effect, na.rm = T)
        ey_bin = model$coef[[1]] + model$coef[[2]]*et_bin
        ey_cont = model$coef[[1]] + model$coef[[2]]*et_cont
        
        prob_bin <-vapply(1:m, prob_rank_givenEffect, 1, et = ey_bin,
                          ey = ey_bin, nrep = nrep, m0 = m0, m1 = m1)
        prob_cont <-vapply(1:m, prob_rank_givenEffect, 1, et = abs(ey_cont),
                           ey = abs(ey_cont), nrep = nrep, m0 = m0, m1 = m1)

        w_bin <- weight_binary(alpha = alpha, et = et_bin, m = m, m1 = m1, tail = 1,
                               delInterval = .0001, prob = prob_bin)
        w_cont = weight_continuous(alpha= alpha, et = et_cont, m = m, tail = 1,
                                   delInterval=.0001 , prob = prob_cont)
        
        ihw_fwer <- ihw(dat$pval, dat$filter, alpha = alpha, adjustment_type = "bonferroni")
        
        pro_bin = sum(odered.pvalue <= alpha*w_bin/m, na.rm = T)
        pro_cont = sum(odered.pvalue <= alpha*w_cont/m, na.rm = T)
        bon = sum(pval <= alpha/m, na.rm = T)
        IHW <- rejections(ihw_fwer)
        
        return(c(pro_bin, pro_cont, bon, IHW))
      }

    Type_I_error_per_simu <- sapply(1:simu, Type_I_error)
    return(apply(Type_I_error_per_simu, 1, mean))
  }



alphaVal = seq(.01, .1, .02)
typeIerror_mat = parSapply(cl, alphaVal, simu_Type_I_error, m = 10000, simu = 1000)



stopCluster(cl)


#------------------------------------------------------------
# save the workspace to the file .RData in the cwd
save.image("smu_Type_I_error.RData")

