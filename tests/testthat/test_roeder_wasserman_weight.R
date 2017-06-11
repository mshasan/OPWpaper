
#-------------------------------------------------------------------------------

roeder_wasserman_weight <- function(pvalue, gamma = .05, alpha = .05,
                                    group = 5L, tail = 1L, c_interval = .01)
    {
        m = length(pvalue)
        # ordered pvalues, thus ordered tests------------
        rankedtests <- qnorm(pvalue/tail, lower.tail = FALSE)

        groupSize <- m/group
        testGroup <- rep(1:group, each = groupSize)

        testMeans <- as.vector(tapply(rankedtests, testGroup, mean))
        testSd <- as.vector(tapply(rankedtests, testGroup, sd))

        pi_hat <- testMeans*testMeans/(testMeans*testMeans + testSd*testSd - 1)
        effect_hat <- testMeans/pi_hat
        effect_hat[pi_hat <= 1/groupSize] <- 0

        if(sum(effect_hat, na.rm = TRUE) == 0){
            norm_wgt <- rep(1, m)
        } else {
            c <- seq(-10, 10, c_interval)
            wgtSum_by_c <- sapply(c, weightSum_by_c, m, gamma = .05,  alpha = .05,
                                     group = 5L, tail = 1L, effect_hat)

            c_out <- c[min(abs(wgtSum_by_c - m)) == abs(wgtSum_by_c - m)]
            c_out <- ifelse(length(c_out) > 1, -1, c_out)
            weight.out <- tail*(m/alpha)*pnorm((effect_hat/2 + c_out/effect_hat),
                                                               lower.tail = FALSE)

            wgt_smooth_cOut <- (1 - gamma)*weight.out + gamma*sum(weight.out)/group
            wgt_per_test_cOut <- rep(wgt_smooth_cOut, each = groupSize)
            sumWeight <- sum(wgt_per_test_cOut, na.rm = TRUE)

            norm_wgt <- if(sumWeight == 0) {
                            rep(1, m)
                        } else {
                            wgt_per_test_cOut/sumWeight*m
                        }
        }

        return(norm_wgt)
    }

