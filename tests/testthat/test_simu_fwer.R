
#===============================================================================

simu_fwer <- function(s, m, alphaVec)
    {
    fwer_per_rep <- function(alpha)
        {
            pval <- runif(m)
            pval_filter <- runif(m)
            test = qnorm(pval, lower.tail = FALSE)
            filter = qnorm(pval_filter, lower.tail = FALSE)

            dat = tibble(test, pval, filter)

            OD = dat[order(dat$filter, decreasing=TRUE), ]
            odered.pvalue = OD$pval

            nullprop = qvalue(pval)$pi0
            m0 = ceiling(m*nullprop)
            m1 = m - m0

            model = lm(filter ~ test)

            test_effect <- if(m1 == 0) {0
                           } else {sort(test, decreasing = TRUE)[1:m1]}

            et_bin = median(test_effect, na.rm = TRUE)
            et_cont = mean(test_effect, na.rm = TRUE)

            ey_bin = model$coef[[1]] + model$coef[[2]]*et_bin
            ey_cont = model$coef[[1]] + model$coef[[2]]*et_cont

            prob_bin <-sapply(1:m, prob_rank_givenEffect, et = ey_bin,
                              ey = ey_bin,m0 = m0, m1 = m1)
            prob_cont <-sapply(1:m, prob_rank_givenEffect, et = ey_cont,
                               ey = ey_cont, m0 = m0, m1 = m1)

            w_bin <- weight_binary(alpha = alpha, et = et_bin, m = m, m1 = m1,
                            tail = 1, delInterval = .0001, ranksProb = prob_bin)
            w_cont = weight_continuous(alpha = alpha, et = et_cont, m = m,
                            tail = 1, delInterval = .0001 , ranksProb = prob_cont)

            ihw_fwer <- ihw(dat$pval, dat$filter, alpha = alpha,
                                            adjustment_type = "bonferroni")

            bon = sum(pval <= alpha/m, na.rm = TRUE)
            pro_bin = sum(odered.pvalue <= alpha*w_bin/m, na.rm = TRUE)
            pro_cont = sum(odered.pvalue <= alpha*w_cont/m, na.rm = TRUE)
            IHW <- rejections(ihw_fwer)

            return(c(bon, pro_bin, pro_cont, IHW))
        }

        fwer_per_rep_mat = sapply(alphaVec, fwer_per_rep)
        return(fwer_per_rep_mat)
    }







