
#
#===============================================================================

fwerPowerFdrPower <- function(i, simu, null, corr = 0, cv = 0, alpha = .05,
                                        groupSize = 100, effectType = c("continuous", "binary"),
                                        filterEffectVec, datWeightByNull)
    {
        W = datWeightByNull[ , i]
        m = length(W)
        weight_pro <- if(sum(W)==0){rep(1, m)} else {W/sum(W)*m}
        ey <- filterEffectVec[i]
        m0 <- ceiling(m*null)
        m1 <- m - m0

        if(effectType == "continuous"){
            xf <- as.vector(runif_by_mean(n = m, mean = ey))
        } else {
            xf <- rep(ey, m)
        }

        xt <- if(cv == 0){xf} else {rnorm(m, ey, cv*ey)}
        Sigma <- matrix(corr, groupSize, groupSize) + diag(groupSize)*(1 - corr)


        fwerPowerFdrPower_simu <- function(s)
            {
                H <- rbinom(m, 1, 1 - null)
                ef <- H*xf
                et <- H*xt
                mGrp = m/groupSize

                test <- if(corr == 0) {rnorm(m, et, 1)
                } else {as.vector(sapply(1:mGrp, test_by_block, eVec = et,
                                         groupSize = groupSize, Sigma = Sigma))}

                filter <- if(corr == 0) {rnorm(m, ef, 1)
                } else {as.vector(sapply(1:mGrp, test_by_block, eVec = ef,
                                         groupSize = groupSize, Sigma = Sigma))}

                pval <- pnorm(test, lower.tail = FALSE)

                dat = tibble(test, pval, et, filter)
                OD = dat[order(dat$filter, decreasing = TRUE), ]

                weight_rdw <- roeder_wasserman_weight(pvalue = OD$pval, alpha = alpha)
                ihw_fwer <- ihw(OD$pval, OD$filter, alpha = alpha, adjustment_type = "bonferroni")
                ihw_fdr <-  ihw(OD$pval, OD$filter, alpha = alpha, adjustment_type = "BH")

                rej_pro <- OD$pval <= alpha*weight_pro/m
                rej_bon <- OD$pval <= alpha/m
                rej_rdw <- OD$pval <= alpha*weight_rdw/m
                rej_ihwFwer <- adj_pvalues(ihw_fwer) <= alpha

                n_null <- max(1, sum(OD$et == 0, na.rm = TRUE))
                n_alt <-  max(1, sum(OD$et != 0, na.rm = TRUE))

                FWER_pro <- sum(rej_pro[OD$et == 0])
                FWER_bon <- sum(rej_bon[OD$et == 0])
                FWER_rdw <- sum(rej_rdw[OD$et == 0])
                FWER_ihw <- sum(rej_ihwFwer[OD$et == 0])

                POWER_pro <- sum(rej_pro[OD$et != 0])/n_alt
                POWER_bon <- sum(rej_bon[OD$et != 0])/n_alt
                POWER_rdw <- sum(rej_rdw[OD$et != 0])/n_alt
                POWER_ihw <- sum(rej_ihwFwer[OD$et != 0])/n_alt

                adjPval_pro <- p.adjust(OD$pval/weight_pro, method="BH")
                adjPval_bon <- p.adjust(OD$pval, method="BH")
                adjPval_rdw <- p.adjust(OD$pval/weight_rdw, method="BH")
                adjPval_ihw <- adj_pvalues(ihw_fdr)

                FDR_pro <- sum(adjPval_pro[OD$et == 0] <= alpha)/max(1, sum(adjPval_pro <= alpha))
                FDR_bh  <- sum(adjPval_bon[OD$et == 0] <= alpha)/max(1, sum(adjPval_bon <= alpha))
                FDR_rdw <- sum(adjPval_rdw[OD$et == 0] <= alpha)/max(1, sum(adjPval_rdw <= alpha))
                FDR_ihw <- sum(adjPval_ihw[OD$et == 0] <= alpha)/max(1, rejections(ihw_fdr))

                FDR_POWER_pro <- sum(adjPval_pro[OD$et != 0] <= alpha)/n_alt
                FDR_POWER_bh  <- sum(adjPval_bon[OD$et != 0] <= alpha)/n_alt
                FDR_POWER_rdw <- sum(adjPval_rdw[OD$et != 0] <= alpha)/n_alt
                FDR_POWER_ihw <- sum(adjPval_ihw[OD$et != 0] <= alpha)/n_alt

                return(c(FWER_pro, FWER_bon, FWER_rdw, FWER_ihw,
                         POWER_pro, POWER_bon, POWER_rdw, POWER_ihw,
                         FDR_pro, FDR_bh, FDR_rdw, FDR_ihw, FDR_POWER_pro,
                         FDR_POWER_bh, FDR_POWER_rdw, FDR_POWER_ihw))
            }

        fwerPowerFdrPower_bysimu <- sapply(1:simu, fwerPowerFdrPower_simu)
        fwerPowerFdrPower <- apply(fwerPowerFdrPower_bysimu, 1, mean, na.rm=TRUE)

        return(fwerPowerFdrPower)

    }





