
#===============================================================================

probRel_filterVstest_effect <- function(r, rho, H0, ed, m0, m1, n_ey = 100)
    {
        mean_ey = rho*ed
        sd_ey = sqrt(1 - rho*rho)
        ey_val = rnorm(n_ey, mean_ey, sd_ey)

        prob_condition_ey <- function(ey)
            {
                et <- ifelse(H0 == 0, 0, ey)
                probs_per_ey = prob_rank_givenEffect(k = r, et = et, ey = ey,
                                                     m0 = m0, m1 = m1)
                return(probs_per_ey)
            }

        prob_per_r = mean(sapply(ey_val,  prob_condition_ey))
        return(prob_per_r)
    }


