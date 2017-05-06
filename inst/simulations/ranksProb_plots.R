#================================start of parrallelcomputing====================
# this data are generated using parallel cmputing system
# null test size = c(20, 50, 75, 90, 99)
sampleSize = 1000 #( use atleast 1,000,000)
effect <- c(0, 1, 2)

# continuous cases==============================================================
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

#===========================end of parallel computing===========================


# this code is to load saved workspace from parallel computing
load("U:/Documents/My Research (UGA)/Multiple Hypoetheses/Article-1/simu_prob_rank_givenEffect.RDATA")



#=======================probability plot for the supplementry materials=========

# Function to plots from the parallel computing outputs=========================
# probability plols for the supplementry materials

ranksProb_plots <- function(m0, effectType = c("binary", "continuous"))
{
    effType <- ifelse(effectType == "binary", "_bin", "_cont")
    prow <- list()
    for(j in 1:3)
    {
        e.one <- j-1
        # make text into varriable name
        probData_com <- eval(parse(text=paste0("prob_",m0,"_",e.one, effType)))
        g <- list()
        for(i in 1:3)
        {
            probData <- probData_com[[i]]
            colnames(probData) <- c("ranks", "SH0","SH1","EH0","EH1","AH0","AH1")
            dat <- melt(probData, id.var = "ranks")

            ey_val <- ifelse(effectType == "binary", i-1,
                             ifelse(i==1, 0, paste0("U(",i-2,", ", i-1,")")))
            eySm <- ifelse(effectType == "binary", "ey = ", "ey ~ ")

            g[[i]] = ggplot(dat, aes(x = ranks, y = value, group = variable,
                                     colour = variable)) +
                geom_line(aes(linetype = variable), size = 1.5) +
                labs(x = "ranks", y = "p(rank | effect)",
                     title = paste0(eySm, ey_val, ", e.one = ", e.one)) +
                theme(legend.title = element_blank())
        }

        # extract the legend from one of the plots
        legend <- get_legend(g[[1]] + theme(legend.direction="horizontal",
                                            legend.position="bottom"))

        # arrange the three plots in a single row
        prow[[j]] <- plot_grid(g[[1]] + theme(legend.position="none"),
                               g[[2]] + theme(legend.position="none"),
                               g[[3]] + theme(legend.position="none"),
                               align = 'hv', nrow = 1)
    }

    effType2 <- ifelse(effectType == "binary", "Binary: m0 = ", "Continuous: m0 = ")
    plots = grid.arrange(prow[[1]], prow[[2]], prow[[3]], legend, ncol=1,
                         heights=c(7,7,7,1), top= paste(effType2, m0, ", m1 = ", 100-m0))
    return(list(plots))
}


# applying plot function========================================================
nullSize <- c(20, 50, 75, 90, 99)
cont_probs_plots <- lapply(nullSize, ranksProb_plots, effectType = "continuous")
bin_probs_plots <- lapply(99, ranksProb_plots, effectType = "binary")

#================================end============================================
