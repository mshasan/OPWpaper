#setwd("C:/Users/mshasan/Google Drive/My R Packages/OPWpaper/inst/simulation_benchmarks")

library(reshape2)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(OPWpaper)


# nice plots---------
filter_test_relationship <- function(m0, m, alt_eff)
{
    m1 = m - m0
    ranks = 1:m

    prob_test0_cor.2 <- sapply(ranks, probRel_filterVstest_effect, rho=.2,
                               H0=0, ed=alt_eff, m0=m0, m1=m1)
    prob_test1_cor.2 <- sapply(ranks, probRel_filterVstest_effect, rho=.2,
                               H0=1, ed=alt_eff, m0=m0, m1=m1)

    prob_test0_cor.5 <- sapply(ranks, probRel_filterVstest_effect, rho=.5,
                               H0=0, ed=alt_eff, m0=m0,m1=m1)
    prob_test1_cor.5 <- sapply(ranks, probRel_filterVstest_effect, rho=.5,
                               H0=1, ed=alt_eff, m0=m0,m1=m1)

    prob_test0_cor.8 <- sapply(ranks, probRel_filterVstest_effect, rho=.8,
                               H0=0, ed=alt_eff,m0=m0,m1=m1)
    prob_test1_cor.8 <- sapply(ranks, probRel_filterVstest_effect, rho=.8,
                               H0=1, ed=alt_eff,m0=m0,m1=m1)

    prob0 <- sapply(ranks, prob_rank_givenEffect, et=0, ey=alt_eff,
                    m0=m0, m1=m1)
    prob1 <- sapply(ranks, prob_rank_givenEffect, et=alt_eff, ey=alt_eff,
                    m0=m0, m1=m1)


    datRelaion1 <- data.frame(ranks, prob0, prob1, prob_test0_cor.2, prob_test1_cor.2)
    colnames(datRelaion1) <- c("ranks", "FH0","FH1","TH0","TH1")
    datRelaion1_melt1 <- melt(datRelaion1, id.var="ranks")

    datRelaion2 <- data.frame(ranks, prob0, prob1, prob_test0_cor.5, prob_test1_cor.5)
    colnames(datRelaion2) <- c("ranks", "FH0","FH1","TH0","TH1")
    datRelaion1_melt2 <- melt(datRelaion2, id.var="ranks")

    datRelaion3 <- data.frame(ranks, prob0, prob1, prob_test0_cor.8, prob_test1_cor.8)
    colnames(datRelaion3) <- c("ranks", "FH0","FH1","TH0","TH1")
    datRelaion1_melt3 <- melt(datRelaion3, id.var="ranks")


    p_.2 <- ggplot(datRelaion1_melt1, aes(x = ranks, y = value, group = variable,
                                          colour = variable)) +
        geom_line(aes(linetype = variable), size = 1.5) +
        labs(x = "ranks", y = "p(rank | effect)", title = "cor = .2") +
        theme(legend.position="none")

    p_.5 <- ggplot(datRelaion1_melt2, aes(x = ranks, y = value, group = variable,
                                          colour = variable)) +
        geom_line(aes(linetype = variable), size = 1.5) +
        labs(x = "ranks", y = "p(rank | effect)", title = "cor = .5") +
        theme(legend.title = element_blank(), legend.position="bottom")

    p_.8 <- ggplot(datRelaion1_melt3, aes(x = ranks, y = value, group = variable,
                                          colour = variable)) +
        geom_line(aes(linetype = variable), size = 1.5) +
        labs(x = "ranks", y = "p(rank | effect)", title = "cor = .8") +
        theme(legend.position="none")

    # extract the legend from one of the plots
    legend_rel <- get_legend(p_.5 + theme(legend.direction="horizontal",
                                          legend.position="bottom"))
    p_.5 = p_.5 + theme(legend.position="none")

    # arrange the plots
    plots <- grid.arrange(p_.2, p_.5, p_.8,ggplot(NULL),legend_rel,ggplot(NULL), ncol=3,
                heights=c(7,1), top= paste("et = ", alt_eff, ", m0 = ", m0, ", m1 = ", m1))

    return(plots)
}


plots100<- filter_test_relationship(m0 = 100,m = 100, alt_eff = 0)
plots50 <- filter_test_relationship(m0 = 50, m = 100, alt_eff = 2)
plots90 <- filter_test_relationship(m0 = 90, m = 100, alt_eff = 2)






#------------------------------------------------------------
# save the workspace to the file .RData in the cwd
save.image("smu_relation_filterTest.RData")







