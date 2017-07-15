#' @title Plots to compare ranks probabilities
#'
#' @description \code{OPWeight} package proposed methods to compute the probabilities
#' of the rank of test given the effect size. This funciton uses the simulation results
#' obtained from the three approahes: 1) simulation, 2) exact formula, and
#' 3) normal approximation to plot the rank probabilities.
#'
#' @param m0 number of true null tests
#' @param effectType type of effect sizes; c("continuous", "binary")
#'
#' @details
#' To use this function to plot, one needs to load "simu_prob_rank_givenEffect.RDATA"
#' and apply for different \code{m0}.
#'
#' @author Mohamad S. Hasan, \email{shakilmohamad7@gmail.com}
#' @export
#'
#'
#' @return \code{Data} A data frame containing seven columns; ranks and null and
#' alternative probabilities of the test from the three approaches
#'
#' @references Hasan and Schliekelman (2017)
#'
#' @examples
#' # just an example
#' # load("simu_prob_rank_givenEffect.RDATA")
#' # nullSize <- 90
#' # lapply(nullSize, ranksProb_compare_plots, effectType = "continuous")
#'
#===============================================================================
# Function to plots from the parallel computing outputs---------
ranksProb_compare_plots <- function(m0, effectType = c("binary", "continuous"))
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

            g[[i]] = ggplot(dat, aes_string(x = "ranks", y = "value", group = "variable",
                                     colour = "variable")) +
                geom_line(aes_string(linetype = "variable"), size = 1.5) +
                labs(x = "Ranks", y = "p(rank | effect)", size = 20,
                     subtitle = paste0(eySm, ey_val, ", e.one = ", e.one)) +
                theme(legend.title = element_blank(),
                      axis.title.x = element_text(size = rel(.7)),
                      axis.title.y = element_text(size = rel(.7)))
        }

        # extract the legend from one of the plots
        legend <- get_legend(g[[1]] + theme(legend.direction="horizontal",
                                            legend.position="bottom"))

        # arrange the three plots in a single row
        prow[[j]] <- plot_grid(g[[1]] + theme(legend.position="none"),
                               g[[2]] + theme(legend.position="none"),
                               g[[3]] + theme(legend.position="none"),
                               align = 'hv', nrow = 1, labels = letters[(3*j-3+1):(3*j)])
    }

    effType2 <- ifelse(effectType == "binary", "Binary: m0 = ", "Continuous: m0 = ")
    pp = plot_grid(prow[[1]], prow[[2]], prow[[3]], nrow = 3, align = 'hv')
    # now add the title
    title <- ggdraw() + draw_label(paste(effType2, m0, ", m1 = ", 100-m0))
    plots = plot_grid(title, pp, legend, ncol = 1, rel_heights=c(0.1, 1, .1))

    return(list(plots))
}
