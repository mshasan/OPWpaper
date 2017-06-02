#' @title Plot ranks' probabilities and the corresponding weights
#'
#' @description \code{OPWeight} package proposed methods to compute the probabilities
#' of the rank of test given the effect size and the corresponding weights.
#' This funciton uses the simulation results of the probabilities and the weights to plot.
#'
#' @param ey_index column numbers of the null proportions of a specific effect
#' @param null_index column numbers of the effects of specific null proportion
#' @param m total number of tests
#' @param ey mean filter effect size used in \code{ey_index}
#' @param null proportion of the null tests used in \code{null_index}
#' @param prob data of the ranks probability
#' @param weight data of the weight
#'
#' @details
#' To use this function to plot, one needs to load two data sets generated earlier
#' 1) ranksProb_byEffect_m10000.csv, and 2) Weight_byEffect_cont_m10000.csv then
#' the function.
#'
#' @author Mohamad S. Hasan, \email{mshasan@uga.edu}
#' @export
#'
#'
#' @return figure of four plots
#'
#' @references Hasan and Schliekelman (2017)
#'
#' @examples
#' # just an example
#' # m = 10000
#' # ranksProb <- read.csv("ranksProb_byEffect_m10000.csv", h=T)
#' # ranksWeight <- read.csv("Weight_byEffect_cont_m10000.csv",h=T)
#' # ranksWeight <- t(t(ranksWeight)/colSums(ranksWeight)*m)
#' # prob_weight_plots(ey_index = c(7, 17, 27, 37), null_index = c(26, 27, 28),
#' #         m = 10000, ey = 2, null = 50, prob = ranksProb, weight = ranksWeight)
#'
#==============================================================================

prob_weight_plots <- function(ey_index, null_index, m, ey, null, prob, weight)
    {
        ranks = 1:m

        prob_by_null = data.frame(ranks, prob[ , ey_index])
        names(prob_by_null) <- c("ranks", "20%", "50%", "75%", "90%")
        prob_by_null_melt <- melt(prob_by_null, id.var = "ranks",
                                  variable.name = "null prop.")

        prob_plot = ggplot(prob_by_null_melt, aes(x = prob_by_null_melt$ranks,
                    y = prob_by_null_melt$value, group = prob_by_null_melt$`null prop.`,
                    colour = prob_by_null_melt$`null prop.`)) +
            geom_line(aes(linetype = prob_by_null_melt$`null prop.`), size = 1.5) +
            labs(x = "Ranks", y = "p(rank | effect)", title = paste0("et = ey = ", ey))+
            theme(legend.position="none")


        prob_by_effect <- data.frame(ranks, prob[ , null_index])
        colnames(prob_by_effect) <- c("ranks", "0.6", "2.0", "5.0")
        prob_by_effect_melt <- melt(prob_by_effect, id.var = "ranks",
                                    variable.name = "effect size")
        prob_plot_by_effect = ggplot(prob_by_effect_melt, aes(x = prob_by_effect_melt$ranks,
                    y = prob_by_effect_melt$value, group = prob_by_effect_melt$`effect size`,
                        colour = prob_by_effect_melt$`effect size`)) +
            geom_line(aes(linetype = prob_by_effect_melt$`effect size`), size = 1.5) +
            labs(x = "Ranks", y = "p(rank | effect)", title = paste0("null = ", null, "%"))+
            theme(legend.position="none")


        weight_by_null = data.frame(ranks, weight[ , ey_index])
        names(weight_by_null) <- c("ranks", "20%", "50%", "75%", "90%")
        weight_by_null_melt <- melt(weight_by_null, id.var = "ranks",
                                    variable.name = "null prop.")
        weight_plot = ggplot(weight_by_null_melt, aes(x = weight_by_null_melt$ranks,
                      y = weight_by_null_melt$value, group = weight_by_null_melt$`null prop.`,
                      colour = weight_by_null_melt$`null prop.`)) +
            geom_line(size = 1.5) +
            labs(x = "Ranks", y = "log(weight)", color = "null prop.")+
            theme(legend.direction = "horizontal", legend.position = "bottom")



        weight_by_effect <- data.frame(ranks, weight[ , null_index])
        colnames(weight_by_effect) <- c("ranks", "0.6", "2.0", "5.0")
        weight_by_effect_melt <- melt(weight_by_effect, id.var = "ranks",
                                      variable.name = "effect size")
        weight_plot_by_effect = ggplot(weight_by_effect_melt, aes(x = weight_by_effect_melt$ranks,
                                y = weight_by_effect_melt$value,
                                group = weight_by_effect_melt$`effect size`,
                                colour = weight_by_effect_melt$`effect size`)) +
            geom_line(size = 1.5) +
            labs(x = "Ranks", y = "log(weight)", color = "effect size (ey)")+
            theme(legend.direction = "horizontal", legend.position = "bottom")


        pp <- plot_grid(prob_plot, prob_plot_by_effect, weight_plot, weight_plot_by_effect,
                        nrow = 2, labels = c("a", "b"), align = 'v')
        title <- ggdraw() + draw_label("Continuous: probability and weight vs. rank")
        plots = plot_grid(title, pp, ncol = 1, rel_heights=c(.1, 1))

        return(plots)

}


