#' @title Plot ranks' probabilities and the corresponding weights
#'
#' @description \code{OPWeight} package proposed a method to compute the ranks
#' probabilities of the covariates and the corresponding weights given the
#' test-effect size. This funciton uses the simulation results of the
#' probabilities and the weights to plot.
#'
#' @param ey_index Integer vector, column numbers of the null proportions of a
#' specific effect
#' @param null_index Integer vector, column numbers of the covariate-effects of
#' a specific null proportion
#' @param m Integer, total number of tests
#' @param ey Numeric, mean covariate-effect size used in \code{ey_index}
#' @param null Numeric, proportion of the null tests used in \code{null_index}
#' @param prob Numerics data of the ranks probabilities
#' @param weight Numeric data of the weights
#'
#' @details
#' To use this function, one needs to load two data sets generated earlier
#' 1) ranksProb_byEffect_m10000.csv, and 2) Weight_byEffect_cont_m10000.csv then
#' use the function.
#'
#' @author Mohamad S. Hasan, \email{shakilmohamad7@gmail.com}
#' @export
#'
#'
#' @return Figure of four plots
#'
#' @references Hasan and Schliekelman (2017)
#'
#' @examples
#' # ranksProb <- read.csv(system.file("simulations/results/ranksProb_byEffect_m10000.csv",
#' #                                   package = "OPWpaper"), h = TRUE)
#' # ranksWeight_cont <- read.csv(system.file("simulations/results/Weight_byEffect_cont_m10000.csv",
#' #                                          package = "OPWpaper"), h = TRUE)
#' # prob_weight_plots(ey_index = c(7, 17, 27, 37), null_index = c(14, 17, 19),
#' # m = 10000, ey = 2, null = 50, prob = ranksProb, weight = ranksWeight_cont)
#'
#==============================================================================

prob_weight_plots <- function(ey_index, null_index, m, ey, null, prob, weight)
{
    ranks = 1:m

    # ranksProb plots for the different null proportions------
    p_null = prob[ , ey_index]
    prob_by_null = data.frame(ranks, cbind(SMA(p_null[,1], n = 100),
                                           SMA(p_null[,2], n = 100),
                                           SMA(p_null[,3], n = 100),
                                           SMA(p_null[,4], n = 100)))
    colnames(prob_by_null) <- c("ranks", "20%", "50%", "75%", "90%")
    prob_by_null_melt <- melt(prob_by_null, id.var = "ranks",
                              variable.name = "nullprop")
    prob_plot = ggplot(prob_by_null_melt, aes_string(x = "ranks", y = "value",
                                group = "nullprop", colour = "nullprop")) +
        geom_line(aes(linetype = "nullprop"), size = 1.5) +
        labs(x = "Ranks", y = "p(rank | effect)",
             title = paste0("et = ey = ", ey)) +
        theme(legend.position = "none")


    # ranksProb plots for the different effect sizes------
    p_ef = prob[ , null_index]
    prob_by_effect <- data.frame(ranks, cbind(SMA(p_ef[,1], n = 100),
                                              SMA(p_ef[,2], n = 100),
                                              SMA(p_ef[,3], n = 100)))
    colnames(prob_by_effect) <- c("ranks", "0.6", "2.0", "5.0")
    prob_by_effect_melt <- melt(prob_by_effect, id.var = "ranks",
                                variable.name = "effectSize")
    prob_plot_by_effect = ggplot(prob_by_effect_melt, aes_string(x = "ranks",
                                y = "value", group = "effectSize",
                                             colour = "effectSize")) +
        geom_line(aes(linetype = "effectSize"), size = 1.5) +
        labs(x = "Ranks", y = "p(rank | effect)",
             title = paste0("null = ", null, "%")) +
        theme(legend.position = "none")


    # weights plots for the different null proportions------
    w_null = weight[ , ey_index]
    weight_by_null = data.frame(ranks, cbind(SMA(w_null[,1], n = 100),
                                             SMA(w_null[,2], n = 100),
                                             SMA(w_null[,3], n = 100),
                                             SMA(w_null[,4], n = 100)))
    names(weight_by_null) <- c("ranks", "20%", "50%", "75%", "90%")
    weight_by_null_melt <- melt(weight_by_null, id.var = "ranks",
                                variable.name = "nullprop")
    weight_plot = ggplot(weight_by_null_melt, aes_string(x = "ranks",
                        y = "value", group = "nullprop", colour = "nullprop")) +
        geom_line(size = 1.5) +
        labs(x = "Ranks", y = "Weight", color = "nullprop") +
        theme(legend.direction = "horizontal", legend.position = "bottom")


    # weights plots for the different effect sizes------
    w_ef = weight[ , null_index]
    weight_by_effect <- data.frame(ranks, cbind(SMA(w_ef[,1], n = 100),
                                                SMA(w_ef[,2], n = 100),
                                                SMA(w_ef[,3], n = 100)))
    colnames(weight_by_effect) <- c("ranks", "0.6", "2.0", "5.0")
    weight_by_effect_melt <- melt(weight_by_effect, id.var = "ranks",
                                  variable.name = "effectSize")
    weight_plot_by_effect = ggplot(weight_by_effect_melt, aes_string(x = "ranks",
                                   y = "value", group = "effectSize",
                                                colour = "effectSize")) +
        geom_line(size = 1.5) +
        labs(x = "Ranks", y = "Weight", color = "effectSize (ey)") +
        theme(legend.direction = "horizontal", legend.position = "bottom")


    pp <- plot_grid(prob_plot, prob_plot_by_effect,
                    weight_plot, weight_plot_by_effect,
                    nrow = 2, labels = c("a", "b"), align = 'v')
    title <- ggdraw() + draw_label("Continuous: probability and weight vs. rank")
    plots = plot_grid(title, pp, ncol = 1, rel_heights=c(.1, 1))

    return(plots)

}



