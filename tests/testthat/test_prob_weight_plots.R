
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
            labs(x = "Ranks", y = "Weight", color = "null prop.")+
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
            labs(x = "Ranks", y = "Weight", color = "effect size (ey)")+
            theme(legend.direction = "horizontal", legend.position = "bottom")


        pp <- plot_grid(prob_plot, prob_plot_by_effect, weight_plot, weight_plot_by_effect,
                        nrow = 2, labels = c("a", "b"), align = 'v')
        title <- ggdraw() + draw_label("Continuous: probability and weight vs. rank")
        plots = plot_grid(title, pp, ncol = 1, rel_heights=c(.1, 1))

        return(plots)

}


