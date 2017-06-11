
#===============================================================================
# function to generate nice plots------------

nice_plots <- function(x_vec, y_matrix, fdr = TRUE, power = TRUE, low_eff_plot = FALSE,
                       null = NULL, cv = NULL, ey = NULL, cor = NULL,
                       figure = c("ranksProb", "nullPropVsPower", "effectVsFPFP", "CV"))
    {
        # configure data sets-------------
        if(figure == "ranksProb"){

            x_axis = "ranks"
            x_lab = "Ranks"
            y_lab = "p(ranks | effect)"
            dat <- data.frame(x_vec, y_matrix)

        } else if(figure == "nullPropVsPower"){

            x_axis <- "nullProp"
            x_lab = "Prop. of null"
            y_lab <- "Power"
            dat <- data.frame(x_vec, y_matrix)

        } else {

            if(fdr == FALSE & power == FALSE){
                row_indx <- 1:4
                y_lab <- "FWER"
            } else if(fdr == FALSE & power == TRUE) {
                row_indx <- 5:8
                y_lab <- "Power"
            } else if(fdr == TRUE & power == FALSE){
                row_indx <- 9:12
                y_lab <- "FDR"
            } else {
                row_indx <- 13:16
                y_lab <- "Power"
            }

            x_axis <- "effectSize"
            x_lab = "Mean filter effect (ey)"
            dat <- data.frame(x_vec, t(y_matrix[row_indx, ]))

       }


        # label the columns--------
        if(figure == "ranksProb"){
            colnames(dat) <- c(x_axis, "FH0","FH1","TH0","TH1")
        } else {
            colnames(dat) <- c(x_axis, "PRO", "BH", "RDW", "IHW")
        }


        # initial plot with melted data-------------
        if(low_eff_plot == FALSE){
            dat_melt <- melt(dat, id.var = x_axis)
            plt <- ggplot(dat_melt, aes_string(x = names(dat_melt)[[1]], y = "value",
                                        group = "variable", col = "variable"))
        } else {
            y_lab <- "log(power)"
            dat <- dat[1:6, ]
            dat[,2:5] <- log(dat[,2:5])
            dat_melt <- melt(dat, id.var = x_axis)
            plt <- ggplot(dat_melt, aes_string(x = names(dat_melt)[[1]], y = "value",
                                        group = "variable", col = "variable"))
        }


        # fixed the tilte of the plot---------------
        if(figure == "ranksProb"){
            titl <- paste0("cor = ", cor)
        } else if(figure == "nullPropVsPower"){
            titl <- paste0("ey = ", ey)
        } else if(figure == "effectVsFPFP"){
            titl <- paste0("null = ", null, "%")
        } else {
            titl <- paste0("cv = ", cv)
        }


        # final plot with titles and labels----------
        plt = plt + geom_line(aes_string(linetype = "variable"), size = 1.5) +
            labs(x = x_lab, y = y_lab, title = if(low_eff_plot == FALSE){titl}) +
            theme(legend.position = "none",
                  axis.title.x = element_text(size = rel(.8)),
                  axis.title.y = element_text(size = rel(.8)))

        return(plt)
    }
