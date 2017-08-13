#' @title Funciton to plot nice plots
#'
#' @description \code{OPWpaper} has stored .RDATA from the simulations. This
#' function will use those simulated data to plot
#'
#' @param x_vec A numeric vector corresponds to the x-axis
#' @param y_matrix A numeric matrix corresponds to the y-axix
#' @param fdr A character vector of ("TRUE" or "FALSE"), determine whether the
#' FDR or FWER will be used, default is FDR.
#' @param power A character vector of ("TRUE" or "FALSE"), determine whether the
#' power will be plotted, default is TRUE
#' @param low_eff_plot A character vector of ("TRUE" or "FALSE"), deteremine
#' whether the power of the low effect sizes will be plotted, default is FALSE
#' @param null Numeric, proportion of the true null if power or FDR/FWER
#' is plotted against the effect sizes
#' @param cv Numeric, coefficient of variation of the test statistics
#' @param ey Numeric, the value of the effect size if power is plotted against
#' the proportion of the true null tests.
#' @param cor Numeric, the correlation coefficient if the figure is for the
#' ranks probability
#' @param figure A character vector of c("ranksProb", "nullPropVsPower",
#' "effectVsFPFP", "CV"), determine the types of figure will be plotted
#'
#' @details
#' \code{OPWeight} package proposed methods to compute the ranks probabilities
#' of the covariate given the test effect size to obtian the optimal power.
#' This function is desigend to plot the power curves under different scenerios.
#' Note that, we alreday simulated power and FDR/FWER for the different scenerios
#' and stored in the packages \code{OPWpaper} as .RDATA. This function will only
#' be able to use those data sets or data with the similar formats.
#'
#'
#' @author Mohamad S. Hasan, \email{shakilmohamad7@gmail.com}
#' @export
#'
#'
#' @return
#' A plot of multiple curves
#'
#' @references Hasan and Schliekelman (2017)
#'
#' @examples
#' # examples from the previously stored .RDATA
#' # plot of power against the effect sizes
#'
#' load(system.file("simulations/results", "simu_fwerPowerFdrPower_cont.RDATA",
#'                 package = "OPWpaper"), envir = environment())
#' ey_vec <- c(seq(0, 1, .2), 2, 3, 5, 8)
#' p_.5_eq_power <- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower2e1,
#'                                 null = 50, figure = "effectVsFPFP")
#'
#' # p_.9_eq_power <- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower4e1,
#' #                                null = 90, figure = "effectVsFPFP")
#' # p_.99_eq_power<- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower5e1,
#' #                                null = 99, figure = "effectVsFPFP")
#'
#' # p_.5_low_ef_eq_power <- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower2e1,
#' #                       null = 50, low_eff_plot = TRUE, figure = "effectVsFPFP")
#' # p_.9_low_ef_eq_power <- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower4e1,
#' #                       null = 90, low_eff_plot = TRUE, figure = "effectVsFPFP")
#' # p_.99_low_ef_eq_power<- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower5e1,
#' #                       null = 99, low_eff_plot = TRUE, figure = "effectVsFPFP")
#'
#' # p_eq_power = plot_grid(p_.5_eq_power, p_.9_eq_power, p_.99_eq_power,
#' #                    p_.5_low_ef_eq_power, p_.9_low_ef_eq_power, p_.99_low_ef_eq_power,
#' #                    ncol = 3, labels = letters[1:3], align = 'hv')
#' # title <- ggdraw() + draw_label("Power: et = ey")
#' # plot_grid(title, p_eq_power, legend, ncol = 1, rel_heights=c(.1, 1, .1))
#'
#' # plot of power against the true propotion of the null
#' # mat_ef.6 <- rbind(FwerPowerFdrPower1f1[13:16, 4], FwerPowerFdrPower2f1[13:16, 4],
#' # FwerPowerFdrPower3f1[13:16, 4], FwerPowerFdrPower4f1[13:16, 4],
#' # FwerPowerFdrPower5f1[13:16, 4])
#' # p_ef.6 <- nice_plots(x_vec = nullProp, y_matrix = mat_ef.6, fdr = TRUE,
#' # power = TRUE, ey = 0.6, figure = "nullPropVsPower")

#'
#===============================================================================
# internal parameters:-----
# x_axis = variable name of the x-axis
# x_lab = label of the x-axis
# y_lab = label of the y_axix
# dat = a data matrix of x-axis and y-axis values
#
# prob0 = probability of rank of a null test by simulation
# prob1 = probability of rank of an alternative test by simulation
# prob0, prob1 =  rank may generate missing valuse because of the large effcet size,
# therefore, to make a matplot equal vectors size are needed. This procedure
# will replace the missing value to make equal sized vector
# probability of rank of a null test
#
# prob0_exact = probability of rank of a null test by exact mehtod
# prob1_exact = probability of rank of an alternative test by exact mehtod
# do not compute probability for a large number of tests
#
# prob0_approx = probability of rank of a null test by normal approximaiton
# prob1_exact = probability of rank of an alternative test by normal approximaiton
#
#===============================================================================
# function to generate nice plots------------

nice_plots <- function(x_vec, y_matrix, fdr = TRUE, power = TRUE,
            low_eff_plot = FALSE, null = NULL, cv = NULL, ey = NULL, cor = NULL,
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
            x_lab = bquote(pi[0])
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
            x_lab = expression(E(tau[i]))
            dat <- data.frame(x_vec, t(y_matrix[row_indx, ]))

       }


        # label the columns--------
        if(figure == "ranksProb"){
            colnames(dat) <- c(x_axis, "CH0","CH1","TH0","TH1")
        } else {
            colnames(dat) <- c(x_axis, "CRW", "BH", "RDW", "IHW")
        }


        # initial plot with melted data-------------
        if(low_eff_plot == FALSE){
            dat_melt <- melt(dat, id.var = x_axis)
            plt <- ggplot(dat_melt, aes_string(x = names(dat_melt)[[1]],
                             y = "value", group = "variable", col = "variable"))
        } else {
            y_lab <- "log10(power)"
            dat <- dat[1:6, ]
            dat[,2:5] <- log10(dat[,2:5])
            dat_melt <- melt(dat, id.var = x_axis)
            plt <- ggplot(dat_melt, aes_string(x = names(dat_melt)[[1]],
                             y = "value", group = "variable", col = "variable"))
        }


        # fixed the tilte of the plot---------------
        if(figure == "ranksProb"){
            titl <- bquote(rho == .(cor))
        } else if(figure == "nullPropVsPower"){
            titl <- bquote(tau[y]==.(ey))
        } else if(figure == "effectVsFPFP"){
            titl <- bquote(pi[0] == .(null)*"%")
        } else {
            titl <- paste0("cv = ", cv)
        }


        # final plot with titles and labels----------
        plt = plt + geom_line(aes_string(linetype = "variable"), size = 1.5) +
           labs(x = x_lab, y = y_lab, title = if(low_eff_plot == FALSE){titl}) +
            theme(legend.position = "none",
                  axis.title.x = element_text(size = rel(1)),
                  axis.title.y = element_text(size = rel(1)))
        return(plt)
    }

