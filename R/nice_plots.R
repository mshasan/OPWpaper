#' @title Funciton to plot nice plots
#'
#' @description \code{OPWpaper} has stored .RDATA from the simulation. This
#' function will use those simulated data to plots
#'
#' @param x_vec a vector of values corresponds to x-axis
#' @param y_matrix a matrix of values correspond to y-axix for multiple plots
#' @param fdr determine whether FDR or FWER will be used, default is FDR
#' @param power determine whether power will be plotted, default is TRUE
#' @param low_eff_plot deteremine whether power of the low effect sizes will be plotted, default is FALSE
#' @param null the proportion of the true null if power or FDR/FWER is plotted against the effect sizes
#' @param cv the value of the coefficient of variation if power is plotted against effect sizes
#' @param ey the value of the effect size if power is plotted against the proportion of the true null
#' @param figure = types of figure will be plotted c("nullPropVsPower", "effectVsFPFP")
#'
#' @details
#' \code{OPWeight} package proposed methods to compute the probabilities
#' of the rank of test given the effect size to obtian the optimal power.
#' This function is desigend to plot the power curves under different scenerio.
#' Note that, we alreday simulated power and FDR/FWER for the different scenerios
#' and stored in the packages *OPWpaper* as .RDATA. This function will only
#' be able to use those data sets or data with the similar formats.
#'
#'
#' @author Mohamad S. Hasan, \email{mshasan@uga.edu}
#' @export
#'
#'
#' @return \code{Data}
#' A plot of multiple curves
#'
#' @references Hasan and Schliekelman (2017)
#'
#' @examples
#' # only just examples from the previously stored .RDATA
#' # plot of power against the effect sizes
#' # p_.5_eq_power <- nice_plots(x_vec = ey_vec, y_matrix = FwerPowerFdrPower2e1,
#' #                                null = 50, figure = "effectVsFPFP")
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
# inpout:----------------
# x_vec = a vector of values corresponds to x-axis
# y_matrix = a matrix of values correspond to y-axix for multiple plots
# fdr = determine whether FDR or FWER will be used, default is FDR
# power = determine whether power will be plotted, default is TRUE
# low_eff_plot = deteremine whether power of the low effect sizes will be plotted, default is FALSE
# null = the proportion of the true null if power or FDR/FWER is plotted against the effect sizes
# cv = the value of the coefficient of variation if power is plotted against effect sizes
# ey = the value of the effect size if power is plotted against the proportion of the true null
# figure = types of figure will be plotted c("nullPropVsPower", "effectVsFPFP")
#
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
# output:---------------
# plots
#
#===============================================================================
# function to generate nice plots------------

nice_plots <- function(x_vec, y_matrix, fdr = TRUE, power = TRUE,
                       low_eff_plot = FALSE, null = NULL, cv = NULL, ey = NULL,
                       figure = c("nullPropVsPower", "effectVsFPFP"))
    {
        if(figure == "nullPropVsPower"){

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


        colnames(dat) <- c(x_axis, "PRO", "BH", "RDW", "IHW")

        if(low_eff_plot == FALSE){
            dat_melt <- melt(dat, id.var = x_axis)
            plt <- ggplot(dat_melt, aes(x = dat_melt[,1], y = dat_melt$value,
                                        group = dat_melt$variable, col = dat_melt$variable))
        } else {
            y_lab <- "log(power)"
            dat <- dat[1:6, ]
            dat_melt <- melt(dat, id.var = x_axis)
            plt <- ggplot(dat_melt, aes(x = dat_melt[,1], y = log(dat_melt$value),
                                        group = dat_melt$variable, col = dat_melt$variable))
        }

        if(is.null(cv) & figure != "nullPropVsPower"){
            titl <- paste0("null = ", null, "%")
        } else if(is.null(cv) & figure == "nullPropVsPower"){
            titl <- paste0("ey = ", ey)
        } else {
            titl <- paste0("cv = ", cv)
        }


        plt = plt + geom_line(aes(linetype = dat_melt$variable), size = 1.5) +
            labs(x = x_lab, y = y_lab, title = if(low_eff_plot == FALSE){titl}) +
            theme(legend.position = "none",
                  axis.title.x = element_text(size = rel(.8)),
                  axis.title.y = element_text(size = rel(.8)))

        return(plt)
    }
