## ----setup, echo = FALSE-------------------------------------------------
knitr::opts_chunk$set(fig.width = 8.5, fig.height = 7)
knitr::opts_chunk$set(tidy = FALSE, cache = TRUE, autodep = TRUE)

## ----loadLib, message=FALSE, warning=FALSE-------------------------------
library(OPWeight)       # library for the proposed method
library(OPWpaper)       
library(ggplot2)
library(reshape2)       # library for the melt function
library(cowplot)        # plot_grid function
library(dplyr)          # for %>%

## ----bot_data------------------------------------------------------------
bottomly <- system.file("real_data_examples/results", package = "OPWpaper")
setwd(bottomly)
load("bottomly_data_example.RDATA")

## ----summary_bot---------------------------------------------------------
# summary statistics of the data
#------------------------------------
p_bot = plot_grid(hist_test, hist_filter, hist_pval, pval_filter, p_ecdf, qqplot,
                   ncol = 3, labels = letters[1:6], align = 'hv')
title_bot <- ggdraw() + draw_label("Bottomly: data summary")
plot_grid(title_bot, p_bot, ncol = 1, rel_heights=c(.1, 1))

## ----result_bot----------------------------------------------------------
# from IHW paper-------------
alpha = seq(.05, .1, length = 5)
bottomly_file <- system.file("real_data_examples/result_files", "RNAseq_benchmark.Rds", package = "IHWpaper")
rnaseq_data <- readRDS(file=bottomly_file)
panel_a_data <- group_by(rnaseq_data$alpha_df, alpha) %>% summarize(BH = max(bh_rejections), IHW=max(rejections))


# rejected test from different methods-------------
rej_mat_bot_FDR <- data.frame(n_rej_bin, n_rej_cont, panel_a_data)
colnames(rej_mat_bot_FDR) <- c("PRO_bin","PRO_cont", "alpha", "BH","IHW")
rej_mat_bot_FDR2 <- melt(rej_mat_bot_FDR, id.var = "alpha")

p_fdr_bot <- ggplot(rej_mat_bot_FDR2, aes(x = alpha, y = value, group = variable,
                                            colour = variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = expression(bold(paste("Nominal ",alpha))),
	 y = "discoveries", title = "Bottomly data") +
    theme(legend.title = element_blank(), legend.position = "bottom")

p_fdr_bot

## ----prot_data-----------------------------------------------------------
proteomics <- system.file("real_data_examples/results", package = "OPWpaper")
setwd(proteomics)
load("proteomics_data_example.RDATA")

## ----summary_prot--------------------------------------------------------
# summary statistics of the data
#------------------------------------
p_prot = plot_grid(hist_test, hist_filter, hist_pval, pval_filter, p_ecdf, qqplot,
                  ncol = 3, labels = letters[1:6], align = 'hv')
title_prot <- ggdraw() + draw_label("Proteomics: data summary")
plot_grid(title_prot, p_prot, ncol = 1, rel_heights=c(.1, 1))

## ----result_prot---------------------------------------------------------
# from IHW paper-------------
alpha = seq(.05, .1, length = 5)
proteomics_file <- system.file("real_data_examples/result_files", "proteomics_benchmark.Rds", package = "IHWpaper")
proteomics_data <- readRDS(file=proteomics_file)
panel_c_data <- group_by(proteomics_data$alpha_df, alpha) %>%
    summarize(BH = max(bh_rejections), IHW=max(rejections))


# rejected test from different methods-------------
rej_mat_prot_FDR <- data.frame(n_rej_bin, n_rej_cont, panel_c_data)
colnames(rej_mat_prot_FDR) <- c("PRO_bin","PRO_cont", "alpha", "BH","IHW")
rej_mat_prot_FDR2 <- melt(rej_mat_prot_FDR, id.var = "alpha")

p_fdr_prot <- ggplot(rej_mat_prot_FDR2, aes(x = alpha, y = value, group = variable,
                                            colour = variable)) +
    geom_line(aes(linetype = variable), size = 1.5) +
    labs(x = expression(bold(paste("Nominal ",alpha))),
	   y = "Discoveries", title = "Proteomics data") +
    theme(legend.title = element_blank(), legend.position="bottom")

p_fdr_prot

