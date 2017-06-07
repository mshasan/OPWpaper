setwd("U:/Documents/My Research (UGA)/Multiple Hypoetheses/Article-1")

library(IHW)
library(IHWpaper)
library(OPWeight)   # proposed library
library(dplyr)      # group_by function
library(reshape2)   # melt
library(ggplot2)    # ggplot
library(cowplot)    # nice plots
library(MASS)       # box-cox transformation
library(gridExtra)  # grid.arrange fucntion


# data processing
#-------------------------
proteomics_file <- system.file("extdata/real_data","science_signaling.csv", package = "IHWpaper")
proteomics_df <- read.csv(proteomics_file, stringsAsFactors = F)
# pvalues were adjusted by BH method so rewrite to obtain orginal pvlaues
proteomics_df$pvalue <- rank(proteomics_df$p1, ties.method="first")*proteomics_df$p1/nrow(proteomics_df)
proteomics_df$test = qnorm(proteomics_df$pvalue, lower.tail = F)
names(proteomics_df)


pval = proteomics_df$pvalue
test = qnorm(pval, lower.tail = FALSE)
test[which(!is.finite(test))] <- NA
filter = proteomics_df$X..peptides


Data <- tibble(test, pval, filter)	# data of filter covariate and pvlaues


bc2 <- boxcox(filter ~ test)
trans2 <- bc2$x[which.max(bc2$y)]
model_prot <- lm(filter^trans2 ~ test)


# summary statistics of the data
#------------------------------------
barlines <- "#1F3552"

hist_test <- ggplot(Data, aes(x = Data$test)) +
        geom_histogram(aes(y = ..density..), binwidth = 1,
	  colour = barlines, fill = "#4271AE") +
		labs(x = "Test statistics")

hist_filter <- ggplot(Data, aes(x = Data$filter)) +
        geom_histogram(aes(y = ..density..),
	  colour = barlines, fill = "#4274AE") +
		labs(x = "Filter statistics")

hist_pval <- ggplot(Data, aes(x = Data$pval)) +
        geom_histogram(aes(y = ..density..),
	  colour = barlines, fill = "#4281AE")+
		labs(x = "P-values")

pval_filter <- ggplot(Data, aes(x = rank(-Data$filter), y = -log10(pval))) +
		geom_point()+
		labs(x = "Ranks of filters", y = "-log(pvalue)")

p_ecdf <- ggplot(Data, aes(x = pval)) +
			stat_ecdf(geom = "step")+
			labs(x = "P-values", title="empirical cumulative distribution")+
			theme(plot.title = element_text(size = rel(.7)))



qqplot.data <- function (vec) # argument: vector of numbers
{
  # following four lines from base R's qqline()
  y <- quantile(vec[!is.na(vec)], c(0.25, 0.75))
  x <- qnorm(c(0.25, 0.75))
  slope <- diff(y)/diff(x)
  int <- y[1L] - slope * x[1L]

  d <- data.frame(resids = vec)

  ggplot(d, aes(sample = resids)) + stat_qq() +
	geom_abline(slope = slope, intercept = int, col="red")+
	labs(x = "Normal quantiles", y = "Fitted values",
	title = expression(paste("Model: ", filter^(-1.414), " ~ ", beta[0] + beta[1]*test)))+
	theme(plot.title = element_text(size = rel(.7)))

}

qqplot <- qqplot.data(model_prot$fit)

p_prot = plot_grid(hist_test, hist_filter, hist_pval, pval_filter, p_ecdf, qqplot,
                  ncol = 3, labels = letters[1:6], align = 'hv')
title_prot <- ggdraw() + draw_label("Proteomics: data summary")
plot_grid(title_prot, p_prot, ncol = 1, rel_heights=c(.1, 1))




# applying proposed methods------------
n_rej_cont <- NULL
n_rej_bin <- NULL
for(alphaVal in seq(.05, .1, length = 5))
    {
    set.seed(123)
    res_cont = opw(pvalue = pval, filter = filter, effectType = "continuous",
                    alpha = alphaVal, method = "BH")
    res_bin = opw(pvalue = pval, filter = filter, effectType = "binary",
              alpha = alphaVal, method = "BH")
    n_rej_cont <- c(n_rej_cont, res_cont$rejections)
    n_rej_bin <- c(n_rej_bin, res_bin$rejections)
    }


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

# save the results-------------------
save.image("proteomics_data_example.RDATA")
