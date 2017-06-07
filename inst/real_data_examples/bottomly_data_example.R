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
library(DESeq2)     # RNA-seq analysis


# data processing
#-------------------------
bottomly_count_table <- read.table("bottomly_count_table.txt",h=T)
bottomly_phenodata <- read.table("bottomly_phenodata.txt",h=T)
countData <- as.matrix(bottomly_count_table[,-1])		# counts
condition <- factor(bottomly_phenodata[,3])				# strain as factor
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), ~ condition)
dds <- DESeq(dds)
bottomly <- as.data.frame(results(dds))
colnames(bottomly)


pval = bottomly$pvalue
test = qnorm(pval/2, lower.tail = FALSE)
test[which(!is.finite(test))] <- NA
filter = bottomly$baseMean + .0001

# Data <- data.frame(pvals, filters)
# write.csv(Data, "Data.csv")

Data <- tibble(test, pval, filter)	# data of filter covariate and pvlaues

# fite box-cox regression
#--------------------------------
bc <- boxcox(filter ~ test)
trans <- bc$x[which.max(bc$y)]
model_bot <- lm(filter^trans ~ test)


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
		labs(x = "Ranks of filters", y = "-log(pvalue)")+
		scale_x_continuous(limits = c(0, 25000), breaks=seq(0, 25000, 10000))

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
	title = expression(paste("Model: ", filter^.061, " ~ ", beta[0] + beta[1]*test)))+
	theme(plot.title = element_text(size = rel(.7)))

}

qqplot <- qqplot.data(model_bot$fit)

p_bot = plot_grid(hist_test, hist_filter, hist_pval, pval_filter, p_ecdf, qqplot,
                   ncol = 3, labels = letters[1:6], align = 'hv')
title_bot <- ggdraw() + draw_label("Bottomly: data summary")
plot_grid(title_bot, p_bot, ncol = 1, rel_heights=c(.1, 1))





# applying proposed methods------------
n_rej_cont <- NULL
n_rej_bin <- NULL
for(alphaVal in seq(.05, .1, length = 5))
{
    set.seed(123)

    res_cont = opw(pvalue = pval, filter = filter, effectType = "continuous",
                   alpha = alphaVal, tail=2, method = "BH")
    res_bin = opw(pvalue = pval, filter = filter, effectType = "binary",
                  alpha = alphaVal, tail=2, method = "BH")
    n_rej_cont <- c(n_rej_cont, res_cont$rejections)
    n_rej_bin <- c(n_rej_bin, res_bin$rejections)
}


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

# save the results-------------------
save.image("bottomly_data_example.RDATA")

