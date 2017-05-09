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


tests = bottomly$stat
pvals = bottomly$pvalue
filters = bottomly$baseMean

Data <- tibble(test=tests, pval=pvals, filter=filters)	# data of filter covariate and pvlaues

# summary statistics plots of the data
#---------------------------------------------
barlines <- "#1F3552"

hist_test <- ggplot(Data, aes(x = Data$test)) +
    geom_histogram(aes(y = ..density..), binwidth = 1,
                   colour = barlines, fill = "#4271AE") +
    labs(x = "test statistics")

hist_pval <- ggplot(Data, aes(x = Data$pval)) +
    geom_histogram(aes(y = ..density..),
                   colour = barlines, fill = "#4281AE")+
    labs(x = "pvalues")

hist_filter <- ggplot(Data, aes(x = Data$filter)) +
    geom_histogram(aes(y = ..density..),
                   colour = barlines, fill = "#4274AE") +
    labs(x = "filter statistics")

test_filter <- ggplot(Data, aes(x = Data$test, y = Data$filter)) +
    geom_point() + labs(x = "test statstics", y = "filter statistics")
#scale_x_continuous(limits = c(0, 25000), breaks=seq(0, 25000, 10000))

pval_filter <- ggplot(Data, aes(x = rank(-Data$filter), y = -log10(pval))) +
    geom_point()+
    labs(x = "ranks of filters", y = "-log(pvalue)")+
    scale_x_continuous(limits = c(0, 25000), breaks=seq(0, 25000, 10000))

p_ecdf <- ggplot(Data, aes(x = pval)) +
    stat_ecdf(geom = "step")+
    labs(x = "pvalues", title="empirical cumulative distribution")+
    theme(plot.title = element_text(size = rel(.7)))


sum_plots <- grid.arrange(hist_test,  hist_pval, hist_filter, test_filter , pval_filter, p_ecdf,
                          ncol = 3, heights = c(7, 7), top = "bottomly: data summary")
sum_plots


# applying proposed methods------------
n_rej_cont <- NULL
n_rej_bin <- NULL
for(alphaVal in seq(.05, .1, length = 5))
{
    res_cont = opw(pvalue = pvals, filter = filters, test = tests, effectType = "continuous",
                   alpha = alphaVal, tail=2, method = "BH")
    res_bin = opw(pvalue = pvals, filter = filters, test = tests, effectType = "binary",
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
    labs(x = "alpha", y = "discoveries", title = "bottomly data") +
    theme(legend.title = element_blank(), legend.position = "bottom")

p_fdr_bot

# save the results-------------------
save.image("bottomly_data_example.RDATA")

