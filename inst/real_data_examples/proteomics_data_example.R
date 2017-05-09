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


tests = proteomics_df$test
pvals = proteomics_df$pvalue
filters = proteomics_df$X..peptides

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
    labs(x = "ranks of filters", y = "-log(pvalue)")
#scale_x_continuous(limits = c(0, 25000), breaks=seq(0, 25000, 10000))

p_ecdf <- ggplot(Data, aes(x = pval)) +
    stat_ecdf(geom = "step")+
    labs(x = "pvalues", title="empirical cumulative distribution")+
    theme(plot.title = element_text(size = rel(.7)))


sum_plots <- grid.arrange(hist_test,  hist_pval, hist_filter, test_filter , pval_filter, p_ecdf,
             ncol = 3, heights = c(7, 7), top = "Airway: data summary")
sum_plots


# applying proposed methods------------
n_rej_cont <- NULL
n_rej_bin <- NULL
for(alphaVal in seq(.05, .1, length = 5))
    {
    res_cont = opw(pvalue = pvals, filter = filters, test = tests, effectType = "continuous",
                    alpha = alphaVal, method = "BH")
    res_bin = opw(pvalue = pvals, filter = filters, test = tests, effectType = "binary",
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
    labs(x = "alpha", y = "discoveries", title = "Proteomics data") +
    theme(legend.title = element_blank(), legend.position="bottom")

p_fdr_prot

# save the results-------------------
save.image("proteomics_data_example.RDATA")
