# Ai: this program analyse the microarray part of the gene expression analysis.
# Input: microarray files (coldata (.tsv), countdata(.tsv), config parameters(.yml)).
# Optional Input: filter_gene.txt
# Output: gene expression analysis for unfilter +/- filtered genes (.tsv)
# and 2 plots (mean microarray_mean_variance_trend(.png) and volcano_plot(.png)

MICRO_ARRAY_MEAN_VARIANCE_TREND_IMAGE_FILE = "mean_variance_trend_plot_UNFILTERED_microarray.png"
MICRO_ARRAY_VOLCANO_IMAGE_FILE = "volcano_plot_UNFILTERED_microarray.png"
#FILTERED_MICRO_ARRAY_MEAN_VARIANCE_TREND_IMAGE_FILE = "filtered_microarray_mean_variance_trend.png"
FILTERED_MICRO_ARRAY_VOLCANO_IMAGE_FILE = "volcano_plot_FILTERED_microarray.png"

suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(limma)))
suppressMessages(suppressWarnings(library(yaml)))
suppressMessages(suppressWarnings(library(ggrepel)))


args = commandArgs(trailingOnly = TRUE)

user_directory <- paste("user_files/", gsub("[][]","",args[1]), "/", sep="")
counts_filepath <- paste(user_directory, "counts.tsv", sep="")
coldata_filepath <- paste(user_directory, "coldata.tsv", sep="")
config_filepath <- paste(user_directory, "config.yml", sep="")
filter_filepath <- paste(user_directory, "filter.txt", sep="")

mean_variance_trend = function(fit, filename){
  micro_array_mean_variance_trend_path <- paste(user_directory, filename)
  png(filename = "micro_array_mean_variance_trend_path", width = 756, height = 756)
  plot<- plotSA(fit, xlab="Average log-expression", ylab="log2(sigma)", zero.weights=FALSE, pch=16, cex=0.2)
  dev.off()
}

volcano_plot = function(fit, filename){
  micro_array_volcano_path <- paste(user_directory, filename)
  de <- fit

  de$differential_expression <- "Not significance"
  de$differential_expression[de$l2fc > 0.6 & de$pval < 0.05] <- "UP"
  de$differential_expression[de$l2fc < -0.6 & de$pval < 0.05] <- "DOWN"
  de$delabel <- NA
  de$delabel[de$differential_expression != "NO"] <- de$symbol[de$differential_expression != "Not significance"]

  plot <- ggplot(data=de, aes(x=l2fc, y=-log10(pval), col=differential_expression, label=delabel)) +
    geom_point() +
    theme_minimal() +
    scale_color_manual(values=c("blue", "black", "red")) +
    geom_vline(xintercept=c(-0.6, 0.6), col="red") +
    geom_hline(yintercept=-log10(0.05), col="red")
  ggsave(micro_array_volcano_path,  width = 20, height = 20, units = "cm")
}


counts_df <- read_tsv(counts_filepath, col_types=cols())
coldata_df <- read_tsv(coldata_filepath, col_types=cols())

config_yml <- read_yaml(config_filepath)

min_expr <- config_yml$min_expr #log2(50) = 5.6438
min_prop <- config_yml$min_prop #0.25
padj_thresh <- config_yml$padj_thresh #0.05
adj_method <- config_yml$adj_method #"BH"
condition_col <- config_yml$condition #"condition"
contrast_level <- config_yml$contrast_level #"endometriosis"
reference_level <- config_yml$reference_level #"normal"
use_qual_weights <- config_yml$use_qual_weights #TRUE

coldata_df <- coldata_df %>%
    mutate(condition = factor(condition, levels = c(reference_level, contrast_level)))

filt_counts_df <- counts_df %>%
    filter(rowSums(.[-1] > min_expr) / (ncol(.) - 1) >= min_prop)
filt_expr <- filt_counts_df %>%
    column_to_rownames("symbol") %>%
    as.matrix()

design <- model.matrix(~ condition, data = coldata_df)
rownames(design) <- coldata_df$sample_name

all(colnames(filt_expr) == rownames(design))

if (use_qual_weights) {
    qual_weights <- arrayWeights(filt_expr, design = design)
} else {
    qual_weights <- NULL
}

lm_fit <- lmFit(filt_expr, design = design, weights = qual_weights)
bayes_fit <- eBayes(lm_fit)
mean_variance_trend(bayes_fit, MICRO_ARRAY_MEAN_VARIANCE_TREND_IMAGE_FILE)
bayes_fit$coefficients %>% colnames()

fit_de_res_df <- topTable(bayes_fit, coef = paste(condition_col, contrast_level, sep=""), number = nrow(filt_counts_df), adjust.method = adj_method, p.value = padj_thresh) %>%
    rename(lfc = logFC, ave_expr = AveExpr, pval = P.Value, padj = adj.P.Val) %>%
    as_tibble(rownames = "symbol")

colnames(fit_de_res_df) <- c("symbol", "l2fc", "base_avg", "test_stat", "pval", "padj", "B")

write_tsv(fit_de_res_df, paste(user_directory,"output.tsv", sep=""))
volcano_plot(fit_de_res_df, MICRO_ARRAY_VOLCANO_IMAGE_FILE)

if (file.info(filter_filepath)$size != 0) {
    filter_list <- scan(filter_filepath, what="character")
    filtered_df <- fit_de_res_df[fit_de_res_df$symbol %in% filter_list,]
    write_tsv(filtered_df, paste(user_directory, "filter_output.tsv", sep=""))
    volcano_plot(filtered_df, FILTERED_MICRO_ARRAY_VOLCANO_IMAGE_FILE)
}
