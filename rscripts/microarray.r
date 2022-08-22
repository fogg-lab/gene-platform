# Ai: this program analyse the microarray part of the gene expression analysis.
# Input: microarray files (coldata (.tsv), countdata(.tsv), config parameters(.yml)).
# Optional Input: filter_gene.txt
# Output: gene expression analysis for unfiltered +/- filtered genes (.tsv)
# and 2 plots - mean variance trend and volcano plot (.png)

suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(limma)))
suppressMessages(suppressWarnings(library(yaml)))
suppressMessages(suppressWarnings(library(ggrepel)))

args = commandArgs(trailingOnly = TRUE)

# Get paths
user_directory <- args[1]
counts_filepath <- file.path(user_directory, "counts.tsv")
coldata_filepath <- file.path(user_directory, "coldata.tsv")
config_filepath <- file.path(user_directory, "config.yml")
filter_filepath <- file.path(user_directory, "filter.txt")

mean_variance_trend = function(fit, filename){
  micro_array_mean_variance_trend_path <- file.path(user_directory, filename)
  png(micro_array_mean_variance_trend_path)
  plot<- plotSA(fit, xlab="Average log-expression", ylab="log2(sigma)", zero.weights=FALSE, pch=16, cex=0.2)
  dev.off()
}

volcano_plot = function(fit, filename){
  micro_array_volcano_path <- file.path(user_directory, filename)
  de <- fit

  de$differential_expression <- "Not sig."
  de$differential_expression[de$l2fc > 0.6 & de$pval < 0.05] <- "Up"
  de$differential_expression[de$l2fc < -0.6 & de$pval < 0.05] <- "Down"
  de$delabel <- NA
  de$delabel[de$differential_expression != "Not sig."] <- de$symbol[de$differential_expression != "Not sig."]

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

min_expr <- config_yml$min_expr
min_prop <- config_yml$min_prop
padj_thresh <- config_yml$padj_thresh
adj_method <- config_yml$adj_method
condition_col <- config_yml$condition
contrast_level <- config_yml$contrast_level
reference_level <- config_yml$reference_level
use_qual_weights <- config_yml$use_qual_weights

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
mean_variance_trend(bayes_fit, "plot_mean_variance_microarray_unfiltered.png")

bayes_fit$coefficients %>% colnames()

fit_de_res_df <- topTable(bayes_fit, coef = paste0(condition_col, contrast_level), number = nrow(filt_counts_df), adjust.method = adj_method, p.value = padj_thresh) %>%
    rename(lfc = logFC, ave_expr = AveExpr, pval = P.Value, padj = adj.P.Val) %>%
    as_tibble(rownames = "symbol")

colnames(fit_de_res_df) <- c("symbol", "l2fc", "base_avg", "test_stat", "pval", "padj", "B")

write_tsv(fit_de_res_df, file.path(user_directory,"output.tsv"))
volcano_plot(fit_de_res_df, "plot_volcano_microarray_unfiltered.png")

if (file.info(filter_filepath)$size != 0) {

    filter_list <- scan(filter_filepath, what="character")
    filtered_df <- fit_de_res_df[fit_de_res_df$symbol %in% filter_list,]
    write_tsv(filtered_df, file.path(user_directory, "filter_output.tsv"))
    volcano_plot(filtered_df, "plot_volcano_microarray_filtered.png")

    #rerun a part of the analysis to make the mean_variance_trend for the filtered gene list
    filtered_count_df <- counts_df[counts_df$symbol %in% filter_list,]
    filtered_counts_df <- filtered_count_df %>%
    filter(rowSums(.[-1] > min_expr) / (ncol(.) - 1) >= min_prop)

    filtered_expr <- filtered_counts_df %>%
    column_to_rownames("symbol") %>%
    as.matrix()

    design <- model.matrix(~ condition, data = coldata_df)
    rownames(design) <- coldata_df$sample_name
    all(colnames(filtered_expr) == rownames(design))

    if (use_qual_weights) {
        qual_weights <- arrayWeights(filtered_expr, design = design)
    } else {
        qual_weights <- NULL
    }
    filtered_lm_fit <- lmFit(filtered_expr, design = design, weights = qual_weights)
    filtered_bayes_fit <- eBayes(filtered_lm_fit)
    mean_variance_trend(filtered_bayes_fit, "plot_mean_variance_microarray_filtered.png")
}
