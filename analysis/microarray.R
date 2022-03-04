suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(limma)))
suppressMessages(suppressWarnings(library(yaml)))

args = commandArgs(trailingOnly = TRUE)

user_directory <- paste("user_files/", gsub("[][]","",args[1]), "/", sep="")
counts_filepath <- paste(user_directory, "counts.tsv", sep="")
coldata_filepath <- paste(user_directory, "coldata.tsv", sep="")
config_filepath <- paste(user_directory, "config.yml", sep="")
filter_filepath <- paste(user_directory, "filter.txt", sep="")

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

bayes_fit$coefficients %>% colnames()

fit_de_res_df <- topTable(bayes_fit, coef = paste(condition_col, contrast_level, sep=""), number = nrow(filt_counts_df), adjust.method = adj_method, p.value = padj_thresh) %>%
    rename(lfc = logFC, ave_expr = AveExpr, pval = P.Value, padj = adj.P.Val) %>%
    as_tibble(rownames = "symbol")

colnames(fit_de_res_df) <- c("symbol", "l2fc", "base_avg", "test_stat", "pval", "padj", "B")

write_tsv(fit_de_res_df, paste(user_directory,"output.tsv", sep=""))

if (file.info(filter_filepath)$size != 0) {
    filter_list <- scan(filter_filepath, what="character")
    filtered_df <- fit_de_res_df[fit_de_res_df$symbol %in% filter_list,]
    write_tsv(filtered_df, paste(user_directory, "filter_output.tsv", sep=""))
}