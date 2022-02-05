library(tidyverse)
library(DESeq2)
library(BiocParallel)
library(parallel)
library(yaml)

args = commandArgs(trailingOnly = TRUE)

write("rnaseq script running", stderr())

user_filepath <- paste("user_files/", gsub("[][]","",args[1]), "/", sep="")
counts_filepath <- paste(user_filepath, "counts.tsv", sep="")
coldata_filepath <- paste(user_filepath, "coldata.tsv", sep="")
config_filepath <- paste(user_filepath, "config.yml", sep="")



n_cores <- detectCores() - 2
BiocParallel::register(MulticoreParam(n_cores))

counts_df <- read_tsv(counts_filepath, col_types=cols())
coldata_df <- read_tsv(coldata_filepath, col_types=cols())

config_yml <- read_yaml(config_filepath)

# Parameters gathered from config file: WORK IN PROGRESS
min_expr <- config_yml$min_expr #0.0
min_prop <- config_yml$min_prop #0.33
padj_thresh <- config_yml$padj_thresh #0.05
adj_method <- config_yml$adj_method #"BH"
condition_col <- config_yml$condition #"condition"
contrast_level <- config_yml$contrast_level #"tumor"
reference_level <- config_yml$reference_level #"healthy"

counts_df <- counts_df %>%
    select(-Entrez_Gene_Id) %>%
    dplyr::rename(symbol = Hugo_Symbol) %>%
    mutate_if(is.numeric, round, 0)
coldata_df <- coldata_df %>%
    mutate(condition = factor(condition, levels = c(reference_level, contrast_level)))

filt_counts_df <- counts_df %>%
    filter(rowSums(.[-1] > min_expr) / (ncol(.) - 1) >= min_prop)
filt_expr <- filt_counts_df %>%
    column_to_rownames("symbol") %>%
    as.matrix()

# filt_expr

coldata <- coldata_df %>%
    column_to_rownames("sample_name")

# coldata

dds <- DESeqDataSetFromMatrix(
    countData = filt_expr,
    colData = coldata,
    design = ~ condition
)

dds_seq <- DESeq(dds, parallel = TRUE)
fit_res <- results(
    dds_seq,
    contrast = c(condition_col, contrast_level, reference_level),
    pAdjustMethod = adj_method,
    alpha = padj_thresh,
    parallel = TRUE
)

dge_res_df <- as_tibble(fit_res, rownames = "symbol")

colnames(dge_res_df) <- c("symbol", "base_avg", "l2fc", "l2fc_se", "test_stat", "pval", "padj")

write("script done running, outputting tsv", stderr())

write_tsv(dge_res_df, paste(user_filepath,"output.tsv", sep=""))