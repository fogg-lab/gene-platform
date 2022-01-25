library(tidyverse)
library(DESeq2)
library(BiocParallel)

counts_filepath = "data/rna_seq_ucec_counts.tsv"
coldata_filepath = "data/rna_seq_ucec_coldata.tsv"

n_cores <- detectCores() - 2
BiocParallel::register(MulticoreParam(n_cores))

counts_df <- read_tsv(counts_filepath)
coldata_df <- read_tsv(coldata_filepath)

counts_df <- counts_df %>%
    select(-Entrez_Gene_Id) %>%
    dplyr::rename(symbol = Hugo_Symbol) %>%
    mutate_if(is.numeric, round, 0)
coldata_df <- coldata_df %>%
    mutate(condition = factor(condition, levels = c("healthy", "tumor")))

min_expr <- 0
min_prop <- 1/3
padj_thresh <- 0.05

filt_counts_df <- counts_df %>%
    filter(rowSums(.[-1] > min_expr) / (ncol(.) - 1) >= min_prop)
filt_expr <- filt_counts_df %>%
    column_to_rownames("symbol") %>%
    as.matrix()

filt_expr

coldata <- coldata_df %>%
    column_to_rownames("sample_name")

coldata

dds <- DESeqDataSetFromMatrix(
    countData = filt_expr,
    colData = coldata,
    design = ~ condition
)

dds_seq <- DESeq(dds, parallel = TRUE)
fit_res <- results(
    dds_seq,
    contrast = c("condition", "tumor", "healthy"),
    pAdjustMethod = "BH",
    alpha = padj_thresh,
    parallel = TRUE
)

dge_res_df <- as_tibble(fit_res, rownames = "symbol")
dge_res_df