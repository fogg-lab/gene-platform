library(tidyverse)
library(limma)

counts_filepath = "data/microarray_endometriosis_counts.tsv"
coldata_filepath = "data/microarray_endometriosis_coldata.tsv"

counts_df <- read_tsv(counts_filepath)
coldata_df <- read_tsv(coldata_filepath)

coldata_df <- coldata_df %>%
    mutate(condition = factor(condition, levels = c("normal", "endometriosis")))

min_expr <- log2(50)
min_prop <- 0.25

filt_counts_df <- counts_df %>%
    filter(rowSums(.[-1] > min_expr) / (ncol(.) - 1) >= min_prop)
filt_expr <- filt_counts_df %>%
    column_to_rownames("symbol") %>%
    as.matrix()

design <- model.matrix(~ condition, data = coldata_df)
rownames(design) <- coldata_df$sample_name

all(colnames(filt_expr) == rownames(design))

qual_weights <- arrayWeights(filt_expr, design = design)

lm_fit <- lmFit(filt_expr, design = design, weights = qual_weights)
bayes_fit <- eBayes(lm_fit)

bayes_fit$coefficients %>% colnames()

fit_de_res_df <- topTable(bayes_fit, coef = "conditionendometriosis", number = nrow(filt_counts_df), adjust.method = "BH") %>%
    rename(lfc = logFC, ave_expr = AveExpr, pval = P.Value, padj = adj.P.Val) %>%
    as_tibble(rownames = "symbol")

fit_de_res_df
