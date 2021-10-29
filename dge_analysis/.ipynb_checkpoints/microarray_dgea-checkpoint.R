library(tidyverse)
library(limma)
# library(WGCNA)
library(qvalue)

# Custom package
library(rutils)


dirs <- rutils::get_dev_directories("../dev_paths.txt")
# vals <- rutils::load_proj_vals("../proj/proj_vals.csv")
lnames <- load("../proj/proj_vals.Rdata")
fold_change <- (2 ^ (vals$deg_lfc_thresh)) %>% str_replace("\\.", "_")

save_results <- FALSE

full_df <- read_tsv(paste0(dirs$data_dir, "/polished_GEO_data/full_bc_GEO_ref_", vals$ref_batch, ".tsv"))
matrisome_df <- load_matrisome_df(paste0(dirs$data_dir, "/matrisome/matrisome_hs_masterlist.tsv")) %>%
    select(gene_symbol, division, category)


dfs <- rutils::split_data(full_df, 1:8, c(1, 9:ncol(full_df)), filter = TRUE, min_expr = vals$min_expr, min_prop = vals$min_prop)
design_df <- dfs$design %>%
    mutate(condition = factor(condition, levels = c("normal", "endometriosis")))
filt_rma_df <- dfs$expr

filt_rma <- filt_rma_df %>%
    column_to_rownames("symbol") %>%
    as.matrix()


phases <- c("proliferative", "secretory_early", "secretory_mid")

for (i in seq_len(length(phases))) {
    phase_grp <- phases[i]
    
    phase_grp_design_df <- design_df %>%
        filter(phase == phase_grp)
    phase_grp_filt_rma <- filt_rma[, phase_grp_design_df$sample_name]

    phase_grp_design <- model.matrix(~ condition, data = phase_grp_design_df)
    rownames(phase_grp_design) <- phase_grp_design_df$sample_name
    stopifnot(all(colnames(phase_grp_filt_rma) == rownames(phase_grp_design)))

    qual_weights <- arrayWeights(phase_grp_filt_rma, design = phase_grp_design)

    phase_grp_lm_fit <- lmFit(phase_grp_filt_rma, design = phase_grp_design, weights = qual_weights)
    phase_grp_bayes_fit <- eBayes(phase_grp_lm_fit)

    phase_grp_de_df <- topTable(phase_grp_bayes_fit, coef = "conditionendometriosis", number = nrow(phase_grp_filt_rma), adjust.method = "BH") %>%
        rename(lfc = logFC, ave_expr = AveExpr, pval = P.Value, padj = adj.P.Val) %>%
        as_tibble(rownames = "symbol") %>%
        mutate(qval = qvalue(pval)$qvalues)

    phase_grp_sig_de_df <- phase_grp_de_df %>%
        filter(abs(lfc) > vals$deg_lfc_thresh, padj < vals$deg_fdr_thresh)

    if (save_results) {
        phase_grp_de_df %>% write_tsv(paste0(dirs$analysis_dir, "/deg/", phase_grp, "_all_deg_res_fc_", fold_change, ".tsv"))
    }

}
