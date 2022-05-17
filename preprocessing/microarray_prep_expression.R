library(tidyverse)
library(affy)
library(affyPLM)
library(affycoretools)
library(hgu133plus2.db)
library(limma)

# Custom package
library(rutils)


dirs <- rutils::get_dev_directories("../dev_paths.txt")
dsets <- c("GSE4888", "GSE6364", "GSE7305", "GSE51981")
dset_suffix_patterns <- c("_[A-Z0-9]*.CEL.gz", ".CEL.gz", ".CEL.gz", "_[A-Z0-9]*.CEL.gz")


cel_files <- list()
data <- list()
data_rma <- list()
data_rma_dfs <- list()


for (idx in seq_len(length(dsets))) {
    cel_files[[dsets[idx]]] <- list.files(paste0(dirs$data_dir, "/raw_GEO_data/untarred/", dsets[idx]), pattern = "*.CEL.gz", full.names = TRUE)
    data[[dsets[idx]]] <- ReadAffy(filenames = cel_files[[dsets[idx]]])
    data_rma[[dsets[idx]]] <- rma(data[[dsets[idx]]])
    # Replace RMA data with annotated RMA data
    data_rma[[dsets[idx]]] <- annotateEset(data_rma[[dsets[idx]]], x = hgu133plus2.db)
    
    # Organize data
    rma_df <- exprs(data_rma[[dsets[idx]]]) %>%
        as_tibble(rownames = "probeid")
    probe2symbol_df <- featureData(data_rma[[dsets[idx]]])@data %>%
        as_tibble() %>%
        rename_with(~ tolower(.))
    # Cols 1-4 are annotation information
    rma_joined_df <- probe2symbol_df %>%
        inner_join(rma_df, by = "probeid")
    # # Average reps among same gene symbol
    # data_rma_dfs[[dsets[idx]]] <- limma::avereps(rma_joined_df[5:ncol(rma_joined_df)], ID = rma_joined_df$symbol) %>%
    #     as_tibble(rownames = "symbol") %>%
    #     rename_with(~ str_remove(., dset_suffix_patterns[idx])) %>%
    #     # Some symbols may be NA -- these are uninformative & can't be used
    #     filter(!is.na(symbol))

    # Collapse gene expression in sample to highest probe expression
    data_rma_dfs[[dsets[idx]]] <- rma_joined_df %>%
        dplyr::select(-c(probeid, entrezid, genename)) %>%
        filter(!is.na(symbol)) %>%
        rename_with(~ str_remove(., dset_suffix_patterns[idx])) %>%
        # Transform to long-form
        gather(sample, value, -symbol) %>%
        group_by(symbol, sample) %>%
        # collapse to max expression
        summarize(expr_max = max(value)) %>%
        # Transform back to wide form
        spread(sample, expr_max) %>%
        ungroup()

    data_rma_dfs[[dsets[idx]]] %>% write_tsv(file = paste0(dirs$data_dir, "/preprocessed_GEO_data/", dsets[idx], ".tsv"))
}

