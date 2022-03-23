library(tidyverse)
library(BiocParallel)
library(parallel)
library(sva)

#Carson's custom library
library(rutils)


##############
#Loading Data#
##############

dsets <- c("GSE4888", "GSE6364", "GSE7305","GSE51981")
batches <- as.list(setNames(seq(1,length(dsets)),dsets))

dirs <- rutils::get_dev_directories("../dev_path.txt")
lnames <- load("../proj/proj_vals.Rdata")

matrisome_df <- rutils::load_matrisome_df(paste0(dirs$data_dir, "/matrisome/matrisome_hs_masterlist.tsv")) %>%
    select(gene_symbol, division, category)
rma_dfs <- list()
clinical_dfs <- list()

for (ds in dsets) {
    rma_dfs[[ds]] <- read_tsv(paste0(dirs$data_dir, "/preprocessed_GEO_data/", ds, ".tsv"))
    clinical_dfs[[ds]] <- read_tsv(paste0(dirs$data_dir, "/preprocessed_GEO_clinical_data/", ds, ".tsv"))
}


# Handle making data uniform
clinical_dfs$GSE4888 <- clinical_dfs$GSE4888 %>%
    filter(condition == "normal") %>%
    rename(sample_name = geo_accession)
rma_dfs$GSE4888 <- rma_dfs$GSE4888 %>%
    select(one_of(c("symbol", clinical_dfs$GSE4888$sample_name)))

clinical_dfs$GSE6364 <- clinical_dfs$GSE6364 %>%
    rename(sample_name = geo_accession)

clinical_dfs$GSE51981 <- clinical_dfs$GSE51981 %>%
    filter(!(tissue == "unknown" | phase == "ambiguous_unknown" | condition == "unspecified_pathology")) %>%
    mutate(title = as.character(title)) %>%
    rename(sample_name = geo_accession)
rma_dfs$GSE51981 <- rma_dfs$GSE51981 %>%
    select(one_of(c("symbol", clinical_dfs$GSE51981$sample_name)))

all_clinical_df <- bind_rows(clinical_dfs$GSE4888, clinical_dfs$GSE6364, clinical_dfs$GSE51981) %>%
    rowwise() %>%
    mutate(batch = as.double(batches[[series]])) %>%
    ungroup() %>%
    mutate(
        condition = factor(condition, levels = c("normal", "endometriosis")),
        phase = factor(phase, levels = c("proliferative", "secretory_early", "secretory_mid", "secretory_late"))
    )


all_rma_df <- bind_cols(rma_dfs$GSE4888, rma_dfs$GSE6364[-1], rma_dfs$GSE51981[-1])

rma_expr <- as.matrix(all_rma_df[-1])
rownames(rma_expr) <- all_rma_df$symbol
model_m <- model.matrix(~ condition, data = all_clinical_df)
batch <- all_clinical_df$batch


#Batch Correction

ref_batch <- match("GSE4888", dsets)

bc_rma_expr <- ComBat(rma_expr, batch = batch, mod = model_m, BPPARAM = mc_param, ref.batch = ref_batch)

bc_rma_expr_df <- bc_rma_expr %>%
    as_tibble(rownames = "symbol")

bc_rma_expr_df