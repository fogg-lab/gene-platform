library(tidyverse)
library(BiocParallel)
library(parallel)
library(sva)
library(dplyr)
library(ggplot2)
library(magrittr)


#make variablity 
#list of para
#condition, 
#add a filter for the condition (normal and endometrosis)

dsets <- c("GSE4888", "GSE6364", "GSE7305","GSE51981")
batch <- as.list(setNames(seq(1,length(dsets)),dsets))

rma_dfs <- list()
clinical_dfs <- list()


for (ds in dsets) 
{
    rma_dfs[[ds]] <- read_tsv(paste0("./batch_correction_example_data/preprocessed_GEO_data/", ds, ".tsv"))
    clinical_dfs[[ds]] <- read_tsv(paste0("./batch_correction_example_data/preprocessed_GEO_clinical_data/", ds, ".tsv"))
}

metadata <- readr::read_tsv("./batch_correction_example_data/preprocessed_GEO_clinical_data/GSE4888.tsv")

# Handle making data uniform
# this is specific to Carson's anaylsis 
# either make variablity or make user do it
# take in any user file and then standardize the file 

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
    mutate(batch = as.double(batch[[series]])) %>%
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

#make variablity 
ref_batch <- match("GSE4888", dsets)

bc_rma_expr <- ComBat(rma_expr, batch = batch, mod = model_m, ref.batch = ref_batch)

#file exploer to select gene 
#_bc to select bc

bc_rma_expr_df <- bc_rma_expr %>%
    as_tibble(rownames = "symbol")

bc_rma_expr_df

drops <- c("symbol")

bc_rma_expr[, !(names(bc_rma_expr) %in% drops)]


bc_rma_expr_df

pca_matrix <- bc_rma_expr_df %>%
	column_to_rownames("symbol") %>%
	as.matrix() %>%
	t()

pca <- prcomp(pca_matrix)

#print("jer")

pca_summary <- summary(pca)					
print(pca_summary)

print("meta data")
print(metadata)

