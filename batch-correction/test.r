library(tidyverse)
library(BiocParallel)
library(parallel)
library(sva)
library(dplyr)

dsets <- c("GSE4888")
batches <- as.list(setNames(seq(1,length(dsets)),dsets))

rma_dfs <- list()
clinical_dfs <- list()

for (ds in dsets) {
    rma_dfs[[ds]] <- read_tsv(paste0("./batch_correction_example_data/preprocessed_GEO_data/", ds, ".tsv"))
    clinical_dfs[[ds]] <- read_tsv(paste0("./batch_correction_example_data/preprocessed_GEO_clinical_data/", ds, ".tsv"))
}

# Handle making data uniform
clinical_dfs$GSE4888 <- clinical_dfs$GSE4888 %>%
    filter(condition == "normal") #%>%
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

