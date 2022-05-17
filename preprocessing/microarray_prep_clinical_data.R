library(tidyverse)
library(GEOquery)

# Custom package
library(rutils)

dirs <- rutils::get_dev_directories("../dev_paths.txt")
Sys.setenv("VROOM_CONNECTION_SIZE" = as.integer(1e7))

keeper_cols <- list(
    GSE4888 = c("geo_accession", "title", "characteristics_ch1", "characteristics_ch1.1", "source_name_ch1"),
    GSE6364 = c("geo_accession", "title", "characteristics_ch1", "source_name_ch1"),
    GSE7305 = c("geo_accession", "title", "characteristics_ch1", "source_name_ch1", "description"),
    GSE51981 = c("geo_accession", "title", "characteristics_ch1", "characteristics_ch1.1", "characteristics_ch1.2", "source_name_ch1")
)

prepGSE4888 <- function(df) {
    phase_dict <- list(    
        "proliferative" = "proliferative",
        "ambiguous" = "ambiguous_unknown",
        "early secretory" = "secretory_early",
        "mid secretory" = "secretory_mid",
        "late secretory" = "secretory_late"
    )
    condition_dict <- list(
        "nml" = "normal"
    )
    df %>%
        rename(tissue = source_name_ch1, phase = characteristics_ch1, condition = characteristics_ch1.1) %>%
        mutate(
            condition = str_extract(tolower(condition), "leiomyomata|nml|adenomyosis|pain|prolapse|ovarian cyst"),
            phase = str_extract(tolower(phase), "proliferative|early secretory|mid secretory|late secretory|ambiguous"),
            endometriosis_stage = "none"
        ) %>%
        rowwise() %>%
        mutate(
            phase = phase_dict[[phase]],
            condition = ifelse(condition == "nml", condition_dict[[condition]], condition)
        ) %>%
        ungroup() %>%
        select(geo_accession, series, title, condition, endometriosis_stage, phase, tissue)
}


prepGSE6364 <- function(df) {
    phase_dict <- list(    
        "proliferative" = "proliferative",
        "early secretory" = "secretory_early",
        "mid secretory" = "secretory_mid",
        "late secretory" = "secretory_late"
    )
    df %>%
        mutate(
            phase = str_extract(tolower(characteristics_ch1), "proliferative|early secretory|mid secretory|late secretory"),
            condition = str_extract(tolower(characteristics_ch1), "endometriosis|normal"),
            tissue = str_extract(tolower(source_name_ch1), "uterus"),
            endometriosis_stage = ifelse(condition == "endometriosis", "moderate_severe", "none")
        ) %>%
        rowwise() %>%
        mutate(phase = phase_dict[[phase]]) %>%
        ungroup() %>%
        select(geo_accession, series, title, condition, endometriosis_stage, phase, tissue)
}


prepGSE7305 <- function(df) {
    phase_dict <- list(
        "Menstrual phase - Follicular" = "menstrual_follicular",
        "Menstrual phase - Luteal" = "menstrual_luteal"
    )

    condition_dict <- list(
        "disease" = "endometriosis",
        "normal" = "normal"
    )
    df %>%
        rename(condition = characteristics_ch1) %>%
        mutate(
            condition = str_extract(tolower(condition), "disease|normal")
        ) %>%
        rowwise() %>%
        mutate(
            condition = condition_dict[[condition]],
            phase = phase_dict[[description]]
        ) %>%
        ungroup() %>%
        mutate(
            endometriosis_stage = ifelse(condition == "normal", "none", "unknown"),
            tissue = ifelse(condition == "normal", "uterus", "ovary")
        ) %>%
        select(geo_accession, series, title, condition, endometriosis_stage, phase, tissue)
}


prepGSE51981 <- function(df) {
    phase_dict <- list(
        "tissue: Early Secretory Endometrial tissue" = "secretory_early",
        "tissue: Mid-Secretory Endometrial tissue" = "secretory_mid",
        "tissue: Late Secretory Endometrial tissue" = "secretory_late",
        "tissue: Proliferative Endometrial tissue" = "proliferative",
        "tissue: Unknown" = "ambiguous_unknown"
    )

    condition_dict <- list(
        "endometriosis severity: Moderate/Severe" = "endometriosis",
        "endometriosis severity: Minimal/Mild" = "endometriosis",
        "presence or absence of uterine/pelvic pathology: No Uterine Pelvic Pathology" = "normal",
        "presence or absence of uterine/pelvic pathology: Uterine Pelvic Pathology" = "unspecified_pathology"
    )

    endometriosis_stage_dict <- list(
        "endometriosis severity: Moderate/Severe" = "moderate_severe",
        "endometriosis severity: Minimal/Mild" = "minimal_mild",
        "presence or absence of uterine/pelvic pathology: No Uterine Pelvic Pathology" = "none",
        "presence or absence of uterine/pelvic pathology: Uterine Pelvic Pathology" = "none"
    )
    df %>%
        rowwise() %>%
        mutate(
            condition = condition_dict[[characteristics_ch1.2]],
            phase = phase_dict[[characteristics_ch1]],
            endometriosis_stage = endometriosis_stage_dict[[characteristics_ch1.2]]
        ) %>%
        ungroup() %>%
        mutate(
            tissue = str_extract(tolower(characteristics_ch1), "endometrial|unknown"),
            tissue = ifelse(tissue == "endometrial", "uterus", tissue)
        ) %>%
        select(geo_accession, series, title, condition, endometriosis_stage, phase, tissue)
}


prepGSE <- function(df, series) {
    df <- df %>% mutate(series = series)
    if (series == "GSE4888") {
        prepGSE4888(df)
    } else if (series == "GSE6364") {
        prepGSE6364(df)
    } else if (series == "GSE7305") {
        prepGSE7305(df)
    } else if (series == "GSE51981") {
        prepGSE51981(df)
    }
}

series_matrix_files <- list.files(paste0(dirs$data_dir, "/raw_GEO_clinical_data/"), pattern = "*.txt.gz", full.names = TRUE)
series_ids <- purrr::map(series_matrix_files, ~ str_extract(., "GSE[0-9]+")) %>% unlist()
series_matrix_files <- as.list(series_matrix_files)
names(series_matrix_files) <- series_ids
gses <- list()
filt_gse_dfs <- list()

# gse <- getGEO(filename = series_matrix_files$GSE4888, getGPL = FALSE)
for (si in series_ids) {
    gses[[si]] <- getGEO(filename = series_matrix_files[[si]], getGPL = FALSE)
    gse <- gses[[si]]
    keepers <- keeper_cols[[si]]
    gse_df <- gse@phenoData@data[keepers] %>% as_tibble()
    filt_gse_dfs[[si]] <- prepGSE(gse_df, si)
    filt_gse_dfs[[si]] %>% write_tsv(file = paste0(dirs$data_dir, "/preprocessed_GEO_clinical_data/", si, ".tsv"))
}
