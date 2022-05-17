library(tidyverse)

#' @export
load_matrisome_df <- function(matrisome_list_file) {
    matrisome_df <- readr::read_tsv(matrisome_list_file, quote = "")
    colnames(matrisome_df) <- purrr::map(sub(" ", "_", colnames(matrisome_df)), tolower)
    matrisome_df <- dplyr::select(matrisome_df, gene_symbol, everything()) %>%
        dplyr::filter(division != "Retired")    # Ignore "Retired" matrisome genes
    return(matrisome_df)
}


#' @export
transpose_df <- function(df, future_colnames_col, previous_colnames_col) {
    temp_df <- as.data.frame(df)
    rownames(temp_df) <- df[[future_colnames_col]]
    temp_df <- temp_df %>% dplyr::select(-(!!future_colnames_col))
    t(temp_df) %>% as_tibble(rownames = previous_colnames_col)
}


#' @export
split_data <- function(df, design_idx, gene_idx, future_colnames_col = "sample_name", previous_colnames_col = "symbol", filter = FALSE, min_expr = NULL, min_prop = NULL, expr_samp_idx = 1) {
    design_df <- df %>%
        select(all_of(design_idx))
    expr_df <- df %>%
        select(all_of(gene_idx)) %>%
        rutils::transpose_df(future_colnames_col, previous_colnames_col)
    if (filter) {
        expr_df <- expr_df %>%
            filter(rowSums(.[-expr_samp_idx] > min_expr) / (ncol(.) - 1) >= min_prop)
    }
    list(design = design_df, expr = expr_df)
}


#' @export
to_one_hot <- function(df, col) {
    one_hot <- model.matrix(
        as.formula(paste0("~ ", col, " - 1")),    # We do not want the intercept
        model.frame(~ ., df[col], na.action = na.pass)
    )
    # Don't want white space
    colnames(one_hot) <- gsub(" ", "_", colnames(one_hot))
    # model.matrix() will prepend original column name to each one-hot column
    colnames(one_hot) <- gsub(col, paste0(col, "_"), colnames(one_hot))
    return(tibble::as_tibble(one_hot))
}
