library(tidyverse)
library(limma)
library(stringr)
library(BiocGenerics)
library(tibble)
library(config)

dbSafeNames = function(names) {
  names = gsub('[^a-z0-9]+','_',str_to_lower(names))
  names = make.names(names, unique=TRUE, allow_=TRUE)
  names = gsub('.','_',names, fixed=TRUE)
  names
}

checkFilePath = function(path){
    return(file.exists(path)}

tryReadTXTFile= function(file, name){
    tryRead <- try(read.table(file, sep = "", header = FALSE), silent = TRUE)
    if (class(tryRead) != "try-error") {
        name <- read.table(file, sep = "", header = FALSE)
    } else {
        return("File doesn't exist, please check")
    }
}


tryReadTSVFile = function(file, name){
    tryRead <- try(read_tsv(file), silent = TRUE)
    if (class(tryRead) != "try-error") {
        name <- read.table(file, sep = "", header = FALSE)
    } else {
        return("File doesn't exist, please check")
    }
}

microarray_counts_path = "C:\Users\trungn\PycharmProjects\DGEAP\validation\microarray_endometriosis_counts.tsv"
microarray_col_path = "C:\Users\trungn\PycharmProjects\DGEAP\validation\microarray_endometriosis_coldata.tsv"
input_config_path = "C:\Users\trungn\PycharmProjects\DGEAP\validation\config.yml"
filter_gene_path = "C:\Users\trungn\PycharmProjects\DGEAP\validation\filter_gene.txt"

if (checkFilePath(microarray_counts_path) &
    checkFilePath(microarray_col_path) &
    config::is_active("microarray_analysis")
    ){
        tryReadTSVFile(microarray_counts_path, counts_df)
        tryReadTSVFile(microarray_col_path, coldata_df)
        microarray_config <- config::get(
            config = Sys.getenv("R_CONFIG_ACTIVE", "microarray_analysis"),
            file = Sys.getenv("R_CONFIG_FILE", input_config_path ),
            use_parent = TRUE)

   # get the parameters from user
    min_expr <- microarray_config$min_exp
    min_prop <- microarray_config$min_prop
   # padj_thresh <- microarray_config$padj_thresh
    condition <- microarray_config$condition
    adj_method <- microarray_config$adj_method
    use_weight <- microarray_config$use_equal_weights

    contrast_level <- microarray_config$contrast_level
    reference_level <- microarray_config$reference_level

    if (nrow(counts_df) > 1 & nrow(coldata_df) > 1 & nrow(user_input_df) > 1){

        coldata_df <- coldata_df %>%
        mutate(condition = factor(condition, levels = c(contrast_level, reference_level)))

        coef_name = gsub(" ", "", paste(condition, reference_level))
        filt_counts_df <- counts_df %>%
            dplyr::filter(rowSums(.[-1] > min_expr) / (ncol(.) - 1) >= min_prop)

        filt_expr <- filt_counts_df %>%
            column_to_rownames("symbol") %>%
            as.matrix()


        design <- model.matrix(~ condition, data = coldata_df)
        rownames(design) <- coldata_df$sample_name


        design
        all(colnames(filt_expr) == rownames(design))

        if (use_weight){
            qual_weights <- arrayWeights(filt_expr, design = design)
        }else{
            qual_weights <- 0
        }
        lm_fit <- lmFit(filt_expr, design = design, weights = qual_weights)
        bayes_fit <- eBayes(lm_fit)

        bayes_fit$coefficients %>% colnames()

        fit_de_res_df_unfiltered_genes <- topTable(bayes_fit, coef = coef_name, number = nrow(filt_counts_df), adjust.method = adj_method) %>%
            rename(lfc = logFC, ave_expr = AveExpr, pval = P.Value, padj = adj.P.Val) %>%
            as_tibble(rownames = "symbol")

        # save the data table output as tsv
        write.table(fit_de_res_df_unfiltered_genes,
            file='microrray_result_unfiltered.tsv', quote=FALSE, sep='\t', col.names = NA)
    }
    if (checkFilePath(filter_gene_path) & tryReadTXTFile(filter_gene_path, filter_gene_df)
        & nrow(filter_gene_df) > 1)){
        gene_vec <- filter_gene_df[['filter_genes']]
        modified_counts_df <- counts_df[counts_df$symbol %in% gene_vec,]
        modified_filt_counts_df <- modified_counts_df %>%
            dplyr::filter(rowSums(.[-1] > min_expr) / (ncol(.) - 1) >= min_prop)

        modified_filt_expr <- modified_filt_counts_df %>%
            column_to_rownames("symbol") %>%
            as.matrix()
        print('==')


        design <- model.matrix(~ condition, data = coldata_df)
        rownames(design) <- coldata_df$sample_name

        design
        all(colnames(modified_filt_expr) == rownames(design))
        modified_qual_weights <- arrayWeights(modified_filt_expr, design = design)
        modified_lm_fit <- lmFit(modified_filt_expr, design = design, weights = modified_qual_weights)
        modified_bayes_fit <- eBayes(modified_lm_fit)
        modified_bayes_fit$coefficients %>% colnames()
        print("====")
        fit_de_res_df_filtered_genes <- topTable(modified_bayes_fit, coef = coef_name, number = nrow(modified_filt_counts_df), adjust.method = "BH") %>%
        rename(lfc = logFC, ave_expr = AveExpr, pval = P.Value, padj = adj.P.Val) %>%
        as_tibble(rownames = "symbol")
        print("======")
        colnames(fit_de_res_df_filtered_genes) = dbSafeNames(colnames(fit_de_res_df_filtered_genes))

        # save the data table output as tsv
        write.table(fit_de_res_df_filtered_genes,
            file='microrray_result_filtered.tsv', quote=FALSE, sep='\t', col.names = NA)
    }

}




















