# Name: Ai Vu
# Source: based on the provided microrray analysis pipeline in R
# Description: Given the input files: microrray_coldata(.tsv), micro_array_count(.tsv),
# parameters (yaml), +/- filter_gene(.txt). Output the results of microarray analysis (.tsv)
# for filter and unfiltered genes.
# To use the program, specify the path to input files in lines 40-43

# Loading libraries required for microarray analysis.
# Suppress Messages/Warnings are used to hide long output from Bioconductor packages
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(limma)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(BiocGenerics)))
suppressMessages(suppressWarnings(library(tibble)))
suppressMessages(suppressWarnings(library(yaml)))

# check if the filepath exists
checkFilePath = function(path){
    return(file.exists(path))}

#try to read the .txt files
tryReadTXTFile= function(file){
    tryRead <- try(read.table(file, sep = "", header = FALSE), silent = TRUE)
    if (class(tryRead) != "try-error") {
        name <- read.table(file, sep = "", header = FALSE)
    }else{
        message("File doesn't exist, please check")
    }

    return(name)
}

# try to read the .tsv file
tryReadTSVFile = function(file){
    tryRead <- try(read_tsv(file), silent = TRUE)
    if (class(tryRead) != "try-error") {
        name <- read_tsv(file)
    } else {
        message("File doesn't exist, please check")
    }
    return(name)
}

# Command arguments stored in args variable, used to get Session ID from Flask app
args = commandArgs(trailingOnly = TRUE)
print(args)

# Grabs directory where User session is located, creates variables for each file required for analysis
user_directory <- paste("user_files/", gsub("[][]","",args[1]), "/", sep="")
counts_filepath <- paste(user_directory, "counts.tsv", sep="")
coldata_filepath <- paste(user_directory, "coldata.tsv", sep="")
config_filepath <- paste(user_directory, "config.yml", sep="")
filter_filepath <- paste(user_directory, "filter.txt", sep="")
print(user_directory)


if (checkFilePath(counts_filepath) &
    checkFilePath(coldata_filepath) &
    checkFilePath(config_filepath)){
    #read files
    counts_df <- tryReadTSVFile(counts_filepath)
    coldata_df <- tryReadTSVFile(coldata_filepath)
    input_config <- read_yaml(config_filepath)

    # get the parameters from user
    min_expr = log(input_config$min_expr,2)
    min_prop <- input_config$min_prop
    condition <- input_config$condition
    adj_method <- input_config$adj_method
    use_qual_weights <- input_config$use_qual_weights
    contrast_level <- input_config$contrast_level
    reference_level <- input_config$reference_level
    # check if there is contents in files
    if (nrow(counts_df) > 1 & nrow(coldata_df) > 1){

        # run the main analysis for microrray
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

        if (isTRUE(use_qual_weights)){
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

        #if there is filter-gene file: create another output for unfiltered gene
        if (checkFilePath(filter_filepath)){
            filter_gene_df <- tryReadTXTFile(filter_filepath)
            if(nrow(filter_gene_df) > 1){
                gene_vec <- filter_gene_df[['V1']]
                modified_counts_df <- counts_df[counts_df$symbol %in% gene_vec,]

                modified_filt_counts_df <- modified_counts_df %>%
                    dplyr::filter(rowSums(.[-1] > min_expr) / (ncol(.) - 1) >= min_prop)

                modified_filt_expr <- modified_filt_counts_df %>%
                    column_to_rownames("symbol") %>%
                    as.matrix()

                design <- model.matrix(~ condition, data = coldata_df)
                rownames(design) <- coldata_df$sample_name

                design
                all(colnames(modified_filt_expr) == rownames(design))
                print(design)

                if (isTRUE(use_qual_weights)){
                    modified_qual_weights <- arrayWeights(modified_filt_expr, design = design)
                }else{
                    modified_qual_weights <- 0
                }
                modified_lm_fit <- lmFit(modified_filt_expr, design = design, weights = modified_qual_weights)
                modified_bayes_fit <- eBayes(modified_lm_fit)
                modified_bayes_fit$coefficients %>% colnames()

                fit_de_res_df_filtered_genes <- topTable(modified_bayes_fit, coef = coef_name, number = nrow(modified_filt_counts_df), adjust.method = adj_method) %>%
                rename(lfc = logFC, ave_expr = AveExpr, pval = P.Value, padj = adj.P.Val) %>%
                as_tibble(rownames = "symbol")

                # save the filtered output as tsv
                write.table(fit_de_res_df_filtered_genes,
                    file=paste(user_directory, "filter_output.tsv", sep=""), quote=FALSE, sep='\t', col.names = NA)
                
                # save the unfiltered output as tsv
                write.table(fit_de_res_df_unfiltered_genes,
                    file=paste(user_directory, "output.tsv", sep=""), quote=FALSE, sep='\t', col.names = NA)
            }
        }
    }

}