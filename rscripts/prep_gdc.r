### Imports ###
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(TCGAbiolinks)))
suppressMessages(suppressWarnings(library(HDF5Array)))
suppressMessages(suppressWarnings(library(SummarizedExperiment)))

###############################################################################
# Utility functions
###############################################################################
rna_seq_query <- function(p) {
    return(GDCquery(
        project = p,
        data.category = "Transcriptome Profiling",
        data.type = "Gene Expression Quantification",
        workflow.type = "STAR - Counts"
    ))
}

prep_and_save_count_data <- function(rses, label_field, dest_dir, dest_subdir) {
    id_symbol_map <- as_tibble(rowData(rses[[1]]))
    
    # Does the matrix data directory exist? If not, create it
    dir.create(paste0(dest_dir, "/", dest_subdir))
    
    for (name in names(rses)) {
        counts_df <- assays(rses[[name]])[["STAR - Counts"]] %>%
            as_tibble(rownames = "ensembl_gene_id") %>%
            inner_join(id_symbol_map, by = "ensembl_gene_id") %>%
            dplyr::select(ensembl_gene_id, external_gene_name, everything()) %>%
            dplyr::select(-original_ensembl_gene_id)
        
        condition_labels <- rses[[name]][[label_field]]
        coldata_df <- as_tibble(colnames(counts_df)[-c(1:2)]) %>%
            dplyr::rename(sample_name = value) %>%
            mutate(condition = condition_labels) %>%
            mutate(project = name)
        
        write_tsv(counts_df, path = paste0(dest_dir, "/", dest_subdir, "/", n, "_counts.tsv"))
        write_tsv(coldata_df, path = paste0(dest_dir, "/", dest_subdir, "/", n, "_coldata.tsv"))
    }
}

load_RSE_objects <- function(RSE_objects_dest_dir, projects) {
    project_prefixes <- paste0(projects, "_RNA_")
    rses <- vector(mode='list', length=length(project_prefixes))
    for (i in seq_along(project_prefixes)) {
        se <- loadHDF5SummarizedExperiment(dir=RSE_objects_dest_dir, prefix=project_prefixes[[i]])
        rses[[i]] <- se
    }
    names(rses) <- projects
    return(rses)
}

###############################################################################
# Main script
###############################################################################

# Parse command line arguments
args = commandArgs(trailingOnly = TRUE)

data_dir <- args[1]

# Read in the project names
projects <- vector(mode='list', length=length(args)-1)
for(i in seq_along(args[-1])) {
    projects[[i]] <- toupper(args[[i+1]])
}

# Data paths
project_paths <- unlist(map(projects, function(prj) paste0(data_dir, "/", prj)))
biolinks_dir <- "tcga_biolinks_downloads"
TCGA_dest_dir <- paste0(data_dir, "/", biolinks_dir)
RSE_objects_dest_dir <- paste0(data_dir, "/saved_RSE_objects")

### Query the GDC API for the data for each project ###
for (i in seq_along(projects)) {

    project_saved <- FALSE
    query_failed <- FALSE
    tries_left <- 5

    while(!project_saved & tries_left>0) {

        q <- NULL
        out = tryCatch({
            q <- rna_seq_query(projects[i])
            query_failed <- FALSE
        }, warning = function(w) {
            # pass
        }, error = function(e) {
            #message(e)
            query_failed <- TRUE
        }, finally = { })

        if(!query_failed) {
            tryCatch({
                GDCdownload(q, method = "api", directory = TCGA_dest_dir, files.per.chunk = 10)
            }, warning= function(w) {
                # pass
            }, error = function(e) {
                #message(e)
                query_failed <- TRUE
            }, finally = { })
        }

        data = NULL
        if(!query_failed) {
            data = tryCatch({
                GDCprepare(q, directory = TCGA_dest_dir)
            }, warning= function(w) {
                # pass
            }, error = function(e) {
                #message(e)
                query_failed <- TRUE
            }, finally = { })
        }

        if(!query_failed) {
            tryCatch({
                saveHDF5SummarizedExperiment(data, dir = RSE_objects_dest_dir, prefix = paste0(projects[i], "_RNA_"), replace=TRUE)
                project_saved <- TRUE
            }, warning= function(w) {
                # pass
            }, error = function(e) {
                #message(e)
                query_failed <- TRUE
                project_saved <- FALSE
            }, finally = { })
        }

        tries_left <- tries_left - 1
        if (tries_left == 0 & !project_saved) {
            message(paste0("Could not download data for project ", projects[i]))
        }
    }
}

# Create output directory if none exists
dest_dir = data_dir
dest_subdir = "TCGA_RNA_matrix_count_data"
dir.create(paste0(dest_dir, "/", dest_subdir))

# Load RangedSummarizedExperiment objects
rses <- load_RSE_objects(RSE_objects_dest_dir, projects)

for (rse_index in seq_along(rses)){
    
    project_name <- names(rses)[[rse_index]]
    
    message(paste0("Getting counts and coldata for project: ", project_name))
    
    ###
    # Get count matrix
    ###
    message("Parsing count matrix...")
    
    id_symbol_map <- as_tibble(rowData(rses[[rse_index]]))
    
    # Create dataframe for counts
    counts_df <- assays(rses[[project_name]]) %>%
    as_tibble(rownames = "gene_id") %>%
    inner_join(id_symbol_map, by = "gene_id") %>%
    dplyr::select(gene_id, gene_name, everything())
    
    # Get sample names
    sample_names <- as.vector(rownames(colData(rses[[project_name]])))
    
    # Filter out extraneous columns
    sample_names <- as.vector(rownames(colData(rses[[project_name]])))
    cols_to_keep <- c("gene_name", sample_names)
    
    counts_df <- counts_df[,cols_to_keep]
    counts_df <- dplyr::rename(counts_df, Hugo_Symbol = gene_name)
    
    ###
    # Get coldata
    ###

    # Get condition labels - could be under "definition" or "tissue_type"
    condition_labels <- rses[[project_name]][["definition"]]
    if (is.null(condition_labels)) {
        condition_labels <- rses[[project_name]][["tissue_type"]]
    }

    # Create dataframe for coldata
    message("Creating dataframe for coldata...")
    coldata_df <- as_tibble(sample_names) %>%
        dplyr::rename(sample_name = value) %>%
        mutate(condition = condition_labels) %>%

    coldata_df = subset(coldata_df, select = -project_name)

    # Write counts and coldata dataframes to tsv files
    message("Writing counts and coldata files...")
    write_tsv(counts_df, path = paste0(dest_dir, "/", dest_subdir, "/", project_name, "_counts.tsv"))
    write_tsv(coldata_df, path = paste0(dest_dir, "/", dest_subdir, "/", project_name, "_coldata.tsv"))
    
    message("Done.")
}
