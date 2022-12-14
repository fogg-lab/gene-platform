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

#############
# Main script
#############

# Arguments
args = commandArgs(trailingOnly = TRUE)
data_dir <- args[1]
output_dir <- args[2]
projects <- toupper(args[c(-1,-2)])

# Data paths
biolinks_dir <- "tcga_biolinks_downloads"
TCGA_dest_dir <- file.path(data_dir, biolinks_dir)
RSE_objects_dest_dir <- file.path(data_dir, "saved_RSE_objects")
dir.create(RSE_objects_dest_dir, showWarnings = FALSE)

### Query the GDC API for the data for each project ###
for (i in seq_along(projects)) {

    dir.create(TCGA_dest_dir, showWarnings = FALSE)

    if (file.exists(file.path(RSE_objects_dest_dir, paste0(projects[[i]], "_RNA_assays.h5")))) {
        unlink(TCGA_dest_dir, recursive = TRUE)
        next
    }

    project_saved <- FALSE
    query_failed <- FALSE
    tries_left <- 5

    while(!project_saved & tries_left>0) {

        q <- NULL
        out = tryCatch({
            q <- rna_seq_query(projects[i])
            query_failed <- FALSE
        }, warning = function(w) {
            message(w)
        }, error = function(e) {
            message(e)
            query_failed <- TRUE
        }, finally = { })

        if(!query_failed) {
            tryCatch({
                GDCdownload(q, method = "api", directory = TCGA_dest_dir, files.per.chunk = 10)
            }, warning= function(w) {
                message(w)
            }, error = function(e) {
                message(e)
                query_failed <- TRUE
            }, finally = { })
        }

        data = NULL
        if(!query_failed) {
            message(paste0("Preparing data for ", projects[i]))
            data = tryCatch({
                GDCprepare(q, directory = TCGA_dest_dir)
            }, warning= function(w) {
                message(w)
            }, error = function(e) {
                message(e)
                query_failed <- TRUE
            }, finally = { })
            message(paste0("Data prepared for ", projects[i]))
        }

        if(!query_failed) {
            message(paste0("Saving RSE object for project ", projects[i]))
            tryCatch({
                saveHDF5SummarizedExperiment(data, dir = RSE_objects_dest_dir, prefix = paste0(projects[i], "_RNA_"), replace=TRUE, as.sparse=TRUE)
                project_saved <- TRUE
            }, warning= function(w) {
                message(w)
            }, error = function(e) {
                message(e)
                query_failed <- TRUE
                project_saved <- FALSE
            }, finally = { })
            message(paste0("Saved RSE object for project ", projects[i]))
        }

        tries_left <- tries_left - 1
        if (tries_left == 0 & !project_saved) {
            message(paste0("Could not download data for project ", projects[i]))
        }
    }

    # Remove the downloaded files (only need the RSE objects)
    unlink(TCGA_dest_dir, recursive = TRUE)
}

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
    message("Gathering sample information...")

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
        mutate(project = project_name)

    # Write counts and coldata dataframes to tsv files
    message("Writing counts and coldata files...")
    write_tsv(counts_df, path = paste0(output_dir, "/", tolower(project_name), "_counts_processed.tsv"))
    write_tsv(coldata_df, path = paste0(output_dir, "/", tolower(project_name), "_coldata_processed.tsv"))

    message("Done.")
}
