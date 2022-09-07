suppressMessages(suppressWarnings(library(GEOquery)))
suppressMessages(suppressWarnings(library(Biobase)))
suppressMessages(suppressWarnings(library(AnnotationDbi)))
suppressMessages(suppressWarnings(library(limma)))
suppressMessages(suppressWarnings(library(tidyverse)))

# Arguments
args = commandArgs(trailingOnly = TRUE)
data_dir <- args[1]
output_dir <- args[2]
gses <- toupper(args[c(-1,-2)])

saveExpression <- function(expr, series, platform) {
    output_file = paste0(user_dir, "/", series, "_counts_unmapped.tsv")
    write.table(expr, file=output_file, sep="\t", quote=F, row.names=FALSE)
}

saveColdata <- function(coldata, series) {
    output_file = paste0(user_dir, "/", series, "_coldata_processed.tsv")
    write.table(coldata, file=output_file, sep="\t", quote=F, row.names = FALSE)
}

for (series in gses) {
    out = tryCatch({
        # Get GEO series and corresponding platform
        gset <- getGEO(series, GSEMatrix=TRUE, getGPL=TRUE)
        if (length(gset) > 1) {
            idx <- grep(, attr(gset, "names"))
        } else idx <- 1
        gset <- gset[[idx]]
        platform <- gset@annotation

        # Get expression matrix
        expr <- exprs(gset)
        expr <- as.data.frame(expr)
        expr <- data.frame(probe = rownames(expr), expr, row.names = NULL)

        # Get coldata
        coldata <- gset@phenoData@data
        match_pattern <- paste0("channel_count|",
                                "last_update|",
                                "submission|",
                                "status|",
                                "contact|",
                                "supplementary|",
                                "country|",
                                "treatment|",
                                "label|",
                                "tax|",
                                "hyb|",
                                "scan|",
                                "data|",
                                "molecule|",
                                "extract|",
                                "protocol|",
                                "description")
        coldata <- select(coldata, -matches(match_pattern))
        coldata <- as.data.frame(coldata)
        coldata <- data.frame(sample_name = rownames(coldata), coldata, row.names = NULL)

        # Save expression matrix and coldata
        saveExpression(expr, series, platform)
        saveColdata(coldata, series)
    }, error = function(e) {
        print(paste('error:', e))
    })
}
