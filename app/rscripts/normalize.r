suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(limma)))
suppressMessages(suppressWarnings(library(edgeR)))
suppressMessages(suppressWarnings(library(DESeq2)))

#command line arguments
args = commandArgs(trailingOnly = TRUE)

input_dir <- args[1]
output_dir <- args[2]
normalization_method <- args[3]
counts_path <- file.path(input_dir, "counts.tsv")
coldata_path <- file.path(input_dir, "coldata.tsv")

counts <- readr::read_tsv(counts_path, show_col_types = FALSE)
metadata <- readr::read_tsv(coldata_path, show_col_types = FALSE)

genes <- counts[,1, drop=FALSE]

# filter out unnecessary columns
cols = colnames(counts)
extraneous_cols = c("hugo_symbol", "entrez_gene_id", "symbol", "gene")
for (x in 1:length(cols)) {
    if(tolower(cols[x]) %in% extraneous_cols) {
        new_counts = counts[, -which(names(counts) == cols[x])]
        counts = new_counts
    }
}

expr_mat <- as.matrix(counts)

# initialize empty matrix as a placeholder for the normalized counts
normalized_counts <- matrix(, nrow = 1, ncol = 1)

if(normalization_method == "mrn") { # median of ratios
    dds <- DESeqDataSetFromMatrix(countData = expr_mat, colData=metadata, design = ~condition)
    dds <- estimateSizeFactors(dds)
    normalized_counts <- counts(dds, normalized=TRUE)
} else if (normalization_method == "tmm") { # trimmed mean of means
    dge_list <- DGEList(counts=expr_mat, samples=metadata)
    dge_list <- calcNormFactors(dge_list, method="TMM")
    normalized_counts <- cpm(dge_list)
}

rownames(normalized_counts) <- NULL

normalized_counts <- as.data.frame(normalized_counts)
final_counts <- do.call(cbind, list(genes, normalized_counts))

output_file = file.path(output_dir, "counts_normalized.tsv")
write.table(final_counts, file=output_file, sep="\t", quote=F, row.names=FALSE)
