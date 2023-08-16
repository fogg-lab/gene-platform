suppressMessages(suppressWarnings(library(readr)))
suppressMessages(suppressWarnings(library(sva)))
suppressMessages(suppressWarnings(library(dplyr)))

args = commandArgs(trailingOnly = TRUE)

# Paths
expression_file = args[1]
clinical_file = args[2]
output_dir = args[3]
data_type = args[4]

# Read in the sample data and expression
expression <- read_tsv(expression_file)
metadata <- read_tsv(clinical_file)

# Prep batch correction
expression_matrix <- as.matrix(expression[-1])
rownames(expression_matrix) <- expression$symbol
model_m <- model.matrix(~ condition, data = metadata)
batch <- metadata$batch

# Get symbol column for adding to output later
symbol <- expression$symbol
symbols <- data.frame(symbol)

# Batch correction
if (data_type == "rnaseq") {
  bc_expression_matrix <- ComBat_seq(expression_matrix, batch = batch, group = NULL, covar_mod=model_m)
} else { # microarray
  bc_expression_matrix <- ComBat(expression_matrix, batch = batch, mod = model_m, ref.batch = 1)
}

# Add symbol column and write the dataframe to a file
bc_expression_df <- data.frame(bc_expression_matrix)
result <- bind_cols(symbols, bc_expression_df)
write_tsv(result, file.path(output_dir, "batch_corrected_counts.tsv"))
