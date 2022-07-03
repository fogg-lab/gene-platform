library(tidyverse)
library(BiocParallel)
library(parallel)
library(sva)
library(dplyr)
library(ggplot2)
library(magrittr)

#command line arguments
args = commandArgs(trailingOnly = TRUE)

countspath = args[1]
coldatapath = args[2]
user_directory = args[3]

counts <- readr::read_tsv(countspath)
metadata <- readr::read_tsv(coldatapath)

#creating correct data structures for batch correction
rma_expr <- as.matrix(counts[-1])
rownames(rma_expr) <- counts$Hugo_Symbol
model_m <- model.matrix(~ condition, data = metadata)
batch <- metadata$batch

#getting label for output
labels <- counts$Hugo_Symbol
label <- data.frame(labels)

# batch correction
bc_rma_expr <- ComBat_seq(rma_expr, batch = batch, group = NULL, covar_mod=model_m)

#editing and combining files
final <- data.frame(bc_rma_expr)
final2 <- bind_cols(label,final)

#writing out
write_tsv(final2,paste0(user_directory, "counts_bc.tsv"))
