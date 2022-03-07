# James McGuire
# Ai: Added visualization plot. Please notice that the rna_seq data needs to preprocess before
# generating the plot. Preprocess data haven't done yet.

RNA_SEQ_MEAN_DIFF_IMAGE_FILE = "rna_mean_difference.png"
RNA-SEQ_VOLCANO_IMAGE_FILE = "rna_seq_volcano.png"
FILTERED_RNA_SEQ_MEAN_DIFF_IMAGE_FILE = "filtered_rna_mean_difference.png"
FILTERED_RNA_SEQ_VOLCANO_IMAGE_FILE = "filtered_rna_seq_volcano.png"

# write("Loading libraries: tidyverse, DESeq2, BiocParallel, parallel, yaml", stderr())
# Loading libraries required for RNA Sequence Analysis.
# Suppress Messages/Warnings are used to hide long output from Bioconductor packages
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(DESeq2)))
suppressMessages(suppressWarnings(library(BiocParallel)))
suppressMessages(suppressWarnings(library(parallel)))
suppressMessages(suppressWarnings(library(yaml)))

# write("Libraries loaded, preforming analysis...", stderr())

# Command arguments stored in args variable, used to get Session ID from Flask app
args = commandArgs(trailingOnly = TRUE)

# Grabs directory where User session is located, creates variables for each file required for analysis
user_directory <- paste("user_files/", gsub("[][]","",args[1]), "/", sep="")
counts_filepath <- paste(user_directory, "counts.tsv", sep="")
coldata_filepath <- paste(user_directory, "coldata.tsv", sep="")
config_filepath <- paste(user_directory, "config.yml", sep="")
filter_filepath <- paste(user_directory, "filter.txt", sep="")

mean_difference = function(fit, name){
  mean_diff_rna_seq_path <- paste(user_directory, name)
  diff_express = data.frame(fit)
  colnames(diff_express) <- c("symbol", "baseMean", "log2FoldChange", "l2fc_se", "test_stat", "pval", "padj")
  options(ggrepel.max.overlaps = Inf)
  ggpubr::ggmaplot(diff_express, main = expression("Group 1" %->% "Group 2"),
                   fdr = 0.05, fc = 2, size = 0.4,
                   palette = c("#B31B21", "#1465AC", "darkgray"),
                   genenames =  as.vector(diff_express$symbol),
                   legend = "top", top = 20,
                   font.label = c("bold", 11), label.rectangle = TRUE,
                   font.legend = "bold",
                   font.main = "bold",
                   select.top.method = c("padj", "fc"),
                   xlab = "Log2 mean expression",
                   ylab = "Log2 fold change",
                   ggtheme = ggplot2::theme_minimal()
                   )
  ggsave(mean_diff_rna_seq_path,  width = 20, height = 20, units = "cm")
}

volcano_plot = function(fit, name){
  rna_volcano_path <- paste(user_directory, name)
  de <- fit
  de$diffexpressed <- "NO"
  de$diffexpressed[de$l2fc > 0.6 & de$pval < 0.05] <- "UP"
  de$diffexpressed[de$l2fc < -0.6 & de$pval < 0.05] <- "DOWN"
  de$delabel <- NA
  de$delabel[de$diffexpressed != "NO"] <- de$symbol[de$diffexpressed != "NO"]
  plot <- ggplot(data=de, aes(x=l2fc, y=-log10(pval), col=diffexpressed, label=delabel)) +
    geom_point() +
    theme_minimal() +
    scale_color_manual(values=c("blue", "black", "red")) +
    geom_vline(xintercept=c(-0.6, 0.6), col="red") +
    geom_hline(yintercept=-log10(0.05), col="red")
  ggsave(rna_volcano_path, width = 20, height = 20, units = "cm")
}

n_cores <- detectCores() - 2
BiocParallel::register(MulticoreParam(n_cores))

counts_df <- read_tsv(counts_filepath, col_types=cols())
coldata_df <- read_tsv(coldata_filepath, col_types=cols())

config_yml <- read_yaml(config_filepath)

# Parameters gathered from config file: WORK IN PROGRESS
min_expr <- config_yml$min_expr #0.0
min_prop <- config_yml$min_prop #0.33
padj_thresh <- config_yml$padj_thresh #0.05
adj_method <- config_yml$adj_method #"BH"
condition_col <- config_yml$condition #"condition"
contrast_level <- config_yml$contrast_level #"tumor"
reference_level <- config_yml$reference_level #"healthy"

counts_df <- counts_df %>%
    select(-Entrez_Gene_Id) %>%
    dplyr::rename(symbol = Hugo_Symbol) %>%
    mutate_if(is.numeric, round, 0)
coldata_df <- coldata_df %>%
    mutate(condition = factor(condition, levels = c(reference_level, contrast_level)))

filt_counts_df <- counts_df %>%
    filter(rowSums(.[-1] > min_expr) / (ncol(.) - 1) >= min_prop)
filt_expr <- filt_counts_df %>%
    column_to_rownames("symbol") %>%
    as.matrix()

# filt_expr

coldata <- coldata_df %>%
    column_to_rownames("sample_name")

# coldata

dds <- DESeqDataSetFromMatrix(
    countData = filt_expr,
    colData = coldata,
    design = ~ condition
)

dds_seq <- DESeq(dds, parallel = TRUE)
fit_res <- results(
    dds_seq,
    contrast = c(condition_col, contrast_level, reference_level),
    pAdjustMethod = adj_method,
    alpha = padj_thresh,
    parallel = TRUE
)

dge_res_df <- as_tibble(fit_res, rownames = "symbol")

colnames(dge_res_df) <- c("symbol", "base_avg", "l2fc", "l2fc_se", "test_stat", "pval", "padj")

# write("Analysis complete, writing output files", stderr())

write_tsv(dge_res_df, paste(user_directory,"output.tsv", sep=""))
volcano_plot(dge_res_df, RNA-SEQ_VOLCANO_IMAGE_FILE)
mean_difference(dge_res_df, RNA_SEQ_MEAN_DIFF_IMAGE_FILE)

if (file.info(filter_filepath)$size != 0) {
    filter_list <- scan(filter_filepath, what="character")
    filtered_df <- dge_res_df[dge_res_df$symbol %in% filter_list,]
    write_tsv(filtered_df, paste(user_directory, "filter_output.tsv", sep=""))
    volcano_plot(filtered_df, FILTERED_RNA_SEQ_VOLCANO_IMAGE_FILE)
    mean_difference(dge_res_df, FILTERED_RNA_SEQ_MEAN_DIFF_IMAGE_FILE)
}