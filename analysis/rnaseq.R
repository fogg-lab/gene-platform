# James McGuire
# Ai: Added visualization plot.

# Loading libraries required for RNA Sequence Analysis.
# Suppress Messages/Warnings are used to hide long output from Bioconductor packages
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(DESeq2)))
suppressMessages(suppressWarnings(library(BiocParallel)))
suppressMessages(suppressWarnings(library(parallel)))
suppressMessages(suppressWarnings(library(yaml)))
suppressMessages(suppressWarnings(library(ggpubr)))
suppressMessages(suppressWarnings(library(ggrepel)))

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
  differential_expression = data.frame(fit)
  colnames(differential_expression) <- c("symbol", "baseMean", "log2FoldChange", "l2fc_se", "test_stat", "pval", "padj")
  options(ggrepel.max.overlaps = Inf)
  ggpubr::ggmaplot(differential_expression, main = expression("Group 1" %->% "Group 2"),
                   fdr = 0.05, fc = 2, size = 0.4,
                   palette = c("#B31B21", "#1465AC", "darkgray"),
                   genenames =  as.vector(differential_expression$symbol),
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

  de$differential_expression <- "Not sig."
  de$differential_expression[de$l2fc > 0.6 & de$pval < 0.05] <- "Up"
  de$differential_expression[de$l2fc < -0.6 & de$pval < 0.05] <- "Down"
  de$delabel <- NA
  de$delabel[de$differential_expression != "Not sig."] <- de$symbol[de$differential_expression != "Not sig."]

  plot <- ggplot(data=de, aes(x=l2fc, y=-log10(pval), col=differential_expression, label=delabel)) +
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

min_expr <- config_yml$min_expr
min_prop <- config_yml$min_prop
padj_thresh <- config_yml$padj_thresh
adj_method <- config_yml$adj_method
condition_col <- config_yml$condition
contrast_level <- config_yml$contrast_level
reference_level <- config_yml$reference_level

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

write("Analysis complete, writing output files", stderr())

write_tsv(dge_res_df, paste(user_directory,"output.tsv", sep=""))
volcano_plot(dge_res_df, "volcano_rnaseq_unfiltered.png")
mean_difference(dge_res_df, "mean_difference_rnaseq_unfiltered.png")

if (file.info(filter_filepath)$size != 0) {
    filter_list <- scan(filter_filepath, what="character")
    filtered_df <- dge_res_df[dge_res_df$symbol %in% filter_list,]
    write_tsv(filtered_df, paste(user_directory, "filter_output.tsv", sep=""))
    volcano_plot(filtered_df, "volcano_rnaseq_filtered.png")
    mean_difference(dge_res_df, "mean_difference_rnaseq_filtered.png")
}
