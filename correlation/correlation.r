# How to use script:
# Rscript correlation.R [path_to_counts_file] [output_path] [correlation_method]
# Supported correlation methods are pearson and spearman
# Example:
#   Rscript correlation.R ./example_data/counts.tsv ./example_output pearson

# load packages
suppressWarnings(library(tidyverse, warn.conflicts = FALSE))
suppressWarnings(library(corrr, warn.conflicts = FALSE))
suppressWarnings(library(ggplot2, warn.conflicts = FALSE))
suppressWarnings(library(corrplot, warn.conflicts = FALSE))

args = commandArgs(trailingOnly = TRUE)

data_path <- args[1]
output_dir <- args[2]
corr_method <- args[3]

output_file <- paste0(output_dir, "/", corr_method, ".pdf")

# read the data
cts <- read_tsv(data_path, col_types=cols())

# filter out unnecessary columns
cols = colnames(cts)
extraneous_cols = c("hugo_symbol", "entrez_gene_id", "symbol", "gene")
for (x in 1:length(cols)) {
    if(tolower(cols[x]) %in% extraneous_cols) {
        new_cts = cts[, -which(names(cts) == cols[x])]
        cts = new_cts
    }
}

num_samples = ncol(cts)

# decide how large the plot should be based on the number of samples
num_samples = ncol(cts)
plot_size = num_samples %/% 2.5

corrplot2 <- function(data,
                      method = "spearman",
                      sig.level = 0.05,
                      order = "original",
                      diag = FALSE,
                      type = "upper",
                      tl.srt = 90,
                      number.font = 1.0,
                      number.cex = 0.65,
                      mar = c(0, 0, 1, 0)) {
    data_incomplete <- data
    data <- data[complete.cases(data), ]
    mat <- cor(data, method = method)
    cor.mtest <- function(mat, method) {
        mat <- as.matrix(mat)
        n <- ncol(mat)
        p.mat <- matrix(NA, n, n)
        diag(p.mat) <- 0
        for (i in 1:(n - 1)) {
          for (j in (i + 1):n) {
            tmp <- cor.test(mat[, i], mat[, j], method = method, exact=FALSE)
            p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
          }
        }
        colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
        p.mat
    }
    p.mat <- cor.mtest(data, method = method)
    col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
    corrplot(mat,
        method = "color", col = col(200), number.font = number.font,
        mar = mar, number.cex = number.cex,
        type = type, order = order,
        addCoef.col = "black", # add correlation coefficient
        tl.col = "black", tl.srt = tl.srt, # rotation of text labels
        # combine with significance level
        p.mat = p.mat, sig.level = sig.level, insig = "blank",
        # hide correlation coefficients on the diagonal
        diag = diag,
        #shrink text labels
        tl.cex=0.65,
        col.lim=c(0,1),
        title=method,
    )
}

pdf(file = output_file, width=plot_size, height=plot_size)

corrplot2(
  data = cts,
  method = corr_method,
  sig.level = 0.5,
  order = "original",
  diag = TRUE,
  type = "full",
  tl.srt = 75
)
