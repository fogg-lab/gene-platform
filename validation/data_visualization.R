# Name: Ai Vu
# Description:
# source: https://rpkgs.datanovia.com/ggpubr/reference/diff_express.html


library(ggrepel)
library(ggpubr)

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

filename = "C:/Users/trungn/PycharmProjects/DGEAP/validation/microrray_result_unfiltered.tsv"

# volcano plot: displays statistical significance (-log10 P value)
# versus magnitude of change (log2 fold change)
volcano_plot = function(fit){
    #read file
    de <- fit
    # add a column of NAs
    de$diffexpressed <- "NO"

    de$diffexpressed[de$lfc > 0.6 & de$pval < 0.05] <- "UP"

    de$diffexpressed[de$lfc < -0.6 & de$pval < 0.05] <- "DOWN"

    # Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
    de$delabel <- NA
    de$delabel[de$diffexpressed != "NO"] <- de$symbol[de$diffexpressed != "NO"]
    png("C:/Users/trungn/PycharmProjects/DGEAP/validation/volcano_micrrarray_unfiltered.png")
    ggplot(data=de, aes(x=lfc, y=-log10(pval), col=diffexpressed, label=delabel)) +
        geom_point() +
        theme_minimal() +
#         geom_text_repel(box.padding = 0.5, max.overlaps = 10) +
        scale_color_manual(values=c("blue", "black", "red")) +
        geom_vline(xintercept=c(-0.6, 0.6), col="red") +
        geom_hline(yintercept=-log10(0.05), col="red")
    dev.off()
}

#mean difference (MD) plot displays log2 fold change versus average log2 expression values
#data frame containing the columns baseMean, log2FoldChange, and padj. Rows are genes

mean_difference = function(fit){
    diff_express = fit
    colnames(diff_express) <- c("symbol", "baseMean", "log2FoldChange", "l2fc_se", "test_stat", "pval", "padj")
    # Add rectangle around labels
    #fdr: Accepted false discovery rate for considering genes as differentially expressed.
    #fc: the fold change threshold. Only genes with a fold change >= fc and padj <= fdr are
    # considered as significantly differentially expressed.
   ggpubr::ggmaplot(diff_express, main = expression("Group 1" %->% "Group 2"),
   fdr = 0.05, fc = 2, size = 0.4,
   palette = c("#B31B21", "#1465AC", "darkgray"),
   genenames =  as.vector(diff_express$symbol)
   legend = "top", top = 20,
   font.label = c("bold", 11), label.rectangle = TRUE,
   font.legend = "bold",
   font.main = "bold",
   select.top.method = c("padj", "fc"),
   xlab = "Log2 mean expression",
   ylab = "Log2 fold change",
   ggtheme = ggplot2::theme_minimal())

}


# check the mean-variance relationship of the expression data, after fitting a linear model
# red line is mean-variance trend approximation
# blue line is constant variance approximation
#save the plot as png file
mean_variance_trend = function(fit){
    png("C:/Users/trungn/PycharmProjects/DGEAP/validation/microarray_mean_variance_trend.png")
    plot <- plotSA(fit, xlab="Average log-expression", ylab="log2(sigma)", zero.weights=FALSE, pch=16, cex=0.2)
    dev.off()
}




