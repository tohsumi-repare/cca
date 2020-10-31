#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(ComplexHeatmap)

args = commandArgs(trailingOnly = TRUE)

ess <- read.delim(args[1], header = TRUE)
ess <- ess[c(1)]
noness <- read.delim(args[2], header = TRUE)
noness <- noness[c(1)]
data <- read.delim(args[3], header = TRUE)
nogenes <- data[-c(1,2)]
datacor <- cor(nogenes)
essdata <- data %>% filter(Gene %in% ess[,1])
nogenes_ess <- essdata[-c(1,2)]
nonessdata <- data %>% filter(Gene %in% noness[,1])
nogenes_noness <- nonessdata[-c(1,2)]
num_guides <- NROW(nogenes)
outhead <- args[4]

norm <- FALSE
ess_noness_line <- 30
line_str <- "The horizontal line is at 30 reads coverage."
if (length(args) >= 5 && args[5] == "norm") {
  norm <- TRUE
  ess_noness_line <- 10
  line_str <- "The horizontal line is at 10 reads coverage."
}


sumdata <- data.frame(value=apply(nogenes,2,sum))
sumdata$key <- rownames(sumdata)
fileout <- paste(outhead, ".Total_reads_per_sample.png", sep = "")
png(fileout, height = 1000, width = 1000)
ggplot(data=sumdata, aes(x=key, y=value, fill=key)) + geom_bar(colour="black", stat="identity") + theme(axis.text.x = element_text(angle = 90), legend.position = "none") + geom_hline(yintercept = 200*num_guides, color = "red") + ylab("Total number of reads") + ggtitle("The horizontal red line is 200X of the number of guides in the library.") + geom_boxplot(width = 0.1, fill = "white")
dev.off()


fileout <- paste(outhead, ".Distribution_of_total_reads_per_sample.png", sep = "")
png(fileout, height = 1000, width = 1000)
ggplot(stack(nogenes), aes(x=ind, y=values, fill=ind)) + geom_violin() + scale_y_log10() + theme(axis.text.x = element_text(angle = 90), legend.position = "none") + geom_hline(yintercept = 30, color = "red") + ggtitle("Distribution of the reads per guide per sample.\nThe horizontal red line is the minimum 30 reads coverage for a given guide.") + geom_boxplot(width = 0.1, fill = "white")
dev.off()


fileout <- paste(outhead, ".Dist_noness_reads_per_sample.png", sep = "")
noness_title <- paste("Distribution of the reads per guide per sample for nonessential genes.", "\n", line_str, sep = "")
png(fileout, height = 1000, width = 1000)
ggplot(stack(nogenes_noness), aes(x=ind, y=values, fill=ind)) + geom_violin() + scale_y_log10() + theme(axis.text.x = element_text(angle = 90), legend.position = "none") + ggtitle(noness_title) + geom_hline(yintercept = ess_noness_line, color = "red") + geom_boxplot(width = 0.1, fill = "white")
dev.off()


fileout <- paste(outhead, ".Dist_ess_reads_per_sample.png", sep = "")
ess_title <- paste("Distribution of the reads per guide per sample for essential genes.", "\n", line_str, sep = "")
png(fileout, height = 1000, width = 1000)
ggplot(stack(nogenes_ess), aes(x=ind, y=values, fill=ind)) + geom_violin() + scale_y_log10() + theme(axis.text.x = element_text(angle = 90), legend.position = "none") + ggtitle(ess_title) + geom_hline(yintercept = ess_noness_line, color = "red") + geom_boxplot(width = 0.1, fill = "white")
dev.off()


fileout <- paste(outhead, ".Clustering_of_samples.png", sep = "")
png(fileout, height = 1000, width = 1000)
draw(Heatmap(datacor, width = unit(12, "cm"), height = unit(12, "cm"), row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8), heatmap_legend_param = list(title = "Pearson correlation", titleposition = "topcenter", direction = "horizontal")), heatmap_legend_side = "top")
dev.off()

