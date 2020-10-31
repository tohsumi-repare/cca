#!/usr/bin/env Rscript

library(ComplexHeatmap)

args = commandArgs(trailingOnly = TRUE)

data <- read.delim(args[1], header = TRUE)
data <- data[-c(1,2)]
# Even though the data isn't normally distributed, we want to use Pearson
# because we want to see if the columns are _linearly_ correlated.   
# Non-parametric correlations correlate monotonicity, 
# e.g. if x=(1:100) and y=exp(x), then spearman cor = 1, but Pearson = 0.25.
# For a normalized table, we want Pearson.
res <- cor(data)
out_res <- round(res, 2)
fileout <- paste(args[2], ".tsv", sep = "")
write.table(out_res, fileout, quote = FALSE, sep = "\t", col.names = NA)
fileout <- paste(args[2], ".png", sep = "")
hw = 1000
png(fileout, height = hw, width = hw)
draw(Heatmap(res, width = unit(12, "cm"), height = unit(12, "cm"), row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8), heatmap_legend_param = list(title = "Pearson correlation", titleposition = "topcenter", direction = "horizontal")), heatmap_legend_side = "top")
dev.off()
