#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 2 && length(args) != 3) {
  stop("This script need arguments 1. input file, 2. output head, and optional 3. Number of top genes to plot the distribution of separately - default 300.")
}

library(ggplot2)
library(dplyr)

output_head <- args[2]
data <- as.data.frame(read.delim(args[1], header = TRUE))
data <- data[data$Test.median.depletion >= -1.0, ]
data <- mutate(data, Core.Ess.Or.Noness = ifelse(Core.Ess.Or.Noness == "", "Other", as.character(Core.Ess.Or.Noness)))

topN <- 300
if (length(args) == 2) {
  topN <- as.numeric(args[3])
}
topN <- min(topN, nrow(data))

ess <- data[data$Core.Ess.Or.Noness == "Essential", ]
png(paste(output_head, ".depletion_ess_med.png", sep = ""), width = 1000, height = 725)
ggplot(ess, aes(x = Test.median.depletion, fill = "Test")) + geom_density(alpha = 0.3) + geom_density(data = ess, aes(x = Ctrl.median.depletion, fill = "Control"), alpha=0.3) +  xlab("Median depletion") + xlim(-1, 1)
dev.off()
ess <- ess[order(-ess$Test.median.depletion), ]
ess$test_ind <- seq(1, dim(ess)[1], 1)
ess$test_pct <- 100.0 - 100.0*(ess$test_ind/dim(ess)[1])
ess <- ess[order(-ess$Ctrl.median.depletion), ]
ess$ctrl_ind <- seq(1, dim(ess)[1], 1)
ess$ctrl_pct <- 100.0 - 100.0*(ess$ctrl_ind/dim(ess)[1])
noness <- data[data$Core.Ess.Or.Noness == "Nonessential", ]
png(paste(output_head, ".depletion_noness_med.png", sep = ""), width = 1000, height = 725)
ggplot(noness, aes(x = Test.median.depletion, fill = "Test")) + geom_density(alpha = 0.3) + geom_density(data = noness, aes(x = Ctrl.median.depletion, fill = "Control"), alpha=0.3) +  xlab("Median depletion") + xlim(-1, 1)
dev.off()
noness <- noness[order(-noness$Test.median.depletion), ]
noness$test_ind <- seq(1, dim(noness)[1], 1)
noness$test_pct <- 100.0 - 100.0*(noness$test_ind/dim(noness)[1])
noness <- noness[order(-noness$Ctrl.median.depletion), ]
noness$ctrl_ind <- seq(1, dim(noness)[1], 1)
noness$ctrl_pct <- 100.0 - 100.0*(noness$ctrl_ind/dim(noness)[1])
other <- data[data$Core.Ess.Or.Noness == "Other", ]
png(paste(output_head, ".depletion_other_med.png", sep = ""), width = 1000, height = 725) 
ggplot(other, aes(x = Test.median.depletion, fill = "Test")) + geom_density(alpha = 0.3) + geom_density(data = other, aes(x = Ctrl.median.depletion, fill = "Control"), alpha=0.3) +  xlab("Median depletion") + xlim(-1, 1)
dev.off()
other <- other[order(-other$Test.median.depletion), ]
other$test_ind <- seq(1, dim(other)[1], 1)
other$test_pct <- 100.0 - 100.0*(other$test_ind/dim(other)[1])
other <- other[order(-other$Ctrl.median.depletion), ]
other$ctrl_ind <- seq(1, dim(other)[1], 1)
other$ctrl_pct <- 100.0 - 100.0*(other$ctrl_ind/dim(other)[1])
top_hits <- data[data$Rank <= topN, ]
png(paste(output_head, ".depletion_top_", as.character(topN), "_hits_med.png", sep = ""), width = 1000, height = 725)
ggplot(top_hits, aes(x = Test.median.depletion, fill = "Test")) + geom_density(alpha = 0.3) + geom_density(data = top_hits, aes(x = Ctrl.median.depletion, fill = "Control"), alpha=0.3) +  xlab("Median depletion") + xlim(-1, 1)
dev.off()

data2 <- mutate(data, Core.Ess.Or.Noness = ifelse((Core.Ess.Or.Noness == "Other") & (Rank <= topN), paste("Other top ", as.character(topN), " hit", sep = ""), as.character(Core.Ess.Or.Noness)))
png(paste(output_head, ".test.depletion_med.png", sep = ""), width = 1000, height = 725)
ggplot(data2, aes(x = Test.median.depletion, fill = Core.Ess.Or.Noness)) + geom_density(alpha = 0.3)
dev.off()
png(paste(output_head, ".control.depletion_med.png", sep = ""), width = 1000, height = 725)
ggplot(data2, aes(x = Ctrl.median.depletion, fill = Core.Ess.Or.Noness)) + geom_density(alpha = 0.3)
dev.off()

data3 <- mutate(data, Core.Ess.Or.Noness = ifelse(Rank <= topN, paste("All top ", as.character(topN), " hit", sep = ""), as.character(Core.Ess.Or.Noness)))
png(paste(output_head, ".test.depletion_all_top_hits_as_group_med.png", sep = ""), width = 1000, height = 725)
ggplot(data3, aes(x = Test.median.depletion, fill = Core.Ess.Or.Noness)) + geom_density(alpha = 0.3)
dev.off()
png(paste(output_head, ".control.depletion_all_top_hits_as_group_med.png", sep = ""), width = 1000, height = 725)
ggplot(data3, aes(x = Ctrl.median.depletion, fill = Core.Ess.Or.Noness)) + geom_density(alpha = 0.3)
dev.off()

