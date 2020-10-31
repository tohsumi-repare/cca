#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2 || length(args) > 3) {
  stop("This script need arguments 1. input file, 2. output head, and 3. if it is the CCA score we are doing (0), non-parametric Z score (1), median depletion (2), or Q3 depletion (3).")
}

data <- as.data.frame(read.delim(args[1], header = TRUE))
output_head <- args[2]
score_choice <- as.numeric(args[3])

library(PRROC)

sub_data <- NULL
if (score_choice == 0) {
  # CCA Score - only analyze positive scores
  sub_data <- data[c("Score", "CCA.score.recall", "CCA.score.precision")]
  names(sub_data) <- c("Score", "Recall", "Precision")
  sub_data <- sub_data[sub_data$Score >= -10, ]
} else if (score_choice == 1) {
  # Z-score
  sub_data <- data[c("Negative.non.parameteric.Z.score", "Z.score.recall", "Z.score.precision")]
  names(sub_data) <- c("Score", "Recall", "Precision")
} else if (score_choice == 2) {
  # Median depletion
  sub_data <- data[c("Test.median.depletion", "Median.test.depletion.recall", "Median.test.depletion.precision")]
  names(sub_data) <- c("Depletion", "Recall", "Precision")
} else if (score_choice == 3) {
  # Q3 depletion
  sub_data <- data[c("Test.Q3", "Q3.test.depletion.recall", "Q3.test.depletion.precision")]
  names(sub_data) <- c("Depletion", "Recall", "Precision")
}
sub_data$FDR <- 1 - sub_data$Precision

if ((score_choice == 0) || (score_choice == 1)) {
  png(paste(output_head, ".score_vs_recall.png", sep = ""), width = 1000, height = 725)
  plot(x=sub_data$Score, y=sub_data$Recall, xlab="Score", ylab="Recall", ylim=c(0, 1))
} else if ((score_choice == 2) || (score_choice == 3)) {
  png(paste(output_head, ".depletion_vs_recall.png", sep = ""), width = 1000, height = 725)
  plot(x=sub_data$Depletion, y=sub_data$Recall, xlab="Depletion", ylab="Recall", ylim=c(0,1))
}
dev.off()

if ((score_choice == 0) || (score_choice == 1)) {
  png(paste(output_head, ".score_vs_FDR.png", sep = ""), width = 1000, height = 725)
  plot(x=sub_data$Score, y=sub_data$FDR, xlab="Score", ylab="FDR (1-Precision)", ylim=c(0,1))
} else if ((score_choice == 2) || (score_choice == 3)) {
  png(paste(output_head, ".depletion_vs_FDR.png", sep = ""), width = 1000, height = 725)
  plot(x=sub_data$Depletion, y=sub_data$FDR, xlab="Depletion", ylab="FDR (1-Precision)", ylim=c(0,1))
}
dev.off()

png(paste(output_head, ".recall_vs_precision.png", sep = ""), width = 1000, height = 725)
plot(x=sub_data$Recall, y=sub_data$Precision, xlab="Recall", ylab="Precision", xlim=c(0,1), ylim=c(0,1))
dev.off()

roc <- roc.curve(scores.class0 = sub_data$Recall, scores.class1 = sub_data$FDR, curve = TRUE)
sink(paste(output_head, ".ROC.tsv", sep = ""))
cat(sprintf("ROC curve\n"))
cat(sprintf("Area under curve:\t%f\n", roc$auc))
cat(sprintf("Curve for scores from\t%f\tto\t%f\n", min(roc$curve[,3]), max(roc$curve[,3])))
sink()
png(paste(output_head, ".ROC.png", sep = ""), width = 1000, height = 725)
plot(roc, xlab="FDR (1-Precison)", ylab="Recall (True positive rate)")
# Old plot below.
# plot(x=sub_data$FDR, y=sub_data$Recall, xlab="FDR (1-Precison)", ylab="Recall (True positive rate)", xlim=c(0,1), ylim=c(0,1))
dev.off()

