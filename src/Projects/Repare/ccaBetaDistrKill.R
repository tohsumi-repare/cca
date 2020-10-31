#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3 && length(args) != 4) {
  stop("This script need arguments 1. input file, 2. output file, 3. plot header, optional 4. top number of genes (default 300)", call.=FALSE)
}

library(fitdistrplus)
library(actuar)
library(ggplot2)

topN <- 300
if (length(args) == 4) {
  topN <- as.numeric(args[4])
}

# Below is needed as a limit of numerical precision.
epsilon <- 1.0e-6

data <- subset(as.data.frame(read.delim(args[1], header = TRUE)),
               select = c("Gene", "Core.Ess.Or.Noness", "Test.Q3"))
colnames(data) <- c("Gene", "Essential", "Depletion")
data$Depletion <- as.numeric(data$Depletion)
min_deplet <- min(data$Depletion)
# We know the max_depletion is 1
data$Transformed <- (2 - data$Depletion)/(3 - min_deplet)

# Fit data on essentials favoring low.
ess_data <- subset(data[data$Essential == "Essential", ], select = "Transformed")
plotfile <- paste(args[3], ".essential.cullen_and_frey.png", sep = "")
png(plotfile, width = 1000, height = 750)
descdist(ess_data$Transformed, boot = 1000)
dev.off()
# While AD2L is probably more accurate, the algorithm cannot often
# evaluate the function at the initial parameters.
# fit_ess <- fitdist(ess_data$Transformed, "beta", method = "mge", gof = "AD2L")
fit_ess <- fitdist(ess_data$Transformed, "beta", method = "mge", gof = "CvM")
plotfile <- paste(args[3], ".essential.denscomp.png", sep = "")
png(plotfile, width = 1000, height = 750)
denscomp(fit_ess)
dev.off()
plotfile <- paste(args[3], ".essential.cdfcomp.png", sep = "")
png(plotfile, width = 1000, height = 750)
cdfcomp(fit_ess)
dev.off()
plotfile <- paste(args[3], ".all_data.density.png", sep = "")
png(plotfile, width = 1000, height = 750)
plot(density(data$Transformed))
dev.off()

shape1 <- fit_ess$estimate[1]
shape2 <- fit_ess$estimate[2]
# Rescale due to numerical instability
min_beta <- pbeta(min(data$Transformed), shape1, shape2)
max_beta <- pbeta(max(data$Transformed), shape1, shape2)
denom <- max_beta - min_beta

ess_kill <- 1 - (pbeta(ess_data$Transformed, shape1, shape2) - min_beta)/denom
max_ess_kill <- max(ess_kill)
min_ess_kill <- min(ess_kill)
sorted_essential_killing <- sort(ess_kill, decreasing = TRUE)
plotfile <- paste(args[3], ".use_this.essential.killing.png", sep = "")
png(plotfile, width = 1000, height = 750)
plot(sorted_essential_killing)
dev.off()
plotfile <- paste(args[3], ".essential.killing_histogram.png", sep = "")
png(plotfile, width = 1000, height = 750)
plot(hist(sorted_essential_killing))
dev.off()

sek <- as.data.frame(sorted_essential_killing)
len_sek <- nrow(sek)
sek$percent_less <- 100*(len_sek - row(sek)[,1])/len_sek
approx_pct <- approxfun(sek$sorted_essential_killing, sek$percent_less)

sek1 <- sek[sek$sorted_essential_killing >= (1.0 - epsilon), ]
sek0 <- sek[sek$sorted_essential_killing <= epsilon, ]
sek01 <- sek[(sek$sorted_essential_killing > epsilon) & 
             (sek$sorted_essential_killing < (1.0 - epsilon)), ]
sek_stat <- data.frame("Statistic" = c("Number of essential genes at 1",
                                       "Percent of essential genes at 1",
                                       "Number of essential genes in the open interval (0, 1)",
                                       "Percent of essential genes in the open interval (0, 1)",
                                       "Number of essential genes at 0",
                                       "Percent of essential genes at 0"),
                       "Value" = c(nrow(sek1),
                                   100.0*nrow(sek1)/len_sek,
                                   nrow(sek01),
                                   100.0*nrow(sek01)/len_sek,
                                   nrow(sek0),
                                   100.0*nrow(sek0)/len_sek))
write.table(sek_stat, 
            paste(args[3], ".essential.killing_statistics.tsv", sep = ""), 
            quote = FALSE, sep = "\t", row.names = FALSE)

data$Killing <- 1 - (pbeta(data$Transformed, shape1, shape2) - min_beta)/denom
out_data = subset(data, select = c("Gene", "Killing"))
for (i in 1:nrow(out_data)) {
  if ((out_data$Killing[i] < max_ess_kill) && 
      (out_data$Killing[i] > min_ess_kill)) {
    out_data$Percent[i] <- approx_pct(out_data$Killing[i])
  } else if (out_data$Killing[i] >= max_ess_kill) {
    out_data$Percent[i] = 100;
  } else {
    out_data$Percent[i] = 0;
  }
}
write.table(out_data, args[2], quote = FALSE, sep = "\t", row.names = FALSE)
