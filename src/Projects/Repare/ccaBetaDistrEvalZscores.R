#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2 || length(args) > 9) {
  stop("This script need arguments 1. input file, 2. output file, 3. optional number of entries to use (default = 3000), 4. optional p-value (default = 0.05), 5. optional whether to output a header (default = FALSE), 6. optional output Jenks natural breaks (default = \"\"), 7. the number of Jenks classes (default = 4), 8. optional how much more than the cutoff to use for the Jenks breaks (default = 0), and 9. optional output file of the classes.", call.=FALSE)
}

topN <- 3000
pval <- 0.05
header <- FALSE
jenks <- ""
jenks_breaks <- 4
add_jenks <- 0
temp_strata <- ""
if (length(args) >= 3) {
  topN <- as.numeric(args[3])
}
if (length(args) >= 4) {
  pval <- as.numeric(args[4])
}
if (length(args) >= 5 && args[5] == "TRUE") {
  header <- TRUE
}
if (length(args) >= 6) {
  jenks <- args[6]
  if (length(args) >= 7) {
    jenks_breaks <- as.numeric(args[7])
  }
  if (length(args) >= 8) {
    add_jenks <- as.numeric(args[8])
  }
  if (length(args) >= 9) {
    temp_strata <- args[9]
  }
}


library(fitdistrplus)
library(actuar)
library(classInt)
library(stringr)



# Below is from but modified to include the plotfile.
# http://cainarchaeology.weebly.com/r-function-for-plotting-jenks-natural-breaks-classification.html
plotJenks <- function(data, n=3, plotfile="jenks.png", temp_strata = "", brks.cex=0.70, top.margin=10, dist=5){ 
  df <- data.frame(sorted.values=sort(data, decreasing=TRUE))
  Jclassif <- classIntervals(df$sorted.values, n, style = "jenks") #requires the 'classInt' package
  print(Jclassif)
  test <- jenks.tests(Jclassif) #requires the 'classInt' package
  print(test)
  df$class <- cut(df$sorted.values, unique(Jclassif$brks), labels=FALSE, include.lowest=TRUE) #the function unique() is used to remove non-unique breaks, should the latter be produced. This is done because the cut() function cannot break the values into classes if non-unique breaks are provided
  print(df)
  if (temp_strata != "") {
    write.table(df, temp_strata, quote = FALSE, sep = '\t')
  }
  if(length(Jclassif$brks)!=length(unique(Jclassif$brks))){
    info <- ("The method has produced non-unique breaks, which have been removed. Please, check '...$classif$brks'")
  } else {info <- ("The method did not produce non-unique breaks.")}
  print(info)
  loop.res <- numeric(nrow(df))
  i <- 1
  repeat{
    i <- i+1
    loop.class <- classIntervals(df$sorted.values, i, style = "jenks")
    loop.test <- jenks.tests(loop.class)
    loop.res[i] <- loop.test[[2]]
    if(loop.res[i]>0.9999){
      break
    }
  }
  max.GoF.brks <- which.max(loop.res)
  print(paste("max.GoF.brks =", max.GoF.brks))
  png(plotfile, width = 1000, height = 725)
  plot(x=df$sorted.values, y=c(1:nrow(df)), type="b", main=paste0("Jenks natural breaks optimization; number of classes: ", n), sub=paste0("Goodness of Fit: ", round(test[[2]],4), ". Max GoF (", round(max(loop.res),4), ") with classes:", max.GoF.brks), ylim =c(0, nrow(df)+top.margin), cex=0.75, cex.main=0.95, cex.sub=0.7, ylab="observation index", xlab="value (increasing order)")
  abline(v=Jclassif$brks, lty=3, col="red")
  text(x=Jclassif$brks, y= max(nrow(df)) + dist, labels=sort(round(Jclassif$brks, 2)), cex=brks.cex, srt=90)
  dev.off()
  results <- list("info"=info, "classif" = Jclassif, "breaks.max.GoF"=max.GoF.brks, "class.data" = df)
  return(results)
}


data <- as.data.frame(read.delim(args[1], header = TRUE))
scores <- as.numeric(data$Negative.non.parameteric.Z.score)
n <- min(topN, length(scores))
top_scores <- scores[1:n]
min_top <- min(top_scores)
max_top <- max(top_scores)
transformed_scores <- (top_scores - min_top)/(max_top - min_top)

cullen_file <- args[2]
if (((str_sub(cullen_file, start= -4) == ".tsv") == TRUE) || 
    ((str_sub(cullen_file, start= -4) == ".txt") == TRUE)) {
  cullen_file <- str_sub(cullen_file, end = -5)
}
cullen_file <- paste(cullen_file, ".cullen_and_frey.png", sep = "")
png(cullen_file, width = 1000, height = 750)
descdist(transformed_scores, boot = 1000)
dev.off()

# While AD2L is probably better more accurate, the algorithm cannot often
# evaluate the function at the initial parameters.
# https://asaip.psu.edu/Articles/beware-the-kolmogorov-smirnov-test
fit_beta <- fitdist(transformed_scores, "beta", method = "mge", gof = "CvM")
shape1 <- fit_beta$estimate[1]
shape2 <- fit_beta$estimate[2]

do_stop = FALSE
stop_at = 0
stop_val = -1
last_val = -1
while (!do_stop && stop_at < n) {
  stop_val <- data$Negative.non.parameteric.Z.score[stop_at+1]
  trans_val <- (stop_val - min_top)/(max_top - min_top)
  stop_pval <- 1 - pbeta(trans_val, shape1, shape2) 
  if (stop_pval > pval) {
    do_stop = TRUE
  } else {
    last_val = stop_val
  }
  stop_at = stop_at + 1
}

out_data = c(shape1, shape2, stop_at, last_val)
if (header) {
  out_head = c("Shape1\tShape2\tMax_Rank\tMin_Z_Score");
  write(out_head, args[2], append = TRUE, sep = "\t")
}
write(out_data, args[2], append = TRUE, sep = "\t")

if (jenks != "") {
  top_jenks <- stop_at + add_jenks
  jenks_scores <- scores[1:top_jenks]
  sink(paste(jenks, ".txt", sep = ""), append = FALSE, type = c("output", "message"))  
  plotJenks(jenks_scores, jenks_breaks, paste(jenks, ".png", sep = ""), temp_strata)
  sink()
}
