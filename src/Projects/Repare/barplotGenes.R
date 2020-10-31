library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
data <- as.data.frame(read.delim(args[1], header = TRUE))
start <- as.numeric(args[3])
end <- as.numeric(args[4])
by <- as.numeric(args[5])
png(args[2], width = 2000, height = 2500);
if (args[6] == "N") {
  ggplot(data, aes(fill = Group, group = Group, y = Foldchange, x = Sample)) +
      geom_bar(position = "dodge", stat = "identity") + facet_grid(. ~ sgRNA) +
      coord_cartesian(ylim=c(start, end)) + 
      scale_y_continuous(breaks=seq(start, end, by)) +
      # ggtitle(paste0(args[3], "\n")) +
      theme(axis.text.x = element_text(angle=90),
            axis.title.x = element_blank(),
            strip.text.x = element_text(angle=90),
            text = element_text(size=30),
            axis.text = element_text(size=30),
            legend.text = element_text(size=30),
            legend.title = element_text(size=30),
            plot.margin = unit(c(300,5,100,50), "points"))
} else {
  ggplot(data, aes(fill = Group, group = Group, y = Depletion, x = Sample)) +
      geom_bar(position = "dodge", stat = "identity") + facet_grid(. ~ sgRNA) +
      coord_cartesian(ylim=c(start, end)) + 
      scale_y_continuous(breaks=seq(start, end, by)) +
      # ggtitle(paste0(args[3], "\n")) +
      theme(axis.text.x = element_text(angle=90),
            axis.title.x = element_blank(),
            strip.text.x = element_text(angle=90),
            text = element_text(size=30),
            axis.text = element_text(size=30),
            legend.text = element_text(size=30),
            legend.title = element_text(size=30),
            plot.margin = unit(c(300,5,100,50), "points"))

}
dev.off()
