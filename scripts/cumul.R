library(Biostrings)
library(ggplot2)
library(gplots)

readPac = function(file) {
  s = readDNAStringSet(file)

  pacReadInfo <- data.frame(do.call('rbind', strsplit(names(s),'/')))
  colnames(pacReadInfo) <- c("run", "ZMW", "Position")
  pacReadInfo$readLen = width(s)
  rm(s)

  return(pacReadInfo)
}

G_089 = readPac("G_089/pacbio_assembly/G_089/filtering/data/filtered_subreads.fasta")
G_089$src= c("G_089");
G_4681 = readPac("G_4681/pacbio_assembly/G_4681/filtering/data/filtered_subreads.fasta")
G_4681$src= c("G_4681")
G_690 = readPac("G_690/pacbio_assembly/G_690/filtering/data/filtered_subreads.fasta")
G_690$src= c("G_690")
G_6914 = readPac("G_6914/pacbio_assembly/G_6914/filtering/data/filtered_subreads.fasta")
G_6914$src= c("G_6914");
G_6920 = readPac("G_6920/pacbio_assembly/G_6920/filtering/data/filtered_subreads.fasta")
G_6920$src= c("G_6920")
G_698 = readPac("G_698/pacbio_assembly/G_698/filtering/data/filtered_subreads.fasta")
G_698$src= c("G_698")

experiment = rbind(G_089,G_4681,G_690,G_6914,G_6920,G_698)
write.csv(experiment, file = "experiment.csv")
# experiment <- read.csv("experiment.csv")
sortedExp <- experiment[order(experiment$readLen, decreasing=TRUE),] 
sortedExp.cumLen <- ddply(sortedExp, .(src), transform, cumLen = cumsum(readLen))

pdf(file = "cumulPlot.pdf", height = 5, width = 7)
ggplot(s, aes(readLen, y=cumLen/1000000, colour=src, group=src)) + geom_line() + scale_x_reverse() + ylab("Cumulative Length in MB")
dev.off()
