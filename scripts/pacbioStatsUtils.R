#' A simple function that leverages Biostrings to generate a bunch of stats from dna FASTA files.
#' 
#' @param filename The input fasta file (not vectorized for now)
#' @param type generic, genome or trinity. Not a lot of differences for now
#' @param output.prefix path+filename prefix to append to output files
#' @param write.stats Boolean, should a csv and pdf of stats be written out
#'
#' @export
library(Biostrings)
library(ggplot2)
library(gplots)

nStats = function(readInfo) {
  orderedReadLen = rev(sort(readInfo$readLen))
  cum.fraction = cumsum(orderedReadLen) / sum(orderedReadLen)
  readsNStats = sapply(as.character(seq(10, 90, by = 10)), function(p) {
    orderedReadLen[which(cum.fraction >= as.numeric(p) / 100)[1]]
  }, simplify = TRUE)
  names(readsNStats) = paste0('N', names(readsNStats), " (bp)")
  return(readsNStats)
}

addStats = function(readInfo) {
  readNStats = nStats(readInfo)
  readInfo.stats = c(
    "Nb. Productive ZMWs"                 = nrow(readInfo),
    "Total Unique Sequences Length (bp)"  = sum(readInfo$readLen),
    readNStats
  )
}

writeReport = function(readInfo, prefix) {
  readInfoStats = addStats(readInfo)

   # Create CSV stats file
  write.csv(data.frame(readInfoStats), file = paste0(prefix, ".csv"))

  xlabel = "Sequence Length (bp, log10 scale)"
  ylabel = "Sequence Count"
  # Create histogram plot of sequence length distribution with N50
  sequence.lengths = ggplot(data.frame(length = readInfo$readLen), aes(x = length))
  sequence.lengths = sequence.lengths +
                     geom_histogram(colour = "darkgrey", fill = "white", binwidth = 0.02) +
                     scale_x_log10(
                     breaks = c(200, 500, 1000, 2000, 5000, 10000, 20000, 50000)) +
                     geom_vline(
                       data = data.frame(Legend = paste0("N50 = ", readInfoStats["N50 (bp)"], " bp"), vals = readInfoStats["N50 (bp)"]),
                       aes(xintercept = vals, shape = Legend),
                       colour = "black",
                       show_guide = TRUE) +
                     xlab(xlabel) +
                     ylab(ylabel) +
                     ggtitle("Sequence Length Distribution") +
                     theme(plot.title = element_text(face = "bold"), legend.title = element_blank())

    # Create PDF file of stats + histogram plot
    pdf(file = paste0(prefix, ".pdf"), height = 5, width = 7)
    textplot(data.frame("Statistic" = format(readInfoStats,
                                             digits = 2,
                                             big.mark = ",",
                                             scientific = FALSE,
                                             drop0trailing = TRUE)))
    suppressWarnings(print(sequence.lengths))
    dev.off()

    # Create JPEG image of histogram plot
    jpeg(file = paste0(prefix, ".jpg"), height = 500, width = 700)
    sequence.lengths = sequence.lengths +
                         theme(plot.title = element_text(size = 16),
                               axis.title = element_text(size = 16),
                               axis.text = element_text(size = 14),
                               legend.text = element_text(size = 16))
    print(sequence.lengths)
    dev.off()
}

pacbioReadStats = function(filename, prefix) {
  # Read sequence
  s = readDNAStringSet(filename)

  pacReadInfo <- data.frame(do.call('rbind', strsplit(names(s),'/')))
  colnames(pacReadInfo) <- c("run", "ZMW", "Position")
  pacReadInfo$readLen = width(s)
  rm(s)
  pacReadInfo$ZMW = as.integer(as.character(pacReadInfo$ZMW))

  sortedPacReadInfo = pacReadInfo[order(pacReadInfo$ZMW, pacReadInfo$readLen, decreasing=TRUE),]

  uniquePacReadInfo.all <- sortedPacReadInfo[!duplicated(sortedPacReadInfo$ZMW),]
  uniquePacReadInfo.3kb <- subset(uniquePacReadInfo.all, readLen >= 3000)
  uniquePacReadInfo.7kb <- subset(uniquePacReadInfo.all, readLen >= 7000)

  writeReport(uniquePacReadInfo.all, paste0(prefix, ".uniquePacbioReadInfo.all"))
  writeReport(uniquePacReadInfo.3kb, paste0(prefix, ".uniquePacbioReadInfo.3kb"))
  writeReport(uniquePacReadInfo.7kb, paste0(prefix, ".uniquePacbioReadInfo.7kb"))

}
