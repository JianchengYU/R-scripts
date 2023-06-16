library(ggseqlogo)
# Modifify data(same width needed)
sample <- " " ## input file name
data300 <- read.table(paste0(sample,".txt"), sep="\t", header = TRUE)
## Extract seq coloum
RawSeq <- data300$Sequence
n <- length(RawSeq)
substrings <- character(n)
## Trim Rawseq to 15bp
for (i in 1:n) {
  substrings[i] <- substr(RawSeq[i], start = 1, stop = 15)
}
## Seqlogo
png(paste0(sample,"_W.png"), width = 1200, height = 300, units = "px")
ggseqlogo(frequency_matrix,method="prob")
dev.off()
