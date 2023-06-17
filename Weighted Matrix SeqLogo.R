library(ggseqlogo)
# Modifify data(same width needed)
sample <- "" ## input file name
data300 <- read.table(paste0(sample,".txt"), sep="\t", header = TRUE)
## Extract seq coloum
RawSeq <- data300$Sequence
n <- length(RawSeq)
substrings <- character(n)
read <- data300$Reads
reads <- as.numeric(read)
## Trim Rawseq to 15bp
for (i in 1:n) {
  substrings[i] <- substr(RawSeq[i], start = 1, stop = 15)
}

## Generate Matix 
letterMatrix <- function(substrings){
  # Ensure kmers are the same length characters 
  seq.len = sapply(substrings, nchar)
  num_pos = seq.len[1]
  if(! all(seq.len == num_pos)) stop('Sequences in alignment must have identical lengths')
  # Construct matrix of letters
  split = unlist( sapply(substrings, function(seq){strsplit(seq, '')}) )
  
  t(matrix(split, seq.len, length(split)/num_pos) )
}
Matrix_Letter <- letterMatrix(substrings)
Matrix_Letter[,1]
#Calculate the weighted matirx 
total_reads <- sum(reads)
frequency_matrix <- matrix(0, nrow = 4, 15 ,dimnames = list(c("A", "T", "C", "G"), NULL))

for (col in 1:15) {
  countA <- 0
  countT <- 0
  countC <- 0
  countG <- 0
  for (i in 1:length(Matrix_Letter[, col])) {
    if (Matrix_Letter[i, col] == "A") {
      countA <- countA + reads[i]
    } else if (Matrix_Letter[i, col] == "T") {
      countT <- countT + reads[i]
    } else if (Matrix_Letter[i, col] == "C") {
      countC <- countC + reads[i]
    } else if (Matrix_Letter[i, col] == "G") {
      countG <- countG + reads[i]
    }
  }
  frequency_matrix[, col] <- c(countA, countT, countC, countG) / total_reads
}

frequency_matrix
## Seqlogo
png(paste0(sample,"_W.png"), width = 1200, height = 300, units = "px")
ggseqlogo(frequency_matrix,method="prob")
dev.off()
