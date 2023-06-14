generate_base_composition_plot <- function(sequences, counts) {
  # Convert sequences to DNAStringSet object
  sequences <- DNAStringSet(sequences)
  
  # Create empty count matrix
  base_counts <- matrix(0, nrow = nchar(sequences[[1]]), ncol = 4)
  colnames(base_counts) <- c("A", "C", "G", "T")
  
  # Update count matrix based on the counts of each sequence
  for (i in seq_along(sequences)) {
    seq_counts <- as.numeric(alphabetFrequency(sequences[i]))
    seq_counts <- seq_counts * counts[i]
    base_counts <- base_counts + matrix(seq_counts, nrow = nchar(sequences[[1]]), ncol = 4)
  }
  
  # Calculate the proportions of each base at each position
  base_proportions <- prop.table(base_counts, margin = 1)
  
  # Create positions vector
  positions <- 1:nchar(sequences[[1]])
  
  # Plot base proportions
  ggplot(data = data.frame(Position = positions), aes(x = Position)) +
    geom_line(aes(y = base_proportions[, "A"], color = "A")) +
    geom_line(aes(y = base_proportions[, "C"], color = "C")) +
    geom_line(aes(y = base_proportions[, "G"], color = "G")) +
    geom_line(aes(y = base_proportions[, "T"], color = "T")) +
    labs(x = "Position", y = "Base Proportion", color = "Base") +
    scale_color_manual(values = c("A" = "blue", "C" = "green", "G" = "red", "T" = "purple")) +
    theme_minimal()
}
