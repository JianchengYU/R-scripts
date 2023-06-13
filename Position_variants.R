library(ggplot2)
library(gridExtra)
# Set the folder path
folder_path <- "./results/vcf_anno/"
# Get the list of tsv files in the folder
tsv_files <- list.files(folder_path, pattern = "\\.tsv$", full.names = TRUE)
# Create an empty list to store the plots
plots <- list()

# Iterate over each tsv file
for (tsv_file in tsv_files) {
  # Read the tsv file
  data <- read.delim(tsv_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  # Extract relevant columns
  var_ids <- data$VARID
  alt_depths <- data$altDepth
  chromosomes <- sub(":.*", "", var_ids)
  positions <- as.numeric(sub(".*:", "", sub("_.*", "", var_ids)))
  mutations <- sub(".*_", "", var_ids)
  # Create a data frame
  df <- data.frame(Chromosome = chromosomes, Position = positions, Mutation = mutations, Depth = alt_depths)
  # Compute the relative position within each chromosome
  df <- df %>%
    group_by(Chromosome) %>%
    mutate(RelativePosition = Position - min(Position) + 1)
  # Create the plot
  plot <- ggplot(df, aes(x = RelativePosition, y = Depth, color = Chromosome)) +
    geom_point() +
    xlab("Relative Position on Chromosome") +
    ylab("Depth") +
    ggtitle(paste("Position vs. Depth on Different Chromosomes -", basename(tsv_file)))
  theme(
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold")
     )
  # Add the plot to the list
  plots[[tsv_file]] <- plot
}

# Combine the plots into a single grid
grid_plot <- do.call(grid.arrange, plots)
# Save the grid plot as a PNG file
output_file <- "./results/Variants Distribution.png"
ggsave(grid_plot, filename = output_file, width = 16, height = 8)
dev.off()
