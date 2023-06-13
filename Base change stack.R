library(ggplot2)
library(gridExtra)
# Set the folder path
folder_path <- "./results/vcf_anno"
# Get the list of tsv files in the folder
tsv_files <- list.files(folder_path, pattern = "\\.tsv$", full.names = TRUE)
# Create an empty list to store the plots
plots <- list()

# Process each TSV file
for (tsv_file in tsv_files) {
  # Read annotated file
  sample_name <- gsub("\\_gatk_filter_anno.tsv$", "", basename(tsv_file))
  data <- read.delim(tsv_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  # Extract mutations
  mutations <- sub(".*_", "", data$VARID)
  # Define valid SNP types
  valid_snp_types <- c("A/T", "A/C", "A/G", "T/A", "T/C", "T/G", "C/A", "C/G", "C/T", "G/A", "G/T", "G/C")
  # Modify mutation names
  mutations[!(mutations %in% valid_snp_types)] <- "indel"
  # Recalculate mutation counts
  mutation_counts <- table(mutations)
  # Create a data frame for plotting
  plot_data <- data.frame(mutation = names(mutation_counts),
                          count = as.numeric(mutation_counts),
                          sample = sample_name)
  # Plot
  bar_chart <- ggplot(plot_data, aes(x = sample, y = count, fill = mutation)) +
    geom_col(position = "stack") +
    labs(title = paste("Mutation Distribution -", basename(sample_name))) +
    scale_colour_viridis_d()
  # Store the plot in the list
  plots[[basename(tsv_file)]] <- bar_chart
}

# Combine the plots into a single grid
combined_plot <- do.call(grid.arrange, c(plots, ncol = 3))
# Save the combined plot to a file
output_file <- "./results/Base Change_combined_stack.png"
ggsave(output_file, combined_plot, width = 12, height = 6)
