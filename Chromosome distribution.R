library(VariantAnnotation)
library(ggplot2)
library(gridExtra)
# Read the VCF files and create a list of plots
folder_path <- "./results/vcf_filter/"
vcf_files <- list.files(folder_path, pattern = "\\.vcf\\.bgz$", full.names = TRUE)
plots <- list()
for (vcf_file in vcf_files) {
  sample_name <- gsub("\\.vcf\\.bgz$", "", basename(vcf_file))
  vcf <- readVcf(vcf_file)
  # Extract chromosome information
  chromosomes <- seqnames(vcf)
  # Count SNPs per chromosome
  snp_counts <- table(chromosomes)
  # Create a data frame for plotting
  snp_data <- data.frame(Chromosome = names(snp_counts), SNPs = as.numeric(snp_counts))
  # Plot the SNP chromosomal distribution
  plot <- ggplot(snp_data, aes(x = Chromosome, y = SNPs)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(x = "Chromosome", y = "Number of SNPs") +
    ggtitle(paste("SNP Chromosomal Distribution -", sample_name))
    theme(
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    legend.text = element_text(size = 10 ),
    legend.title = element_text(size = 10 ),
    axis.title.y = element_text(size = 14 )
  )
  # Store the plot in the list
  plots[[sample_name]] <- plot
}
# Arrange and save the plots as a single file
output_file <- "./results/chromosomal_distribution_combined.png"
png(output_file, width = 1600, height = 800)
grid.arrange(grobs = plots, ncol = 3)  # Adjust the number of columns as needed
dev.off()
