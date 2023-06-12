#Library 
library(VariantAnnotation)
library(ggplot2)
#Function to comparsion
compare_variants <- function(sample_names) {
  ## Create data frame
  all_data <- data.frame(Sample = character(), Variants = numeric())
  ## Iterate each sample name
  for (sample_name in sample_names) {
    ## Define file paths
    before_gatk_file <- paste0("./results/create_vcf/", sample_name, "_gatk.vcf")
    after_gatk_file <- paste0("./results/vcf_filter/", sample_name, "_gatk_filter.vcf.bgz")
    before_bcf_file <- paste0("./results/create_vcf_BCFtool/", sample_name, "_bcf.vcf")
    after_bcf_file <- paste0("./results/vcf_filter_BCFtools/", sample_name, "_bcf_filter.vcf.bgz")
    ## Read VCF files and get variant counts
    before_gatk_count <- length(as(readVcf(before_gatk_file, genome = "Ath"), "VRanges")[, 1])
    after_gatk_count <- length(as(readVcf(after_gatk_file, genome = "Ath"), "VRanges")[, 1])
    before_bcf_count <- length(as(readVcf(before_bcf_file, genome = "Ath"), "VRanges")[, 1])
    after_bcf_count <- length(as(readVcf(after_bcf_file, genome = "Ath"), "VRanges")[, 1])
    ## Create data frame
    data <- data.frame(Sample = c(paste0(sample_name,"_gatk_before"), paste0(sample_name, "_gatk_filtered"),
                                  paste0(sample_name,"_bcf_before" ), paste0(sample_name,"_bcf_filtered" )),
                       Variants = c(before_gatk_count, after_gatk_count,before_bcf_count,after_bcf_count))
    # Append result to the overall result
    all_data <- rbind(all_data, data)
  }
  ## Plot the results
  plot <- ggplot(all_data, aes(x = Sample, y = Variants, fill = Sample)) +
    geom_bar(stat = "identity", width = 0.75, color = "black") +
    labs(x = " ", y = "Variant Number", title = paste0("Variant Count Comparison")) +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    legend.text = element_text(size = 10 ),
    legend.title = element_text(size = 10 ),
    axis.title.y = element_text(size = 14 )
  )
  ggsave(plot, file = "./results/variant_filter_plot.png", width = 1600*length(sample_names), height = 1600, units = "px")
  ##Save data to file
  write.csv(all_data, file = "./results/variant_filter.csv", row.names = FALSE)
  return(list(png_file = "./results/variant_filter_plot.png", csv_file = "./results/variant_filter.csv"))
  dev.off()
}
