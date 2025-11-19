#!/usr/bin/env Rscript
# =============================================================================
# GWAS Manhattan Plot with Gene Annotation - Usage Examples
# =============================================================================

# Set working directory
setwd("~/Documents/Software/manhattan2D/")

# Load the annotated plotting script
cat("Loading gene annotation functions...\n")
source("gwas_manhattan_annotated.R")

# =============================================================================
# Example 1: Basic Usage with Gene Annotation
# =============================================================================

cat("\n", strrep("=", 70), "\n")
cat("Example 1: Basic Gene Annotation\n")
cat(strrep("=", 70), "\n\n")

# Simplest usage - just provide trait name and GTF file
result <- plot_manhattan_with_genes(
  trait_name = "betaine",
  gtf_file = "Homo_sapiens.GRCh38.gtf.gz"  # Or your GTF file path
)

# Output files:
# - betaine.GWAS.annotated.tiff (Manhattan plot with gene labels)
# - betaine.QTL_regions.txt (Table of QTL regions and nearest genes)

# =============================================================================
# Example 2: Customize Significance Threshold
# =============================================================================

cat("\n", strrep("=", 70), "\n")
cat("Example 2: Genome-wide Significance (p < 5e-8)\n")
cat(strrep("=", 70), "\n\n")

result <- plot_manhattan_with_genes(
  trait_name = "betaine",
  gtf_file = "genome.gtf.gz",
  p_threshold = 5e-8,        # Genome-wide significance
  merge_distance = 1e6,      # 1 Mb to merge nearby QTLs
  max_labels = 20            # Label up to 20 QTL regions
)

# =============================================================================
# Example 3: Label Only Top N Most Significant QTLs
# =============================================================================

cat("\n", strrep("=", 70), "\n")
cat("Example 3: Label Top 10 Most Significant QTLs\n")
cat(strrep("=", 70), "\n\n")

result <- plot_manhattan_with_genes(
  trait_name = "betaine",
  gtf_file = "genome.gtf.gz",
  p_threshold = 1e-5,
  label_only_top = 10       # Only label top 10 by significance
)

# =============================================================================
# Example 4: Strict QTL Definition (Smaller Merge Distance)
# =============================================================================

cat("\n", strrep("=", 70), "\n")
cat("Example 4: Strict QTL Definition (500kb merge distance)\n")
cat(strrep("=", 70), "\n\n")

result <- plot_manhattan_with_genes(
  trait_name = "betaine",
  gtf_file = "genome.gtf.gz",
  p_threshold = 1e-6,
  merge_distance = 5e5,     # 500kb - stricter QTL definition
  max_labels = 15
)

# =============================================================================
# Example 5: Output to PDF (Vector Graphics)
# =============================================================================

cat("\n", strrep("=", 70), "\n")
cat("Example 5: Output to PDF\n")
cat(strrep("=", 70), "\n\n")

result <- plot_manhattan_with_genes(
  trait_name = "betaine",
  gtf_file = "genome.gtf.gz",
  output_format = "pdf"      # PDF format for publications
)

# =============================================================================
# Example 6: Without Binning (For Smaller Datasets)
# =============================================================================

cat("\n", strrep("=", 70), "\n")
cat("Example 6: No Binning (for datasets < 1M SNPs)\n")
cat(strrep("=", 70), "\n\n")

result <- plot_manhattan_with_genes(
  trait_name = "betaine",
  gtf_file = "genome.gtf.gz",
  use_binning = FALSE        # Plot all SNPs (slower but exact)
)

# =============================================================================
# Example 7: Access and Export QTL Summary
# =============================================================================

cat("\n", strrep("=", 70), "\n")
cat("Example 7: Working with QTL Summary Table\n")
cat(strrep("=", 70), "\n\n")

result <- plot_manhattan_with_genes(
  trait_name = "betaine",
  gtf_file = "genome.gtf.gz"
)

# Access QTL regions table
qtl_table <- result$qtl_regions

cat("\nQTL Regions Summary:\n")
print(qtl_table)

# Save to Excel-friendly format
library(data.table)
fwrite(qtl_table, "betaine_QTL_summary.csv")

# =============================================================================
# Example 8: Batch Processing Multiple Traits with Gene Annotation
# =============================================================================

cat("\n", strrep("=", 70), "\n")
cat("Example 8: Batch Processing with Gene Annotation\n")
cat(strrep("=", 70), "\n\n")

traits <- c("betaine", "glucose", "cholesterol")
gtf_file <- "Homo_sapiens.GRCh38.gtf.gz"

for(trait in traits) {
  cat("\n### Processing:", trait, "###\n")

  tryCatch({
    result <- plot_manhattan_with_genes(
      trait_name = trait,
      gtf_file = gtf_file,
      p_threshold = 1e-5,
      label_only_top = 15
    )

    cat("SUCCESS:", trait, "\n")

  }, error = function(e) {
    cat("ERROR processing", trait, ":", conditionMessage(e), "\n")
  })
}

# =============================================================================
# Example 9: Parallel Batch Processing
# =============================================================================

if(FALSE) {  # Uncomment to run

  cat("\n", strrep("=", 70), "\n")
  cat("Example 9: Parallel Processing\n")
  cat(strrep("=", 70), "\n\n")

  library(parallel)

  traits <- c("betaine", "glucose", "cholesterol", "BMI", "height")
  gtf_file <- "genome.gtf.gz"

  # Process in parallel using 4 cores
  results <- mclapply(traits, function(trait) {
    plot_manhattan_with_genes(
      trait_name = trait,
      gtf_file = gtf_file,
      p_threshold = 1e-5,
      label_only_top = 10
    )
  }, mc.cores = 4)

  names(results) <- traits

  # Access results for each trait
  for(trait in names(results)) {
    cat("\nQTL regions for", trait, ":\n")
    print(results[[trait]]$qtl_regions)
  }
}

# =============================================================================
# Example 10: Custom Plot Styling
# =============================================================================

if(FALSE) {  # Uncomment to run

  cat("\n", strrep("=", 70), "\n")
  cat("Example 10: Custom Styling\n")
  cat(strrep("=", 70), "\n\n")

  result <- plot_manhattan_with_genes(
    trait_name = "betaine",
    gtf_file = "genome.gtf.gz",
    p_threshold = 1e-5
  )

  # Customize the plot further
  library(ggplot2)

  custom_plot <- result$plot +
    theme(
      plot.title = element_text(size = 16, face = "bold", color = "darkblue"),
      axis.title = element_text(size = 14, face = "bold"),
      panel.background = element_rect(fill = "white"),
      panel.grid.major.y = element_line(color = "grey90")
    ) +
    ggtitle("Custom Title: Betaine GWAS Results")

  # Save custom plot
  ggsave("betaine_custom.png", custom_plot, width = 14, height = 8, dpi = 300)
}

# =============================================================================
# Example 11: Download and Use Public GTF Files
# =============================================================================

if(FALSE) {  # Example code - adapt to your needs

  # For Human (GRCh38):
  # Download from Ensembl:
  # wget http://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz

  # For Mouse (GRCm39):
  # wget http://ftp.ensembl.org/pub/release-109/gtf/mus_musculus/Mus_musculus.GRCm39.109.gtf.gz

  # For Arabidopsis (TAIR10):
  # wget ftp://ftp.ensemblgenomes.org/pub/plants/release-52/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.52.gtf.gz

  # Then use in your analysis:
  result <- plot_manhattan_with_genes(
    trait_name = "my_trait",
    gtf_file = "Homo_sapiens.GRCh38.109.gtf.gz"
  )
}

# =============================================================================
# Complete
# =============================================================================

cat("\n\n", strrep("=", 70), "\n")
cat("All examples completed!\n")
cat("Check your working directory for output files:\n")
cat("  - *.GWAS.annotated.tiff (Manhattan plots with gene labels)\n")
cat("  - *.QTL_regions.txt (QTL summary tables)\n")
cat(strrep("=", 70), "\n")
