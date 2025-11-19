#!/usr/bin/env Rscript
# =============================================================================
# ULTRA-FAST GWAS Manhattan Plot with Rcpp + scattermore
# Expected speedup: 10-100x for large datasets (10M+ SNPs)
# =============================================================================

library(data.table)
library(ggplot2)
library(Rcpp)

# Compile and load C++ functions
cat("Loading C++ acceleration module...\n")
sourceCpp("gwas_rcpp_utils.cpp")

# Theme
theme_manhattan <- theme(
  plot.title = element_text(size = rel(1.3), vjust = 2, hjust = 0.5),
  legend.title = element_text(size = 15, face = "bold"),
  legend.text = element_text(size = 12),
  axis.title.x = element_text(face = "bold", vjust = -1, size = 16),
  axis.title.y = element_text(face = "bold", vjust = 2, size = 16, angle = 90),
  axis.text = element_text(size = rel(1.1)),
  axis.text.x = element_text(hjust = 0.5, vjust = 0, size = 12),
  axis.text.y = element_text(vjust = 0.5, hjust = 0, size = 14),
  axis.line = element_line(colour = "black"),
  axis.ticks = element_line(colour = 'black'),
  panel.background = element_blank(),
  aspect.ratio = 0.4,
  panel.border = element_rect(colour = "grey", fill = NA, size = 1),
  panel.grid.major = element_line(colour = NA),
  complete = TRUE
)

# =============================================================================
# ULTRA-FAST VERSION - Use this for 1M+ SNPs
# =============================================================================

plot_manhattan_ultrafast <- function(
    trait_name,
    use_binning = TRUE,          # Reduce to ~5000 points (HUGE speedup)
    n_bins = 5000,                # Number of bins for downsampling
    use_smart_sampling = FALSE,   # Alternative: keep sig + sample non-sig
    sample_fraction = 0.1,        # Fraction of non-sig SNPs to keep
    sig_threshold = 0.001,        # Significance threshold for sampling
    use_scattermore = TRUE,       # Use scattermore for ultra-fast plotting
    point_size = 1.5,             # Point size
    output_format = "tiff"        # "tiff", "png", or "pdf"
) {

  cat("\n", strrep("=", 70), "\n")
  cat("ULTRA-FAST GWAS Manhattan Plot\n")
  cat(strrep("=", 70), "\n\n")

  # 1. Read data with data.table (fastest reader)
  cat("1. Reading data...\n")
  t1 <- Sys.time()
  gwas <- fread(paste0(trait_name, ".mlma.gz"), header = TRUE)
  cat(sprintf("   Loaded %s SNPs in %.2f sec\n",
              format(nrow(gwas), big.mark = ","),
              difftime(Sys.time(), t1, units = "secs")))

  # 2. Data processing with C++
  cat("\n2. Processing data with C++...\n")
  t2 <- Sys.time()

  # Calculate cumulative positions (C++)
  gwas$BPcum <- calc_cumulative_pos_cpp(gwas$Chr, gwas$bp)

  # Calculate -log10(p) (C++)
  gwas$log10p <- fast_log10p_cpp(gwas$p)

  # Assign colors
  gwas$Color <- ifelse(gwas$Chr %% 2 == 0, "turquoise", "darkcyan")

  cat(sprintf("   Data processing completed in %.2f sec\n",
              difftime(Sys.time(), t2, units = "secs")))

  # 3. Intelligent downsampling
  original_n <- nrow(gwas)

  if (use_binning) {
    cat("\n3. Binning SNPs by position (C++)...\n")
    t3 <- Sys.time()
    keep <- bin_by_position_cpp(gwas$BPcum, gwas$log10p, n_bins = n_bins)
    gwas <- gwas[keep, ]
    cat(sprintf("   Binning completed in %.2f sec\n",
                difftime(Sys.time(), t3, units = "secs")))

  } else if (use_smart_sampling) {
    cat("\n3. Smart sampling (keep significant + sample non-significant)...\n")
    t3 <- Sys.time()
    keep <- smart_sample_cpp(gwas$p,
                            sig_threshold = sig_threshold,
                            sample_fraction = sample_fraction)
    gwas <- gwas[keep, ]
    cat(sprintf("   Sampling completed in %.2f sec\n",
                difftime(Sys.time(), t3, units = "secs")))
  } else {
    cat("\n3. No downsampling (plotting all SNPs)\n")
  }

  cat(sprintf("   Final dataset: %s SNPs (%.1f%% of original)\n",
              format(nrow(gwas), big.mark = ","),
              100 * nrow(gwas) / original_n))

  # 4. Calculate plot elements with C++
  cat("\n4. Calculating plot elements (C++)...\n")
  t4 <- Sys.time()

  X_axis <- calc_chr_centers_cpp(gwas$Chr, gwas$BPcum)
  rect_data <- calc_backgrounds_cpp(gwas$Chr, gwas$BPcum)
  y_max <- max(gwas$log10p, na.rm = TRUE) * 1.05

  cat(sprintf("   Calculations completed in %.2f sec\n",
              difftime(Sys.time(), t4, units = "secs")))

  # 5. Create plot
  cat("\n5. Creating plot...\n")
  t5 <- Sys.time()

  p <- ggplot(gwas, aes(x = BPcum, y = log10p))

  # Add background rectangles
  if (nrow(rect_data) > 0) {
    p <- p + geom_rect(
      data = rect_data,
      aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = y_max),
      fill = 'grey80',
      alpha = 0.2,
      inherit.aes = FALSE
    )
  }

  # Add points - use scattermore for HUGE speedup
  if (use_scattermore && requireNamespace("scattermore", quietly = TRUE)) {
    cat("   Using scattermore for ultra-fast plotting...\n")
    p <- p + scattermore::geom_scattermore(
      aes(color = Color),
      pointsize = point_size,
      alpha = 0.6,
      pixels = c(3600, 1200)  # High resolution
    )
  } else {
    if (use_scattermore) {
      cat("   scattermore not available, using standard geom_point\n")
      cat("   Install with: install.packages('scattermore')\n")
    }
    p <- p + geom_point(
      aes(color = Color),
      alpha = 0.6,
      size = 0.8,
      show.legend = FALSE
    )
  }

  # Finalize plot
  p <- p +
    scale_color_identity() +
    scale_x_continuous(labels = X_axis$Chr, breaks = X_axis$center) +
    geom_hline(yintercept = -log10(1e-5), linetype = 'dashed', color = "blue") +
    geom_hline(yintercept = -log10(5e-8), color = "red") +
    theme_manhattan +
    theme(
      panel.border = element_blank(),
      axis.line.y = element_line(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
    xlab("Chromosome") +
    ylab("-log10(P-value)")

  cat(sprintf("   Plot created in %.2f sec\n",
              difftime(Sys.time(), t5, units = "secs")))

  # 6. Save plot
  cat("\n6. Saving plot...\n")
  t6 <- Sys.time()

  output_file <- paste0(trait_name, ".GWAS.ultrafast.", output_format)

  if (output_format == "tiff") {
    tiff(output_file, width = 12, height = 8, units = "in", res = 300, compression = "lzw")
  } else if (output_format == "png") {
    png(output_file, width = 12, height = 8, units = "in", res = 300)
  } else if (output_format == "pdf") {
    pdf(output_file, width = 12, height = 8)
  }

  print(p)
  dev.off()

  cat(sprintf("   Plot saved to '%s' in %.2f sec\n",
              output_file,
              difftime(Sys.time(), t6, units = "secs")))

  # Summary
  total_time <- difftime(Sys.time(), t1, units = "secs")
  cat("\n", strrep("=", 70), "\n")
  cat(sprintf("TOTAL TIME: %.2f seconds\n", total_time))
  cat(strrep("=", 70), "\n\n")

  return(p)
}

# =============================================================================
# FAST VERSION - Simpler, still much faster than original
# =============================================================================

plot_manhattan_fast <- function(trait_name, output_format = "tiff") {

  cat("FAST GWAS Manhattan Plot (Rcpp accelerated)\n\n")

  # Read data
  cat("Reading data...\n")
  gwas <- fread(paste0(trait_name, ".mlma.gz"), header = TRUE)
  cat(sprintf("Loaded %s SNPs\n", format(nrow(gwas), big.mark = ",")))

  # Process with C++
  cat("Processing with C++...\n")
  gwas$BPcum <- calc_cumulative_pos_cpp(gwas$Chr, gwas$bp)
  gwas$log10p <- fast_log10p_cpp(gwas$p)
  gwas$Color <- ifelse(gwas$Chr %% 2 == 0, "turquoise", "darkcyan")

  # Calculate plot elements
  X_axis <- calc_chr_centers_cpp(gwas$Chr, gwas$BPcum)
  rect_data <- calc_backgrounds_cpp(gwas$Chr, gwas$BPcum)
  y_max <- max(gwas$log10p, na.rm = TRUE) * 1.05

  # Create plot
  cat("Creating plot...\n")
  p <- ggplot(gwas, aes(x = BPcum, y = log10p)) +
    geom_rect(
      data = rect_data,
      aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = y_max),
      fill = 'grey80',
      alpha = 0.2,
      inherit.aes = FALSE
    ) +
    geom_point(aes(color = Color), alpha = 0.6, size = 0.8, show.legend = FALSE) +
    scale_color_identity() +
    scale_x_continuous(labels = X_axis$Chr, breaks = X_axis$center) +
    geom_hline(yintercept = -log10(1e-5), linetype = 'dashed', color = "blue") +
    geom_hline(yintercept = -log10(5e-8), color = "red") +
    theme_manhattan +
    theme(
      panel.border = element_blank(),
      axis.line.y = element_line(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
    xlab("Chromosome") +
    ylab("-log10(P-value)")

  # Save
  output_file <- paste0(trait_name, ".GWAS.fast.", output_format)

  if (output_format == "tiff") {
    tiff(output_file, width = 12, height = 8, units = "in", res = 300, compression = "lzw")
  } else if (output_format == "png") {
    png(output_file, width = 12, height = 8, units = "in", res = 300)
  } else {
    pdf(output_file, width = 12, height = 8)
  }

  print(p)
  dev.off()

  cat(sprintf("Saved to: %s\n", output_file))

  return(p)
}

# =============================================================================
# USAGE EXAMPLES
# =============================================================================

if (FALSE) {  # Don't run automatically, just examples

  # Example 1: Ultra-fast with binning (RECOMMENDED for 1M+ SNPs)
  plot_manhattan_ultrafast("betaine",
                          use_binning = TRUE,
                          n_bins = 5000,
                          use_scattermore = TRUE)

  # Example 2: Ultra-fast with smart sampling
  plot_manhattan_ultrafast("betaine",
                          use_binning = FALSE,
                          use_smart_sampling = TRUE,
                          sample_fraction = 0.1,
                          sig_threshold = 0.001)

  # Example 3: Fast version (all SNPs, but C++ accelerated)
  plot_manhattan_fast("betaine")

  # Example 4: Batch processing
  traits <- c("betaine", "glucose", "cholesterol")
  for (trait in traits) {
    plot_manhattan_ultrafast(trait)
  }

  # Example 5: Parallel processing
  library(parallel)
  traits <- c("betaine", "glucose", "cholesterol")
  mclapply(traits, plot_manhattan_ultrafast, mc.cores = 4)
}
