#!/usr/bin/env Rscript
# =============================================================================
# GWAS Manhattan Plot with Gene Annotation
# Automatically labels nearest genes to lead SNPs in QTL regions
# =============================================================================

library(data.table)
library(ggplot2)
library(ggrepel)
library(Rcpp)

# Load C++ functions (make sure gwas_rcpp_utils.cpp is compiled)
if(!exists("calc_cumulative_pos_cpp")) {
  cat("Loading C++ acceleration module...\n")
  sourceCpp("gwas_rcpp_utils.cpp")
}

# =============================================================================
# Helper Functions
# =============================================================================

#' Parse GTF file and extract gene information
#'
#' @param gtf_file Path to GTF file (can be .gtf or .gtf.gz)
#' @param feature_type Feature type to extract (default: "gene")
#' @return data.table with gene information
#' @export
parse_gtf <- function(gtf_file, feature_type = "gene") {

  cat("Reading GTF file:", gtf_file, "\n")

  # Read GTF file
  if(grepl("\\.gz$", gtf_file)) {
    gtf_lines <- readLines(gzfile(gtf_file))
  } else {
    gtf_lines <- readLines(gtf_file)
  }

  # Filter comment lines
  gtf_lines <- gtf_lines[!grepl("^#", gtf_lines)]

  cat("Parsing", format(length(gtf_lines), big.mark = ","), "lines...\n")

  # Parse GTF format
  genes <- data.table()

  for(line in gtf_lines) {
    fields <- strsplit(line, "\t")[[1]]

    if(length(fields) < 9) next
    if(fields[3] != feature_type) next

    # Extract chromosome (remove "chr" prefix if present)
    chr <- fields[1]
    chr <- gsub("^chr", "", chr)
    chr <- gsub("^Chr", "", chr)

    # Skip if not numeric chromosome
    if(!grepl("^[0-9]+$", chr)) next

    chr_num <- as.integer(chr)
    start <- as.numeric(fields[4])
    end <- as.numeric(fields[5])

    # Extract gene name from attributes
    attr_str <- fields[9]
    gene_name <- NA

    # Try to extract gene_name or gene_id
    if(grepl('gene_name "([^"]+)"', attr_str)) {
      gene_name <- sub('.*gene_name "([^"]+)".*', '\\1', attr_str)
    } else if(grepl('gene_id "([^"]+)"', attr_str)) {
      gene_name <- sub('.*gene_id "([^"]+)".*', '\\1', attr_str)
    }

    if(!is.na(gene_name)) {
      genes <- rbind(genes, data.table(
        Chr = chr_num,
        start = start,
        end = end,
        gene_name = gene_name
      ))
    }
  }

  cat("Extracted", format(nrow(genes), big.mark = ","), "genes\n")

  return(genes)
}

#' Faster GTF parsing using fread
#'
#' @param gtf_file Path to GTF file
#' @param feature_type Feature type to extract (default: "gene")
#' @return data.table with gene information
#' @export
parse_gtf_fast <- function(gtf_file, feature_type = "gene") {

  cat("Reading GTF file (fast mode):", gtf_file, "\n")

  # Read GTF with fread (much faster)
  gtf <- fread(gtf_file, sep = "\t", header = FALSE, skip = "#",
               select = c(1, 3, 4, 5, 9))

  setnames(gtf, c("chr", "feature", "start", "end", "attributes"))

  # Filter for genes only
  gtf <- gtf[feature == feature_type]

  cat("Found", format(nrow(gtf), big.mark = ","), "genes\n")

  # Clean chromosome names
  gtf[, chr := gsub("^chr", "", chr)]
  gtf[, chr := gsub("^Chr", "", chr)]

  # Keep only numeric chromosomes
  gtf <- gtf[grepl("^[0-9]+$", chr)]
  gtf[, Chr := as.integer(chr)]

  # Extract gene names
  gtf[, gene_name := NA_character_]

  # Extract gene_name
  gtf[grepl('gene_name "', attributes),
      gene_name := sub('.*gene_name "([^"]+)".*', '\\1', attributes)]

  # Fallback to gene_id
  gtf[is.na(gene_name) & grepl('gene_id "', attributes),
      gene_name := sub('.*gene_id "([^"]+)".*', '\\1', attributes)]

  # Clean up and return
  result <- gtf[!is.na(gene_name), .(Chr, start, end, gene_name)]

  cat("Extracted", format(nrow(result), big.mark = ","), "genes\n")

  return(result)
}

# =============================================================================
# Main Plotting Function with Gene Annotation
# =============================================================================

#' Plot GWAS Manhattan plot with gene annotations
#'
#' @param trait_name Trait name (file prefix for .mlma.gz)
#' @param gtf_file Path to GTF annotation file
#' @param p_threshold Significance threshold for QTL identification (default: 1e-5)
#' @param merge_distance Distance to merge nearby QTLs in bp (default: 1Mb)
#' @param max_labels Maximum number of QTL regions to label (default: 20)
#' @param label_only_top Label only top N regions by significance (default: NULL = all)
#' @param use_binning Use binning for speed (default: TRUE)
#' @param n_bins Number of bins (default: 5000)
#' @param use_scattermore Use scattermore for plotting (default: TRUE)
#' @param output_format Output format ("tiff", "png", "pdf")
#' @export
plot_manhattan_with_genes <- function(
    trait_name,
    gtf_file,
    p_threshold = 1e-5,
    merge_distance = 1e6,
    max_labels = 20,
    label_only_top = NULL,
    use_binning = TRUE,
    n_bins = 5000,
    use_scattermore = TRUE,
    output_format = "tiff"
) {

  cat("\n", strrep("=", 70), "\n")
  cat("GWAS Manhattan Plot with Gene Annotation\n")
  cat(strrep("=", 70), "\n\n")

  # 1. Read GWAS data
  cat("1. Reading GWAS data...\n")
  t1 <- Sys.time()
  gwas <- fread(paste0(trait_name, ".mlma.gz"), header = TRUE)
  cat(sprintf("   Loaded %s SNPs in %.2f sec\n",
              format(nrow(gwas), big.mark = ","),
              difftime(Sys.time(), t1, units = "secs")))

  # Store original data before any filtering
  gwas_original <- copy(gwas)

  # 2. Process GWAS data
  cat("\n2. Processing GWAS data with C++...\n")
  t2 <- Sys.time()

  gwas$BPcum <- calc_cumulative_pos_cpp(gwas$Chr, gwas$bp)
  gwas$log10p <- fast_log10p_cpp(gwas$p)
  gwas$Color <- ifelse(gwas$Chr %% 2 == 0, "turquoise", "darkcyan")
  gwas$original_idx <- 1:nrow(gwas)  # Track original indices

  cat(sprintf("   Processing completed in %.2f sec\n",
              difftime(Sys.time(), t2, units = "secs")))

  # 3. Identify QTL regions (before binning!)
  cat("\n3. Identifying QTL regions...\n")
  t3 <- Sys.time()

  qtl_regions <- identify_qtl_regions_cpp(
    gwas_original$Chr,
    gwas_original$bp,
    gwas_original$p,
    p_threshold = p_threshold,
    merge_distance = merge_distance
  )

  cat(sprintf("   Identified %d QTL regions in %.2f sec\n",
              nrow(qtl_regions),
              difftime(Sys.time(), t3, units = "secs")))

  # 4. Read GTF annotations
  cat("\n4. Reading GTF annotations...\n")
  t4 <- Sys.time()

  genes <- parse_gtf_fast(gtf_file, feature_type = "gene")

  cat(sprintf("   GTF parsing completed in %.2f sec\n",
              difftime(Sys.time(), t4, units = "secs")))

  # 5. Find nearest genes for lead SNPs
  cat("\n5. Finding nearest genes for lead SNPs...\n")
  t5 <- Sys.time()

  if(nrow(qtl_regions) > 0) {
    # Get lead SNP positions
    lead_snps <- gwas_original[qtl_regions$lead_snp_idx + 1]  # R is 1-indexed

    # Find nearest genes using C++
    nearest <- find_nearest_genes_cpp(
      lead_snps$Chr,
      lead_snps$bp,
      genes$Chr,
      genes$start,
      genes$end,
      genes$gene_name
    )

    # Add gene annotations to QTL regions
    qtl_regions$gene_name <- nearest$nearest_gene
    qtl_regions$gene_distance <- nearest$distance

    # Add cumulative positions for labeling
    qtl_regions$BPcum <- lead_snps$BPcum
    qtl_regions$log10p <- lead_snps$log10p

    # Filter out regions without gene annotations
    qtl_regions <- qtl_regions[!is.na(gene_name)]

    # Limit number of labels
    if(!is.null(label_only_top) && nrow(qtl_regions) > label_only_top) {
      qtl_regions <- qtl_regions[order(lead_snp_p)][1:label_only_top]
    } else if(nrow(qtl_regions) > max_labels) {
      qtl_regions <- qtl_regions[order(lead_snp_p)][1:max_labels]
    }

    cat(sprintf("   Found nearest genes in %.2f sec\n",
                difftime(Sys.time(), t5, units = "secs")))
    cat(sprintf("   Labeling %d QTL regions\n", nrow(qtl_regions)))
  } else {
    cat("   No QTL regions found!\n")
  }

  # 6. Bin data if requested
  original_n <- nrow(gwas)

  if(use_binning) {
    cat("\n6. Binning SNPs...\n")
    t6 <- Sys.time()
    keep <- bin_by_position_cpp(gwas$BPcum, gwas$log10p, n_bins = n_bins)
    gwas <- gwas[keep, ]
    cat(sprintf("   Binning completed in %.2f sec\n",
                difftime(Sys.time(), t6, units = "secs")))
    cat(sprintf("   Final dataset: %s SNPs (%.1f%% of original)\n",
                format(nrow(gwas), big.mark = ","),
                100 * nrow(gwas) / original_n))
  }

  # 7. Calculate plot elements
  cat("\n7. Calculating plot elements...\n")

  X_axis <- calc_chr_centers_cpp(gwas$Chr, gwas$BPcum)
  rect_data <- calc_backgrounds_cpp(gwas$Chr, gwas$BPcum)
  y_max <- max(gwas$log10p, na.rm = TRUE) * 1.05

  # 8. Create plot
  cat("\n8. Creating plot with gene annotations...\n")
  t8 <- Sys.time()

  p <- ggplot(gwas, aes(x = BPcum, y = log10p))

  # Add background rectangles
  if(nrow(rect_data) > 0) {
    p <- p + geom_rect(
      data = rect_data,
      aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = y_max),
      fill = 'grey80',
      alpha = 0.2,
      inherit.aes = FALSE
    )
  }

  # Add points
  if(use_scattermore && requireNamespace("scattermore", quietly = TRUE)) {
    p <- p + scattermore::geom_scattermore(
      aes(color = Color),
      pointsize = 1.5,
      alpha = 0.6,
      pixels = c(3600, 1200)
    )
  } else {
    p <- p + geom_point(
      aes(color = Color),
      alpha = 0.6,
      size = 0.8,
      show.legend = FALSE
    )
  }

  # Add gene labels if we have QTL regions
  if(exists("qtl_regions") && nrow(qtl_regions) > 0) {
    cat("   Adding", nrow(qtl_regions), "gene labels...\n")

    p <- p + geom_text_repel(
      data = qtl_regions,
      aes(x = BPcum, y = log10p, label = gene_name),
      size = 3.5,
      fontface = "italic",
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "grey50",
      segment.size = 0.3,
      max.overlaps = Inf,
      min.segment.length = 0
    )
  }

  # Finalize plot
  p <- p +
    scale_color_identity() +
    scale_x_continuous(labels = X_axis$Chr, breaks = X_axis$center) +
    geom_hline(yintercept = -log10(1e-5), linetype = 'dashed', color = "blue") +
    geom_hline(yintercept = -log10(5e-8), color = "red") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
      axis.title.x = element_text(face = "bold", size = 14),
      axis.title.y = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 11),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.border = element_rect(colour = "grey", fill = NA, size = 1)
    ) +
    xlab("Chromosome") +
    ylab("-log10(P-value)") +
    ggtitle(paste0("GWAS Manhattan Plot - ", trait_name))

  cat(sprintf("   Plot created in %.2f sec\n",
              difftime(Sys.time(), t8, units = "secs")))

  # 9. Save plot
  cat("\n9. Saving plot...\n")
  t9 <- Sys.time()

  output_file <- paste0(trait_name, ".GWAS.annotated.", output_format)

  if(output_format == "tiff") {
    tiff(output_file, width = 14, height = 8, units = "in", res = 300, compression = "lzw")
  } else if(output_format == "png") {
    png(output_file, width = 14, height = 8, units = "in", res = 300)
  } else if(output_format == "pdf") {
    pdf(output_file, width = 14, height = 8)
  }

  print(p)
  dev.off()

  cat(sprintf("   Plot saved to '%s' in %.2f sec\n",
              output_file,
              difftime(Sys.time(), t9, units = "secs")))

  # 10. Save QTL summary table
  if(exists("qtl_regions") && nrow(qtl_regions) > 0) {
    qtl_summary_file <- paste0(trait_name, ".QTL_regions.txt")

    qtl_summary <- qtl_regions[, .(
      Chr,
      start,
      end,
      width_kb = (end - start) / 1000,
      lead_snp_p,
      gene_name,
      gene_distance_kb = gene_distance / 1000
    )]

    fwrite(qtl_summary, qtl_summary_file, sep = "\t")

    cat(sprintf("\n   QTL summary saved to '%s'\n", qtl_summary_file))

    # Print summary
    cat("\n", strrep("-", 70), "\n")
    cat("QTL Summary:\n")
    cat(strrep("-", 70), "\n")
    print(qtl_summary)
  }

  # Summary
  total_time <- difftime(Sys.time(), t1, units = "secs")
  cat("\n", strrep("=", 70), "\n")
  cat(sprintf("TOTAL TIME: %.2f seconds\n", total_time))
  cat(strrep("=", 70), "\n\n")

  return(list(plot = p, qtl_regions = qtl_regions))
}

# =============================================================================
# Usage Examples (for reference)
# =============================================================================

if(FALSE) {

  # Example 1: Basic usage with gene annotation
  result <- plot_manhattan_with_genes(
    trait_name = "betaine",
    gtf_file = "Homo_sapiens.GRCh38.gtf.gz",
    p_threshold = 1e-5,
    merge_distance = 1e6
  )

  # Example 2: Label only top 10 most significant QTLs
  result <- plot_manhattan_with_genes(
    trait_name = "betaine",
    gtf_file = "genome.gtf.gz",
    label_only_top = 10,
    p_threshold = 5e-8  # Genome-wide significance
  )

  # Example 3: More stringent QTL merging
  result <- plot_manhattan_with_genes(
    trait_name = "betaine",
    gtf_file = "genome.gtf.gz",
    p_threshold = 1e-6,
    merge_distance = 5e5,  # 500kb
    max_labels = 15
  )

  # Access results
  print(result$qtl_regions)  # QTL summary table
  print(result$plot)          # ggplot object
}
