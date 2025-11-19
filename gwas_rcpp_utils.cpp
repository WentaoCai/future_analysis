// GWAS Manhattan Plot - Rcpp Acceleration Module
// This provides 10-100x speedup for data processing

#include <Rcpp.h>
#include <algorithm>
#include <unordered_map>
#include <cmath>

using namespace Rcpp;

//' Fast cumulative position calculation
//' @param chr Chromosome numbers
//' @param bp Base pair positions
//' @export
// [[Rcpp::export]]
NumericVector calc_cumulative_pos_cpp(IntegerVector chr, NumericVector bp) {
  int n = chr.size();
  NumericVector BPcum(n);

  // Calculate max bp per chromosome
  std::unordered_map<int, double> chr_max;
  for(int i = 0; i < n; i++) {
    if(chr_max.find(chr[i]) == chr_max.end() || bp[i] > chr_max[chr[i]]) {
      chr_max[chr[i]] = bp[i];
    }
  }

  // Calculate offsets
  std::unordered_map<int, double> chr_offset;
  double offset = 0.0;

  // Sort chromosomes to ensure consistent order
  std::vector<int> unique_chrs;
  for(auto& p : chr_max) {
    unique_chrs.push_back(p.first);
  }
  std::sort(unique_chrs.begin(), unique_chrs.end());

  for(int c : unique_chrs) {
    chr_offset[c] = offset;
    double chr_len = chr_max[c];
    offset += std::ceil(chr_len / 10.0) * 10.0;
  }

  // Calculate cumulative positions
  for(int i = 0; i < n; i++) {
    BPcum[i] = bp[i] + chr_offset[chr[i]];
  }

  return BPcum;
}

//' Fast -log10(p) calculation
//' @param p P-values
//' @export
// [[Rcpp::export]]
NumericVector fast_log10p_cpp(NumericVector p) {
  int n = p.size();
  NumericVector result(n);

  for(int i = 0; i < n; i++) {
    if(NumericVector::is_na(p[i]) || p[i] <= 0) {
      result[i] = NA_REAL;
    } else {
      result[i] = -std::log10(p[i]);
    }
  }

  return result;
}

//' Intelligent SNP sampling - keep significant + sample non-significant
//' @param p P-values
//' @param sig_threshold Significance threshold (default 0.001)
//' @param sample_fraction Fraction of non-sig SNPs to keep (default 0.05)
//' @export
// [[Rcpp::export]]
LogicalVector smart_sample_cpp(NumericVector p,
                               double sig_threshold = 0.001,
                               double sample_fraction = 0.05) {
  int n = p.size();
  LogicalVector keep(n);

  std::vector<int> nonsig_indices;
  int n_sig = 0;

  // Identify significant and non-significant SNPs
  for(int i = 0; i < n; i++) {
    if(!NumericVector::is_na(p[i]) && p[i] < sig_threshold) {
      keep[i] = true;
      n_sig++;
    } else {
      keep[i] = false;
      nonsig_indices.push_back(i);
    }
  }

  // Randomly sample non-significant SNPs
  int n_to_sample = std::max(1, (int)(nonsig_indices.size() * sample_fraction));

  if(nonsig_indices.size() > 0) {
    // Shuffle indices
    for(int i = nonsig_indices.size() - 1; i > 0; i--) {
      int j = rand() % (i + 1);
      std::swap(nonsig_indices[i], nonsig_indices[j]);
    }

    // Keep first n_to_sample
    for(int i = 0; i < n_to_sample && i < nonsig_indices.size(); i++) {
      keep[nonsig_indices[i]] = true;
    }
  }

  Rcpp::Rcout << "Kept " << n_sig << " significant + "
              << n_to_sample << " sampled = "
              << (n_sig + n_to_sample) << " / " << n << " SNPs ("
              << (100.0 * (n_sig + n_to_sample) / n) << "%)\n";

  return keep;
}

//' Bin SNPs by position - keep most significant in each bin
//' This is CRITICAL for speed with 10M+ SNPs
//' @param BPcum Cumulative positions
//' @param log10p -log10(p) values
//' @param n_bins Number of bins (default 5000 for high-res plot)
//' @export
// [[Rcpp::export]]
LogicalVector bin_by_position_cpp(NumericVector BPcum,
                                  NumericVector log10p,
                                  int n_bins = 5000) {
  int n = BPcum.size();
  LogicalVector keep(n, false);

  if(n == 0) return keep;

  // Find min/max positions
  double min_pos = BPcum[0];
  double max_pos = BPcum[0];
  for(int i = 1; i < n; i++) {
    if(BPcum[i] < min_pos) min_pos = BPcum[i];
    if(BPcum[i] > max_pos) max_pos = BPcum[i];
  }

  double bin_width = (max_pos - min_pos) / n_bins;
  if(bin_width <= 0) {
    std::fill(keep.begin(), keep.end(), true);
    return keep;
  }

  // For each bin, track index of SNP with highest -log10(p)
  std::unordered_map<int, int> bin_best_idx;
  std::unordered_map<int, double> bin_best_value;

  for(int i = 0; i < n; i++) {
    int bin = (int)((BPcum[i] - min_pos) / bin_width);
    if(bin >= n_bins) bin = n_bins - 1;

    double val = NumericVector::is_na(log10p[i]) ? 0.0 : log10p[i];

    if(bin_best_idx.find(bin) == bin_best_idx.end() || val > bin_best_value[bin]) {
      bin_best_idx[bin] = i;
      bin_best_value[bin] = val;
    }
  }

  // Mark best SNPs from each bin
  for(auto& kv : bin_best_idx) {
    keep[kv.second] = true;
  }

  int n_kept = 0;
  for(int i = 0; i < n; i++) {
    if(keep[i]) n_kept++;
  }

  Rcpp::Rcout << "Binning: " << n << " SNPs -> " << n_kept
              << " SNPs (" << (100.0 * n_kept / n) << "%)\n";

  return keep;
}

//' Calculate chromosome centers for axis labels
//' @param chr Chromosome numbers
//' @param BPcum Cumulative positions
//' @export
// [[Rcpp::export]]
DataFrame calc_chr_centers_cpp(IntegerVector chr, NumericVector BPcum) {
  std::unordered_map<int, double> chr_min;
  std::unordered_map<int, double> chr_max;

  for(int i = 0; i < chr.size(); i++) {
    int c = chr[i];
    double pos = BPcum[i];

    if(chr_min.find(c) == chr_min.end()) {
      chr_min[c] = pos;
      chr_max[c] = pos;
    } else {
      if(pos < chr_min[c]) chr_min[c] = pos;
      if(pos > chr_max[c]) chr_max[c] = pos;
    }
  }

  std::vector<int> chrs;
  std::vector<double> centers;

  for(auto& p : chr_min) {
    chrs.push_back(p.first);
  }
  std::sort(chrs.begin(), chrs.end());

  for(int c : chrs) {
    centers.push_back((chr_min[c] + chr_max[c]) / 2.0);
  }

  return DataFrame::create(
    Named("Chr") = chrs,
    Named("center") = centers
  );
}

//' Calculate background rectangles
//' @param chr Chromosome numbers
//' @param BPcum Cumulative positions
//' @export
// [[Rcpp::export]]
DataFrame calc_backgrounds_cpp(IntegerVector chr, NumericVector BPcum) {
  std::unordered_map<int, double> chr_min;
  std::unordered_map<int, double> chr_max;

  for(int i = 0; i < chr.size(); i++) {
    int c = chr[i];
    double pos = BPcum[i];

    if(chr_min.find(c) == chr_min.end()) {
      chr_min[c] = pos;
      chr_max[c] = pos;
    } else {
      if(pos < chr_min[c]) chr_min[c] = pos;
      if(pos > chr_max[c]) chr_max[c] = pos;
    }
  }

  std::vector<int> chrs;
  std::vector<double> xmins, xmaxs;

  for(auto& p : chr_min) {
    if(p.first % 2 == 1) {  // Only odd chromosomes
      chrs.push_back(p.first);
      xmins.push_back(p.second);
      xmaxs.push_back(chr_max[p.first]);
    }
  }

  return DataFrame::create(
    Named("Chr") = chrs,
    Named("xmin") = xmins,
    Named("xmax") = xmaxs
  );
}

// =============================================================================
// Gene Annotation Functions
// =============================================================================

//' Find nearest gene for each SNP (C++ optimized)
//' @param snp_chr SNP chromosome numbers
//' @param snp_bp SNP base pair positions
//' @param gene_chr Gene chromosome numbers
//' @param gene_start Gene start positions
//' @param gene_end Gene end positions
//' @param gene_names Gene names
//' @export
// [[Rcpp::export]]
DataFrame find_nearest_genes_cpp(IntegerVector snp_chr,
                                 NumericVector snp_bp,
                                 IntegerVector gene_chr,
                                 NumericVector gene_start,
                                 NumericVector gene_end,
                                 CharacterVector gene_names) {
  int n_snps = snp_chr.size();
  int n_genes = gene_chr.size();

  CharacterVector nearest_genes(n_snps);
  NumericVector distances(n_snps);

  // Build chromosome-based index for genes
  std::unordered_map<int, std::vector<int>> chr_to_genes;
  for(int i = 0; i < n_genes; i++) {
    chr_to_genes[gene_chr[i]].push_back(i);
  }

  // For each SNP, find nearest gene on same chromosome
  for(int i = 0; i < n_snps; i++) {
    int chr = snp_chr[i];
    double pos = snp_bp[i];

    if(chr_to_genes.find(chr) == chr_to_genes.end()) {
      nearest_genes[i] = NA_STRING;
      distances[i] = NA_REAL;
      continue;
    }

    double min_dist = R_PosInf;
    int nearest_idx = -1;

    // Check all genes on this chromosome
    for(int gene_idx : chr_to_genes[chr]) {
      double dist;

      // Calculate distance to gene
      if(pos < gene_start[gene_idx]) {
        // SNP is upstream
        dist = gene_start[gene_idx] - pos;
      } else if(pos > gene_end[gene_idx]) {
        // SNP is downstream
        dist = pos - gene_end[gene_idx];
      } else {
        // SNP is within gene
        dist = 0.0;
      }

      if(dist < min_dist) {
        min_dist = dist;
        nearest_idx = gene_idx;
      }
    }

    if(nearest_idx >= 0) {
      nearest_genes[i] = gene_names[nearest_idx];
      distances[i] = min_dist;
    } else {
      nearest_genes[i] = NA_STRING;
      distances[i] = NA_REAL;
    }
  }

  return DataFrame::create(
    Named("nearest_gene") = nearest_genes,
    Named("distance") = distances
  );
}

//' Identify QTL regions (clusters of significant SNPs)
//' @param chr Chromosome numbers
//' @param bp Base pair positions
//' @param p P-values
//' @param p_threshold Significance threshold (default 1e-5)
//' @param merge_distance Max distance to merge nearby QTLs (default 1Mb)
//' @export
// [[Rcpp::export]]
DataFrame identify_qtl_regions_cpp(IntegerVector chr,
                                   NumericVector bp,
                                   NumericVector p,
                                   double p_threshold = 1e-5,
                                   double merge_distance = 1e6) {
  int n = chr.size();

  // Find significant SNPs
  std::vector<int> sig_indices;
  for(int i = 0; i < n; i++) {
    if(!NumericVector::is_na(p[i]) && p[i] < p_threshold) {
      sig_indices.push_back(i);
    }
  }

  if(sig_indices.empty()) {
    // No significant SNPs
    return DataFrame::create(
      Named("Chr") = IntegerVector(),
      Named("start") = NumericVector(),
      Named("end") = NumericVector(),
      Named("lead_snp_idx") = IntegerVector(),
      Named("lead_snp_p") = NumericVector()
    );
  }

  // Sort by chromosome and position
  std::sort(sig_indices.begin(), sig_indices.end(),
            [&](int a, int b) {
              if(chr[a] != chr[b]) return chr[a] < chr[b];
              return bp[a] < bp[b];
            });

  // Merge nearby significant SNPs into QTL regions
  std::vector<int> qtl_chr;
  std::vector<double> qtl_start, qtl_end;
  std::vector<int> lead_snp_idx;
  std::vector<double> lead_snp_p;

  int current_chr = chr[sig_indices[0]];
  double current_start = bp[sig_indices[0]];
  double current_end = bp[sig_indices[0]];
  int current_lead_idx = sig_indices[0];
  double current_min_p = p[sig_indices[0]];

  for(size_t i = 1; i < sig_indices.size(); i++) {
    int idx = sig_indices[i];

    // Check if this SNP should be merged with current region
    if(chr[idx] == current_chr && bp[idx] - current_end <= merge_distance) {
      // Extend current region
      current_end = bp[idx];

      // Update lead SNP if this one is more significant
      if(p[idx] < current_min_p) {
        current_lead_idx = idx;
        current_min_p = p[idx];
      }
    } else {
      // Save current region
      qtl_chr.push_back(current_chr);
      qtl_start.push_back(current_start);
      qtl_end.push_back(current_end);
      lead_snp_idx.push_back(current_lead_idx);
      lead_snp_p.push_back(current_min_p);

      // Start new region
      current_chr = chr[idx];
      current_start = bp[idx];
      current_end = bp[idx];
      current_lead_idx = idx;
      current_min_p = p[idx];
    }
  }

  // Save last region
  qtl_chr.push_back(current_chr);
  qtl_start.push_back(current_start);
  qtl_end.push_back(current_end);
  lead_snp_idx.push_back(current_lead_idx);
  lead_snp_p.push_back(current_min_p);

  Rcpp::Rcout << "Identified " << qtl_chr.size() << " QTL regions\n";

  return DataFrame::create(
    Named("Chr") = wrap(qtl_chr),
    Named("start") = wrap(qtl_start),
    Named("end") = wrap(qtl_end),
    Named("lead_snp_idx") = wrap(lead_snp_idx),
    Named("lead_snp_p") = wrap(lead_snp_p)
  );
}
