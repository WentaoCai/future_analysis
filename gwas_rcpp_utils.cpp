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
