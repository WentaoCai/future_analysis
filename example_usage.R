#!/usr/bin/env Rscript
# =============================================================================
# GWAS Manhattan Plot - Ultra-Fast Version Usage Examples
# =============================================================================

# 设置工作目录（包含你的 .mlma.gz 文件）
setwd("~/Documents/Software/manhattan2D/")

# 加载ultra-fast脚本（会自动编译C++代码）
cat("Loading ultra-fast plotting functions...\n")
source("gwas_manhattan_ultrafast.R")

# =============================================================================
# 示例 1: 最简单的用法 - Ultra-Fast版本（推荐）
# =============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("Example 1: Ultra-Fast Version (Recommended)\n")
cat(strrep("=", 70), "\n\n")

# 只需一行代码！
plot_manhattan_ultrafast("betaine")

# 输出: betaine.GWAS.ultrafast.tiff
# 对于1000万SNPs，预期时间：5-10秒

# =============================================================================
# 示例 2: 自定义参数
# =============================================================================

cat("\n\n")
cat(strrep("=", 70), "\n")
cat("Example 2: Custom Parameters\n")
cat(strrep("=", 70), "\n\n")

plot_manhattan_ultrafast(
  trait_name = "betaine",
  use_binning = TRUE,              # 使用分箱（推荐）
  n_bins = 5000,                   # 分箱数量
  use_scattermore = TRUE,          # 使用快速绘图
  point_size = 2.0,                # 点大小
  output_format = "png"            # PNG格式（比TIFF快）
)

# =============================================================================
# 示例 3: 智能采样（保留所有显著SNP）
# =============================================================================

cat("\n\n")
cat(strrep("=", 70), "\n")
cat("Example 3: Smart Sampling (Keep All Significant SNPs)\n")
cat(strrep("=", 70), "\n\n")

plot_manhattan_ultrafast(
  trait_name = "betaine",
  use_binning = FALSE,              # 不使用分箱
  use_smart_sampling = TRUE,        # 使用智能采样
  sample_fraction = 0.05,           # 保留5%非显著SNP
  sig_threshold = 1e-4              # p < 0.0001 为显著
)

# =============================================================================
# 示例 4: Fast版本（C++加速但不分箱）
# =============================================================================

cat("\n\n")
cat(strrep("=", 70), "\n")
cat("Example 4: Fast Version (No Binning)\n")
cat(strrep("=", 70), "\n\n")

# 适合中等规模数据（<100万SNPs）
plot_manhattan_fast("betaine", output_format = "pdf")

# =============================================================================
# 示例 5: 批量处理多个性状
# =============================================================================

cat("\n\n")
cat(strrep("=", 70), "\n")
cat("Example 5: Batch Processing\n")
cat(strrep("=", 70), "\n\n")

# 定义性状列表
traits <- c("betaine", "glucose", "cholesterol")

# 串行处理
for (trait in traits) {
  cat("\n### Processing:", trait, "###\n")
  tryCatch({
    plot_manhattan_ultrafast(trait)
  }, error = function(e) {
    cat("Error processing", trait, ":", conditionMessage(e), "\n")
  })
}

# =============================================================================
# 示例 6: 并行处理（最快）
# =============================================================================

if (FALSE) {  # 取消注释以运行

  cat("\n\n")
  cat(strrep("=", 70), "\n")
  cat("Example 6: Parallel Processing\n")
  cat(strrep("=", 70), "\n\n")

  library(parallel)

  traits <- c("betaine", "glucose", "cholesterol", "BMI", "height")

  # 使用4个CPU核心并行处理
  mclapply(traits, function(trait) {
    plot_manhattan_ultrafast(trait)
  }, mc.cores = 4)

  cat("All traits processed!\n")
}

# =============================================================================
# 示例 7: 性能对比测试
# =============================================================================

if (FALSE) {  # 取消注释以运行

  cat("\n\n")
  cat(strrep("=", 70), "\n")
  cat("Example 7: Performance Comparison\n")
  cat(strrep("=", 70), "\n\n")

  library(microbenchmark)

  # 对比不同方法的速度
  benchmark <- microbenchmark(
    ultrafast_binning = {
      plot_manhattan_ultrafast("betaine", use_binning = TRUE)
    },
    ultrafast_sampling = {
      plot_manhattan_ultrafast("betaine",
                              use_binning = FALSE,
                              use_smart_sampling = TRUE)
    },
    fast_no_binning = {
      plot_manhattan_fast("betaine")
    },
    times = 3
  )

  print(benchmark)
  boxplot(benchmark)
}

# =============================================================================
# 示例 8: 生成模拟数据进行测试
# =============================================================================

if (FALSE) {  # 取消注释以运行

  cat("\n\n")
  cat(strrep("=", 70), "\n")
  cat("Example 8: Generate Simulated Data\n")
  cat(strrep("=", 70), "\n\n")

  library(data.table)

  # 生成1000万个SNPs的模拟数据
  n_snps <- 10000000
  n_chr <- 22

  cat("Generating", format(n_snps, big.mark = ","), "simulated SNPs...\n")

  simulated_gwas <- data.table(
    Chr = rep(1:n_chr, each = n_snps / n_chr),
    SNP = paste0("rs", 1:n_snps),
    bp = rep(seq(1, 1e8, length.out = n_snps / n_chr), n_chr),
    A1 = "A",
    A2 = "G",
    Freq = runif(n_snps, 0.1, 0.9),
    b = rnorm(n_snps, 0, 0.01),
    se = abs(rnorm(n_snps, 0.005, 0.001)),
    p = 10^(-runif(n_snps, 0, 8))  # p值从1到10^-8
  )

  # 添加一些显著的峰
  sig_indices <- sample(1:n_snps, 50)
  simulated_gwas$p[sig_indices] <- 10^(-runif(50, 8, 15))

  # 保存
  output_file <- "simulated_test.mlma.gz"
  cat("Saving to", output_file, "...\n")
  fwrite(simulated_gwas, output_file, compress = "gzip")

  cat("Simulated data saved!\n")

  # 绘制
  cat("\nPlotting simulated data...\n")
  plot_manhattan_ultrafast("simulated_test")

  cat("\nDone! Check simulated_test.GWAS.ultrafast.tiff\n")
}

# =============================================================================
# 完成
# =============================================================================

cat("\n\n")
cat(strrep("=", 70), "\n")
cat("All examples completed!\n")
cat("Check the output files in your working directory.\n")
cat(strrep("=", 70), "\n")
