# =============================================================================
# GWAS 曼哈顿图性能对比脚本
# =============================================================================

library(ggplot2)
library(data.table)
library(microbenchmark)

# 设置工作目录
setwd("~/Documents/Software/manhattan2D/")

# 主题设置（两个版本共用）
theme_ydz0.3 <- theme(
  plot.title = element_text(size = rel(1.3), vjust = 2, hjust = 0.5, lineheight = 0.8),
  legend.title = element_text(size = 15, face = "bold"),
  legend.text = element_text(size = 12),
  axis.title.x = element_text(face = "bold", vjust = -1, size = 16),
  axis.title.y = element_text(face = "bold", vjust = 2, size = 16, angle = 90),
  axis.text = element_text(size = rel(1.1)),
  axis.text.x = element_text(hjust = 0.5, vjust = 0, size = 12),
  axis.text.y = element_text(vjust = 0.5, hjust = 0, size = 14),
  axis.line = element_line(colour = "black"),
  axis.ticks = element_line(colour = 'black'),
  strip.text = element_text(size = rel(1.3)),
  panel.background = element_blank(),
  aspect.ratio = 0.4,
  panel.border = element_rect(colour = "grey", fill = NA, size = 1),
  panel.grid.major = element_line(colour = NA),
  complete = TRUE
)

# =============================================================================
# 原始版本（使用 tidyverse）
# =============================================================================

plot_manhattan_original <- function(trait_name) {
  library(tidyverse)

  gwasResults <- fread(paste0(trait_name, ".mlma.gz"), head = TRUE)

  chr_len <- gwasResults %>%
    group_by(Chr) %>%
    summarise(chr_len = max(bp))

  chr_pos <- chr_len %>%
    mutate(total = cumsum(chr_len/10)*10 - chr_len) %>%
    select(-chr_len)

  Snp_pos <- chr_pos %>%
    left_join(gwasResults, ., by = "Chr") %>%
    arrange(Chr, bp) %>%
    mutate(BPcum = bp + total)

  Snp_pos$Color <- ifelse(Snp_pos$Chr %% 2 == 0, "turquoise", "darkcyan")

  X_axis <- Snp_pos %>%
    group_by(Chr) %>%
    summarize(center = (max(BPcum) + min(BPcum)) / 2)

  p <- ggplot(Snp_pos, aes(x = BPcum, y = -log10(p)))

  region <- as.data.frame(chr_pos)
  length <- as.data.frame(chr_len)
  anno <- merge(region, length, by.x = "Chr")
  anno$end1 <- anno$total + anno$chr_len

  # 性能瓶颈：循环调用 annotate()
  for (i in 1:(length(anno$Chr))) {
    if (i %% 2 == 1) {
      p <- p + annotate("rect",
                        xmin = anno$total[i],
                        xmax = anno$end1[i],
                        ymin = 0,
                        ymax = as.data.frame(layer_scales(p)$y$range$range)[2,],
                        alpha = .2,
                        fill = 'grey80')
    }
  }

  p <- p +
    geom_point(color = Snp_pos$Color, alpha = 0.6, size = 0.8) +
    scale_x_continuous(label = X_axis$Chr, breaks = X_axis$center) +
    theme_ydz0.3 +
    geom_hline(yintercept = -log10(0.00001), linetype = 'dashed', color = "blue") +
    geom_hline(yintercept = -log10(0.00000005), color = "red") +
    theme(
      panel.border = element_blank(),
      axis.line.y = element_line(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
    xlab("Position") +
    ylab("-log10(P-value)")

  return(p)
}

# =============================================================================
# 优化版本（使用 data.table + geom_rect）
# =============================================================================

plot_manhattan_optimized <- function(trait_name) {

  gwasResults <- fread(paste0(trait_name, ".mlma.gz"), header = TRUE)
  setDT(gwasResults)

  # 使用 data.table 进行数据处理
  chr_len <- gwasResults[, .(chr_len = max(bp)), by = Chr]

  chr_pos <- chr_len[order(Chr)]
  chr_pos[, total := cumsum(chr_len / 10) * 10 - chr_len]
  chr_pos[, chr_len := NULL]

  Snp_pos <- gwasResults[chr_pos, on = "Chr"]
  setorder(Snp_pos, Chr, bp)
  Snp_pos[, BPcum := bp + total]
  Snp_pos[, Color := ifelse(Chr %% 2 == 0, "turquoise", "darkcyan")]
  Snp_pos[, log10p := -log10(p)]  # 预计算

  X_axis <- Snp_pos[, .(center = (max(BPcum) + min(BPcum)) / 2), by = Chr]

  # 准备背景矩形数据
  anno <- merge(chr_pos, chr_len, by = "Chr")
  anno[, end1 := total + chr_len]
  rect_data <- anno[Chr %% 2 == 1]

  # 预计算 y 轴范围
  y_max <- max(Snp_pos$log10p) * 1.05

  # 使用 geom_rect 一次性绘制所有背景
  p <- ggplot(Snp_pos, aes(x = BPcum, y = log10p)) +
    geom_rect(
      data = rect_data,
      aes(xmin = total, xmax = end1, ymin = 0, ymax = y_max),
      fill = 'grey80',
      alpha = 0.2,
      inherit.aes = FALSE
    ) +
    geom_point(aes(color = Color), alpha = 0.6, size = 0.8, show.legend = FALSE) +
    scale_color_identity() +
    scale_x_continuous(labels = X_axis$Chr, breaks = X_axis$center) +
    geom_hline(yintercept = -log10(0.00001), linetype = 'dashed', color = "blue") +
    geom_hline(yintercept = -log10(0.00000005), color = "red") +
    theme_ydz0.3 +
    theme(
      panel.border = element_blank(),
      axis.line.y = element_line(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
    xlab("Position") +
    ylab("-log10(P-value)")

  return(p)
}

# =============================================================================
# 性能测试
# =============================================================================

cat("开始性能对比测试...\n\n")

trait_name <- "betaine"

cat("提示：如果数据文件不存在，可以生成模拟数据进行测试\n\n")

# 检查文件是否存在
if (!file.exists(paste0(trait_name, ".mlma.gz"))) {
  cat("数据文件不存在，生成模拟数据...\n")

  # 生成模拟 GWAS 数据
  set.seed(123)
  n_snps <- 1000000  # 100万个SNP
  n_chr <- 22

  simulated_data <- data.table(
    Chr = rep(1:n_chr, each = n_snps / n_chr),
    SNP = paste0("rs", 1:n_snps),
    bp = rep(seq(1, 100000000, length.out = n_snps / n_chr), n_chr),
    p = 10^(-runif(n_snps, 0, 8))  # p值在 1 到 10^-8 之间
  )

  fwrite(simulated_data, paste0(trait_name, ".mlma.gz"), compress = "gzip")
  cat("模拟数据已生成\n\n")
}

# 性能测试
cat("=" %+% strrep("=", 70) %+% "\n")
cat("运行性能基准测试（每个版本运行3次）\n")
cat("=" %+% strrep("=", 70) %+% "\n\n")

benchmark_result <- microbenchmark(
  original = {
    p1 <- plot_manhattan_original(trait_name)
  },
  optimized = {
    p2 <- plot_manhattan_optimized(trait_name)
  },
  times = 3
)

print(benchmark_result)

cat("\n\n")
cat("=" %+% strrep("=", 70) %+% "\n")
cat("性能提升分析\n")
cat("=" %+% strrep("=", 70) %+% "\n\n")

summary_result <- summary(benchmark_result)
speedup <- summary_result$median[1] / summary_result$median[2]

cat(sprintf("原始版本中位数时间: %.2f 秒\n", summary_result$median[1] / 1000))
cat(sprintf("优化版本中位数时间: %.2f 秒\n", summary_result$median[2] / 1000))
cat(sprintf("性能提升: %.2fx\n", speedup))

# 绘制对比图
library(ggplot2)

comparison_df <- data.frame(
  Version = rep(c("Original", "Optimized"), each = 3),
  Time = c(benchmark_result$time[benchmark_result$expr == "original"] / 1e9,
           benchmark_result$time[benchmark_result$expr == "optimized"] / 1e9)
)

p_comparison <- ggplot(comparison_df, aes(x = Version, y = Time, fill = Version)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.5) +
  labs(
    title = "GWAS Manhattan Plot Performance Comparison",
    x = "Version",
    y = "Time (seconds)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("performance_comparison.png", p_comparison, width = 8, height = 6, dpi = 300)

cat("\n性能对比图已保存至: performance_comparison.png\n")

# 保存结果
saveRDS(benchmark_result, "benchmark_results.rds")
cat("基准测试结果已保存至: benchmark_results.rds\n")
