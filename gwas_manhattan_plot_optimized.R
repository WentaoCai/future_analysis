library(ggplot2)
library(data.table)
library(ggrepel)

# 设置工作目录
setwd("~/Documents/Software/manhattan2D/")

# 主题设置
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
# 优化后的曼哈顿图绘制函数
# =============================================================================

plot_manhattan_optimized <- function(trait_name) {

  # 1. 使用 data.table 读取数据（保持原样，已经很快）
  gwasResults <- fread(paste0(trait_name, ".mlma.gz"), header = TRUE)

  # 2. 使用 data.table 进行数据处理（比 tidyverse 更快）
  setDT(gwasResults)

  # 计算每条染色体的长度
  chr_len <- gwasResults[, .(chr_len = max(bp)), by = Chr]

  # 计算每条染色体的起始位置
  chr_pos <- chr_len[order(Chr)]
  chr_pos[, total := cumsum(chr_len / 10) * 10 - chr_len]
  chr_pos[, chr_len := NULL]

  # 计算累积位置和颜色
  Snp_pos <- gwasResults[chr_pos, on = "Chr"]
  setorder(Snp_pos, Chr, bp)
  Snp_pos[, BPcum := bp + total]
  Snp_pos[, Color := ifelse(Chr %% 2 == 0, "turquoise", "darkcyan")]
  Snp_pos[, log10p := -log10(p)]  # 预计算，避免重复计算

  # 计算 X 轴中心点
  X_axis <- Snp_pos[, .(center = (max(BPcum) + min(BPcum)) / 2), by = Chr]

  # 3. 【关键优化】预先准备背景矩形数据框，避免循环
  anno <- merge(chr_pos, chr_len, by = "Chr")
  anno[, end1 := total + chr_len]

  # 只保留奇数染色体的背景（偶数染色体用白色）
  rect_data <- anno[Chr %% 2 == 1]

  # 4. 预计算 y 轴范围（避免在循环中调用 layer_scales）
  y_max <- max(Snp_pos$log10p) * 1.05  # 稍微增加一点上边界

  # 5. 使用 geom_rect 一次性绘制所有背景，而不是循环调用 annotate
  p <- ggplot(Snp_pos, aes(x = BPcum, y = log10p)) +
    # 先添加背景矩形
    geom_rect(
      data = rect_data,
      aes(xmin = total, xmax = end1, ymin = 0, ymax = y_max),
      fill = 'grey80',
      alpha = 0.2,
      inherit.aes = FALSE  # 不继承主图层的 aes
    ) +
    # 添加点（使用预计算的颜色列）
    geom_point(aes(color = Color), alpha = 0.6, size = 0.8, show.legend = FALSE) +
    scale_color_identity() +  # 直接使用颜色值
    scale_x_continuous(labels = X_axis$Chr, breaks = X_axis$center) +
    # 添加阈值线
    geom_hline(yintercept = -log10(0.00001), linetype = 'dashed', color = "blue") +
    geom_hline(yintercept = -log10(0.00000005), color = "red") +
    # 应用主题
    theme_ydz0.3 +
    theme(
      panel.border = element_blank(),
      axis.line.y = element_line(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
    xlab("Position") +
    ylab("-log10(P-value)")

  # 6. 可选：添加标签（如果需要）
  # Name_mark <- Snp_pos[p <= 0.00000001]
  # if (nrow(Name_mark) > 0) {
  #   p <- p + geom_text_repel(
  #     data = Name_mark,
  #     aes(x = BPcum, y = log10p, label = SNP),
  #     size = 3
  #   )
  # }

  # 7. 保存图片
  tiff(
    paste0(trait_name, ".GWAS.optimized.tiff"),
    width = 12,
    height = 8,
    units = "in",
    res = 300,
    compression = "none"
  )
  print(p)
  dev.off()

  return(p)
}

# =============================================================================
# 使用示例
# =============================================================================

j <- "betaine"
plot_manhattan_optimized(j)

# 如果需要批量处理多个性状
# Traits <- c("betaine", "trait2", "trait3")
# for (trait in Traits) {
#   plot_manhattan_optimized(trait)
# }
