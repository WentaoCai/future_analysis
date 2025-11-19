# 🧬 GWAS曼哈顿图基因注释指南

## 功能概述

这个增强版本可以自动：

1. **识别QTL区域**（显著性SNP聚集的区域）
2. **找到每个QTL区域最显著的lead SNP**
3. **查找距离lead SNP最近的基因**
4. **在曼哈顿图上标注基因名**

---

## 🚀 快速开始（3步）

### 步骤1：准备GTF注释文件

下载你的物种的GTF文件：

**人类 (GRCh38):**
```bash
wget http://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz
```

**小鼠 (GRCm39):**
```bash
wget http://ftp.ensembl.org/pub/release-109/gtf/mus_musculus/Mus_musculus.GRCm39.109.gtf.gz
```

**拟南芥 (TAIR10):**
```bash
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-52/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.52.gtf.gz
```

### 步骤2：加载脚本

```r
setwd("~/Documents/Software/manhattan2D/")
source("gwas_manhattan_annotated.R")
```

### 步骤3：绘制带基因注释的曼哈顿图

```r
result <- plot_manhattan_with_genes(
  trait_name = "betaine",
  gtf_file = "Homo_sapiens.GRCh38.109.gtf.gz"
)
```

**输出文件：**
- `betaine.GWAS.annotated.tiff` - 带基因标注的曼哈顿图
- `betaine.QTL_regions.txt` - QTL区域汇总表

---

## 📊 核心功能

### 1. QTL区域识别

**算法：**

1. 识别所有显著性SNP（p < 阈值）
2. 将相邻的显著性SNP合并成QTL区域
3. 对于每个QTL区域，找到最显著的"lead SNP"

**参数控制：**

```r
result <- plot_manhattan_with_genes(
  trait_name = "betaine",
  gtf_file = "genome.gtf.gz",
  p_threshold = 1e-5,      # 显著性阈值
  merge_distance = 1e6     # 合并距离（1 Mb）
)
```

**示例：**

假设你有以下显著性SNP：

```
Chr1: 1,000,000 bp (p=1e-7)
Chr1: 1,200,000 bp (p=1e-8)  ← Lead SNP
Chr1: 1,500,000 bp (p=1e-6)
Chr1: 5,000,000 bp (p=1e-7)  ← 新QTL（距离>1Mb）
```

会被识别为2个QTL区域：
- **QTL1**: Chr1:1,000,000-1,500,000，lead SNP在1,200,000
- **QTL2**: Chr1:5,000,000-5,000,000，lead SNP在5,000,000

### 2. 最近基因查找

对于每个lead SNP，计算到所有基因的距离：

- **SNP在基因内**：距离 = 0
- **SNP在基因上游**：距离 = gene_start - SNP_position
- **SNP在基因下游**：距离 = SNP_position - gene_end

选择距离最小的基因进行标注。

### 3. 基因标注

使用 `ggrepel` 智能标注基因名，避免重叠。

---

## ⚙️ 参数详解

### `plot_manhattan_with_genes()` 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `trait_name` | **必需** | 性状名称（文件前缀） |
| `gtf_file` | **必需** | GTF注释文件路径 |
| `p_threshold` | `1e-5` | QTL识别的显著性阈值 |
| `merge_distance` | `1e6` (1Mb) | 合并邻近QTL的最大距离 |
| `max_labels` | `20` | 最多标注的QTL数量 |
| `label_only_top` | `NULL` | 只标注最显著的N个QTL |
| `use_binning` | `TRUE` | 是否使用分箱加速 |
| `n_bins` | `5000` | 分箱数量 |
| `use_scattermore` | `TRUE` | 是否使用scattermore快速绘图 |
| `output_format` | `"tiff"` | 输出格式（tiff/png/pdf） |

---

## 📝 使用场景

### 场景1：基本用法

```r
result <- plot_manhattan_with_genes(
  trait_name = "betaine",
  gtf_file = "genome.gtf.gz"
)
```

使用默认参数，适合大多数情况。

### 场景2：全基因组显著性（p < 5e-8）

```r
result <- plot_manhattan_with_genes(
  trait_name = "betaine",
  gtf_file = "genome.gtf.gz",
  p_threshold = 5e-8,      # GWAS通用阈值
  merge_distance = 1e6
)
```

### 场景3：只标注Top 10最显著的QTL

```r
result <- plot_manhattan_with_genes(
  trait_name = "betaine",
  gtf_file = "genome.gtf.gz",
  label_only_top = 10      # 只标注前10个
)
```

适合有很多显著性区域，但只想关注最重要的。

### 场景4：严格的QTL定义（500kb）

```r
result <- plot_manhattan_with_genes(
  trait_name = "betaine",
  gtf_file = "genome.gtf.gz",
  p_threshold = 1e-6,
  merge_distance = 5e5     # 500kb
)
```

更保守的QTL识别，避免将远距离的SNP合并。

### 场景5：PDF输出（用于发表）

```r
result <- plot_manhattan_with_genes(
  trait_name = "betaine",
  gtf_file = "genome.gtf.gz",
  output_format = "pdf"    # 矢量图，可无限缩放
)
```

### 场景6：批量处理多个性状

```r
traits <- c("betaine", "glucose", "cholesterol")
gtf_file <- "Homo_sapiens.GRCh38.109.gtf.gz"

for(trait in traits) {
  result <- plot_manhattan_with_genes(
    trait_name = trait,
    gtf_file = gtf_file,
    label_only_top = 15
  )
}
```

### 场景7：并行批量处理

```r
library(parallel)

traits <- c("betaine", "glucose", "cholesterol", "BMI", "height")

results <- mclapply(traits, function(trait) {
  plot_manhattan_with_genes(
    trait_name = trait,
    gtf_file = "genome.gtf.gz"
  )
}, mc.cores = 4)
```

---

## 📂 输出文件

### 1. 曼哈顿图（TIFF/PNG/PDF）

**文件名：** `{trait_name}.GWAS.annotated.{format}`

包含：
- 所有SNP的散点图
- 显著性阈值线（蓝色虚线：p=1e-5，红色实线：p=5e-8）
- QTL区域lead SNP的基因标注

### 2. QTL汇总表（TXT）

**文件名：** `{trait_name}.QTL_regions.txt`

**列说明：**

| 列名 | 说明 |
|------|------|
| `Chr` | 染色体编号 |
| `start` | QTL区域起始位置（bp） |
| `end` | QTL区域终止位置（bp） |
| `width_kb` | QTL区域宽度（kb） |
| `lead_snp_p` | Lead SNP的p值 |
| `gene_name` | 最近的基因名 |
| `gene_distance_kb` | 到基因的距离（kb） |

**示例：**

```
Chr  start      end        width_kb  lead_snp_p  gene_name  gene_distance_kb
1    1000000    1500000    500       1.23e-8     APOE       0.0
2    5000000    5200000    200       4.56e-7     BRCA2      15.3
5    10500000   10800000   300       7.89e-9     TP53       2.1
```

---

## 🔬 GTF文件格式

### 标准GTF格式

```
Chr1  HAVANA  gene  11869  14409  .  +  .  gene_id "ENSG00000223972"; gene_name "DDX11L1";
Chr1  HAVANA  gene  14404  29570  .  -  .  gene_id "ENSG00000227232"; gene_name "WASH7P";
```

### 关键字段

- **列1**: 染色体（如Chr1, 1, chrI）
- **列3**: 特征类型（gene, transcript, exon等）
- **列4-5**: 起始和终止位置
- **列9**: 属性（包含gene_name或gene_id）

### 染色体命名规范

脚本会自动处理以下格式：
- `Chr1` → `1`
- `chr1` → `1`
- `1` → `1`

**注意：** 只保留数字染色体（1-22, 1-19等），跳过X, Y, MT等。

---

## 💡 最佳实践

### 1. 选择合适的 `p_threshold`

| 应用场景 | 推荐阈值 | 说明 |
|---------|---------|------|
| 探索性分析 | `1e-5` | 更宽松，发现更多候选区域 |
| GWAS标准 | `5e-8` | 全基因组显著性 |
| 严格筛选 | `1e-10` | 只关注极显著的信号 |

### 2. 选择合适的 `merge_distance`

| 应用场景 | 推荐值 | 说明 |
|---------|--------|------|
| LD衰减快（如植物） | `100kb` | 更严格的QTL定义 |
| 人类GWAS | `500kb-1Mb` | 标准LD范围 |
| 宽松定义 | `2Mb` | 合并更多候选区域 |

### 3. 控制标注数量

```r
# 方法1：限制最大标注数
plot_manhattan_with_genes(..., max_labels = 20)

# 方法2：只标注最显著的
plot_manhattan_with_genes(..., label_only_top = 10)

# 方法3：两者结合（推荐）
plot_manhattan_with_genes(..., label_only_top = 15, max_labels = 20)
```

### 4. 性能优化

对于大数据集（1000万+ SNPs）：

```r
result <- plot_manhattan_with_genes(
  trait_name = "betaine",
  gtf_file = "genome.gtf.gz",
  use_binning = TRUE,         # 开启分箱
  n_bins = 5000,              # 高分辨率
  use_scattermore = TRUE      # 快速绘图
)
```

预期时间：10-15秒

---

## 🎨 自定义图形

### 方法1：修改返回的ggplot对象

```r
result <- plot_manhattan_with_genes(
  trait_name = "betaine",
  gtf_file = "genome.gtf.gz"
)

# 自定义
library(ggplot2)

custom_plot <- result$plot +
  theme(
    plot.title = element_text(size = 18, face = "bold", color = "darkblue"),
    axis.title = element_text(size = 16),
    panel.background = element_rect(fill = "white")
  ) +
  ggtitle("My Custom Title")

# 保存
ggsave("custom_plot.png", custom_plot, width = 14, height = 8, dpi = 300)
```

### 方法2：修改源代码

编辑 `gwas_manhattan_annotated.R`，找到绘图部分自定义：

```r
p <- p +
  geom_text_repel(
    data = qtl_regions,
    aes(x = BPcum, y = log10p, label = gene_name),
    size = 4,                    # 修改字体大小
    fontface = "bold.italic",    # 修改字体样式
    color = "red",               # 修改颜色
    box.padding = 1.0,           # 增加标签间距
    ...
  )
```

---

## 🧪 使用结果

### 访问QTL区域表

```r
result <- plot_manhattan_with_genes(...)

# 查看QTL表
qtl_table <- result$qtl_regions
print(qtl_table)

# 导出为CSV
library(data.table)
fwrite(qtl_table, "my_qtls.csv")
```

### 筛选特定QTL

```r
# 筛选特定染色体
chr1_qtls <- qtl_table[Chr == 1]

# 筛选极显著的
highly_sig <- qtl_table[lead_snp_p < 1e-10]

# 筛选基因内的QTL（距离=0）
in_gene <- qtl_table[gene_distance == 0]

# 排序
top_qtls <- qtl_table[order(lead_snp_p)][1:10]
```

---

## 🐛 故障排除

### 问题1：GTF文件读取失败

**错误：**
```
Error in fread(gtf_file) : ...
```

**解决方案：**

1. 检查文件是否存在：
   ```r
   file.exists("genome.gtf.gz")
   ```

2. 确认GTF格式正确（使用标准Ensembl格式）

3. 尝试手动解压：
   ```bash
   gunzip genome.gtf.gz
   ```
   然后使用 `.gtf` 文件。

### 问题2：没有找到QTL区域

**可能原因：**

1. `p_threshold` 太严格
   ```r
   # 尝试更宽松的阈值
   plot_manhattan_with_genes(..., p_threshold = 1e-4)
   ```

2. 数据中没有显著性信号
   ```r
   # 检查p值分布
   gwas <- fread("betaine.mlma.gz")
   summary(gwas$p)
   sum(gwas$p < 1e-5)  # 有多少显著SNP？
   ```

### 问题3：基因名为NA

**可能原因：**

1. GTF文件没有该染色体的基因
2. SNP距离所有基因太远

**解决方案：**

检查GTF文件覆盖的染色体：
```r
genes <- parse_gtf_fast("genome.gtf.gz")
table(genes$Chr)
```

### 问题4：标签重叠

**解决方案：**

1. 减少标注数量：
   ```r
   plot_manhattan_with_genes(..., label_only_top = 5)
   ```

2. 增加图片宽度：
   ```r
   # 修改脚本中的保存部分
   tiff(..., width = 20, height = 10)  # 更宽的图
   ```

3. 调整 `ggrepel` 参数（修改脚本）：
   ```r
   geom_text_repel(
     ...,
     box.padding = 1.0,      # 增加
     max.overlaps = 20       # 增加
   )
   ```

---

## 📊 性能基准

| 数据规模 | GTF大小 | 处理时间 |
|---------|---------|---------|
| 100万 SNPs | 50MB | ~5秒 |
| 500万 SNPs | 50MB | ~8秒 |
| 1000万 SNPs | 50MB | ~12秒 |
| 1000万 SNPs | 200MB | ~15秒 |

**注：** 使用 `use_binning=TRUE` 和 `use_scattermore=TRUE`

---

## 🎓 引用和致谢

如果你在发表的文章中使用了这个工具，建议引用：

- **Rcpp**: Eddelbuettel & François (2011)
- **data.table**: Dowle & Srinivasan (2021)
- **ggrepel**: Slowikowski (2021)
- **ggplot2**: Wickham (2016)

---

## 📚 延伸阅读

- [Ensembl GTF 文件格式](https://www.ensembl.org/info/website/upload/gff.html)
- [GWAS Catalog](https://www.ebi.ac.uk/gwas/)
- [ggrepel 文档](https://ggrepel.slowkow.com/)

---

## 💬 常见问题（FAQ）

### Q1: 支持非人类物种吗？

**A:** 支持！只需提供对应物种的GTF文件：
- 小鼠、大鼠：Ensembl
- 拟南芥、水稻：Ensembl Plants
- 斑马鱼、果蝇：Ensembl Metazoa

### Q2: 可以使用GFF3格式吗？

**A:** GTF和GFF3格式类似但不完全相同。建议：
1. 使用GTF格式（推荐）
2. 或将GFF3转换为GTF：
   ```bash
   gffread input.gff3 -T -o output.gtf
   ```

### Q3: 如何标注特定基因？

**A:** 修改脚本，手动添加标注：
```r
# 在生成plot后
specific_genes <- data.frame(
  gene_name = c("APOE", "BRCA2"),
  Chr = c(19, 13),
  position = c(45411941, 32315474)
)

# 添加到图中
p <- p + geom_text(data = specific_genes, ...)
```

### Q4: 性能还是慢怎么办？

**A:**
1. 确保使用了分箱：`use_binning = TRUE`
2. 减少分箱数：`n_bins = 2000`
3. 预先筛选GTF（只保留目标染色体）
4. 使用SSD硬盘

---

**最后更新：** 2025-11-19

祝你的GWAS分析顺利！🎉
