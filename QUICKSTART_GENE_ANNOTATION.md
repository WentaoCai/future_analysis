# 🧬 基因注释快速开始

## TL;DR - 3步搞定

```r
# 1. 下载GTF文件（只需一次）
# wget http://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz

# 2. 加载脚本
source("gwas_manhattan_annotated.R")

# 3. 绘制带基因注释的曼哈顿图
result <- plot_manhattan_with_genes(
  trait_name = "betaine",
  gtf_file = "Homo_sapiens.GRCh38.109.gtf.gz"
)
```

**输出：**
- `betaine.GWAS.annotated.tiff` - 曼哈顿图（带基因标注）
- `betaine.QTL_regions.txt` - QTL区域汇总表

---

## 🎯 功能说明

这个工具会自动：

1. ✅ **识别QTL区域**（显著性SNP聚集区域）
2. ✅ **找到每个QTL的lead SNP**（最显著的SNP）
3. ✅ **查找最近的基因**
4. ✅ **在图上标注基因名**

**示例输出：**

```
找到15个QTL区域：

Chr  Position    P-value    Gene      Distance
1    1,234,567   1.2e-8     APOE      0 kb (在基因内)
2    5,678,901   4.5e-7     BRCA2     15 kb
5    10,234,567  7.8e-9     TP53      2 kb
...
```

---

## 📦 准备GTF文件

### 人类 (GRCh38)

```bash
wget http://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz
```

大小：~60 MB（压缩），包含约6万个基因

### 小鼠 (GRCm39)

```bash
wget http://ftp.ensembl.org/pub/release-109/gtf/mus_musculus/Mus_musculus.GRCm39.109.gtf.gz
```

### 拟南芥 (TAIR10)

```bash
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-52/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.52.gtf.gz
```

### 其他物种

访问 [Ensembl FTP](http://ftp.ensembl.org/pub/) 下载你的物种的GTF文件。

---

## 💻 基本用法

### 最简单的用法

```r
source("gwas_manhattan_annotated.R")

result <- plot_manhattan_with_genes(
  trait_name = "betaine",
  gtf_file = "genome.gtf.gz"
)
```

### 自定义参数

```r
result <- plot_manhattan_with_genes(
  trait_name = "betaine",
  gtf_file = "genome.gtf.gz",
  p_threshold = 5e-8,      # 全基因组显著性
  merge_distance = 1e6,    # 1 Mb合并距离
  label_only_top = 10      # 只标注Top 10
)
```

### 批量处理

```r
traits <- c("betaine", "glucose", "cholesterol")

for(trait in traits) {
  plot_manhattan_with_genes(
    trait_name = trait,
    gtf_file = "genome.gtf.gz"
  )
}
```

---

## 🎨 输出示例

### 曼哈顿图

![示例图](https://via.placeholder.com/800x400?text=Manhattan+Plot+with+Gene+Labels)

- 蓝色虚线：p = 1e-5
- 红色实线：p = 5e-8
- 斜体文字：基因名（自动标注，避免重叠）

### QTL汇总表

```
Chr  start      end        width_kb  lead_snp_p  gene_name  gene_distance_kb
1    1000000    1500000    500       1.23e-8     APOE       0.0
2    5000000    5200000    200       4.56e-7     BRCA2      15.3
5    10500000   10800000   300       7.89e-9     TP53       2.1
```

---

## ⚙️ 常用参数组合

### 场景1：探索性分析（发现更多候选）

```r
plot_manhattan_with_genes(
  trait_name = "betaine",
  gtf_file = "genome.gtf.gz",
  p_threshold = 1e-5       # 宽松阈值
)
```

### 场景2：发表级图片（严格筛选）

```r
plot_manhattan_with_genes(
  trait_name = "betaine",
  gtf_file = "genome.gtf.gz",
  p_threshold = 5e-8,      # GWAS标准
  label_only_top = 15,     # 只标注Top 15
  output_format = "pdf"    # 矢量图
)
```

### 场景3：关注最重要的信号

```r
plot_manhattan_with_genes(
  trait_name = "betaine",
  gtf_file = "genome.gtf.gz",
  label_only_top = 5,      # 只看Top 5
  p_threshold = 1e-6
)
```

### 场景4：高密度LD区域（人类）

```r
plot_manhattan_with_genes(
  trait_name = "betaine",
  gtf_file = "genome.gtf.gz",
  merge_distance = 5e5     # 500 kb合并
)
```

---

## 📊 参数速查表

| 参数 | 默认值 | 推荐值 | 说明 |
|------|--------|--------|------|
| `p_threshold` | `1e-5` | `5e-8` (GWAS), `1e-5` (探索) | 显著性阈值 |
| `merge_distance` | `1e6` | `5e5-2e6` | 合并QTL距离（bp） |
| `label_only_top` | `NULL` | `10-20` | 只标注最显著的N个 |
| `max_labels` | `20` | `15-30` | 最大标注数 |
| `output_format` | `"tiff"` | `"pdf"` (发表), `"png"` (预览) | 输出格式 |

---

## 🔍 检查结果

### 查看QTL表

```r
result <- plot_manhattan_with_genes(...)

# 查看QTL汇总
print(result$qtl_regions)

# 保存为CSV
library(data.table)
fwrite(result$qtl_regions, "qtl_summary.csv")
```

### 筛选特定QTL

```r
# 只看第1号染色体的QTL
chr1 <- result$qtl_regions[Chr == 1]

# 只看极显著的（p < 1e-10）
highly_sig <- result$qtl_regions[lead_snp_p < 1e-10]

# 只看基因内的QTL（距离=0）
in_gene <- result$qtl_regions[gene_distance == 0]
```

---

## 🚀 性能

| 数据规模 | 处理时间 |
|---------|---------|
| 100万 SNPs | ~5秒 |
| 500万 SNPs | ~8秒 |
| 1000万 SNPs | ~12秒 |

*使用默认参数（binning + scattermore）*

---

## 🐛 常见问题

### Q: 找不到GTF文件？

```bash
# 检查文件是否存在
ls -lh *.gtf.gz

# 如果没有，重新下载
wget http://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz
```

### Q: 没有发现QTL区域？

```r
# 1. 检查是否有显著性SNP
gwas <- fread("betaine.mlma.gz")
sum(gwas$p < 1e-5)  # 有多少个？

# 2. 降低阈值
plot_manhattan_with_genes(..., p_threshold = 1e-4)
```

### Q: 标签重叠太多？

```r
# 方法1：减少标注数
plot_manhattan_with_genes(..., label_only_top = 5)

# 方法2：增加图片宽度（修改脚本）
tiff(..., width = 20, height = 10)
```

### Q: 基因名显示为NA？

```r
# 检查GTF文件是否覆盖该染色体
genes <- parse_gtf_fast("genome.gtf.gz")
table(genes$Chr)
```

---

## 📚 完整文档

详细说明请参考：

- **完整指南**: `GENE_ANNOTATION_GUIDE.md`
- **使用示例**: `example_gene_annotation.R`
- **性能优化**: `RCPP_INSTALLATION_GUIDE.md`

---

## 🎉 完整示例

```r
#!/usr/bin/env Rscript

# 设置工作目录
setwd("~/Documents/Software/manhattan2D/")

# 加载脚本
source("gwas_manhattan_annotated.R")

# 单个性状
result <- plot_manhattan_with_genes(
  trait_name = "betaine",
  gtf_file = "Homo_sapiens.GRCh38.109.gtf.gz",
  p_threshold = 5e-8,
  label_only_top = 15,
  output_format = "pdf"
)

# 查看结果
print(result$qtl_regions)

# 批量处理
traits <- c("betaine", "glucose", "cholesterol")

for(trait in traits) {
  plot_manhattan_with_genes(
    trait_name = trait,
    gtf_file = "Homo_sapiens.GRCh38.109.gtf.gz"
  )
}

cat("All done! 🎉\n")
```

---

**开始使用吧！** 任何问题查看 `GENE_ANNOTATION_GUIDE.md` 📖
