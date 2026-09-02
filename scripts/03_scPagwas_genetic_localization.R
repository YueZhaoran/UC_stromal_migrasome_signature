##################################################
# Script: 03_scPagwas_genetic_localization.R
#
# Purpose:
# Analysis workflow for the study:
# Integrating genetic risk and single-cell transcriptomics
# identifies a stromal migrasome-related signature in ulcerative colitis
#
# This script was prepared for reproducible research.
# Please configure input/output paths before execution.
##################################################
# Input files should be placed under data/ and resources/.
# Output files should be written under results/.
##################################################

rm(list = ls())
gc()

# ============================ 环境准备 ============================
# 设置项目目录
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
project_dir  <- dirname(rstudioapi::getActiveDocumentContext()$path)
cat("工作目录：", getwd(), "/n")

# 加载必要的R包
library(Seurat)
library(gwasvcf)
library(VariantAnnotation)
library(scPagwas)
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpubr) 
library(patchwork)
library(ggsci)
library(rtracklayer)

library(pdftools)
#读取GWAS原始VCF数据，提取字段

vcf <- readVcf("ieu-a-32.vcf.gz")
vcf_tbl <- vcf_to_granges(vcf) %>% dplyr::as_tibble()

head(vcf_tbl)


gwas_raw <- vcf_tbl[,c("seqnames", "start", "ID", "SE", "ES")]
gwas_raw$maf <- vcf_tbl$AF

head(gwas_raw)
colnames(gwas_raw) <- c("chrom", "pos", "rsid", "se", "beta", "maf")
write.table(gwas_raw, "GWAS_summ.txt", quote = F, sep = " ", row.names = F, col.names = T)


#基因组版本坐标转换

#读取 GWAS 文件和 Chain 文件
gwas_summ <- fread("GWAS_summ.txt")
chain_file <- import.chain("hg19ToHg38.over.chain")

#转为基因组区间对象
gr <- GRanges(
  seqnames = paste0("chr", gwas_summ$chrom), 
  ranges = IRanges(start = gwas_summ$pos, end = gwas_summ$pos)
)

#执行 LiftOver 转换
gr_hg38_list <- liftOver(gr, chain_file)

success_mask <- lengths(gr_hg38_list) > 0 
um(success_mask)

# 提取原始数据中转换成功的那些行
gwas_intersect <- gwas_summ[success_mask, ]


gr_hg38 <- unlist(gr_hg38_list)

#把旧的 hg19 坐标直接替换为新的 hg38 坐标
gwas_intersect$chrom <- gsub("^chr", "", as.character(seqnames(gr_hg38)))
gwas_intersect$pos <- start(gr_hg38)

# 检查最终行数
row(gwas_intersect)

#保存最终文件
write.csv(gwas_intersect, "GWAS_summ5.csv", quote = F, row.names = F)


#处理好的单细胞数据
v5_file <- '../01singleCell/05_scData_annotated.rds'


#运行 scPagwas 核心算法
Pagwas <- scPagwas_main(
  Pagwas = NULL,
  gwas_data = "GWAS_summ5.csv",
  Single_data = v5_file,
  output.prefix = "af",
  output.dirs = "scPagwas",
  Pathway_list = Genes_by_pathway_kegg,
  n.cores = 15,
  assay = "RNA",
  singlecell = TRUE, 
  iters_singlecell = 100,      # 做100次单细胞水平置换检验，计算显著性
  celltype = TRUE,             # 开启细胞类型特异性分析
  iters_celltype = 200,
  block_annotation = block_annotation,
  chrom_ld = chrom_ld
)
saveRDS(Pagwas,"0.Pagwas.rds")



str(Pagwas@misc$bootstrap_results)
head(Pagwas@misc$bootstrap_results)
colnames(Pagwas@misc$bootstrap_results)

# ==========================================
# 结果可视化
# ==========================================
# Bootstrap 置换检验统计图
Pagwas <- readRDS("0.Pagwas.rds")

Bootstrap_P_Barplot(
  p_results = Pagwas@misc$bootstrap_results$bp_value[-1],
  p_names = rownames(Pagwas@misc$bootstrap_results)[-1],
  figurenames = paste0("./scPagwas/","3.Bootstrap_P_Barplot.pdf"),
  do_plot = T, title = "Bootstrap results"
)
pdf_convert("./scPagwas/3.Bootstrap_P_Barplot.pdf", format = "png", 
            dpi = 300, filenames = "./scPagwas/3.Bootstrap_P_Barplot.png")



Bootstrap_estimate_Plot(
  bootstrap_results = Pagwas@misc$bootstrap_results,
  figurenames = paste0("./scPagwas/","3.estimateplot.pdf"),
  do_plot = T
)
p_record <- recordPlot()

# 保存为正确的 PDF
pdf("./scPagwas/3.estimateplot.pdf", width = 8, height = 5)
replayPlot(p_record)
dev.off()

# 保存为 PNG
png("./scPagwas/3.estimateplot.png", width = 8, height = 5, units = "in", res = 300)
replayPlot(p_record)
dev.off()



#TRS（性状相关得分）在细胞上的分布
#目标：比较"目标"和其他所有细胞的评分差异
select_cells <- "Stromal"
compare_cells <- levels(Idents(Pagwas))[!levels(Idents(Pagwas)) %in% select_cells]
my_comparisons <- lapply(compare_cells, function(x) c(select_cells, x))

afcolor <- ggsci::pal_npg()(7)
p1 <- DimPlot(Pagwas) + ggtitle("Cell type") + scale_color_manual(values=afcolor) + NoLegend() 
p2 <- FeaturePlot(Pagwas, "scPagwas.TRS.Score1", cols = c("#C3DBD9", "#990000"))
p3 <- VlnPlot(Pagwas, features = "scPagwas.TRS.Score1", pt.size = 0) +
  stat_compare_means(comparisons = my_comparisons) + ylim(0, 1.5) +
  scale_fill_manual(values=afcolor) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

p1 + p2 + p3 + plot_layout(ncol = 3)
ggsave("./scPagwas/1.Cell_TRS.pdf", width = 15, height = 5)
ggsave("./scPagwas/1.Cell_TRS.png", width = 15, height = 5)

#细胞类型显著性柱状图与森林图
df <- read.csv("./scPagwas/af_Merged_celltype_pvalue.csv")
df <- df %>% select(celltype, pvalue) %>%
  mutate(logP = -log10(pvalue), Significant = ifelse(pvalue < 0.05, "Yes", "No")) %>%
  arrange(logP)

p <- ggplot(df, aes(x = reorder(celltype, logP), y = logP, fill = Significant)) +
  geom_col(width = 0.7) + coord_flip() +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "red") +
  geom_text(aes(label = signif(pvalue, 3)), hjust = -0.1) +
  labs(x = "Cell type", y = expression(-log[10](P)), title = "scPagwas cell-type enrichment") +
  scale_fill_manual(values = c("Yes" = "#D55E00", "No" = "#999999")) + theme_bw()
ggsave("./scPagwas/scPagwas_celltype_barplot.pdf", p, width = 8, height = 5)
ggsave("./scPagwas/scPagwas_celltype_barplot.png", p, width = 8, height = 5)

# ==========================================
# 针对特定目标细胞，提取核心相关基因
# ==========================================
key_celltype <- "Stromal"
cells_use <- rownames(Pagwas@meta.data)[Pagwas@meta.data$external_annotation == key_celltype]

# 提取gPAS（基因活性得分）和细胞表达矩阵
gpas_sub <- Pagwas@meta.data[cells_use, "scPagwas.gPAS.score"]
expr_mat <- GetAssayData(Pagwas, assay = "RNA", layer = "data") # 提取归一化后的表达量
data_sub <- expr_mat[, cells_use, drop = FALSE]

# 计算Pearson相关系数（PCC）：基因表达量 vs 疾病关联分值
PCC_macro <- scGet_PCC(scPagwas.gPAS.score = gpas_sub, data_mat = data_sub)
write.csv(PCC_macro, "./scPagwas/Stromal_gene_PCC.csv")



# 提取出正/负相关性最强的基因（筛选出PCC > 0.05 和 < -0.05的）
pcc <- read.csv("./scPagwas/Stromal_gene_PCC.csv", row.names = 1)
pcc$gene <- rownames(pcc)
top_pos <- pcc[pcc$PCC > 0.05, ] %>% arrange(desc(PCC))
top_neg <- pcc[pcc$PCC < -0.05, ] %>% arrange(PCC)
top_500 <- rbind(top_pos %>% mutate(direction="Positive"), 
                 top_neg %>% mutate(direction="Negative"))
write.csv(top_500, "Stromal_scpaGWAS.csv", row.names = FALSE)



# 读取 PCC 数据
nk_pcc <- read.csv("./scPagwas/Stromal_gene_PCC.csv", row.names = 1)

# 添加细胞类型列
nk_pcc$celltype <- "Stromal"



# 提取 Top 30 基因画棒棒糖图

plot_pcc_lollipop <- function(df, cell_name, color_code) {
  # 按 PCC 绝对值排序，取前 30 个（正负各取最强的）
  df_top <- df %>%
    mutate(genes = rownames(.)) %>%
    arrange(desc(abs(PCC))) %>%
    head(30) %>%
    mutate(genes = factor(genes, levels = genes[order(PCC)])) # 按 PCC 排序因子
  
  p <- ggplot(df_top, aes(x = genes, y = PCC, color = PCC)) +
    geom_segment(aes(x = genes, xend = genes, y = 0, yend = PCC), size = 1.2) +
    geom_point(size = 3) +
    coord_flip() +
    scale_color_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0) +
    labs(x = "Gene", y = "Pearson Correlation (PCC)", 
         title = paste0("Top 30 PCC Genes in ", cell_name)) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 10, color = "black"),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  return(p)
}

# 生成 Stromal 棒棒糖图
p_nk <- plot_pcc_lollipop(nk_pcc, "Stromal", "#D55E00")
ggsave("./scPagwas/Stromal_gene_PCC_Top30.pdf", p_nk, width = 8, height = 6)
ggsave("./scPagwas/Stromal_gene_PCC_Top30.png", p_nk, width = 8, height = 6)



# 使用 scPagwas 官方函数绘图

# -------------- Stromal ----------------
heritability_cor_scatterplot(
  gene_heri_cor = nk_pcc,      # 数据框
  topn_genes_label = 10,       # 标注 PCC 绝对值 Top 10 的基因
  color_low = "#035397",       # 负相关（蓝色）
  color_high = "#F32424",      # 正相关（红色）
  color_mid = "white",         # 中间色
  text_size = 3,               # 基因标签字体大小
  do_plot = TRUE,
  max.overlaps = 20,
  width = 7,
  height = 7
)
ggsave("./scPagwas/Stromal_heritability_cor_scatter.pdf", width = 7, height = 7)
ggsave("./scPagwas/Stromal_heritability_cor_scatter.png", width = 7, height = 7)
