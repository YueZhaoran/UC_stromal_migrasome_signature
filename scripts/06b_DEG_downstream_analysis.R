##################################################
# Script: 06b_DEG_downstream_analysis.R
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

library(VennDiagram)
library(ggplot2)
library(dplyr)
library(ggvenn)

# 整理三个基因列表

# 读取两份 PCC 数据
Stromal_pcc <- read.csv("../02scPagwas/scPagwas/Stromal_gene_PCC.csv", row.names = 1)


# 提取“UC致病基因”（阈值：|PCC| >= 0.1）
uc_genes <- unique(c(
  rownames(Stromal_pcc)[abs(Stromal_pcc$PCC) >= 0.1]
))
print(paste0("提取出 UC致病基因 数量：", length(uc_genes)))

# 读取“细胞特异性迁移体相关基因”
cell_mrgs <- read.csv("../03Scissor/cell_specific_MRGs.csv")
cell_mrgs_genes <- cell_mrgs$gene
print(paste0("细胞特异性 MRGs 数量：", length(cell_mrgs_genes)))

# 读取“Bulk DEGs” (你刚跑出来的)
bulk_degs <- read.csv("../04limma/bulk_DEGs_GSE87466.csv")
bulk_degs_genes <- bulk_degs$gene
print(paste0("Bulk DEGs 数量：", length(bulk_degs_genes)))



# 取交集并绘制韦恩图


# 计算三者的交集
intersection_genes <- Reduce(intersect, list(uc_genes, cell_mrgs_genes, bulk_degs_genes))
print(paste0("最终锁定候选基因数量：", length(intersection_genes)))



write.csv(data.frame(gene = intersection_genes), "Venn_intersection_genes.csv", row.names = FALSE)

# 用 VennDiagram 画图
datalist <- list(
  DEG  = bulk_degs_genes,   # Bulk 差异基因
  MRGs = cell_mrgs_genes,   # Scissor基因
  scRNA = uc_genes          # scPagwas UC 致病基因
)

p <- ggvenn(
  datalist,
  columns = c('DEG', 'MRGs', 'scRNA'),       # 指定要画哪些集合
  fill_color = c("#ffb2b2", "turquoise", '#D6E7A3'),
  show_percentage = TRUE,                    # 显示交集基因数量占该集合的百分比
  stroke_alpha = 0.5,
  stroke_size = 0.5,                         # 交集处白边的大小
  stroke_color = "white",
  stroke_linetype = "solid",
  text_size = 5,
  set_name_color = c("#ffb2b2", "turquoise", '#D6E7A3'),
  text_color = 'black'
)


ggsave('01.hub_gene.pdf', p, width = 8, height = 8)
ggsave('01.hub_gene.png', p, width = 8, height = 8)


