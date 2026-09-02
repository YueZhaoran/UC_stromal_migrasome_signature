##################################################
# Script: 11_GSEA_analysis.R
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
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
project_dir  <- dirname(rstudioapi::getActiveDocumentContext()$path)
cat("工作目录：", getwd(), "/n")



library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
library(tidyverse)
library(GseaVis)
library(psych)
library(reshape2)

# ============================ 1. 读取数据 ============================
# 关键基因
key_gene <- read.csv("../05DEGs/hub_gene.csv" , check.names = FALSE)

# 训练集表达矩阵（行=基因，列=样本）
dat <- read.csv("../00bulkdata/dat.GSE87466.csv", row.names = 1, check.names = FALSE)

# 训练集分组信息
group <- read.csv("../00bulkdata/group.GSE87466.csv", stringsAsFactors = FALSE)
colnames(group) <- c("sample", "group")

# 只保留训练集中属于组内的样本
dat <- dat[, colnames(dat) %in% group$sample]

# ============================ 2. 整理表达矩阵 ============================
# 基因名去重
dat$temp_symbol <- gsub("-.*", "", rownames(dat))

# 计算每个基因在所有样本中的最大值，按最大值降序排列
max_expr <- apply(dat[, -ncol(dat)], 1, max) 
dat <- dat[order(dat$temp_symbol, max_expr, decreasing = TRUE), ]

# 去重
dat <- dat[!duplicated(dat$temp_symbol), ]


rownames(dat) <- dat$temp_symbol
dat <- dat[, -which(names(dat) == "temp_symbol")]

# ============================ 3. 循环每个关键基因做 GSEA ============================
gmt_file <- "h.all.v2026.1.Hs.symbols .gmt"


for (i in 1:nrow(key_gene)) {
  gene <- key_gene$symbol[i]
  cat("正在处理基因：", gene, "/n")
  
  # 排除当前基因，计算其余基因与它的相关性
  train_data <- dat[!grepl(gene, rownames(dat)), ] %>% t() %>% as.data.frame()
  key_exp <- dat[gene, ] %>% t() %>% as.data.frame()
  
  # 计算 Spearman 相关性
  d <- corr.test(key_exp, train_data, use = "complete", method = "spearman")
  r <- data.frame(t(d$r))
  colnames(r) <- "cor_coefficient"
  r$symbol <- rownames(r)
  r <- r[order(r$cor_coefficient, decreasing = TRUE), ]
  r <- dplyr::select(r, -symbol)
  r <- na.omit(r)
  
  # 创建基因排序列表（向量）
  geneList <- r$cor_coefficient
  names(geneList) <- rownames(r)
  
  # GSEA 分析
  gmt <- read.gmt(gmt_file)
  set.seed(1)
  res_GSEA <- GSEA(geneList, TERM2GENE = gmt, pvalueCutoff = 1, eps = 0)
  
  # 提取结果
  sortGESA <- data.frame(res_GSEA)
  sortGESA <- sortGESA[order(sortGESA$pvalue), ]
  sortGESA <- sortGESA[sortGESA$pvalue < 0.05, ]
  sortGESA <- subset(sortGESA, abs(sortGESA$NES) >= 1)
  
  if (nrow(sortGESA) > 5) {
    paths <- rownames(sortGESA)[1:5]
  } else {
    paths <- sortGESA$ID
  }
  
  write.csv(sortGESA, file = paste0(i, "_", gene, "_gsea_res.csv"))
  
  if (length(paths) > 0) {
    p <- gseaNb(object = res_GSEA,
                geneSetID = paths,
                subPlot = 2,
                termWidth = 45,
                addPval = TRUE,
                rmHt = FALSE,
                pvalX = 1,
                pvalY = 1.2,
                newGsea = FALSE,
                curveCol = c('#FFDAB9', '#00FFFF', 'pink', 'gray', '#90EE90'))#+
      #ggtitle(paste(gene, "GSEA")) + 
      #theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    pdf(file = paste0(i, "_", gene, "_GSEA.pdf"), width = 8, height = 6)
    print(p)
    dev.off()
    
    png(file = paste0(i, "_", gene, "_GSEA.png"), width = 8, height = 6, units = "in", res = 600)
    print(p)
    dev.off()
  } else {
    cat("基因", gene, "没有显著富集通路。/n")
  }
}