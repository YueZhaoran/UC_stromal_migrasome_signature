##################################################
# Script: 04_migrasome_activity_scoring.R
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
# 设置项目目录（RStudio 中运行自动定位；Rscript 运行时依赖外部工作目录）
if (rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
cat("工作目录：", getwd(), "/n")

library(GSVA)
library(readxl)
library(ggplot2)
library(ggpubr)
library(dplyr)
set.seed(798)

# ============================ 1. 读取数据 ============================
# 读取训练集表达矩阵
dat_bulk <- read.csv("../00bulkdata/dat.GSE87466.csv", row.names = 1, check.names = FALSE)
# 读取分组信息
group_info <- read.csv("../00bulkdata/group.GSE87466.csv", stringsAsFactors = FALSE)

table(group_info$group)

# 迁移体相关基因集
datasets <- read_excel("../migrasome-genes.xls")
mrg_list <- list(MRGs = datasets$genes)

# ============================ 2. ssGSEA 计算每个样本的 MRGs 评分 ============================
ssgsea_param <- ssgseaParam(expr = as.matrix(dat_bulk),
                            geneSets = mrg_list)
mrg_scores <- gsva(ssgsea_param, verbose = FALSE)

# 整理成数据框并合并分组信息
mrg_df <- data.frame(sample = colnames(mrg_scores),
                     MRGs_score = as.numeric(mrg_scores[1, ]))
mrg_df <- merge(mrg_df, group_info, by = "sample")

# 保存评分结果
write.csv(mrg_df, "MRGs_ssGSEA_scores.csv", row.names = FALSE)

# ============================ 3. UC组和对照组间MRGs评分差异检验 ============================
wilcox_res <- wilcox.test(MRGs_score ~ group, data = mrg_df)
cat("Wilcoxon秩和检验 P值：", wilcox_res$p.value, "/n")
cat("UC组评分中位数：", median(mrg_df$MRGs_score[mrg_df$group == "UC"]), "/n")
cat("Normal组评分中位数：", median(mrg_df$MRGs_score[mrg_df$group == "Normal"]), "/n")

# ============================ 4. 绘制MRGs评分差异箱线图（Figure 2A） ============================
p2A <- ggplot(mrg_df, aes(x = group, y = MRGs_score, fill = group)) +
  stat_boxplot(geom = "errorbar", width = 0.4) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = c("Normal" = "#4682B4", "UC" = "#CD3700")) +
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     comparisons = list(c("Normal", "UC"))) +
  labs(x = "", y = "MRGs ssGSEA score",
       title = "Figure 2A: MRGs score (UC vs Normal)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none",
        panel.grid = element_blank())

# 保存
ggsave("Figure2A_MRGs_score_boxplot.pdf", p2A, width = 6, height = 5)
ggsave("Figure2A_MRGs_score_boxplot.png", p2A, width = 6, height = 5, dpi = 600)
