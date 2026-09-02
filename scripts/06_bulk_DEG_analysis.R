##################################################
# Script: 06_bulk_DEG_analysis.R
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


library(limma)
library(ggplot2)
library(dplyr)
library(pheatmap)


# 读取 Bulk 数据
dat_bulk <- read.csv("../00bulkdata/dat.GSE87466.csv", row.names = 1, check.names = FALSE)
# 读取分组信息
group_info <- read.csv("../00bulkdata/group.GSE87466.csv", stringsAsFactors = FALSE)



# limma 差异表达分析
group_vec <- factor(group_info$group, levels = c("Normal", "UC"))
design <- model.matrix(~ 0 + group_vec)
colnames(design) <- levels(group_vec)

fit <- lmFit(dat_bulk, design)
cont.matrix <- makeContrasts(UC - Normal, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

all_DEGs <- topTable(fit2, coef = 1, n = Inf, sort.by = "P")

# 筛选差异基因并保存
bulk_DEGs <- all_DEGs %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1) %>%
  rownames_to_column("gene")

write.csv(bulk_DEGs, "bulk_DEGs_GSE87466.csv", row.names = FALSE)

# 绘制火山图
df_vol <- all_DEGs %>%
  mutate(group = case_when(
    logFC > 1 & adj.P.Val < 0.05 ~ "Up",
    logFC < -1 & adj.P.Val < 0.05 ~ "Down",
    TRUE ~ "Not Sig"
  ))

p_vol <- ggplot(df_vol, aes(x = logFC, y = -log10(adj.P.Val), color = group)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c("Up" = "#D55E00", "Down" = "#0072B2", "Not Sig" = "#999999")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
  theme_bw() +
  labs(x = expression(log[2]("Fold Change")), y = expression(-log[10]("Adj. P-value")),
       title = "Bulk DEGs Volcano Plot")


ggsave("bulk_DEGs_Volcano.pdf", p_vol, width = 7, height = 6)
ggsave("bulk_DEGs_Volcano.png", p_vol, width = 7, height = 6)



# 绘制热图（选取前50个显著差异基因）
top_genes <- head(rownames(all_DEGs[order(all_DEGs$P.Value), ]), 50)
mat <- dat_bulk[top_genes, ]

annotation_col <- data.frame(Group = group_info$group)
rownames(annotation_col) <- colnames(mat)

pheatmap(mat, scale = "row", show_rownames = FALSE, show_colnames = FALSE,
         annotation_col = annotation_col, cluster_cols = FALSE,main = "Top 50 Differentially Expressed Genes",
         filename = "bulk_DEGs_Heatmap.pdf", width = 8, height = 6)

pheatmap(mat, scale = "row", show_rownames = FALSE, show_colnames = FALSE,
         annotation_col = annotation_col, cluster_cols = FALSE,main = "Top 50 Differentially Expressed Genes",
         filename = "bulk_DEGs_Heatmap.png", width = 8, height = 6)



