##################################################
# Script: 12_microenvironment_analysis.R
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



# 加载必要的包
library(xCell)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(psych)

# ============================ 1. 读取数据 ============================
# 读取训练集表达矩阵（行=基因，列=样本）
train_data <- read.csv("../00bulkdata/dat.GSE87466.csv", row.names = 1, check.names = FALSE)

# 读取分组信息
train_group <- read.csv("../00bulkdata/group.GSE87466.csv", stringsAsFactors = FALSE)
colnames(train_group) <- c("sample", "group")

# 读取关键基因列表（从机器学习的交集文件中读取）
hub_gene <- read.csv("../05DEGs/intersection_gene.csv", stringsAsFactors = FALSE)
key_genes <- hub_gene$symbol

# 确保样本顺序一致
train_data <- train_data[, train_group$sample]

# ============================ 2. xCell 免疫细胞浸润打分 ============================
cat("正在使用 xCell 计算 64 种免疫细胞评分，请稍候.../n")
# 运行 xCell（注意：xCellAnalysis 的输入要求行为基因，列为样本）
xcell_res <- xCellAnalysis(train_data, rnaseq = FALSE)

# xCell 返回的结果是 行=细胞类型，列=样本，需要转置为 行=样本，列=细胞类型
xcell_score <- as.data.frame(t(xcell_res))
xcell_score$sample <- rownames(xcell_score)

# 合并分组信息
xcell_score <- merge(xcell_score, train_group, by = "sample")
rownames(xcell_score) <- xcell_score$sample
xcell_score <- xcell_score %>% dplyr::select(-sample)

# 保存原始得分
write.csv(xcell_score, "01.xCell_scores.csv")

# ============================ 3. 绘制免疫细胞浸润热图（图 6B） ============================
# 按分组排序样本，让热图分组更清晰
xcell_score_plot <- xcell_score[order(xcell_score$group), ]
group_order <- xcell_score_plot$group
xcell_score_plot <- xcell_score_plot %>% dplyr::select(-group)

# 准备注释信息
annotation_col <- data.frame(Group = group_order)
rownames(annotation_col) <- rownames(xcell_score_plot)

# 设置热图颜色（行标准化后使用以0为中心的发散色带）
ann_colors <- list(Group = c(Normal = "#45a9b8", UC = "#f76a56"))
color_key <- c("#2166AC", "#67A9CF", "white", "#EF8A62", "#B2182B")

# 绘制并保存热图（按行做z-score标准化，突出各细胞类型在样本间的相对高低）
p_heatmap <- pheatmap(
  t(xcell_score_plot),  # 转置，使行为细胞类型，列为样本
  scale = "row",        # 按行（细胞类型）z-score标准化，解决原始得分整体偏低导致全图偏蓝的问题
  color = colorRampPalette(color_key)(50),
  border_color = "darkgrey",
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  clustering_method = "ward.D2",
  show_rownames = TRUE,
  show_colnames = FALSE,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  main = "Immune cell infiltration (row z-score)",
  filename = "01.Immune_Cell_Heatmap.pdf",
  width = 8,
  height = 8
)
png("01.Immune_Cell_Heatmap.png", width = 8, height = 8, units = "in", res = 600)
print(p_heatmap)
dev.off()






# ============================ 4. 差异分析（筛选差异免疫细胞） ============================
# 把长数据整理成画图用的格式
plot_data <- xcell_score %>%
  pivot_longer(
    cols = -group,
    names_to = "ImmuneCell",
    values_to = "Score"
  )

# 对每种免疫细胞做 Wilcoxon 秩和检验
wilcox_res <- plot_data %>%
  group_by(ImmuneCell) %>%
  rstatix::wilcox_test(Score ~ group) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance("p")

# 筛选出 P < 0.05 的差异显著细胞
DE_cells <- wilcox_res %>% filter(p < 0.05)

# 保存差异分析结果
write.csv(wilcox_res, "02.wilcox_res.csv", row.names = FALSE)
write.csv(DE_cells, "03.DE_cells.csv", row.names = FALSE)

# ============================ 5. 绘制差异免疫细胞小提琴图（图 6C） ============================
# 仅保留差异显著的细胞用于绘图
plot_data_diff <- plot_data %>%
  filter(ImmuneCell %in% DE_cells$ImmuneCell)

# 小提琴图绘制（这里使用箱线图+抖动点，更加直观）
p_violin <- ggplot(plot_data_diff, aes(x = ImmuneCell, y = Score, fill = group)) +
  stat_boxplot(geom = "errorbar", width = 0.5, position = position_dodge(0.9)) +
  geom_boxplot(width = 0.5, position = position_dodge(0.9), outlier.shape = NA) +
  scale_fill_manual(values = c("Normal" = "#4682B4", "UC" = "#CD3700"), name = "Group") +
  labs(title = "Differences in Immune Cell Infiltration", x = "", y = "xCell Score") +
  stat_compare_means(aes(group = group), label = "p.signif", method = "wilcox.test") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 12),
    legend.position = "top",
    panel.grid = element_blank()
  )

ggsave("02.Immune_Cell_Diff_Boxplot.pdf", p_violin, width = 12, height = 8)
ggsave("02.Immune_Cell_Diff_Boxplot.png", p_violin, width = 12, height = 8, dpi = 600)

# ============================ 6. 关键基因与差异免疫细胞的相关性（图 6D） ============================
if (!"sample" %in% colnames(xcell_score)) {
  xcell_score$sample <- rownames(xcell_score)
}

# 提取关键基因表达矩阵
expr_genes <- t(train_data[key_genes, ]) %>% as.data.frame()
colnames(expr_genes) <- key_genes
expr_genes$sample <- rownames(expr_genes)

# 提取差异免疫细胞的评分（现在 sample 列已经存在，不会报错了）
cell_scores <- xcell_score[, c("sample", DE_cells$ImmuneCell)]

# 合并表达量和细胞评分
cor_merge <- merge(expr_genes, cell_scores, by = "sample")

# 定义基因列和免疫细胞列的名称
gene_cells <- key_genes
immune_cells <- DE_cells$ImmuneCell

# 【关键修正】强制将免疫细胞列全部转为数值型（防止字符型报错）
cor_merge[, immune_cells] <- lapply(cor_merge[, immune_cells], as.numeric)

# ============================ 2. 计算 Spearman 相关性 ============================
cor_res <- data.frame()
for (gene in gene_cells) {
  for (cell in immune_cells) {
    test <- cor.test(cor_merge[[gene]], cor_merge[[cell]], method = "spearman")
    cor_res <- rbind(cor_res, data.frame(
      Gene = gene,
      Cell = cell,
      Cor = test$estimate,
      P = test$p.value
    ))
  }
}

# 保存相关性表格
write.csv(cor_res, "04.Gene_Cell_Correlation.csv", row.names = FALSE)

# ============================ 3. 添加显著性标记（星号） ============================
cor_res <- cor_res %>%
  mutate(text = case_when(
    P <= 0.001 ~ "***",
    P <= 0.01 ~ "**",
    P <= 0.05 ~ "*",
    TRUE ~ ""
  ))

# ============================ 4. 绘制相关性热图（图 6D） ============================
p_cor <- ggplot(cor_res, aes(x = Gene, y = Cell, fill = Cor)) +
  geom_tile(color = "grey", size = 0.5) +
  scale_fill_gradient2(low = "#035397", mid = "white", high = "#F32424", midpoint = 0) +
  geom_text(aes(label = text), color = "black", size = 4) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 10),
    legend.title = element_text(size = 12, face = "bold")
  ) +
  labs(fill = "Spearman Correlation")

# 保存图片
ggsave("04.correlation_biomarker.pdf", p_cor, width = 8, height = 6)
ggsave("04.correlation_biomarker.png", p_cor, width = 8, height = 6, dpi = 600)
