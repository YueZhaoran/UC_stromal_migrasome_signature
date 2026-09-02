##################################################
# Script: 05_Scissor_analysis.R
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

library(Seurat)
library(GSVA)
library(Scissor)
library(data.table)
library(dplyr)
library(tibble)
library(readxl)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggrepel)
library(ggsci)
set.seed(798) 

# 读取 Bulk 数据
# 读取表达矩阵
dat_bulk <- read.csv("../00bulkdata/dat.GSE87466.csv", row.names = 1, check.names = FALSE)
# 读取分组信息
group_info <- read.csv("../00bulkdata/group.GSE87466.csv", stringsAsFactors = FALSE)
# 读取临床/表型信息
pd_info <- read.csv("../00bulkdata/pdGSE87466.csv", stringsAsFactors = FALSE)


table(group_info$group)

#RGs 基因集
datasets <- read_excel("../migrasome-genes.xls")
mrg_list <- list(MRGs = datasets$genes) 




#ssGSEA 计算每个样本的 MRGs 评分
#参数对象
ssgsea_param <- ssgseaParam(expr = as.matrix(dat_bulk), 
                            geneSets = mrg_list)

#运行gsva
mrg_scores <- gsva(ssgsea_param, verbose = FALSE)


# mrg_scores 是一个 1 行 x 样本数 的矩阵，我们需要转置并转为数据框
mrg_df <- as.data.frame(t(mrg_scores))
colnames(mrg_df) <- "MRGs_score"



#只提取 UC 患者的样本，划分高/低表型
uc_samples <- group_info$sample[group_info$group == "UC"]
mrg_uc <- mrg_df[uc_samples, , drop = FALSE]

# 取中位数分割高/低
median_score <- median(mrg_uc$MRGs_score)
mrg_uc$phenotype <- ifelse(mrg_uc$MRGs_score > median_score, 1, 0) 
# 1代表高MRGs，0代表低MRGs

# 整理成 Scissor 需要的格式 (一个包含 0/1 的命名向量)
phenotype <- mrg_uc$phenotype
names(phenotype) <- rownames(mrg_uc)


# 打印一下看看高/低各有多少样本
table(phenotype)


#Scissor 核心代码
pheno_names <- names(phenotype)
bulk_names <- colnames(dat_bulk)


if(!all(pheno_names == bulk_names)){
  #对齐表达矩阵
  dat_bulk <- dat_bulk[, pheno_names, drop = FALSE]
}

print(paste0("Bulk样本数: ", ncol(dat_bulk), "； 表型数: ", length(phenotype)))



#读取单细胞数据，提取表达矩阵
scData <- readRDS("../01singleCell/05_scData_annotated.rds")


#运行 Scissor 算法 原版Scissor包不适配V5 这里用了修改过源码的版本 增加了Scissor_V5函数
Scissor_res <- Scissor_V5(bulk_dataset = dat_bulk, 
                         sc_dataset = scData,
                         phenotype = phenotype,
                         tag = c("Low_MRGs", "High_MRGs"),
                         alpha = 0.05,
                         family = "binomial",
                         Save_file = "Scissor_output.RData") 

save(Scissor_res, file = "Scissor_res_result.RData")
# 提取结果
scissor_pos <- Scissor_res$Scissor_pos
scissor_neg <- Scissor_res$Scissor_neg

cat("Scissor+:", length(scissor_pos), "/n")
cat("Scissor-:", length(scissor_neg), "/n")

# 创建新列，默认所有细胞都是 Background
scData$Scissor_label <- "Background"

# 将 Scissor+ 细胞的 ID 对应位置打上标签
scData$Scissor_label[scissor_pos] <- "Scissor+"

# 将 Scissor- 细胞的 ID 对应位置打上标签
scData$Scissor_label[scissor_neg] <- "Scissor-"

# 查看三类细胞的分布情况（完美对接你方案里的 Figure 2B 和 2C）
table(scData$Scissor_label)


Idents(scData) <- "Scissor_label"


# 差异分析：Scissor+ 细胞 vs 所有 Background 细胞
markers_scissor <- FindMarkers(scData, 
                               ident.1 = "Scissor+", 
                               ident.2 = "Background", 
                               min.pct = 0.1,     # 基因在至少 10% 的细胞中表达
                               logfc.threshold = 0.5)  # 最小 Log2FC 阈值


# 方案要求：|log2FC| > 1 且 Padj < 0.05
cell_spec_MRGs <- markers_scissor %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 1) %>%
  rownames_to_column("gene")

#细胞特异性迁移体相关基因
nrow(cell_spec_MRGs)

write.csv(cell_spec_MRGs, "cell_specific_MRGs.csv", row.names = FALSE)


# 准备数据
scData$celltype <- scData$external_annotation

# 堆叠柱状图：各细胞类型的Scissor比例 
scissor_summary <- scData@meta.data %>%
  filter(!is.na(external_annotation)) %>%
  group_by(external_annotation, Scissor_label) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(external_annotation) %>%
  mutate(percent = count / sum(count) * 100)

p_stack <- ggplot(scissor_summary, aes(x = external_annotation, y = percent, fill = Scissor_label)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Scissor+" = "#E41A1C",
                               "Scissor-" = "#377EB8",
                               "Background" = "grey85")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Percentage (%)") +
  xlab("Cell Type") +
  ggtitle("Scissor Cell Distribution by Cell Type")

print(p_stack)
ggsave("Figure_Scissor_StackedBar_byCellType.pdf", p_stack, width = 10, height = 6, dpi = 300)
ggsave("Figure_Scissor_StackedBar_byCellType.png", p_stack, width = 10, height = 6, dpi = 300)


# Scissor 标签在 UMAP 上的分布图

p2C <- DimPlot(scData, group.by = "Scissor_label") + 
  scale_color_manual(values = c("Scissor+" = "#E64B35", "Scissor-" = "#377EB8", "Background" = "#999999")) +
  labs(title = "Figure 2C: Scissor label UMAP") +
  theme(plot.title = element_text(hjust = 0.5))

# 保存
ggsave("Figure2C_Scissor_UMAP.pdf", p2C, width = 7, height = 6)
ggsave("Figure2C_Scissor_UMAP.png", p2C, width = 7, height = 6)


combined_plot <- p_stack + p2C + plot_layout(ncol = 2, widths = c(1, 1))
ggsave("Figure2B_2C_Combined.pdf", combined_plot, width = 14, height = 6)
ggsave("Figure2B_2C_Combined.png", combined_plot, width = 14, height = 6)






# 差异分析结果
# 将行名（基因名）变成一列
df_volcano <- markers_scissor %>%
  rownames_to_column("gene") %>%
  mutate(
    # 定义阈值
    logFC = avg_log2FC,
    Pval = p_val_adj,
    
    # p 值截断：0 值替换为机器最小正数，防止 -log10(0) = Inf 撑大 Y 轴
    # 同时设置上限 1e-300，超出该值的点统一按 1e-300 绘制
    Pval_cap = pmax(Pval, 1e-300),
    Pval_cap = ifelse(Pval_cap == 0, .Machine$double.xmin, Pval_cap),
    
    # 根据阈值打标签（用于区分颜色）
    group = case_when(
      Pval < 0.05 & logFC > 1 ~ "Up",
      Pval < 0.05 & logFC < -1 ~ "Down",
      TRUE ~ "Not Sig"
    )
  )

# 计算每个方向最显著的 Top 5 基因，用于图中标注名字
top_up <- df_volcano %>% filter(group == "Up") %>% arrange(Pval) %>% head(5)
top_down <- df_volcano %>% filter(group == "Down") %>% arrange(Pval) %>% head(5)
top_genes <- bind_rows(top_up, top_down)


# 绘制火山图（使用截断后的 Pval_cap，并限制 Y 轴上限避免被极值撑散）
p2D <- ggplot(df_volcano, aes(x = logFC, y = -log10(Pval_cap), color = group)) +
  geom_point(alpha = 0.8, size = 1.5) + # 散点
  scale_color_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8", "Not Sig" = "#999999")) +
  
  # 添加垂直和水平阈值辅助线
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
  
  # 标注 Top 基因的名称
  geom_text_repel(data = top_genes, aes(label = gene), 
                  size = 3.5, box.padding = 0.5, max.overlaps = Inf) +
  
  # 限制 Y 轴显示范围，防止少数极显著基因把图拉得太长
  coord_cartesian(ylim = c(0, 330)) +
  
  # 美化主题
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = expression(log[2]("Fold Change")), 
       y = expression(-log[10]("Adjusted P-value")), 
       title = "Figure 2D: Cell-specific MRGs (Scissor+ vs Background)") +
  theme(legend.position = "right", 
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))


ggsave("Figure2D_Volcano_cell_specific_MRGs.pdf", p2D, width = 8, height = 6)
ggsave("Figure2D_Volcano_cell_specific_MRGs.png", p2D, width = 8, height = 6)















