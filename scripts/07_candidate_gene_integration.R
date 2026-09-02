##################################################
# Script: 07_candidate_gene_integration.R
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

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
cat("工作目录：", getwd(), "/n")

library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(GSVA)
library(msigdbr)

# ============================ 读取数据 ============================
scData <- readRDS(file.path("..", "01singleCell", "05_scData_annotated.rds"))
cat("原始数据：", ncol(scData), "个细胞，", nrow(scData), "个基因/n")



# ============================ 按需基因映射 ============================
# 先获取分析要用到的通路基因集
msig <- msigdbr(species = "Homo sapiens", collection = "C2",
                subcollection = "CP:REACTOME")
pathway_genes <- unique(msig$gene_symbol)

my_genes <- rownames(scData)
shared_genes <- intersect(my_genes, pathway_genes)
match_rate <- length(shared_genes) / length(pathway_genes)
cat(sprintf("基因集匹配率：%.2f (%d/%d)/n", match_rate, length(shared_genes), length(pathway_genes)))

# 设定阈值：若匹配率低于 0.5，则尝试对未匹配的基因做别名映射
if (match_rate < 0.5) {
  
  unfound_genes <- setdiff(my_genes, pathway_genes)
  # 查询 ALIAS -> SYMBOL
  alias_map <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db,
                                                      keys = unfound_genes,
                                                      columns = "SYMBOL",
                                                      keytype = "ALIAS"))
  alias_map <- alias_map %>% filter(!is.na(SYMBOL)) %>% distinct(ALIAS, .keep_all = TRUE)
  
  # 构建映射字典（仅替换能查到的基因）
  if (nrow(alias_map) > 0) {
    id_to_symbol <- setNames(alias_map$SYMBOL, alias_map$ALIAS)
    mappable <- intersect(my_genes, names(id_to_symbol))
    cat(sprintf("通过别名映射纠正了 %d 个基因/n", length(mappable)))
    # 替换表达矩阵的行名
    rownames(scData)[rownames(scData) %in% mappable] <- id_to_symbol[rownames(scData)[rownames(scData) %in% mappable]]
  } else {
    cat("没有找到可映射的别名基因。/n")
  }
  # 映射后重新计算匹配率（可选）
  new_shared <- intersect(rownames(scData), pathway_genes)
  cat(sprintf("映射后匹配率：%.2f (%d/%d)/n", length(new_shared)/length(pathway_genes), length(new_shared), length(pathway_genes)))
} else {
  cat("基因名与通路基因集匹配良好，无需映射。/n")
}

# 处理重复基因：同名基因保留平均表达量最高者
mean_expr <- Matrix::rowMeans(GetAssayData(scData, layer = "counts"))
scData <- scData[order(mean_expr, decreasing = TRUE), ]
scData <- scData[!duplicated(rownames(scData)), ]

cat("清洗后基因数:", nrow(scData), "/n")


# 降采样

scData$cell_type_group <- scData$external_annotation
Idents(scData) <- "cell_type_group"

table(scData$cell_type_group)
#B/Plasma Epithelial    Myeloid     NK/ILC    Stromal     T cell Unassigned 
#9121       1947        786         73        423       3981       1610 

#NK/ILC 数量太少 移除 Unassigned 无意义移除
keep_types <- setdiff(names(table(scData$cell_type_group)),
                      c("NK/ILC", "Unassigned"))

scData <- scData[, scData$cell_type_group %in% keep_types]
scData$cell_type_group <- factor(scData$cell_type_group)   # 转成因子，便于后续使用

Idents(scData) <- "cell_type_group"   # 确保 WhichCells 能按此分组

set.seed(42)
cells_per_type <- 2000


keep_cells <- unlist(lapply(unique(scData$cell_type_group), function(ct) {
  cells <- WhichCells(scData, idents = ct)
  if (length(cells) > cells_per_type) sample(cells, cells_per_type) else cells
}))

scData_down <- scData[, keep_cells]

cat("降采样后细胞数量:/n")
print(table(scData_down$cell_type_group))

#B/Plasma Epithelial    Myeloid    Stromal     T cell 
#2000       1947        786        423       2000 

# ============================ ssGSEA 通路富集（GSVA + Reactome 基因集） ============================
avg_expr <- as.matrix(AverageExpression(scData_down, layer = "data",
                                        group.by = "cell_type_group")$RNA)
cat("平均表达矩阵:", nrow(avg_expr), "基因 x", ncol(avg_expr), "细胞类型/n")

# 利用前面已加载的 msig 对象构建基因集列表
gene_sets <- split(msig$gene_symbol, msig$gs_name)
cat("Reactome 通路基因集数:", length(gene_sets), "/n")



pathway_mat <- gsva(ssgseaParam(avg_expr, gene_sets, minSize = 5, normalize = TRUE),
                    verbose = TRUE)


pathway_df <- data.frame(Name = rownames(pathway_mat), pathway_mat,
                         check.names = FALSE)



# ============================ 计算细胞类型间相对活性 ============================
cell_type_cols <- colnames(pathway_df)[-1]
cat("通路数:", nrow(pathway_df), " 细胞类型:", paste(cell_type_cols, collapse = ", "), "/n")

pathway_df$mean_activity <- rowMeans(pathway_df[, cell_type_cols], na.rm = TRUE)
for (col in cell_type_cols) {
  pathway_df[[paste0("diff_", col)]] <- pathway_df[[col]] - pathway_df$mean_activity
}

write.csv(pathway_df, "Reactome_ssGSEA_pathways.csv", row.names = FALSE)



# ============================ 热图：相对活性（高于均值为红，低于为蓝） ============================
top_n <- 10
diff_cols <- paste0("diff_", cell_type_cols)

unique_pathways <- unique(unlist(lapply(diff_cols, function(dcol) {
  pathway_df %>% arrange(desc(!!sym(dcol))) %>% head(top_n) %>% pull(Name)
})))

plot_data <- pathway_df %>%
  dplyr::filter(Name %in% unique_pathways) %>%
  dplyr::select(Name, all_of(diff_cols)) %>%
  pivot_longer(cols = all_of(diff_cols),
               names_to = "Cell_Type",
               values_to = "Activity_Score")
plot_data$Cell_Type <- gsub("^diff_", "", plot_data$Cell_Type)
plot_data$Name <- factor(plot_data$Name, levels = unique_pathways)

fig_heatmap <- ggplot(plot_data, aes(x = Cell_Type, y = Name, fill = Activity_Score)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, name = "Relative/nActivity"
  ) +
  labs(title = "Relative Pathway Activity Across Cell Types", x = NULL, y = NULL) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8, face = "italic"),
    axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.title = element_text(size = 14, face = "bold"),
    plot.title.position = "plot"
  )

print(fig_heatmap)
ggsave("Figure_Pathway_Heatmap.pdf", fig_heatmap, width = 10, height = 12)
ggsave("Figure_Pathway_Heatmap.png", fig_heatmap, width = 10, height = 12, dpi = 300)

