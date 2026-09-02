##################################################
# Script: 02_single_cell_processing.R
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
# 创建输出子文件夹
dirs <- c("QC", "cell_cycle", "PCA", "singleR", "CellChat",
          "AUCell", "feature_plots", "Milo", "Augur")
for (d in dirs) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# 加载全部所需包
library(Seurat)
library(Matrix)
library(data.table)
library(stringr)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(celldex)
library(SingleR)
library(harmony)
library(clustree)
library(CellChat)
library(scCustomize)
library(paletteer)
library(viridis)
library(AUCell)
library(miloR)
library(Augur)
library(ggpubr)
library(BiocParallel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggsci)

sample_info <- data.frame(
  sample_id = c(
    "GSM6614348", "GSM6614349", "GSM6614350", "GSM6614351", "GSM6614352", "GSM6614353",
    "GSM6614354", "GSM6614355", "GSM6614356", "GSM6614357", "GSM6614358", "GSM6614359"
  ),
  sample_name = c(
    "HC-1", "HC-2", "HC-3", "HC-4", "HC-5", "HC-6",
    "UC-1", "UC-2", "UC-3", "UC-4", "UC-5", "UC-6"
  ),
  group = c(
    rep("HC", 6),
    rep("UC", 6)
  ),
  stringsAsFactors = FALSE
)


data_path <- "data/"


read_one_sample_10x <- function(sample_id, sample_name, grade, group, path) {
  
  # ---- 查找三个文件 ----
  barcode_file <- list.files(
    path = path,
    pattern = paste0(sample_id, ".*barcodes//.tsv//.gz$"),
    full.names = TRUE
  )
  features_file <- list.files(
    path = path,
    pattern = paste0(sample_id, ".*features//.tsv//.gz$"),
    full.names = TRUE
  )
  matrix_file <- list.files(
    path = path,
    pattern = paste0(sample_id, ".*matrix//.mtx//.gz$"),
    full.names = TRUE
  )
  
  # 检查文件完整性
  if (length(barcode_file) == 0) stop(paste("找不到 barcodes 文件:", sample_id))
  if (length(features_file) == 0) stop(paste("找不到 features 文件:", sample_id))
  if (length(matrix_file) == 0) stop(paste("找不到 matrix 文件:", sample_id))
  
  barcode_file <- barcode_file[1]
  features_file <- features_file[1]
  matrix_file <- matrix_file[1]
  
  message("读取样本: ", sample_name)
  
  # ---- 读取数据 ----
  barcodes <- read.table(barcode_file, header = FALSE, stringsAsFactors = FALSE)[, 1]
  features <- read.table(features_file, header = FALSE, stringsAsFactors = FALSE,sep = "/t")
  mat <- readMM(matrix_file)   # 稀疏矩阵（行 = 基因，列 = 细胞）
  
  cat("列数",ncol(features))
  if (ncol(features) >= 2) {
    # 直接提取第二列，并用 make.unique 确保绝对不会重复
    gene_symbols <- make.unique(features[, 2])  
  } else {
    gene_symbols <- features[, 1]
  }
  
  # 赋予行列名
  rownames(mat) <- gene_symbols
  colnames(mat) <- barcodes
  

  # ---- 创建 Seurat 对象 ----
  seu <- CreateSeuratObject(
    counts = mat,
    project = sample_name,
    min.cells = 3,
    min.features = 200
  )
  

  seu$sample_id <- sample_id
  seu$sample_name <- sample_name
  seu$grade <- grade
  seu$group <- group
  seu$orig.ident <- sample_name
  
  
  return(seu)
}


seurat_list <- list()

for (i in seq_len(nrow(sample_info))) {
  seurat_list[[sample_info$sample_name[i]]] <- read_one_sample_10x(
    sample_id = sample_info$sample_id[i],
    sample_name = sample_info$sample_name[i],
    grade = sample_info$grade[i],
    group = sample_info$group[i],
    path = data_path
  )
}


scData <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list),
  project = "GSE214695"
)


# ========== 查看合并结果 ==========
scData
table(scData$sample_name)
table(scData$group)
head(rownames(scData))

saveRDS(scData, "00_scData_raw.rds")






# ============================ 质量控制模块 ============================
# 目的：过滤低质量细胞，保留用于后续分析的高质量数据
setwd(file.path(project_dir, "QC"))
scData <- readRDS(file.path(project_dir, "00_scData_raw.rds"))

# 计算每个细胞的线粒体基因比例和红细胞基因比例
scData[["percent.mt"]] <- PercentageFeatureSet(scData, pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes, rownames(scData))
scData[["percent.HB"]] <- PercentageFeatureSet(scData, features = HB.genes)

# 质控前小提琴图：展示过滤前的指标分布
plot.features <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB")
theme.set2 <- theme(axis.title.x = element_blank())
plots <- lapply(plot.features, function(f) {
  VlnPlot(scData, group.by = "group", pt.size = 0, features = f) + theme.set2 + NoLegend()
})

p <- wrap_plots(plots, nrow = 2)
ggsave("QC_before.pdf", plot = p, width = 14, height = 8)
ggsave("QC_before.png", plot = p, width = 14, height = 8)


# 设置质控阈值并过滤细胞
minGene <- 300; maxGene <- 10000
minUMI  <- 600; pctMT <- 10; pctHB <- 1
scData <- subset(scData,
               subset = nFeature_RNA > minGene & nFeature_RNA < maxGene &
                 nCount_RNA > minUMI & percent.mt < pctMT & percent.HB < pctHB)

# 质控后小提琴图：展示过滤后的指标分布
plots <- lapply(plot.features, function(f) {
  VlnPlot(scData, group.by = "group", pt.size = 0, features = f) + theme.set2 + NoLegend()
})

p <- wrap_plots(plots, nrow = 2)

ggsave("QC_after.pdf", plot = p, width = 14, height = 8)
ggsave("QC_after.png", plot = p, width = 14, height = 8)

saveRDS(scData, file.path(project_dir, "01_scData_qc.rds"))
file.remove(file.path(project_dir, "00_scData_raw.rds"))









# ============================ 细胞周期评分模块 ============================
# 目的：评估细胞周期状态，便于后续回归掉周期效应
setwd(file.path(project_dir, "cell_cycle"))
scData <- readRDS(file.path(project_dir, "01_scData_qc.rds"))

# 标准化并寻找高变基因，以便正确计算周期分数
scData <- NormalizeData(scData) %>%
  FindVariableFeatures() %>%
  ScaleData()

# 获取 G2/M 期和 S 期基因（Seurat 内置）
g2m_genes <- CaseMatch(cc.genes$g2m.genes, rownames(scData))
s_genes   <- CaseMatch(cc.genes$s.genes,   rownames(scData))

# 检查匹配数量
cat("匹配到的 S 期基因数:", length(s_genes), "/n")
cat("匹配到的 G2M 期基因数:", length(g2m_genes), "/n")

# 为每个细胞计算细胞周期评分
scData <- CellCycleScoring(scData, g2m.features = g2m_genes, s.features = s_genes)


s_genes_not_present <- setdiff(s_genes, rownames(scData))
g2m_genes_not_present <- setdiff(g2m_genes, rownames(scData))
cat("S 期未找到的基因数：", length(s_genes_not_present), "/n")
cat("G2M 期未找到的基因数：", length(g2m_genes_not_present), "/n")

# 生成 t-SNE 图并按细胞周期阶段着色
scData <- RunPCA(scData) %>% RunTSNE()

p <- DimPlot(scData, group.by = "Phase", reduction = "tsne")

ggsave("cellcycle_tsne.pdf", plot = p, width = 8, height = 6)
ggsave("cellcycle_tsne.png", plot = p, width = 8, height = 6)

saveRDS(scData, file.path(project_dir, "02_scData_cycle.rds"))
file.remove(file.path(project_dir, "01_scData_qc.rds"))




# ============================ 标准化与高变基因 ============================
# 目的：准备用于降维和聚类的标准化数据
scData <- readRDS(file.path(project_dir, "02_scData_cycle.rds"))
DefaultAssay(scData) <- "RNA"

# 对数标准化
scData <- NormalizeData(scData, normalization.method = "LogNormalize", scale.factor = 10000)
# 识别高变基因（默认 2000 个）
scData <- FindVariableFeatures(scData, selection.method = "vst", nfeatures = 2000)
# 缩放数据并回归细胞周期分数
scData <- ScaleData(scData, vars.to.regress = c("S.Score", "G2M.Score"))

# ============================ PCA、Harmony 与聚类 ============================
# 目的：降维、校正批次效应、聚类并可视化
setwd(file.path(project_dir, "PCA"))

Idents(scData) <- "group"

# 绘制高变基因散点图

top10 <- head(VariableFeatures(scData), 10)
p <- VariableFeaturePlot(scData)
pdf("variable_features.pdf", width = 7, height = 6); print(p); dev.off()
ggsave("variable_features.png", plot = p, width = 7, height = 6)

p <- LabelPoints(plot = VariableFeaturePlot(scData), points = top10, repel = TRUE)
pdf("variable_features_labeled.pdf", width = 7, height = 6); print(p); dev.off()
ggsave("variable_features_labeled.png", plot = p, width = 7, height = 6)

scData <- RunPCA(scData, verbose = FALSE)
p <- DimPlot(scData, reduction = "pca", group.by = "group")
pdf("pca_by_group.pdf", width = 7, height = 6); print(p); dev.off()
ggsave("pca_by_group.png", plot = p, width = 7, height = 6)

p <- VizDimLoadings(scData, dims = 1:4, reduction = "pca", nfeatures = 20)
pdf("pca_loadings.pdf", width = 10, height = 9); print(p); dev.off()
ggsave("pca_loadings.png", plot = p, width = 10, height = 9)

png("pca_heatmap.png", width = 10, height = 9, units = "in", res = 300)
DimHeatmap(scData, dims = 1:4, cells = 500, balanced = TRUE, nfeatures = 30, ncol = 2)
dev.off()

# 保存为 PDF
pdf("pca_heatmap.pdf", width = 10, height = 9)
DimHeatmap(scData, dims = 1:4, cells = 500, balanced = TRUE, nfeatures = 30, ncol = 2)
dev.off()

p <- ElbowPlot(scData, ndims = 50)
pdf("elbow_plot.pdf", width = 7, height = 6); print(p); dev.off()
ggsave("elbow_plot.png", plot = p, width = 7, height = 6)




pcs <- 1:30

pdf("harmony_convergence.pdf", width = 8, height = 6) 
scData <- RunHarmony(scData, group.by.vars = "orig.ident", max_iter = 50,
                   plot_convergence = TRUE)
dev.off()
png("harmony_convergence.png", width = 8, height = 6, units = "in", res = 300)
scData <- RunHarmony(scData, group.by.vars = "orig.ident", max_iter = 50,
                     plot_convergence = TRUE)
dev.off()




seq_res <- seq(0.1, 2, by = 0.1)
scData <- FindNeighbors(scData, dims = pcs)
for (res in seq_res) {
  scData <- FindClusters(scData, resolution = res)
}
p <- clustree(scData, prefix = "RNA_snn_res.") + coord_flip()
ggsave("clustree.png", plot = p, width = 30, height = 14)        # 原有是PNG，现额外保存PDF
ggsave("clustree.pdf", plot = p, width = 30, height = 14)



scData <- FindNeighbors(scData, reduction = "harmony", dims = pcs) %>%
  FindClusters(resolution = 1.1) %>%
  RunUMAP(reduction = "harmony", dims = pcs) %>%
  RunTSNE(reduction = "harmony", dims = pcs)

p <- DimPlot(scData, reduction = "umap", label = TRUE)
pdf("umap_clusters.pdf", width = 7, height = 6); print(p); dev.off()
ggsave("umap_clusters.png", plot = p, width = 7, height = 6)


p <- DimPlot(scData, reduction = "umap", group.by = "group")
pdf("umap_group.pdf", width = 7, height = 6); print(p); dev.off()
ggsave("umap_group.png", plot = p, width = 7, height = 6)


p <- DimPlot(scData, reduction = "tsne", label = TRUE)
pdf("tsne_clusters.pdf", width = 7, height = 6); print(p); dev.off()
ggsave("tsne_clusters.png", plot = p, width = 7, height = 6)


p <- DimPlot(scData, reduction = "tsne", group.by = "group")
pdf("tsne_group.pdf", width = 7, height = 6); print(p); dev.off()
ggsave("tsne_group.png", plot = p, width = 7, height = 6)


saveRDS(scData, file.path(project_dir, "04_scData_clustered.rds"))











# ============================ 细胞注释（SingleR） ============================
#利用参考数据集自动注释细胞类型，并结合标记基因手动调整
setwd(file.path(project_dir, "singleR"))
#scData <- readRDS(file.path(project_dir, "04_scData_clustered.rds"))
DefaultAssay(scData) <- "RNA"


scData <- JoinLayers(scData)
# 使用 Human Primary Cell Atlas 作为参考进行自动注释
hpca.se <- HumanPrimaryCellAtlasData()
pred <- SingleR(test = GetAssayData(scData, layer = "data"), ref = hpca.se,
                labels = hpca.se$label.main)
scData$labels <- pred$labels

# 对比 SingleR 注释与原始聚类结果
pdf("singler_umap.pdf", width = 20, height = 6)
DimPlot(scData, group.by = c("seurat_clusters", "labels"), reduction = "umap")
dev.off()

png("singler_umap.png", width = 20, height = 6, units = "in", res = 300)
DimPlot(scData, group.by = c("seurat_clusters", "labels"), reduction = "umap")
dev.off()



#GSE214695_cell_annotation

#读取单细胞文件作者提供的细胞注释文件 进一步映射为大类
anno_df <- read.table("GSE214695_cell_annotation.csv", header = T, sep = ",", stringsAsFactors = FALSE)
colnames(anno_df) <- c("sample_barcode", "sample_id", "barcode", "annotation", "nanostring_ref")
head(anno_df)



# 定义映射表
mapping <- c(
  "CD4" = "T cell",
  "CD8" = "T cell",
  "DN" = "T cell",
  "gd IEL" = "T cell",
  "MAIT" = "T cell",
  "Tregs" = "T cell",
  "Ribhi T cells" = "T cell",
  "T cells CCL20" = "T cell",
  "MT T cells" = "T cell",
  "Cycling T cells" = "T cell",
  
  "B cell" = "B/Plasma",
  "Naive B cell" = "B/Plasma",
  "Memory B cell" = "B/Plasma",
  "GCB cell" = "B/Plasma",
  "PC immediate early response" = "B/Plasma",
  "PC IgA" = "B/Plasma",
  "PC IgA heat shock" = "B/Plasma",
  "PC IgA IgM" = "B/Plasma",
  "PC IgG" = "B/Plasma",
  
  "MO" = "Myeloid",
  "M1" = "Myeloid",
  "M2" = "Myeloid",
  "IDA macrophage" = "Myeloid",
  "Inflammatory monocytes" = "Myeloid",
  "DCs" = "Myeloid",
  "Mast" = "Myeloid",
  "Eosinophils" = "Myeloid",
  "N1" = "Myeloid",
  "N2" = "Myeloid",
  "N3" = "Myeloid",
  "Cycling myeloid" = "Myeloid",
  
  "Colonocytes" = "Epithelial",
  "Goblet" = "Epithelial",
  "Tuft cells" = "Epithelial",
  "Enteroendocrine" = "Epithelial",
  "BEST4 OTOP2" = "Epithelial",
  "Epithelium Ribhi" = "Epithelial",
  "Secretory progenitor" = "Epithelial",
  "S1" = "Epithelial",
  "S2a" = "Epithelial",
  "S2b" = "Epithelial",
  "S3" = "Epithelial",
  "Cycling TA" = "Epithelial",
  
  "Endothelium" = "Stromal",
  "Fibroblasts" = "Stromal",
  "Inflammatory fibroblasts" = "Stromal",
  "Myofibroblasts" = "Stromal",
  "FRCs" = "Stromal",
  "Pericytes" = "Stromal",
  "Glia" = "Stromal",
  
  "NK" = "NK/ILC",
  "ILC4" = "NK/ILC",
  
  "Cycling cells" = "Unassigned",
  "Unassigned" = "Unassigned",
  
  
  "PC  immediate early response" = "B/Plasma",
  "Naïve B cell" = "B/Plasma",
  "GC B cell" = "B/Plasma",
  "Paneth-like" = "Epithelial",
  "M0" = "Myeloid"
)

#  将大类标签添加到注释表的新列中
anno_df$major_group <- mapping[anno_df$nanostring_ref]

#  检查是否有未匹配上的类别
unmatched <- unique(anno_df$nanostring_ref[is.na(anno_df$major_group)])
if(length(unmatched) > 0) {
  message("以下类别未在映射表中找到，将被设为 Unassigned：")
  print(unmatched)
  anno_df$major_group[is.na(anno_df$major_group)] <- "Unassigned"
}


head(anno_df)




seurat_clean_barcodes <- sub("^[^_]*_", "", rownames(scData@meta.data))

print(head(seurat_clean_barcodes))

anno_map <- setNames(anno_df$major_group, anno_df$barcode)

# 执行匹配
matched_annotations <- anno_map[seurat_clean_barcodes]


cat("未匹配到的细胞数量:", sum(is.na(matched_annotations)), "/n")


matched_annotations[is.na(matched_annotations)] <- "Unassigned"

names(matched_annotations) <- rownames(scData@meta.data)

# 将注释写入 metadata 中，新建一列叫 "external_annotation"
scData$external_annotation <- matched_annotations
Idents(scData) <- "external_annotation"



#
print(table(Idents(scData)))

# 画 UMAP 图验证
pdf("umap_celltype.pdf", width = 7, height = 6)
DimPlot(scData, reduction = "umap", group.by = "external_annotation")
dev.off()

png("umap_celltype.png", width = 10, height = 6, units = "in", res = 300)
DimPlot(scData, reduction = "umap", group.by = "external_annotation")
dev.off()

saveRDS(scData, file.path(project_dir, "05_scData_annotated.rds"))
#file.remove(file.path(project_dir, "04_scData_clustered.rds"))




plan("multisession", workers = 15)

# Marker
scData.markers <- FindAllMarkers(
  scData,
  group.by = "external_annotation",
  only.pos = TRUE,
  logfc.threshold = 0.5,
  min.pct = 0.25,
  min.diff.pct = 0.1,
  verbose = TRUE
)

write.csv(scData.markers, "markers_for_annotation.csv")
plan("sequential")


#验证cluster之间Marker差异

scData.markers <- read.table("markers_for_annotation.csv", header = T, sep = ",", stringsAsFactors = FALSE)

#找top基因
top_markers <- scData.markers %>%
  filter(p_val_adj < 0.05) %>%            # 只保留显著的
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)           # 每个 cluster 取 logFC 最高的 5 个基因


pdf("MajorCellType_Marker_DotPlot.pdf", width = 14, height = 8)
DotPlot(scData, 
        features = unique(top_markers$gene), 
        group.by = "external_annotation") + 
  RotatedAxis() +
  labs(x = "Genes", y = "Cell type")
dev.off()

# 保存PNG
png("MajorCellType_Marker_DotPlot.png", width = 14, height = 8, units = "in", res = 300)
DotPlot(scData, 
        features = unique(top_markers$gene), 
        group.by = "external_annotation") + 
  RotatedAxis() +
  labs(x = "Genes", y = "Cell type")
dev.off()


# ============================ 差异分析与 KEGG 富集 ============================

scData <- readRDS(file.path(project_dir, "05_scData_annotated.rds"))
Idents(scData) <- "external_annotation"
scData$celltype.stim <- paste(scData$external_annotation, scData$group, sep = "_")
Idents(scData) <- "celltype.stim"


plan("multisession", workers = 15)
# 寻找所有组合的标记基因
sce.marker <- FindAllMarkers(scData, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
plan("sequential")

saveRDS(sce.marker, file.path(project_dir, "markers_split.rds"))

# 筛选显著基因并转换 ID
markers_sig <- sce.marker %>% filter(p_val_adj < 0.05)
gid <- bitr(unique(markers_sig$gene), "SYMBOL", "ENTREZID", OrgDb = "org.Hs.eg.db")
colnames(gid)[1] <- "gene"
markers_sig <- merge(markers_sig, gid, by = "gene")
markers_sig <- markers_sig %>%
  separate(cluster, into = c("celltype", "group"), sep = "_", remove = FALSE)


# Proxy settings removed for public repository
# KEGG 富集分析
x <- compareCluster(ENTREZID ~ celltype + group, data = markers_sig, fun = "enrichKEGG")

# 绘制点图，分别按细胞类型和分组展示通路
pdf("KEGG_dotplot_by_celltype.pdf", width = 20, height = 10)
dotplot(x, label_format = 60, x = "group") + facet_grid(~ celltype) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
pdf("KEGG_dotplot_by_group.pdf", width = 20, height = 10)
dotplot(x, label_format = 60, x = "celltype") + facet_grid(~ group) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


png("KEGG_dotplot_by_celltype.png", width = 20, height = 10,units = "in", res = 500)
dotplot(x, label_format = 60, x = "group") + facet_grid(~ celltype) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
png("KEGG_dotplot_by_group.png", width = 20, height = 10,units = "in", res = 500)
dotplot(x, label_format = 60, x = "celltype") + facet_grid(~ group) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



# ============================ 细胞通讯（CellChat） ============================
# 目的：推断细胞类型之间的配体-受体相互作用网络
setwd(file.path(project_dir, "CellChat"))
scData <- readRDS(file.path(project_dir, "05_scData_annotated.rds"))
data.input <- scData@assays$RNA$data
meta <- scData@meta.data

# 创建 CellChat 对象并设置数据库为 Secreted Signaling
cellchat <- createCellChat(object = data.input[, rownames(meta)], meta = meta,
                           group.by = "external_annotation")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "external_annotation")
cellchat@idents <- droplevels(cellchat@idents)

CellChatDB.use <- subsetDB(CellChatDB.human, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use

# 执行通讯概率计算
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# 导出所有配体-受体对
write.csv(subsetCommunication(cellchat), "cellchat_interactions.csv", row.names = FALSE)

# 交互作用热图
pdf("heatmap_interactions.pdf", width = 6, height = 5)
netVisual_heatmap(cellchat)
dev.off()

p <- netVisual_heatmap(cellchat)
pdf("heatmap_interactions.pdf", width = 6, height = 5); print(p); dev.off()
png("heatmap_interactions.png", width = 6, height = 5, units = "in", res = 500); print(p); dev.off()

# 交互数目与权重圈图
groupSize <- as.numeric(table(cellchat@idents))
p <- netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE,
                      label.edge = FALSE, title.name = "Number of interactions")
pdf("circle_number.pdf", width = 7, height = 9);print(p);dev.off()
png("circle_number.png", width = 7, height = 9, units = "in", res = 500); print(p); dev.off()


p <- netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE,
                 label.edge = FALSE, title.name = "Interaction weight/strength")
pdf("circle_weight.pdf", width = 7, height = 9); print(p); dev.off()
png("circle_weight.png", width = 7, height = 9, units = "in", res = 500); print(p); dev.off()

# 气泡图展示主要信号通路
p <- netVisual_bubble(cellchat, remove.isolate = FALSE)
pdf("bubble.pdf", width = 9, height = 12); print(p); dev.off()
png("bubble.png", width = 9, height = 12, units = "in", res = 500); print(p); dev.off()









#====================================END===========================================














# ============================ 基因集打分 ============================
# 目的：评估自定义基因集在各细胞类型中的平均表达活性
setwd(file.path(project_dir, "AUCell"))
scData <- readRDS(file.path(project_dir, "05_scData_annotated.rds"))
DefaultAssay(scData) <- "RNA"

# 定义待评估基因集


gene_set <- list(c("ADAMTS1","FBN1","SPARC","VCAM1"))
# 使用 AddModuleScore 计算每个细胞的基因集活性得分
scData <- AddModuleScore(scData, features = gene_set, ctrl = 100, name = "AUCell_Score")
colnames(scData@meta.data)[colnames(scData@meta.data) == "AUCell_Score1"] <- "AUCell_Score"

# 小提琴图展示各细胞类型的得分分布
AUCellP <- VlnPlot(scData, features = "AUCell_Score", pt.size = 0, adjust = 2,
                  group.by = "external_annotation")
ggsave("AUCell_violin.pdf",
       plot = AUCellP,
       width = 16, height = 8)
ggsave("AUCell_violin.png",
       plot = AUCellP,
       width = 16, height = 8)

saveRDS(scData, file.path(project_dir, "06_scData_AUCell.rds"))

# ============================ 关键基因可视化 ============================
# 在 UMAP 和小提琴图上展示指定基因的表达模式
setwd(file.path(project_dir, "feature_plots"))
scData <- readRDS(file.path(project_dir, "06_scData_AUCell.rds"))
genes_to_plot <- c("ADAMTS1","FBN1","SPARC","VCAM1")


theme.set <- theme(
  axis.title.x = element_blank(),
  axis.title.y = element_text(size = 20, face = "bold"),
  axis.text.x = element_text(size = 14, face = "bold"),
  axis.text.y = element_text(size = 14, face = "bold"),
  legend.text = element_text(size = 16, face = "bold"),
  legend.title = element_text(size = 18, face = "bold"),
  panel.border = element_rect(fill = NA, color = "black", size = 1.5, linetype = "solid")
)

genes_vec <- unlist(genes_to_plot)

plots <- lapply(genes_vec, function(g) {
  FeaturePlot_scCustom(scData, features = g, colors_use = viridis_magma_dark_high) +
    theme.set 
})


final_plot <- wrap_plots(plots, ncol = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right") 

ggsave("feature_umap_viridis.pdf", plot = final_plot, width = 18, height = 9)
ggsave("feature_umap_viridis.png", plot = final_plot, width = 18, height = 9, dpi = 300)


#气泡图===========
expr_mat <- GetAssayData(scData, assay = "RNA", layer = "data")
meta <- scData@meta.data


# 提取 3 个基因的表达，并转置
expr_sub <- t(as.matrix(expr_mat[genes_to_plot, ])) %>% as.data.frame()
expr_sub$cell <- rownames(expr_sub)

# 合并分组和细胞类型信息
plot_data <- merge(expr_sub, meta[, c("external_annotation", "group")], by.x = "cell", by.y = "row.names")
colnames(plot_data) <- c("cell", genes_to_plot, "celltype", "group")

# 计算每个分组中，每种细胞类型里各基因的表达比例和平均表达量
plot_summary <- plot_data %>%
  pivot_longer(cols = all_of(genes_to_plot), names_to = "gene", values_to = "expression") %>%
  group_by(celltype, group, gene) %>%
  summarise(
    avg_exp = mean(expression, na.rm = TRUE),               # 平均表达量（决定颜色）
    pct = sum(expression > 0) / n() * 100,                  # 表达比例（决定气泡大小）
    .groups = "drop"
  )

# 把分组设为因子，确保顺序为 Normal -> UC
plot_summary$group <- factor(plot_summary$group, levels = c("HC", "UC"))

# 绘制气泡图（横坐标=分组，纵坐标=细胞类型，用基因分面）
p <- ggplot(plot_summary, aes(x = group, y = celltype, color = avg_exp, size = pct)) +
  geom_point(alpha = 0.9) +
  scale_color_gradient2(low = "#3288BD", mid = "white", high = "#D53E4F", midpoint = 0) +
  scale_size(range = c(2, 8)) +
  facet_wrap(~gene, ncol = 2) +  # 三个基因并排分面展示
  labs(
    x = "Group",
    y = "Cell Type",
    color = "Avg Expression",
    size = "Percent Expressed (%)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 14, color = "black"),
    axis.text.y = element_text(face = "bold", size = 12, color = "black"),
    axis.title = element_text(face = "bold", size = 16),
    strip.text = element_text(face = "bold", size = 14),
    legend.position = "right"
  )

# 保存图片
ggsave("DotPlot_genes_by_group_and_celltype.pdf", p, width = 12, height = 6)
ggsave("DotPlot_genes_by_group_and_celltype.png", p, width = 12, height = 6, dpi = 300)




# UMAP 特征图（红蓝配色）
plots <- lapply(genes_to_plot, function(g) {
  FeaturePlot_scCustom(scData, features = g,
                       colors_use = colorRampPalette(c("#3288BD","white","#D53E4F"))(50))# +
    #NoAxes()
})
ggsave("feature_umap_rdbu.pdf", plot = wrap_plots(plots, ncol = 2), width = 12, height = 10)
ggsave("feature_umap_rdbu.png", plot = wrap_plots(plots, ncol = 2), width = 12, height = 10,dpi = 300)



# 小提琴图：按细胞类型展示基因表达分布
violin_genes <- VlnPlot(scData, features = genes_to_plot, group.by = "external_annotation",
        pt.size = 0.01) + theme.set + plot_layout(ncol = 2)

ggsave("violin_genes.pdf",plot = violin_genes,width = 16, height = 12)
ggsave("violin_genes.png",plot = violin_genes,width = 16, height = 12,dpi = 300)

p <- RidgePlot(scData, features = genes_to_plot, group.by = "external_annotation") + theme.set + plot_layout(ncol = 2)
ggsave("Ridge_genes.pdf",plot = p,width = 16, height = 12)
ggsave("Ridge_genes.png",plot = p,width = 16, height = 12,dpi = 300)




# 用箱线图展示 AUCell_Score 在不同细胞类型中的差异，并加上组间对比（UC vs HC）
p_box_score <- ggplot(scData@meta.data, aes(x = external_annotation, y = AUCell_Score, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  stat_compare_means(aes(group = group), method = "wilcox.test", label = "p.signif", hide.ns = TRUE) +
  scale_fill_manual(values = c("HC" = "#4682B4", "UC" = "#CD3700")) +
  labs(x = "Cell type", y = "AUCell Score", title = "Key Gene Set Activity by Cell Type") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("AUCell_score_boxplot.pdf", p_box_score, width = 14, height = 6)
ggsave("AUCell_score_boxplot.png", p_box_score, width = 14, height = 6, dpi = 300)

boxplot_raw_data <- scData@meta.data[, c("external_annotation", "AUCell_Score", "group")]
write.csv(boxplot_raw_data, "AUCell_score_boxplot_raw_data.csv", row.names = FALSE)



# 计算 UC 组和 HC 组中，每种细胞类型所占的比例
cell_prop <- as.data.frame(prop.table(table(scData$external_annotation, scData$group), margin = 2))
colnames(cell_prop) <- c("CellType", "Group", "Proportion")

# 画堆叠柱状图，直观对比哪类细胞在 UC 中变多/变少
p_prop <- ggplot(cell_prop, aes(x = Group, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_npg() +  # 使用 Nature 配色
  labs(x = "", y = "Proportion of cells (%)", title = "Cell Composition (UC vs HC)") +
  theme_bw() +
  theme(legend.position = "right",
        axis.text.x = element_text(face = "bold", size = 14),
        axis.text.y = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("Cell_Composition_StackBar.pdf", p_prop, width = 10, height = 8)
ggsave("Cell_Composition_StackBar.png", p_prop, width = 10, height = 8, dpi = 300)


write.csv(cell_prop, "Cell_composition_proportion_raw_data.csv", row.names = FALSE)

# ============================ 额外的统计检验（筛选关键细胞） ============================
# 用 Wilcoxon 检验找出 UC 中比例显著变化的细胞类型


cell_counts <- scData@meta.data %>%
  group_by(orig.ident, group, external_annotation) %>%  # 包含样本ID列
  tally() %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

# 2. 对每种细胞类型，用 Wilcoxon 检验比较 Normal 和 UC 组的比例差异
cell_stat_correct <- cell_counts %>%
  group_by(external_annotation) %>%
  rstatix::wilcox_test(proportion ~ group) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance("p")

# 3. 查看正确的 P 值（现在应该有真正的数字了，比如 0.03, 0.001 等）
print(cell_stat_correct %>% filter(p < 0.05))













# ============================ 细胞通讯 ============================

scData <- readRDS(file.path(project_dir, "05_scData_annotated.rds"))


# 定义分组列表
groups <- c("HC", "UC")
cellchat_list <- list()

# ============================ 循环处理两个分组 ============================
for (grp in groups) {
  cat("正在处理分组：", grp, "/n")
  
  # 1. 提取该分组的细胞和表达矩阵
  cells <- rownames(scData@meta.data)[scData@meta.data$group == grp]
  data.input <- scData@assays$RNA$data[, cells]
  meta <- scData@meta.data[cells, ]
  
  # 2. 创建 CellChat 对象
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "external_annotation")
  cellchat <- addMeta(cellchat, meta = meta)
  cellchat <- setIdent(cellchat, ident.use = "external_annotation")
  
  # 3. 加载数据库（沿用你的 Secreted Signaling 子集）
  CellChatDB.use <- subsetDB(CellChatDB.human, search = "Secreted Signaling")
  cellchat@DB <- CellChatDB.use
  
  # 4. 运行标准计算流程
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  # 5. 保存该组的 CellChat 对象，以备后续对比
  cellchat_list[[grp]] <- cellchat
  
  # 6. 输出该组的通讯结果表格
  write.csv(subsetCommunication(cellchat),
            file = paste0("cellchat_interactions_", grp, ".csv"),
            row.names = FALSE)
  
  # 7. 绘制该组的网络热图和圈图（单独保存）
  pdf(paste0("heatmap_interactions_", grp, ".pdf"), width = 6, height = 5)
  print(netVisual_heatmap(cellchat))
  dev.off()
  
  png(paste0("heatmap_interactions_", grp, ".png"), width = 6, height = 5, units = "in", res = 500)
  print(netVisual_heatmap(cellchat))
  dev.off()
  
  groupSize <- as.numeric(table(cellchat@idents))
  
  p1 <- netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                         weight.scale = TRUE, label.edge = FALSE,
                         title.name = paste0("Number of interactions (", grp, ")"))
  pdf(paste0("circle_number_", grp, ".pdf"), width = 7, height = 9)
  print(p1)
  dev.off()
  png(paste0("circle_number_", grp, ".png"), width = 7, height = 9, units = "in", res = 500)
  print(p1)
  dev.off()
  
  p2 <- netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
                         weight.scale = TRUE, label.edge = FALSE,
                         title.name = paste0("Interaction weight/strength (", grp, ")"))
  pdf(paste0("circle_weight_", grp, ".pdf"), width = 7, height = 9)
  print(p2)
  dev.off()
  png(paste0("circle_weight_", grp, ".png"), width = 7, height = 9, units = "in", res = 500)
  print(p2)
  dev.off()
  
  cat("✅ 完成分组：", grp, "/n/n")
}

# ============================ 两组对比分析（关键） ============================
library(ComplexHeatmap)

# 生成两组热图对象
heatmap_normal <- netVisual_heatmap(cellchat_list[["HC"]])
heatmap_UC <- netVisual_heatmap(cellchat_list[["UC"]])

# 保存合并的热图
pdf("heatmap_compare.pdf", width = 12, height = 5)
draw(heatmap_normal + heatmap_UC, column_title = "Interaction Heatmap (HC vs UC)")
dev.off()

png("heatmap_compare.png", width = 12, height = 5, units = "in", res = 300)
draw(heatmap_normal + heatmap_UC, column_title = "Interaction Heatmap (HC vs UC)")
dev.off()

bubble_normal <- netVisual_bubble(cellchat_list[["HC"]], remove.isolate = FALSE)
bubble_UC <- netVisual_bubble(cellchat_list[["UC"]], remove.isolate = FALSE)
combined_bubble <- bubble_normal + bubble_UC + plot_annotation(title = "Signaling Pathways (HC vs UC)")
ggsave("bubble_compare.pdf", combined_bubble, width = 16, height = 12)
ggsave("bubble_compare.png", combined_bubble, width = 16, height = 12, dpi = 300)



