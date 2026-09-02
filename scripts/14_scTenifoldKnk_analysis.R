##################################################
# Script: 14_scTenifoldKnk_analysis.R
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
project_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
cat("工作目录：", getwd(), "/n")

library(dplyr)
library(Seurat)
library(scTenifoldKnk)
library(ggplot2)
library(ggrepel)
library(forcats)
library(Matrix)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)


# 读取数据
scData <- readRDS(file.path("../", "01singleCell", "05_scData_annotated.rds"))

# 选择关键细胞类型
target_celltype <- "Stromal"
scData <- subset(scData, subset = external_annotation == target_celltype)

# 只保留 UC 组（疾病组）
scData <- subset(scData, subset = group == "UC")

# 处理 Seurat v5 的图层
if ("JoinLayers" %in% ls("package:Seurat")) {
  scData <- JoinLayers(scData, assay = "RNA")
}
DefaultAssay(scData) <- "RNA"

# 标准化、找高变基因
scData <- NormalizeData(scData, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
scData <- FindVariableFeatures(scData, assay = "RNA", selection.method = "vst", nfeatures = 3000, verbose = FALSE)
hvg <- VariableFeatures(scData)
hvg <- unique(hvg)

# 提取 counts 矩阵
countMatrix <- GetAssayData(scData, assay = "RNA", layer = "counts")

# 目标基因

genes_to_ko <- c("ADAMTS1", "FBN1", "SPARC","VCAM1")

# 过滤低表达基因：至少 5 个细胞表达，或为目标基因
keep_genes <- Matrix::rowSums(countMatrix > 0) >= 5 | rownames(countMatrix) %in% genes_to_ko
keep_cells <- Matrix::colSums(countMatrix > 0) >= 200
countMatrix <- countMatrix[keep_genes, keep_cells, drop = FALSE]

# 循环每个基因进行虚拟敲除
for (gene in genes_to_ko) {
  message("/n======================================")
  message("Running virtual knockout for: ", gene)
  
  if (!gene %in% rownames(countMatrix)) {
    warning("Gene ", gene, " not found in count matrix. Skip.")
    next
  }
  
  expr_cells <- sum(countMatrix[gene, ] > 0)
  cat(gene, " expressed in ", expr_cells, " cells/n")
  if (expr_cells < 10) {
    warning("Gene ", gene, " expressed in fewer than 10 cells. Skip.")
    next
  }
  
  # 构建用于 scTenifoldKnk 的基因集（高变基因 + 目标基因）
  hvg_use <- unique(c(hvg, genes_to_ko))
  hvg_use <- intersect(hvg_use, rownames(countMatrix))
  
  safe_min_cells <- max(50, round(ncol(countMatrix) * 0.03))
  hvg_use <- hvg_use[Matrix::rowSums(countMatrix[hvg_use, , drop = FALSE] > 0) >= safe_min_cells | hvg_use %in% genes_to_ko]
  if (!gene %in% hvg_use) hvg_use <- unique(c(hvg_use, gene))
  
  countMatrix_hvg <- countMatrix[hvg_use, , drop = FALSE]
  countMatrix_hvg <- as(countMatrix_hvg, "dgCMatrix")
  
  # 运行 scTenifoldKnk
  set.seed(123)
  result <- tryCatch({
    scTenifoldKnk(countMatrix = countMatrix_hvg, gKO = gene, qc = FALSE)
  }, error = function(e) {
    message("scTenifoldKnk failed for ", gene, ": ", e$message)
    return(NULL)
  })
  
  if (is.null(result)) next
  
  df_KO <- result[["diffRegulation"]]
  if (is.null(df_KO) || nrow(df_KO) == 0) {
    warning("No diffRegulation result for ", gene)
    next
  }
  
  if (!"gene" %in% colnames(df_KO)) df_KO$gene <- rownames(df_KO)
  df_KO$logFC <- ifelse(df_KO$FC > 0, log2(df_KO$FC), NA)
  
  # 保存结果
  write.csv(df_KO, file.path(paste0(gene, "_df_KO.result.csv")), row.names = FALSE)
  save(df_KO, file = file.path(paste0(gene, "_df_KO.rda")))
  
  # 绘制气泡图（火山图风格）
  df_plot <- df_KO %>%
    filter(!is.na(logFC), !is.na(p.value), p.value > 0) %>%
    mutate(neg_log10_pval = -log10(p.value),
           sig = ifelse(p.value < 0.05, "P < 0.05", "NS"))
  
  if (nrow(df_plot) == 0) next
  
  top_label <- df_plot %>% arrange(p.value) %>% head(10)
  
  p1 <- ggplot(df_plot, aes(x = logFC, y = neg_log10_pval)) +
    geom_point(aes(size = abs(logFC), color = logFC), alpha = 0.75) +
    scale_color_gradient2(low = "#2166AC", mid = "grey85", high = "#B2182B", midpoint = 0, name = "log2FC") +
    scale_size_continuous(range = c(1.5, 6), name = "|log2FC|") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    ggrepel::geom_text_repel(data = top_label, aes(label = gene), size = 3.5, max.overlaps = Inf) +
    labs(x = "log2 Fold Change", y = "-log10(P value)", title = paste0(gene, " virtual knockout")) +
    theme_classic(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(file.path(paste0(gene, "_bubble_plot.pdf")), p1, width = 7, height = 5.5)
  
  # Top10 上调柱状图
  df_top10_up <- df_KO %>%
    filter(!is.na(logFC), !is.na(p.value), p.value > 0) %>%
    arrange(desc(logFC)) %>%
    head(10) %>%
    mutate(neg_log10_pval = -log10(p.value), gene = fct_reorder(gene, logFC))
  
  if (nrow(df_top10_up) > 0) {
    p2 <- ggplot(df_top10_up, aes(x = logFC, y = gene)) +
      geom_col(aes(fill = neg_log10_pval), width = 0.7) +
      geom_text(aes(label = format(logFC, scientific = TRUE, digits = 3)), hjust = -0.15, size = 3.5) +
      scale_fill_gradient(low = "#C6DBEF", high = "#08519C", name = "-log10(P)") +
      labs(x = "log2 Fold Change", y = NULL, title = paste0(gene, " KO: Top 10 up-regulated genes")) +
      theme_classic(base_size = 14) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
      expand_limits(x = max(df_top10_up$logFC, na.rm = TRUE) * 1.3)
    ggsave(file.path(paste0(gene, "_TOP10_up_bar_plot.pdf")), p2, width = 7, height = 5)
  }
  
  # Top10 下调柱状图
  df_top10_down <- df_KO %>%
    filter(!is.na(logFC), !is.na(p.value), p.value > 0) %>%
    arrange(logFC) %>%
    head(10) %>%
    mutate(neg_log10_pval = -log10(p.value), gene = fct_reorder(gene, logFC))
  
  if (nrow(df_top10_down) > 0) {
    p3 <- ggplot(df_top10_down, aes(x = logFC, y = gene)) +
      geom_col(aes(fill = neg_log10_pval), width = 0.7) +
      geom_text(aes(label = format(logFC, scientific = TRUE, digits = 3)), hjust = 1.1, size = 3.5, color = "white") +
      scale_fill_gradient(low = "#FEE0D2", high = "#A50F15", name = "-log10(P)") +
      labs(x = "log2 Fold Change", y = NULL, title = paste0(gene, " KO: Top 10 down-regulated genes")) +
      theme_classic(base_size = 14) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
      expand_limits(x = min(df_top10_down$logFC, na.rm = TRUE) * 1.3)
    ggsave(file.path(paste0(gene, "_TOP10_down_bar_plot.pdf")), p3, width = 7, height = 5)
  }
  
  # 保存全部差异基因按P值排序
  df_top_sig <- df_KO %>% filter(!is.na(logFC), !is.na(p.value)) %>% arrange(p.value)
  write.csv(df_top_sig, file.path(paste0(gene, "_KO_all_genes_sorted_by_pvalue.csv")), row.names = FALSE)
  
  # --- 富集分析 ---
  message("Running GO and KEGG enrichment for ", gene)
  
  sig_genes <- df_KO %>% filter(p.value < 0.05, !is.na(logFC), abs(logFC) > 0.5) %>% pull(gene)
  if (length(sig_genes) < 10) {
    warning("Significant genes less than 10, skip enrichment.")
    next
  }
  
  gene_trans <- tryCatch({
    bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) NULL)
  
  if (is.null(gene_trans) || nrow(gene_trans) == 0) next
  
  # GO 富集
  go_res <- enrichGO(gene = gene_trans$ENTREZID, OrgDb = org.Hs.eg.db, ont = "ALL",
                     pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
  if (!is.null(go_res) && nrow(go_res) > 0) {
    write.csv(as.data.frame(go_res), file.path(paste0(gene, "_GO_ALL_results.csv")), row.names = FALSE)
    p_go <- dotplot(go_res, showCategory = 5, split = "ONTOLOGY", title = paste0("GO Enrichment - ", gene, " KO")) +
      facet_grid(ONTOLOGY ~ ., scale = "free") +
      theme_bw() + theme(plot.title = element_text(face = "bold", hjust = 0.5))
    ggsave(file.path(paste0(gene, "_GO_ALL_dotplot.pdf")), p_go, width = 8, height = 8)
  }
  
  # KEGG 富集（使用本地数据避免网络问题）
  kegg_res <- tryCatch({
    enrichKEGG(gene = gene_trans$ENTREZID, organism = 'hsa', pAdjustMethod = "BH",
               pvalueCutoff = 0.99, qvalueCutoff = 0.99, use_internal_data = TRUE)
  }, error = function(e) NULL)
  
  if (!is.null(kegg_res) && nrow(kegg_res) > 0) {
    kegg_res <- setReadable(kegg_res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    kegg_df <- as.data.frame(kegg_res) %>% filter(pvalue < 0.05) %>% arrange(desc(Count))
    if (nrow(kegg_df) > 0) {
      write.csv(kegg_df, file.path(paste0(gene, "_KEGG_results.csv")), row.names = FALSE)
      top_kegg <- kegg_df %>% arrange(p.adjust) %>% head(10) %>%
        mutate(GeneRatio_Num = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2])),
               Description = fct_reorder(Description, Count))
      p_kegg <- ggplot(top_kegg, aes(x = GeneRatio_Num, y = Description)) +
        geom_point(aes(size = Count, color = p.adjust)) +
        scale_color_gradient(low = "red", high = "blue", name = "p.adjust") +
        scale_size(range = c(3, 8), name = "Count") +
        theme_bw() +
        labs(title = paste0("KEGG Enrichment - ", gene, " KO"), x = "GeneRatio", y = NULL) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5))
      ggsave(file.path(paste0(gene, "_KEGG_dotplot.pdf")), p_kegg, width = 8, height = 6)
    }
  }
  message("Finished KO for: ", gene)
}

