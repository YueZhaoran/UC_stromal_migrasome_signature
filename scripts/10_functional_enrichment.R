##################################################
# Script: 10_functional_enrichment.R
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
if (rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
cat("工作目录：", getwd(), "/n")

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
set.seed(798)

# KEGG 在线查询需要网络，设置本地代理
# Proxy settings removed for public repository

# ============================  读取候选基因 ============================
candidate_genes <- read.csv("Venn_intersection_genes.csv", stringsAsFactors = FALSE)$gene
cat("候选基因数：", length(candidate_genes), "/n")

# 基因 ID 转换（SYMBOL -> ENTREZID）
gene_id <- bitr(candidate_genes,
                fromType = "SYMBOL",
                toType   = "ENTREZID",
                OrgDb    = org.Hs.eg.db)
cat("成功转换基因数：", nrow(gene_id), "/n")

# ============================  GO 富集分析 ============================
go_res <- enrichGO(gene          = gene_id$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.2,
                   readable      = TRUE)

cat("GO 显著富集条目数：", nrow(go_res@result[go_res@result$p.adjust < 0.05, ]), "/n")
write.csv(go_res@result, "GO_enrichment_result.csv", row.names = FALSE)

# GO 点图（BP/CC/MF 三个模块各展示 Top5）
p_go <- dotplot(go_res, showCategory = 5, split = "ONTOLOGY") +
  facet_grid(ONTOLOGY ~ ., scales = "free") +
  labs(title = "GO enrichment of candidate genes") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("05.GO_enrichment_dotplot.pdf", p_go, width = 9, height = 8)
ggsave("05.GO_enrichment_dotplot.png", p_go, width = 9, height = 8, dpi = 600)

# ============================  KEGG 富集分析 ============================
kegg_res <- enrichKEGG(gene          = gene_id$ENTREZID,
                       organism      = "hsa",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.2)

# 将 ENTREZID 转回 SYMBOL 便于阅读
kegg_res <- setReadable(kegg_res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

cat("KEGG 显著富集通路数：", nrow(kegg_res@result[kegg_res@result$p.adjust < 0.05, ]), "/n")
write.csv(kegg_res@result, "KEGG_enrichment_result.csv", row.names = FALSE)

# KEGG 点图（展示 Top10 通路）
p_kegg <- dotplot(kegg_res, showCategory = 10) +
  labs(title = "KEGG enrichment of candidate genes") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("06.KEGG_enrichment_dotplot.pdf", p_kegg, width = 9, height = 6)
ggsave("06.KEGG_enrichment_dotplot.png", p_kegg, width = 9, height = 6, dpi = 600)

# ============================ 打印 Top 结果 ============================
cat("/n===== GO 各模块 Top5 =====/n")
print(go_res@result %>%
        filter(p.adjust < 0.05) %>%
        group_by(ONTOLOGY) %>%
        arrange(p.adjust, .by_group = TRUE) %>%
        slice_head(n = 5) %>%
        dplyr::select(ONTOLOGY, Description, Count, p.adjust) %>%
        as.data.frame())

cat("/n===== KEGG Top10 =====/n")
print(kegg_res@result %>%
        arrange(p.adjust) %>%
        head(10) %>%
        dplyr::select(Description, Count, p.adjust) %>%
        as.data.frame())
