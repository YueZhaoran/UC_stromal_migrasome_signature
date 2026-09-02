##################################################
# Script: 01_bulk_transcriptome_preprocessing.R
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

library(GEOquery)

# 工作目录和代理
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
cat("工作目录：", getwd(), "/n")

#网络代理
# Proxy settings removed for public repository

GEO_data <- "GSE87466"


gse_file <- "GSE87466_series_matrix.txt.gz"
gse <- getGEO(filename = gse_file)


expr <- exprs(gse)
expr <- expr[, colSums(is.na(expr)) == 0, drop = FALSE]

# log2 转换判断
qx <- as.numeric(quantile(expr, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
LogC <- (qx[5] > 100) ||
  (qx[6] - qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) {
  expr[expr <= 0] <- 0
  expr <- log2(expr + 1)
  print("log2 transform finished")
} else {
  print("log2 transform not needed")
}

#探针注释清洗
probe_anno <- fData(gse)
# 找到实际列名
id_col <- "ID"
symbol_col <- grep("Gene.*Symbol", colnames(probe_anno), value = TRUE)[1]
if (is.na(symbol_col)) stop("未找到 Gene Symbol 列")

probe_df <- probe_anno[, c(id_col, symbol_col)]
colnames(probe_df) <- c("ID", "Gene_Symbol")


probe_df <- probe_df[probe_df$Gene_Symbol != "" & !is.na(probe_df$Gene_Symbol), ]
# 分割 "///" 取第一个基因符号
split_symbols <- strsplit(as.character(probe_df$Gene_Symbol), "///", fixed = TRUE)
probe_df$symbol <- trimws(sapply(split_symbols, `[`, 1))
# 去除 " --- " 和空格
probe_df <- probe_df[probe_df$symbol != "---" & probe_df$symbol != " --- ", ]
probe_df$symbol <- gsub(" ", "", probe_df$symbol)
# 只保留 ID 和 symbol
probe2symbol <- probe_df[, c("ID", "symbol")]

#合并表达矩阵与注释
#统一 ID 为字符型
probe_id <- as.character(probe2symbol$ID)
expr_ids <- rownames(expr)

# 匹配并去除无注释的探针
match_idx <- match(probe_id, expr_ids)
valid <- !is.na(match_idx)
match_idx <- match_idx[valid]
probe_symbols <- probe2symbol$symbol[valid]

# 提取子矩阵，计算行均值
sub_expr <- expr[match_idx, , drop = FALSE]
row_means <- rowMeans(sub_expr, na.rm = TRUE)

# 按 symbol 和 -mean 排序，保留每个基因第一次出现（最高表达）
ord <- order(probe_symbols, -row_means)
sub_expr <- sub_expr[ord, , drop = FALSE]
probe_symbols <- probe_symbols[ord]

keep <- !duplicated(probe_symbols)
sub_expr <- sub_expr[keep, , drop = FALSE]
rownames(sub_expr) <- probe_symbols[keep]

dat <- sub_expr   # 基因表达矩阵（matrix）

# 过滤蛋白编码基因
protein_gene <- read.table("PCGv50.xls", header = TRUE)  # 假设有表头 gene_name
# 如果文件没有表头，可以用 read.table(..., col.names = c("gene_name"))
dat <- dat[rownames(dat) %in% protein_gene$gene_name, , drop = FALSE]

# 构建分组信息
pd <- pData(gse)
table(pd$title)
table(pd$`disease:ch1`)

colnames(pd)

group <- data.frame(
  sample = pd$geo_accession,
  group = pd$`disease:ch1`,
  stringsAsFactors = FALSE
)


group$group <- ifelse(group$group == "Normal", "Normal", "UC")

# 查看分组
table(group$group)


common_samples <- intersect(colnames(dat), group$sample)
dat <- dat[, common_samples, drop = FALSE]
group <- group[group$sample %in% common_samples, ]


write.csv(dat, file = paste0("dat.", GEO_data, ".csv"))
write.csv(group, file = paste0("group.", GEO_data, ".csv"), row.names = FALSE)
write.csv(pd, file = paste0("pd", GEO_data, ".csv"), row.names = FALSE)










rm(list = ls())
gc()

library(GEOquery)

# 工作目录和代理
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
cat("工作目录：", getwd(), "/n")

#网络代理
# Proxy settings removed for public repository

GEO_data <- "GSE75214"


gse_file <- "GSE75214_series_matrix.txt.gz"
gse <- getGEO(filename = gse_file)


expr <- exprs(gse)
expr <- expr[, colSums(is.na(expr)) == 0, drop = FALSE]

# log2 转换判断
qx <- as.numeric(quantile(expr, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
LogC <- (qx[5] > 100) ||
  (qx[6] - qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) {
  expr[expr <= 0] <- 0
  expr <- log2(expr + 1)
  print("log2 transform finished")
} else {
  print("log2 transform not needed")
}



#探针注释提取（基于 gene_assignment，取第一个基因符号）
probe_anno <- fData(gse)
gene_assign <- as.character(probe_anno$gene_assignment)
gene_assign[is.na(gene_assign)] <- ""

extract_first_symbol <- function(gene_str) {
  if (is.na(gene_str) || gene_str == "" || grepl("^---", gene_str)) return("")
  first_block <- strsplit(gene_str, " /// ")[[1]][1]
  fields <- strsplit(first_block, " // ")[[1]]
  if (length(fields) >= 2) {
    sym <- trimws(fields[2])
    if (sym == "---") sym <- ""
    return(sym)
  } else {
    return("")
  }
}

first_symbols <- sapply(gene_assign, extract_first_symbol, USE.NAMES = FALSE)
probe2symbol <- data.frame(
  ID = as.character(probe_anno$ID),
  symbol = first_symbols,
  stringsAsFactors = FALSE
)
probe2symbol <- probe2symbol[probe2symbol$symbol != "", ]

head(probe2symbol)


#合并表达矩阵与注释
#统一 ID 为字符型
probe_id <- as.character(probe2symbol$ID)
expr_ids <- rownames(expr)

# 匹配并去除无注释的探针
match_idx <- match(probe_id, expr_ids)
valid <- !is.na(match_idx)
match_idx <- match_idx[valid]
probe_symbols <- probe2symbol$symbol[valid]

# 提取子矩阵，计算行均值
sub_expr <- expr[match_idx, , drop = FALSE]
row_means <- rowMeans(sub_expr, na.rm = TRUE)

# 按 symbol 和 -mean 排序，保留每个基因第一次出现（最高表达）
ord <- order(probe_symbols, -row_means)
sub_expr <- sub_expr[ord, , drop = FALSE]
probe_symbols <- probe_symbols[ord]

keep <- !duplicated(probe_symbols)
sub_expr <- sub_expr[keep, , drop = FALSE]
rownames(sub_expr) <- probe_symbols[keep]

dat <- sub_expr   # 基因表达矩阵（matrix）

# 过滤蛋白编码基因
protein_gene <- read.table("PCGv50.xls", header = TRUE)  # 假设有表头 gene_name
# 如果文件没有表头，可以用 read.table(..., col.names = c("gene_name"))
dat <- dat[rownames(dat) %in% protein_gene$gene_name, , drop = FALSE]

# 构建分组信息
pd <- pData(gse)
table(pd$title)
table(pd$`disease:ch1`)

colnames(pd)


keep <- grepl("^(control_colon|UC_colon_active)", pd$title)
group <- data.frame(
  sample = pd$geo_accession[keep],
  title  = pd$title[keep],
  stringsAsFactors = FALSE
)

group$group <- ifelse(grepl("^control_colon", group$title), "Normal", "UC")

# 如果不再需要原始 title 列，可以删掉
group$title <- NULL

# 检查分组结果
table(group$group)

head(group)

common_samples <- intersect(colnames(dat), group$sample)
dat <- dat[, common_samples, drop = FALSE]
group <- group[group$sample %in% common_samples, ]


ncol(dat)

write.csv(dat, file = paste0("dat.", GEO_data, ".csv"))
write.csv(group, file = paste0("group.", GEO_data, ".csv"), row.names = FALSE)
write.csv(pd, file = paste0("pd", GEO_data, ".csv"), row.names = FALSE)







#核验数据

set.seed(123)


selected_genes <- sample(rownames(dat), 2)
cat("随机选取的两个基因：", selected_genes, "/n/n")


final_samples <- colnames(dat)


for (gene in selected_genes) {
  cat("========== 基因：", gene, " ==========/n")
  
  probe_ids <- probe2symbol$ID[probe2symbol$symbol == gene]
  cat("该基因对应的探针数目：", length(probe_ids), "/n")
  cat("探针ID：", paste(probe_ids, collapse = ", "), "/n/n")
  
  
  probes_expr <- expr[probe_ids, final_samples, drop = FALSE]
  
  
  probe_means <- rowMeans(probes_expr, na.rm = TRUE)
  cat("各探针平均表达量：/n")
  print(round(probe_means, 4))
  
  
  best_probe <- names(which.max(probe_means))
  cat("/n最高表达的探针（应被保留）：", best_probe, "/n")
  
  
  values_dat <- dat[gene, ]
  values_probe <- expr[best_probe, final_samples]
  
  
  diff_max <- max(abs(values_dat - values_probe), na.rm = TRUE)
  cat("两向量最大差异：", diff_max, "/n")
  
  if (diff_max < 1e-10) {
    cat("✅ 验证通过：数据完全一致，逻辑正确！/n/n")
  } else {
    cat("❌ 存在差异，请检查数据")
  }
}




