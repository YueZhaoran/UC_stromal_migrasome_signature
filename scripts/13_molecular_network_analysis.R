##################################################
# Script: 13_molecular_network_analysis.R
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
project_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
cat("工作目录：", getwd(), "/n")

# 加载包
library(dplyr)
library(tidyr)
library(multiMiR)
library(httr)

# 关键基因列表
hub_gene <- read.csv("../05DEGs/intersection_gene.csv", stringsAsFactors = FALSE)
key_genes <- hub_gene$symbol
print(paste0("关键基因列表：", paste(key_genes, collapse = ", ")))

# TF-基因调控网络
trrust <- read.table("trrust_rawdata.human.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(trrust) <- c("TF", "Target", "Direction", "PMID")

tf_network <- trrust %>%
  filter(Target %in% key_genes) %>%
  dplyr::select(TF, Target, Direction) %>%
  mutate(Interaction = paste0(Direction, " (", TF, "->", Target, ")"))

write.csv(tf_network, "01.TF_KeyGene_Network.csv", row.names = FALSE)
print(paste0("找到 ", nrow(tf_network), " 对 TF-基因调控关系，已保存。"))

# ============================ mirdb miRNA-mRNA 预测============================
mirna_res <- get_multimir(org = "hsa", 
                          target = key_genes, 
                          table = "mirdb") 

data_df <- mirna_res@data

mirna_df <- data_df %>%
  dplyr::select(mature_mirna_id, target_symbol, score) %>%
  filter(score >= 80) %>%
  distinct()

write.csv(mirna_df, "miRNA_mRNA_Network.csv", row.names = FALSE)
print(paste0("找到 ", nrow(mirna_df), " 对 miRNA-mRNA 调控关系。"))



#lncRNA-miRNA 查询

mirna_df <- read.csv("miRNA_mRNA_Network.csv", stringsAsFactors = FALSE)
unique_mirnas <- unique(mirna_df$mature_mirna_id)
print(paste0("共提取到 ", length(unique_mirnas), " 个唯一 miRNA，开始查询..."))

# ============================ 查询 ENCORI API  ============================


# 1. 在循环外初始化空数据框
all_lncRNA_full <- data.frame()

# 2. 开始循环
for (mir in unique_mirnas) {
  cat("/n===== 正在查询 miRNA:", mir, "=====/n")
  
  api_url <- paste0(
    "https://rnasysu.com/encori/api/miRNATarget/?",
    "assembly=hg38",
    "&geneType=lncRNA",   # 这里限定了就是查询 lncRNA
    "&miRNA=", URLencode(mir, reserved = TRUE),
    "&clipExpNum=0",
    "&degraExpNum=0",
    "&pancancerNum=0",
    "&programNum=0",
    "&program=None",
    "&target=all",
    "&cellType=all"
  )
  
  response <- GET(api_url)
  status <- status_code(response)
  cat("HTTP状态码:", status, "/n")
  
  if (status != 200) { cat("⚠️ 请求失败/n"); next }
  
  raw_text <- content(response, "text", encoding = "UTF-8")
  lines <- strsplit(raw_text, "/n")[[1]]
  clean_lines <- lines[!grepl("^#", lines)]
  
  if (length(clean_lines) < 2) { cat("⚠️ 无数据行返回/n"); next }
  
  clean_text <- paste(clean_lines, collapse = "/n")
  con <- textConnection(clean_text)
  temp_df <- read.table(con, header = TRUE, sep = "/t", fill = TRUE, quote = "", stringsAsFactors = FALSE)
  close(con)
  
  cat("原始行数:", nrow(temp_df), "/n")
  
  # 清洗列名：去除隐藏的 /r
  colnames(temp_df) <- trimws(colnames(temp_df), whitespace = "[ /t/r/n]")
  
  # 拦截 API 错误信息
  if (nrow(temp_df) == 1) {
    first_row_text <- paste(temp_df[1, ], collapse = " ")
    if (grepl("haven't been set correctly|not available|parameter", first_row_text, ignore.case = TRUE)) {
      cat("⚠️ API 返回错误，跳过。/n"); next
    }
  }
  
  # 动态匹配列名
  miRNA_col <- colnames(temp_df)[grep("miRNA", colnames(temp_df), ignore.case = TRUE)][1]
  gene_col <- colnames(temp_df)[grep("geneName|targetGene", colnames(temp_df), ignore.case = TRUE)][1]
  # 【新增】额外匹配 geneType 列，以便验证数据类别
  type_col <- colnames(temp_df)[grep("geneType", colnames(temp_df), ignore.case = TRUE)][1]
  
  if (is.na(miRNA_col) || is.na(gene_col)) { cat("⚠️ 缺少关键列，跳过/n"); next }
  
  # 提取需要保留的列（把 geneType 也加进 any_of 中）
  quality_cols <- intersect(c("clipExpNum", "RBP", "geneType"), colnames(temp_df))
  temp_sub <- temp_df %>%
    dplyr::select(any_of(c(miRNA_col, gene_col, quality_cols)))
  
  # 使用 Base R 精确重命名
  colnames(temp_sub)[colnames(temp_sub) == miRNA_col] <- "miRNA"
  colnames(temp_sub)[colnames(temp_sub) == gene_col] <- "lncRNA"
  
  # 关键修复：ENCORI 返回的 miRNA 列为 MIMAT 编号（如 MIMAT0000100），
  # 而 multiMiR 预测的 miRNA-mRNA 网络使用成熟 miRNA 名字（如 hsa-miR-29b-3p）。
  # 若保留 MIMAT 编号，Cytoscape 会把同一 miRNA 识别为两个不同节点，
  # 导致 lncRNA-miRNA 子网和 miRNA-mRNA 子网无法连通，出现大量孤立成团。
  # 因此将 miRNA 列统一替换为查询时使用的成熟 miRNA 名字。
  temp_sub$miRNA <- mir
  
  # 去重
  temp_sub <- temp_sub %>%
    distinct()
  
  cat("去重后行数:", nrow(temp_sub), "/n")
  
  # 拼接到总表
  all_lncRNA_full <- bind_rows(all_lncRNA_full, temp_sub)
  
  Sys.sleep(0.8) # 延时避封
}

# 运行结束后，查看总数据量
cat("/n收集完毕，总数据行数:", nrow(all_lncRNA_full), "/n")

# ============================ 数据清洗 ============================
# 删除 miRNA 或 lncRNA 为 NA 的行
all_lncRNA_full <- all_lncRNA_full %>%
  filter(!is.na(miRNA) & !is.na(lncRNA) & lncRNA != "")

# 将 clipExpNum 中的 NA 全部转为 0（表示无实验证据）
if ("clipExpNum" %in% colnames(all_lncRNA_full)) {
  all_lncRNA_full$clipExpNum[is.na(all_lncRNA_full$clipExpNum)] <- 0
} else {
  # 如果连 clipExpNum 都没有，就新建一列全为 0
  all_lncRNA_full$clipExpNum <- 0
}

# 如果有 RBP 列，将 NA 转为空字符串
if ("RBP" %in% colnames(all_lncRNA_full)) {
  all_lncRNA_full$RBP[is.na(all_lncRNA_full$RBP)] <- ""
}

# 保存清洗后的完整数据
write.csv(all_lncRNA_full, "lncRNA_miRNA_Network_Full.csv", row.names = FALSE)
print(paste0("清洗后，总共获得 ", nrow(all_lncRNA_full), " 对 lncRNA-miRNA 关系（已保存）。"))

# 新的质量分级策略
lncrna_full <- read.csv("lncRNA_miRNA_Network_Full.csv", stringsAsFactors = FALSE)

cat("/n========== 质量分布统计 ==========/n")
cat("clipExpNum 分布（CLIP-seq 实验证据数）：/n")
print(table(lncrna_full$clipExpNum))

# ---- 分级 ----
# S级：有 CLIP 实验证据（clipExpNum >= 1）
s_level <- lncrna_full %>% filter(clipExpNum >= 1)

# A级：无 CLIP 证据，但在 RBP 列里有记录（说明有蛋白结合证据，也算间接支持）
# 注意：RBP 列记录的是结合的 RNA 结合蛋白，如果非空，说明该 lncRNA 有结合证据
if ("RBP" %in% colnames(lncrna_full)) {
  a_level <- lncrna_full %>% filter(clipExpNum == 0 & RBP != "")
} else {
  a_level <- data.frame()  # 如果没有 RBP 列，A级为空
}

# B级：既无 CLIP 也无 RBP，纯软件预测（保底数据）
b_level <- lncrna_full %>% 
  filter(clipExpNum == 0) %>%
  anti_join(a_level, by = c("miRNA", "lncRNA"))

cat("/n========== 分级统计 ==========/n")
cat("S级（有 CLIP 实验证据）: ", nrow(s_level), " 条/n")
cat("A级（无 CLIP，但有 RBP 结合证据）: ", nrow(a_level), " 条/n")
cat("B级（纯预测，无任何实验证据）: ", nrow(b_level), " 条/n")

# 生成高质量网络
high_quality <- bind_rows(s_level, a_level) %>% distinct()

# 如果高质量数据（S+A）太少（比如少于 10 条），
# 就从 B 级里随机抽取一些来补全（优先选出现频率高的 lncRNA）
if (nrow(high_quality) < 20 & nrow(b_level) > 0) {
  cat("/n高质量关系较少，从 B 级中补充高频 lncRNA.../n")
  
  # 统计 B 级中每个 lncRNA 出现的次数（即被几个 miRNA 靶向）
  b_freq <- b_level %>%
    group_by(lncRNA) %>%
    summarise(freq = n(), .groups = "drop") %>%
    arrange(desc(freq))
  
  # 取出现频率最高的前 30 个 lncRNA 对应的边
  top_b_edges <- b_level %>%
    filter(lncRNA %in% head(b_freq$lncRNA, 30)) %>%
    head(50)  # 最多补充 50 条
  
  high_quality <- bind_rows(high_quality, top_b_edges) %>% distinct()
  cat("补充了 ", nrow(top_b_edges), " 条高频 B 级关系/n")
}


# 生成最终用于 Cytoscape 的边列表（包含类型列）
final_edges <- high_quality %>%
  #head(100) %>%
  dplyr::select(Source = lncRNA, Target = miRNA) %>%
  mutate(
    Interaction = "lncRNA-miRNA",
    Source_Type = "lncRNA",
    Target_Type = "miRNA"
  )

nrow(final_edges)
write.csv(final_edges, "High_Quality_LncRNA_miRNA_for_Cytoscape.csv", row.names = FALSE)




tf_df <- read.csv("01.TF_KeyGene_Network.csv", stringsAsFactors = FALSE)   # 含 TF, Target, Direction
mirna_df <- read.csv("miRNA_mRNA_Network.csv", stringsAsFactors = FALSE)   # 含 mature_mirna_id, target_symbol, score
lncrna_df <- read.csv("High_Quality_LncRNA_miRNA_for_Cytoscape.csv", stringsAsFactors = FALSE)  # 我们刚刚生成的

# ============================ 构建三种边表（统一列名） ============================
# TF -> Gene
tf_edges <- tf_df %>%
  dplyr::select(Source = TF, Target = Target) %>%
  mutate(
    Interaction = "TF-Gene",
    Source_Type = "TF",
    Target_Type = "Gene"
  )

#  miRNA -> Gene
mir_edges <- mirna_df %>%
  dplyr::select(Source = mature_mirna_id, Target = target_symbol) %>%
  mutate(
    Interaction = "miRNA-mRNA",
    Source_Type = "miRNA",
    Target_Type = "Gene"
  )

# lncRNA -> miRNA
# 如果 lncrna_df 已有类型列则直接使用，否则手动添加
if (!("Source_Type" %in% colnames(lncrna_df)) || !("Target_Type" %in% colnames(lncrna_df))) {
  lncrna_edges <- lncrna_df %>%
    dplyr::select(Source, Target) %>%
    mutate(
      Interaction = "lncRNA-miRNA",
      Source_Type = "lncRNA",
      Target_Type = "miRNA"
    )
} else {
  lncrna_edges <- lncrna_df %>%
    dplyr::select(Source, Target, Interaction, Source_Type, Target_Type)
}

# 合并所有边
all_edges <- bind_rows(tf_edges, mir_edges, lncrna_edges)

# 去重（保留第一次出现）
all_edges <- all_edges %>%
  distinct(Source, Target, Interaction, .keep_all = TRUE)



# 生成节点属性表
source_attrs <- all_edges %>% dplyr::select(Node = Source, Type = Source_Type)
target_attrs <- all_edges %>% dplyr::select(Node = Target, Type = Target_Type)

node_attrs <- bind_rows(source_attrs, target_attrs) %>%
  distinct(Node, .keep_all = TRUE) %>%   # 按节点名去重，保留第一次出现的类型
  arrange(Node)

write.csv(node_attrs, "Node_Attributes.csv", row.names = FALSE)
cat("节点属性表已保存，共 ", nrow(node_attrs), " 个节点/n")

# 保存最终边表
# 简洁版（只含三列 Source, Target, Interaction）
clean_edges <- all_edges %>%
  dplyr::select(Source, Target, Interaction)

write.csv(clean_edges, "Complete_Regulatory_Network.csv", row.names = FALSE, quote = FALSE)

# 保存带类型的完整版
write.csv(all_edges, "Complete_Regulatory_Network_with_Type.csv", row.names = FALSE)

