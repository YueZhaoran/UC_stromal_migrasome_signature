##################################################
# Script: 09b_validation_statistics.R
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



library(ggplot2)
library(ggpubr)
library(rstatix)
library(ggh4x)
library(pROC)
library(tidyr)
library(dplyr)

# ============================ 1. 读取数据 ============================
hub_gene <- read.csv("intersection_gene.csv", stringsAsFactors = FALSE)

# 读取训练集 (GSE87466)
train_data <- read.csv("../00bulkdata/dat.GSE87466.csv", row.names = 1, check.names = FALSE)
train_group <- read.csv("../00bulkdata/group.GSE87466.csv", stringsAsFactors = FALSE)
colnames(train_group) <- c("sample", "group")

# 读取验证集 (GSE75214)
ver_data <- read.csv("../00bulkdata/dat.GSE75214.csv", row.names = 1, check.names = FALSE)
ver_group <- read.csv("../00bulkdata/group.GSE75214.csv", stringsAsFactors = FALSE)
colnames(ver_group) <- c("sample", "group")

#训练集处理与绘图
temp1 <- train_data[rownames(train_data) %in% hub_gene$symbol, ]

# 计算均值与方向
disease_sample <- train_group[train_group$group == 'UC', ]
disease_data <- temp1[, colnames(temp1) %in% disease_sample$sample]
disease_mean <- apply(disease_data, 1, mean) %>% as.data.frame()
colnames(disease_mean) <- 'disease_mean'

control_sample <- train_group[train_group$group == 'Normal', ]
control_data <- temp1[, colnames(temp1) %in% control_sample$sample]
control_mean <- apply(control_data, 1, mean) %>% as.data.frame()
colnames(control_mean) <- 'control_mean'

train_mean <- cbind(disease_mean, control_mean)
train_mean$compare_res <- ifelse(train_mean$disease_mean > train_mean$control_mean, 1, 0)
train_mean$symbol <- rownames(train_mean)

# 整理长数据做差异检验
temp1$symbol <- rownames(temp1)
temp1_long <- pivot_longer(temp1, cols = -symbol, names_to = 'sample', values_to = 'expr')
train_hub_exp <- merge(temp1_long, train_group, by = 'sample')
train_hub_exp$group <- as.factor(train_hub_exp$group)

# 统计检验
train_stat.test <- train_hub_exp %>%
  group_by(symbol) %>%
  wilcox_test(expr ~ group) %>%
  adjust_pvalue(method = 'fdr')

train_stat.test$p_label <- ifelse(train_stat.test$p < 0.001, "***",
                                  ifelse(train_stat.test$p < 0.01, "**",
                                         ifelse(train_stat.test$p < 0.05, "*", 'ns')))

# 合并结果
train_finally_res <- merge(train_stat.test, train_mean, by = 'symbol')
train_finally_res <- train_finally_res %>% select(symbol, p, p_label, compare_res)
train_finally_res$p_sig <- ifelse(train_finally_res$p < 0.05, 1, 0)

# 绘图（训练集）
train_stat.test <- train_stat.test %>% add_xy_position(x = 'group', dodge = 0.5, step.increase = 0.1)

train_va <- ggplot(train_hub_exp, aes(x = group, y = expr, color = group)) +
  geom_violin(trim = FALSE, color = "black", aes(fill = group)) +
  stat_boxplot(geom = "errorbar", width = 0.1, position = position_dodge(0.9), color = "black") +
  geom_boxplot(width = 0.4, position = position_dodge(0.9), outlier.shape = NA, color = "black") +
  scale_fill_manual(values = c("#4682B4", "#CD3700"), name = "Group") +
  scale_colour_manual(values = c("#4682B4", "#CD3700")) +
  labs(title = "Expression in GSE87466 (Training)", x = "", y = "Expression") +
  stat_pvalue_manual(train_stat.test, size = 3.2, label = "p_label", face = "bold") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
        axis.text.y = element_text(hjust = 0.5, face = "bold", size = 12),
        axis.title = element_text(size = 16, face = "bold"),
        legend.position = "top",
        panel.grid = element_blank()) +
  facet_wrap(~symbol, scales = "free", nrow = 1)

ggsave("01.train_res.pdf", train_va, width = 10, height = 8)
ggsave("01.train_res.png", train_va, width = 10, height = 8, dpi = 600)


# ============================ 3. 验证集处理与绘图 ============================

temp1_ver <- ver_data[rownames(ver_data) %in% hub_gene$symbol, ]

disease_sample_ver <- ver_group[ver_group$group == 'UC', ]
disease_data_ver <- temp1_ver[, colnames(temp1_ver) %in% disease_sample_ver$sample]
disease_mean_ver <- apply(disease_data_ver, 1, mean) %>% as.data.frame()
colnames(disease_mean_ver) <- 'disease_mean'

control_sample_ver <- ver_group[ver_group$group == 'Normal', ]
control_data_ver <- temp1_ver[, colnames(temp1_ver) %in% control_sample_ver$sample]
control_mean_ver <- apply(control_data_ver, 1, mean) %>% as.data.frame()
colnames(control_mean_ver) <- 'control_mean'

ver_mean <- cbind(disease_mean_ver, control_mean_ver)
ver_mean$compare_res <- ifelse(ver_mean$disease_mean > ver_mean$control_mean, 1, 0)
ver_mean$symbol <- rownames(ver_mean)

temp1_ver$symbol <- rownames(temp1_ver)
temp1_ver_long <- pivot_longer(temp1_ver, cols = -symbol, names_to = 'sample', values_to = 'expr')
ver_hub_exp <- merge(temp1_ver_long, ver_group, by = 'sample')
ver_hub_exp$group <- as.factor(ver_hub_exp$group)

ver_stat.test <- ver_hub_exp %>%
  group_by(symbol) %>%
  wilcox_test(expr ~ group) %>%
  adjust_pvalue(method = 'fdr')

ver_stat.test$p_label <- ifelse(ver_stat.test$p < 0.001, "***",
                                ifelse(ver_stat.test$p < 0.01, "**",
                                       ifelse(ver_stat.test$p < 0.05, "*", 'ns')))

ver_finally_res <- merge(ver_stat.test, ver_mean, by = 'symbol')
ver_finally_res <- ver_finally_res %>% select(symbol, p, p_label, compare_res)
ver_finally_res$p_sig <- ifelse(ver_finally_res$p < 0.05, 1, 0)

# 绘图（验证集）
ver_stat.test <- ver_stat.test %>% add_xy_position(x = 'group', dodge = 0.5, step.increase = 0.1)

ver_va <- ggplot(ver_hub_exp, aes(x = group, y = expr, color = group)) +
  geom_violin(trim = FALSE, color = "black", aes(fill = group)) +
  stat_boxplot(geom = "errorbar", width = 0.1, position = position_dodge(0.9), color = "black") +
  geom_boxplot(width = 0.4, position = position_dodge(0.9), outlier.shape = NA, color = "black") +
  scale_fill_manual(values = c("#4682B4", "#CD3700"), name = "Group") +
  scale_colour_manual(values = c("#4682B4", "#CD3700")) +
  labs(title = "Expression in GSE75214 (Validation)", x = "", y = "Expression") +
  stat_pvalue_manual(ver_stat.test, size = 3.2, label = "p_label", face = "bold") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
        axis.text.y = element_text(hjust = 0.5, face = "bold", size = 12),
        axis.title = element_text(size = 16, face = "bold"),
        legend.position = "top",
        panel.grid = element_blank()) +
  facet_wrap(~symbol, scales = "free", nrow = 1)

ggsave("02.verification_res.pdf", ver_va, width = 10, height = 8)
ggsave("02.verification_res.png", ver_va, width = 10, height = 8, dpi = 600)


# ============================ 4. 筛选最终一致的关键基因 ============================
finally_res <- merge(ver_finally_res, train_finally_res, by = 'symbol')
res_keep <- data.frame()

for (i in 1:nrow(finally_res)) {
  if (finally_res$p_sig.x[i] == finally_res$p_sig.y[i] & 
      finally_res$compare_res.x[i] == finally_res$compare_res.y[i] & 
      finally_res$p_sig.y[i] != 0) {
    res_keep <- rbind(res_keep, data.frame(symbol = finally_res$symbol[i], res = 'TRUE'))
  }
}

if (nrow(res_keep) >= 1) {
  print('通过验证，最终基因列表：')
  print(res_keep$symbol)
  hub_gene_final <- data.frame(symbol = res_keep$symbol)
  write.csv(hub_gene_final, 'hub_gene.csv', row.names = FALSE)
} else {
  print('未通过验证')
}


# ============================ 5. 单基因 ROC 曲线绘制（训练集） ============================
# 准备好用于 ROC 的表达矩阵（转置为 行=样本，列=基因）
biomarker <- read.csv("hub_gene.csv", stringsAsFactors = FALSE)

# 定义绘图颜色
color_pal <- c("#FF0000", "#4169E1", "#FFD700", "#8A2BE2", "#FF4500")

# --- 训练集 ROC ---
biomarker_exp_train <- t(train_data[biomarker$symbol, ])
roc_objects_train <- list()
auc_values_train <- c()

for (gene in biomarker$symbol) {
  if (gene %in% colnames(biomarker_exp_train)) {
    roc_obj <- roc(train_group$group, biomarker_exp_train[, gene], levels = c('Normal', 'UC'), quiet = TRUE)
    if (roc_obj$auc > 0.1) {
      roc_objects_train[[gene]] <- roc_obj
      auc_values_train <- c(auc_values_train, round(roc_obj$auc, 3))
    }
  }
}

# 保存训练集 ROC
pdf("03.GSE87466_ROC_all.pdf", family = "Times", width = 6, height = 6)
plot(roc_objects_train[[1]], col = color_pal[1], main = "GSE87466 (Training)", print.auc = FALSE)
if (length(roc_objects_train) > 1) {
  for (i in 2:length(roc_objects_train)) {
    lines(roc_objects_train[[i]], col = color_pal[i])
  }
}
legend_labels <- paste(names(roc_objects_train), " (AUC=", auc_values_train, ")", sep = "")
legend("bottomright", legend = legend_labels, col = color_pal[1:length(roc_objects_train)], lwd = 2)
dev.off()
png("03.GSE87466_ROC_all.png", width = 6, height = 6, units = 'in', res = 600)
plot(roc_objects_train[[1]], col = color_pal[1], main = "GSE87466 (Training)", print.auc = FALSE)
if (length(roc_objects_train) > 1) {
  for (i in 2:length(roc_objects_train)) {
    lines(roc_objects_train[[i]], col = color_pal[i])
  }
}
legend("bottomright", legend = legend_labels, col = color_pal[1:length(roc_objects_train)], lwd = 2)
dev.off()


# --- 验证集 ROC ---
biomarker_exp_ver <- t(ver_data[biomarker$symbol, ])
roc_objects_ver <- list()
auc_values_ver <- c()

for (gene in biomarker$symbol) {
  if (gene %in% colnames(biomarker_exp_ver)) {
    roc_obj <- roc(ver_group$group, biomarker_exp_ver[, gene], levels = c('Normal', 'UC'), quiet = TRUE)
    if (roc_obj$auc > 0.1) {
      roc_objects_ver[[gene]] <- roc_obj
      auc_values_ver <- c(auc_values_ver, round(roc_obj$auc, 3))
    }
  }
}

# 保存验证集 ROC
pdf("04.GSE75214_ROC_all.pdf", family = "Times", width = 6, height = 6)
plot(roc_objects_ver[[1]], col = color_pal[1], main = "GSE75214 (Validation)", print.auc = FALSE)
if (length(roc_objects_ver) > 1) {
  for (i in 2:length(roc_objects_ver)) {
    lines(roc_objects_ver[[i]], col = color_pal[i])
  }
}
legend_labels_ver <- paste(names(roc_objects_ver), " (AUC=", auc_values_ver, ")", sep = "")
legend("bottomright", legend = legend_labels_ver, col = color_pal[1:length(roc_objects_ver)], lwd = 2)
dev.off()
png("04.GSE75214_ROC_all.png", width = 6, height = 6, units = 'in', res = 600)
plot(roc_objects_ver[[1]], col = color_pal[1], main = "GSE75214 (Validation)", print.auc = FALSE)
if (length(roc_objects_ver) > 1) {
  for (i in 2:length(roc_objects_ver)) {
    lines(roc_objects_ver[[i]], col = color_pal[i])
  }
}
legend("bottomright", legend = legend_labels_ver, col = color_pal[1:length(roc_objects_ver)], lwd = 2)
dev.off()

print("✅ 全部代码运行完毕！图已保存至 ./06_expression_verification/")