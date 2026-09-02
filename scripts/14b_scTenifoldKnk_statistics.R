##################################################
# Script: 14b_scTenifoldKnk_statistics.R
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

# Quick extraction of scTenifoldKnk stats
# Working directory should be configured from project root

for (gene in c("ADAMTS1", "FBN1", "SPARC", "VCAM1")) {
  cat("/n==========", gene, "==========/n")
  
  # KO results
  ko_file <- paste0(gene, "_df_KO.result.csv")
  if (file.exists(ko_file)) {
    ko <- read.csv(ko_file, stringsAsFactors = FALSE)
    sig <- ko[ko$p.value < 0.05 & !is.na(ko$logFC) & abs(ko$logFC) > 0.5, ]
    cat("Significant DRGs (p < 0.05, |logFC| > 0.5):", nrow(sig), "/n")
    
    # Top 10 up and down
    up <- sig[order(sig$logFC, decreasing = TRUE), ]
    down <- sig[order(sig$logFC, decreasing = FALSE), ]
    cat("Top 10 up:", paste(head(up$gene, 10), collapse = ", "), "/n")
    cat("Top 10 down:", paste(head(down$gene, 10), collapse = ", "), "/n")
  }
  
  # GO results
  go_file <- paste0(gene, "_GO_ALL_results.csv")
  if (file.exists(go_file)) {
    go <- read.csv(go_file, stringsAsFactors = FALSE)
    cat("GO terms:", nrow(go), "/n")
    cat("Top 5 GO:", paste(head(go$Description, 5), collapse = "; "), "/n")
  } else {
    cat("GO: not available/n")
  }
  
  # KEGG results
  kegg_file <- paste0(gene, "_KEGG_results.csv")
  if (file.exists(kegg_file)) {
    kegg <- read.csv(kegg_file, stringsAsFactors = FALSE)
    cat("KEGG pathways:", nrow(kegg), "/n")
    cat("Top 5 KEGG:", paste(head(kegg$Description, 5), collapse = "; "), "/n")
  } else {
    cat("KEGG: not available/n")
  }
}
