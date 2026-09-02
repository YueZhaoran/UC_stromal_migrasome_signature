# Computational Analysis Code Description

## Software environment

All analyses were performed in R. The complete software environment, including R version and package versions, is provided in the environment directory.

## Bulk transcriptome analysis

Bulk RNA expression datasets were obtained from GEO accession numbers GSE87466 and GSE75214.
Preprocessing included expression matrix extraction, annotation, probe aggregation, and quality control.

## Single-cell RNA-seq analysis

Single-cell data from GSE214695 were processed using Seurat-based workflows.
Quality control, normalization, dimensional reduction, batch correction, clustering, and cell-type annotation were performed according to the manuscript Methods.

## Genetic-risk localization

GWAS summary statistics from OpenGWAS dataset ieu-a-32 were integrated with single-cell transcriptomic profiles using scPagwas.

## Machine learning analysis

Candidate genes were evaluated using multiple feature-selection algorithms including LASSO regression, SVM-RFE, and random forest.

## Functional and network analyses

GO/KEGG enrichment, GSEA, microenvironment analysis, molecular interaction networks, CellChat-related analyses, and perturbation analyses were performed using the corresponding scripts.

## Reproducibility

Scripts should be executed sequentially after downloading the public datasets and configuring local paths.
