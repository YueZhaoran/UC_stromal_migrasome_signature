# UC_stromal_migrasome_signature

Analysis repository accompanying:

**Integrating genetic risk and single-cell transcriptomics identifies a stromal migrasome-related signature in ulcerative colitis**

## Overview

This repository contains the R workflows used for integration of:
- bulk transcriptomic profiles,
- single-cell RNA sequencing,
- ulcerative colitis GWAS summary statistics,
- stromal cell state analysis,
- machine learning model construction,
- functional and network analyses.

## Workflow

01 Bulk transcriptome preprocessing  
02 Single-cell preprocessing and annotation  
03 scPagwas genetic-risk localization  
04 Migrasome activity scoring  
05 Scissor analysis  
06 Differential expression analysis  
07 Candidate gene integration  
08 Machine learning signature construction  
09 Independent validation  
10 Functional enrichment  
11 GSEA  
12 Microenvironment analysis  
13 Molecular network analysis  
14 scTenifoldKnk perturbation analysis  
15 GraphBAN-based compound prioritization  

## Data sources

GEO:
- GSE87466
- GSE75214
- GSE214695

OpenGWAS:
- ieu-a-32

Raw public datasets are not redistributed.

## Reproducibility

Before running:

1. Install required R packages.
2. Download public datasets.
3. Place small annotation files in resources/.
4. Configure local paths.
5. Run scripts sequentially.

## Repository structure

scripts/      R analysis workflows  
data/         dataset instructions  
resources/    gene lists and annotations  
environment/  software versions  
results/      output descriptions  

