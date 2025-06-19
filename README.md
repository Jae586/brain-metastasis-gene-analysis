# Brain Metastasis Gene Analysis

This repository contains RNA-seq analysis code and results for studying how extracellular matrix (ECM) proteins affect microglia gene expression during brain metastasis in breast cancer. The work supports research from the Snyder Lab on cancer immunosurveillance in the brain.

## Background

Microglia are immune cells in the brain that help suppress cancer through phagocytosis. However, certain ECM proteins in the brain microenvironment may inhibit this function, allowing tumor cells to evade immune responses.

This project investigates transcriptional changes in microglia cultured on various ECM proteins to identify affected gene pathways.

## Objectives

- Analyze RNA-seq data from microglia treated with ECM proteins
- Identify differentially expressed gene (DEG) groups
- Use Gene Set Enrichment Analysis (GSEA) and DAVID to identify enriched pathways
- Compare experimental groups against control coatings (PLL, Matrigel)

## Data & Tools

- **Input Data:** FASTQ RNA-seq files
- **Controls:** PLL, Matrigel
- **ECM Conditions:** Agrin, Pan-Laminin, Laminin-211, Collagen
- **Tools Used:**
  - `STAR` – alignment to reference genome
  - `HTSeq-count` – generate read count matrix
  - `DESeq2` – differential expression analysis
  - `GSEA` – gene set enrichment analysis
  - `DAVID` – functional annotation

## Workflow

1. **Align RNA-seq data** with STAR using Ensembl reference genome
2. **Generate read counts** using HTSeq-count
3. **Identify DEGs** between control and ECM-treated groups using DESeq2
4. **Rank genes** by log2 fold change
5. **Perform enrichment analysis** with GSEA and DAVID
6. **Visualize results** with enrichment plots and pathway reports



