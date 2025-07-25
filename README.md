# Brain Metastasis Gene Expression Analysis

This project analyzes how extracellular matrix (ECM) proteins influence gene expression in human microglia during breast cancer brain metastasis. It performs RNA-seq preprocessing, differential gene expression analysis (DESeq2), and pathway enrichment analysis (GSEA) to identify changes in immune and metabolic signaling.

## Background

Microglia are brain-resident immune cells. Certain ECM proteins have been shown to suppress microglial phagocytosis of tumor cells. This project investigates whether those effects correspond to specific transcriptomic changes.

## Data and Tools

- Input: HTSeq-count outputs merged into `combined_counts.csv`
- Gene ID Mapping: MyGene (Ensembl → HGNC symbols)
- DEG Analysis: pyDESeq2 (Python implementation of DESeq2)
- Enrichment Analysis: GSEA via gseapy (using KEGG 2016)
- Visualization: Bubble plots of top enriched pathways

## Visual Results

The following plots show GSEA results for each ECM treatment. Bubble size indicates gene set coverage; color represents up/downregulation.

<table>
  <tr>
    <td><strong>Agrin</strong></td>
    <td><strong>Pan-Laminin</strong></td>
  </tr>
  <tr>
    <td><img src="figures/Agrin_GSEA.webp" width="400"/></td>
    <td><img src="figures/Pan-Laminin_GSEA.webp" width="400"/></td>
  </tr>
  <tr>
    <td><strong>Laminin-211</strong></td>
    <td><strong>Collagen I</strong></td>
  </tr>
  <tr>
    <td><img src="figures/Laminin-211_GSEA.webp" width="400"/></td>
    <td><img src="figures/CollagenI_GSEA.webp" width="400"/></td>
  </tr>
</table>

## Repository Structure

```
├── data/                       # Raw and processed input data
│   └── combined_counts.csv
├── scripts/                    # Python pipeline
│   └── pipeline.py
├── HTSeq_count_results/                    # Output GSEA summary tables
│   └── 3L_output.csv, etc.
├── figures/                    # Bubble plots per condition
│   └── Agrin_GSEA.webp, etc.
├── docs/                       # Final project report
│   └── final_report.pdf
├── README.md                   # This file
```

## Project Report

Full details on experimental motivation, data handling, and statistical results are available in `docs/final_report.pdf`.

## Visual Results

Enrichment plots for each ECM condition show pathway-level differences in microglial gene expression. Bubble color indicates up/downregulation; size represents the proportion of genes affected.

- Agrin
- Pan-Laminin
- Laminin-211
- Collagen I

Plots are located in the `figures/` folder.

## Reproducibility

To rerun the pipeline:

1. Place your HTSeq-count merged matrix as `combined_counts.csv` in `data/`
2. Run the analysis:
   ```
   python scripts/pipeline.py
   ```

Make sure all dependencies (pyDESeq2, gseapy, pandas, etc.) are installed.

## Authors

- Autumn Davis
- Jasdeep Singh
- Jed Pagcaliwagan
- Vivian White

Conducted as part of CSCI 474: Bioinformatics, Western Washington University, Fall 2024.


