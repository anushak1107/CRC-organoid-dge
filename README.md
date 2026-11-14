# CRC-organoid-dge
End-to-end bulk RNA-seq analysis of colorectal cancer organoids investigating extracellular matrix remodeling using DESeq2 and clusterProfiler

# 🧬 Transcriptomic Profiling of Extracellular Matrix Remodeling in CRC Organoids (GSE349730)

[![R Version](https://img.shields.io/badge/R-%3E%3D_4.0.0-blue.svg)](https://www.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-3.18-green.svg)](https://bioconductor.org/)
[![DESeq2](https://img.shields.io/badge/DESeq2-Differential%20Expression-orange.svg)](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## 📌 Project Overview
This repository contains a reproducible, end-to-end bulk RNA-seq analysis pipeline implemented in **R**. The project investigates how extracellular matrix (ECM) composition alters the transcriptional landscape of **colorectal cancer (CRC) patient-derived organoids (PDOs)** using NCBI GEO dataset **GSE349730**[cite: 2].

Specifically, this workflow evaluates the transcriptional reprogramming of organoids cultured in a **Collagen I-enriched extracellular matrix (`MGCOL`)** compared to **standard Matrigel (`MG`)**, while rigorously controlling for inter-patient baseline heterogeneity using multifactorial statistical modeling[cite: 2].

---

## 🧪 Dataset & Experimental Design
* **GEO Accession:** [GSE349730](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE349730)[cite: 2]
* **Sample Size:** 48 RNA-seq libraries derived from 4 distinct patient backgrounds (`P23`, `P25`, `P30`, `P33`)[cite: 2].
* **Experimental Variables:**
  * **Matrix Composition:** Standard Matrigel (`MG`) vs. Collagen I-Enriched Matrigel (`MGCOL`)[cite: 2]
  * **Treatment Conditions:** Untreated Control vs. 5-Fluorouracil (`5FU`) chemotherapy[cite: 2]

---

## ⚙️ Analytical Workflow

```text
1. Raw Count Matrix (data/349730dge.csv: 61,552 genes × 48 samples)
       │
       ▼
2. Metadata Parsing & Low-Count Filtering (≥10 counts in ≥3 samples)
       │
       ▼
3. Multifactorial DESeq2 Modeling (~ Patient + Matrix | Contrast: MGCOL vs MG)
       │
       ▼
4. Exploratory & DEG Visualization (VST PCA, EnhancedVolcano, pheatmap)
       │
       ▼
5. Functional Enrichment (clusterProfiler & org.Hs.eg.db: GO BP/CC/MF & KEGG)

📊 Key Biological Findings & Visualizations
1. Patient Heterogeneity & Batch Effect Correction (PCA)
Principal Component Analysis revealed that inter-patient baseline variation accounts for the primary axes of global transcriptomic variance (PC1: 46%, PC2: 32%). Incorporating Patient into the DESeq2 design formula (~ Patient + Matrix) effectively isolated the true ECM-driven biological signal.

2. Differential Gene Expression (Volcano & Heatmap)
Comparing MGCOL against MG identified statistically robust differentially expressed genes (padj < 0.05), demonstrating coordinated transcriptional shifts in response to collagen enrichment.

3. Focal Adhesion, ECM Remodeling & Hypoxia Adaptation (GO Enrichment)
Over-representation analysis of upregulated genes revealed significant enrichment in extracellular matrix organization, cell-matrix adhesion, focal adhesion complexes, and response to hypoxia, mirroring in vivo stromal desmoplasia and 3D microenvironmental stiffness.

📂 Repository Structure
scripts/01_preprocessing.R — Programmatic metadata extraction and count matrixfiltering.
scripts/02_deseq2_analysis.R — Multifactorial differential expression, VST normalization, PCA, Volcano, and Heatmap generation.
scripts/03_functional_enrichment.R — Ensembl-to-Entrez/Symbol mapping and GO/KEGG enrichment dotplots.
results/ — Exported differential expression tables and pathway enrichment CSVs.
plots/ — High-resolution publication figures.

🚀 Getting Started

Prerequisites
Install R (>= 4.0.0) and the required CRAN/Bioconductor packages:

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "clusterProfiler", "org.Hs.eg.db", "EnhancedVolcano", "pheatmap"))
install.packages(c("dplyr", "tidyr", "ggplot2", "enrichplot"))

Running the Pipeline
Clone the repository and execute the scripts sequentially from the repository root:

git clone [https://github.com/YOUR-USERNAME/crc-organoid-rnaseq-gse349730.git](https://github.com/YOUR-USERNAME/crc-organoid-rnaseq-gse349730.git)
cd crc-organoid-rnaseq-gse349730
Rscript scripts/01_preprocessing.R
Rscript scripts/02_deseq2_analysis.R
Rscript scripts/03_functional_enrichment.R

📜 **License**
This project is licensed under the MIT License — see the LICENSE file for details.
