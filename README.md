# RNAseq-Data-Analysis_Major-Depression-Disorder
Transcriptomic analysis of ventral white matter in MDD using RNA-seq. Includes DEG analysis (DESeq2), pathway enrichment, PPI networks, drug repurposing via CMap, and molecular docking to identify potential therapeutic targets.

# 🧠 Transcriptomic and Drug Repurposing Analysis in Major Depressive Disorder (MDD)

This repository contains code and documentation for an integrative analysis of RNA-seq data from ventral white matter samples in individuals with Major Depressive Disorder (MDD). The study aims to identify differentially expressed genes (DEGs), functionally annotate them, and explore potential drug repurposing opportunities.

## 📊 Project Overview

- **Dataset:** GSE248260 (NCBI GEO)
- **Tissue:** Ventral white matter (Brodmann Area 47)
- **Samples:** 15 MDD, 9 Control
- **Analysis Pipeline:**
  - RNA-seq preprocessing and normalization (TPM, RPKM, VST)
  - DEG analysis using DESeq2
  - Functional enrichment via DAVID, Enrichr, Reactome
  - Network mapping using STRING and Cytoscape
  - Drug repurposing via Connectivity Map (CLUE)
  - Molecular docking using AutoDock Vina and HDOCK

## 🔬 Key Findings

- Significant DEGs: *SERPINE1*, *FAT2*, *SCML2P2*, etc.
- Enriched pathways: Neuroinflammation, synaptic signaling, oxidative stress
- Top drug candidates: **JWH-015**, **RG-108**, **Bortezomib**
- Promising targets: SIRT1, DNMT1, CB2, ROCK, Proteasome

## 📁 Repository Structure

```
├── data/                # Input RNA-seq and metadata
├── scripts/             # R scripts for DESeq2, plots, enrichment, etc.
├── docking/             # Ligand files, target PDBs, docking output
├── results/             # Output plots, DEGs, pathway results
└── README.md
```

## 🚀 How to Run

1. Clone the repo:
   ```
   git clone https://github.com/yourusername/mdd-transcriptome-analysis.git
   ```
2. Open RStudio and run the scripts in `/scripts/` step-by-step.
3. Use `docking/` to perform ligand-target docking (AutoDock Vina / HDOCK).

## 🧪 Tools & Packages

- `DESeq2`, `ggplot2`, `EnhancedVolcano`, `clusterProfiler`
- `STRINGdb`, `Cytoscape`, `AutoDock Vina`, `HDOCK`
- Connectivity Map (https://clue.io)

## 📚 Citation

If you use this work, please cite:
> Nair, N. et al. (2025). Integrative RNA-seq and Drug Repurposing Analysis for Major Depressive Disorder. *Unpublished project*.


