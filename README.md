# qtl_analysis

[![R](https://img.shields.io/badge/R-%3E%3D4.1-blue.svg)](https://www.r-project.org/) 
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This project performs QTL analysis and Gene Ontology (GO) enrichment in *Drosophila*. The pipeline is fully automated and reproducible, including:

1. **QTLseq analysis** on sequencing data.  
2. **Annotation and GO enrichment** of significant regions.  
3. **Generation of supplementary tables and plots** for publications.

---

## Project Structure
```
project_root/
├── data/
│ ├── raw/ # Raw input files
│ │ ├── REFSEQ_FLYBASE_Dana.txt
│ │ ├── dmel_dana_orthologs.xlsx
│ │ ├── genomic.gtf
│ │ ├── VCF_File_Statistics_Combined.csv
│ │ ├── fast_O.table
│ │ ├── slow_O.table
│ │ ├── fast_P.table
│ │ └── slow_P.table
│ └── processed/ # Intermediate processed files
│ └── sigQTL.csv
├── results/
│ ├── plots/ # Plots (QTL maps, GO enrichment)
│ ├── enrichment/ # GO enrichment results
│ └── supplementary_tables/ # Tables for publication
├── objects/ # R objects / intermediate RDS files
├── scripts/ # All R scripts
│ ├── 01_QTLseq_Analysis.R
│ ├── 02_QTL_Annotation_and_GO_Enrichment.R
│ ├── 03_Supplementary_Tables_and_Reports.R
│ └── run_pipeline.R # Master script to run the pipeline
├── create_project_folders.sh # Shell script to create folder structure
└── README.md
```


---

## Installation

The pipeline requires **R ≥ 4.1** and the following packages:

```r
install.packages(c(
  "vcfR", "QTLseqr", "dplyr", "tidyr", "readr", "readxl", 
  "ggplot2", "ggpubr", "GenomicRanges", "rtracklayer", 
  "AnnotationDbi", "org.Dm.eg.db", "clusterProfiler", 
  "enrichR", "openxlsx", "writexl", "kableExtra", "knitr"
))
```
## Setup

- Clone or download the repository.

- Create folder structure:

```chmod +x create_project_folders.sh
./setup_project.sh
```

- Place all raw input files in `data/raw/`.

- Ensure intermediate files, such as `sigQTL.csv`, are in `data/processed/`.

## Running the Pipeline

From the project root, launch R and run:

```r
source("scripts/run_pipeline.R")
```
This executes the following steps in order:

1. QTLseq analysis (`01_QTLseq_Analysis.R`)

2. QTL annotation and GO enrichment (`02_QTL_Annotation_and_GO_Enrichment.R`)

3. Supplementary tables and reports generation (`03_Supplementary_Tables_and_Reports.R`)

All results will be saved in the `results/` directory.

## Outputs

- Plots: `results/plots/`

  - QTL maps

  - SNP density plots

  - GO enrichment bar and dot plots

- GO enrichment tables: `results/enrichment/`

- Supplementary tables: `results/supplementary_tables/`

  - Excel sheets for gene IDs

  - HTML tables for significant QTLs and VCF statistics

  - Raw sequencing QC tables

## Citation

> **Dissecting Cold Tolerance in *Drosophila ananassae***
> 
> *Yılmaz VM, Kara FT, Grath S*. (2025). Dissecting Cold Tolerance in *Drosophila ananassae*. *bioRxiv*. https://doi.org/10.1101/2025.04.23.650207v1

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

