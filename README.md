# Bulk-Segregant Analysis for Cold Tolerance in _Drosophila ananassae_

[![R](https://img.shields.io/badge/R-%3E%3D4.1-blue.svg)](https://www.r-project.org/) 
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## About

This repository contains the complete workflow for QTL-seq analysis of *Drosophila ananassae*. It includes data processing, QTL identification, gene annotation, GO enrichment analysis, and generation of publication-ready tables and plots. All scripts are modular and reproducible, supporting the analysis presented in the associated manuscript.

---

## Citation

If you use this pipeline, please cite:

> **Dissecting Cold Tolerance in *Drosophila ananassae***
> 
> *Yılmaz VM, Kara FT, Grath S*. (2025). Dissecting Cold Tolerance in *Drosophila ananassae*.
> 
> *bioRxiv*. https://doi.org/10.1101/2025.04.23.650207v1

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
│ ├── plots/ # QTL maps, GO enrichment plots
│ ├── enrichment/ # GO enrichment results
│ └── supplementary_tables/ # Tables for publication
├── objects/ # R objects / intermediate RDS files
├── scripts/ # All R scripts
│ ├── 01_QTLseq_Analysis.R
│ ├── 02_QTL_Annotation_and_GO_Enrichment.R
│ ├── 03_Supplementary_Tables_and_Reports.R
│ ├── _install_packages.R # Helper script to install required R packages
│ ├── _save_objects.R # Helper script to save R objects
│ └── run_pipeline.R # Master script to run the pipeline
├── setup_project.sh # Shell script to create folder structure
└── README.md
```

---

## Installation

1. Clone the repository:
   
```bash
git clone https://github.com/sonjagrath/bsa_paper.git
```
2. Run the `_install_packages.R` script to ensure all required R packages are installed:

```r
source("scripts/_install_packages.R")
```
   
---

## Usage

```r
source("scripts/run_pipeline.R")
```

This will:

- Process raw data files

- Perform QTL-seq and G' analysis

- Annotate significant QTL regions

- Run GO enrichment analysis

- Generate publication-ready plots and tables

- Save all intermediate R objects in `objects/`

## Output

- Plots: `results/plots/`

  - QTL maps

  - SNP density plots

  - GO enrichment bar and dot plots

- GO enrichment tables: `results/enrichment/`

- Supplementary tables: `results/supplementary_tables/`

  - Excel sheets for gene IDs

  - HTML tables for significant QTLs and VCF statistics

  - Raw sequencing QC tables



## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

