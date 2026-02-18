# Changelog 

Vincent Ohlhauser

## [26/01/18]

### Changed

- 01.QTLseq_Analysis.R : Renamed object df -> SNPset
- 01.QTLseq_Analysis.R : Added file.exists check to pools file
- 01.QTLseq_Analysis.R : Re-structured script: Put helper functions at the top
- 01.QTLseq_Analysis.R : Used Veras values for `filterSNPs` and `runQTLseqAnalysis`
- 01.QTLseq_Analysis.R : Added G' plot
- 01.QTLseq_Analysis.R : Exported significant QTL to csv file sigQTL.csv
- 02_QTL_Annotation[...] : Re-structured script: Put helper functions at the top

## [26/02/17]

### Added
- BSA_incl.annotatio.qmd : Collects Veras code including some more annotation into one document

## [26/02/13]

### Changed
- .gitignore : added path data/processed and objects/
- _install_packaged.R : Fixed install of missing packages by letting pak handle it (line 38)
- _install_packaged.R : Removed R objects after script is finished
- 01.QTLseq_Analysis.R : Removed package loading, is handled by _install_packages.R
- 01.QTLseq_Analysis.R : Changed input paths to read_tsv(file.path("data", "raw", file))
- 01.QTLseq_Analysis.R : Changed read_tsv -> read.delim 
- 01.QTLseq_Analysis.R : Added some more linebreaks
- 01.QTLseq_Analysis.R : Added pdf device to save plots at results/plots including timestamp
- 01.QTLseq_Analysis.R : Added titles to plots to make them more informative
- 02_QTL_Annotation[...] : Removed package loading, is handled by _install_packages.R
- 02_QTL_Annotation[...] : Changed input paths to read_tsv(file.path("data", "raw", file))

## [26/02/10]

### Added
- Files from Vera, added at BSA_manuscript for later use : 
  * BSA_manuscript.rmd
  * phenotype_RIL_II.xlsx
  * phenotyping_script.R

### Changed
- _install_packages.R : Re-arranged packages by source (CRAN/Bioconductor/Github)
- _install_packages.R : Switched to package install by `pak::pkg_install()`, this migh require some fine tuning

## [26/02/09]

### Added
- CHANGELOG.md

### Changed
- .gitignore : Excluded bsa_paper.Rproj from Github