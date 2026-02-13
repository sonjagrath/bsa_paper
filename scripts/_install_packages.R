##############################################
# _install_packages.R
# Ensures that all required R packages are installed and loaded
# Used by all analysis scripts in this project
##############################################

required_packages <- c(
  ## CRAN packages
  "tidyverse", # Data handling and visualization
  "ggpubr", # visualization
  "writexl", # Import and export utilities
  "readxl",
  "openxlsx",
  "kableExtra", # Visualization of tables
  "knitr", # Report generation
  "vcfR", # manipulation of variant call format (VCF) data
  "enrichR", # R interface to enrichr databases
  
  # Bioconductor packages 
  "rtracklayer", # R interface to genome browsers
  "AnnotationDbi", # R interface for SQLite annotations
  "GenomicRanges", # Manipulate genomic intervals
  "org.Dm.eg.db", # Drosophila melanogaster annotation
  "clusterProfiler", # Analyse and visualize functional profiles
  
  # GitHub packages
  "bmansfeld/QTLseqr" #  QTL mapping
)

# Install package manager pak for installation handling
options(repos = c(CRAN = "https://cloud.r-project.org"))
if(!requireNamespace("pak", quietly = TRUE)){
  message(paste("Installing pak for library handling"))
  install.packages("pak", dependencies = TRUE)
}

# Install missing dependencies
if (any(pak::pkg_status(required_packages)$diff != "OK")) {
  pak::pkg_install(required_packages)
}

# Load all required packages
required_packages[length(required_packages)] <- "QTLseqr"
invisible(lapply(required_packages, library, character.only = TRUE))

message("✅   All required packages are installed and loaded.")

# remove handling objects
rm(required_packages)

###EOF