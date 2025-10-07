##############################################
# _install_packages.R
# Ensures that all required R packages are installed and loaded
# Used by all analysis scripts in this project
##############################################

required_packages <- c(
  # Core data handling
  "dplyr",
  "tidyr",
  "readr",
  "readxl",

  # Statistical and genomic analysis
  "vcfR",
  "QTLseqr",
  "GenomicRanges",
  "rtracklayer",
  "AnnotationDbi",
  "org.Dm.eg.db",
  "clusterProfiler",
  "enrichR",

  # Visualization
  "ggplot2",
  "ggpubr",

  # Table and export utilities
  "writexl",
  "openxlsx",
  "kableExtra",
  "knitr"
)

# Function to install missing packages
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing missing package:", pkg))
    install.packages(pkg, dependencies = TRUE)
  }
}

# Install missing packages
invisible(lapply(required_packages, install_if_missing))

# Load all required packages
invisible(lapply(required_packages, library, character.only = TRUE))

message("✅ All required packages are installed and loaded.")

