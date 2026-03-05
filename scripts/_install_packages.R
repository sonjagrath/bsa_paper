##############################################
# _install_packages.R
# Ensures that all required R packages are installed and loaded
# Used by all analysis scripts in this project
##############################################

# Check if the required build tools are installed
if(!pkgbuild::has_build_tools()){
  stop("Build tools not fount!\nWindows: Install Rtools\nmacOS: Run xcode-select --install.\n\n")
}

# required packages
required_packages <- c(
  ## CRAN packages
  "tidyverse", # Data handling and visualization
  "ggpubr", # visualization
  "ggh4x",
  "ggtext",
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

# Check which packages are missing
required_pkgs <- required_packages
required_pkgs[length(required_pkgs)] <- "QTLseqr"
missing_pkgs <- required_pkgs[!vapply(required_pkgs,
                                          requireNamespace,
                                          logical(1), quietly = TRUE)]

# Install missing dependencies
if (length(missing_pkgs) > 0) {
  pak::pkg_install(required_packages)
}

# Load all required packages
invisible(lapply(required_pkgs, library, character.only = TRUE))

message("✅   All required packages are installed and loaded.")

# remove handling objects
rm(required_pkgs, required_packages, missing_pkgs)

###EOF