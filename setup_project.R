###############################################################################
# setup_project.R
#
# Description:
#   Automatically creates the full folder structure for the QTLseq + GO 
#   enrichment project starting from the project root directory, including
#   a folder for R scripts.
#
# Usage:
#   Rscript setup_project.sh
###############################################################################

# Define project root
args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grep(file_arg, args)])
PROJECT_ROOT <- dirname(normalizePath(script_path))

setwd(PROJECT_ROOT)

message(paste0("Creating project folder structure in:", PROJECT_ROOT))

# Create main directories
dir.create("data/raw")
dir.create("data/processed")
dir.create("results")
dir.create("results/plots")
dir.create("results/supplementary_tables")
dir.create("results/logs")
dir.create("objects")
dir.create("scripts")
