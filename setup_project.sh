#!/bin/bash

###############################################################################
# setup_project.sh
#
# Description:
#   Automatically creates the full folder structure for the QTLseq + GO 
#   enrichment project starting from the project root directory, including
#   a folder for R scripts.
#
# Usage:
#   chmod +x setup_project.sh
#   ./setup_project.sh
###############################################################################

# Define root directory (current directory)
PROJECT_ROOT=$(pwd)

echo "Creating project folder structure in: $PROJECT_ROOT"

# Create main directories
mkdir -p "$PROJECT_ROOT/data/raw"
mkdir -p "$PROJECT_ROOT/data/processed"
mkdir -p "$PROJECT_ROOT/results/plots"
mkdir -p "$PROJECT_ROOT/results/enrichment"
mkdir -p "$PROJECT_ROOT/results/supplementary_tables"
mkdir -p "$PROJECT_ROOT/objects"
mkdir -p "$PROJECT_ROOT/scripts"

# Optional: create empty README or placeholder files
# touch "$PROJECT_ROOT/README.md"

echo "✅ Folder structure created successfully!"
echo ""
echo "Structure:"
tree -L 3 "$PROJECT_ROOT"

