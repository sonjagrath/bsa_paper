##############################################
# 03_Supplementary_Tables_and_Reports.R
# Creates publication-ready tables and Excel sheets
##############################################

# ---- 1. Load packages ----
source(file.path("scripts", "_install_packages.R"))
source(file.path("scripts", "_save_objects.R"))

library(writexl)
library(readr)
library(kableExtra)
library(knitr)

# ---- 2. Prepare GeneID Excel sheet ----
GeneID_combined_data <- GO_combined_data[, c("Description", "geneID")]
write_xlsx(GeneID_combined_data, file.path("results", "supplementary_tables", "GeneID_combined_data.xlsx"))
save_object(GeneID_combined_data, "GeneID_combined_data")

# ---- 3. Significant QTL table ----
chromosome_lookup <- c(
  "NC_057927.1" = "Chr 2L",
  "NC_057928.1" = "Chr 2R",
  "NC_057929.1" = "Chr 3L",
  "NC_057930.1" = "Chr 3R",
  "NC_057931.1" = "Chr XL",
  "NC_057932.1" = "Chr XR"
)
sigQTL$CHROM <- chromosome_lookup[sigQTL$CHROM]

qtl_table <- sigQTL[, c("CHROM", "start", "end", "length", "nSNPs", "avgDeltaSNP")]

sig_qtl_table <- qtl_table %>%
  kable("html", caption = "Significant QTL Regions") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                full_width = FALSE, position = "center")
save_object(sig_qtl_table, "sig_qtl_table")

# ---- 4. VCF statistics table ----
VCF_File_Statistics_Combined <- read_csv(file.path("data", "VCF_File_Statistics_Combined.csv"))
table_vcf <- kbl(VCF_File_Statistics_Combined, booktabs = TRUE, caption = "VCF File Statistics") %>%
  kable_styling(full_width = FALSE, position = "center", font_size = 12) %>%
  column_spec(1, bold = TRUE)
save_object(table_vcf, "VCF_File_Statistics_Combined_table")

# ---- 5. Raw sequencing QC metrics ----
QC_rawdata <- data.frame(
  Sample = c("slow O", "slow O", "slow O", "fast O", "fast O"),
  Library_Flowcell_Lane = c("EKDN230016689 1A HF3MMDSX7 L1",
                            "EKDN230016689 1A HF37MDSX7 L3",
                            "EKDN230016689 1A HF35WDSX7 L3",
                            "EKDN230016690 1A HF3WTDSX7 L2",
                            "EKDN230016690 1A HF7TMDSX7 L3"),
  Raw_reads = c(14739612, 5796120, 100075174, 89512044, 37727344),
  Effective = c(98.78, 98.83, 98.84, 98.28, 98.43),
  Error = c(0.03, 0.03, 0.03, 0.03, 0.03),
  Q20 = c(95.87, 96.98, 97.46, 96.43, 97.12),
  Q30 = c(89.85, 92.14, 92.98, 91.08, 92.37),
  GC = c(45.14, 45.36, 45.00, 43.67, 43.72)
)

QC_raw_data_table <- kable(QC_rawdata, "html",
                           col.names = c("Sample", "Library Flowcell Lane", "Raw reads",
                                         "Effective (%)", "Error (%)", "Q20 (%)", "Q30 (%)", "GC (%)")) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
  column_spec(1, bold = TRUE) %>%
  add_header_above(c(" " = 1, "Details" = 2, "Quality Metrics" = 5))

save_object(QC_raw_data_table, "QC_raw_data_table")

