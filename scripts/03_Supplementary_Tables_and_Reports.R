##############################################
# 03_Supplementary_Tables_and_Reports.R
# Creates publication-ready tables and Excel sheets
##############################################

# ---- 1. Load packages ----
source(file.path("scripts", "_install_packages.R"))
source(file.path("scripts", "_save_objects.R"))

# Helper function to import newest RDS
import_newest_object <- function(name, dir = "objects/"){
  files <- list.files(dir, pattern = name, full.names = TRUE)
  
  if(length(file)==0){stop(paste0("No object named ", name, "found\naborting..."))}
  
  info <- file.info(files)
  newest <- rownames(info)[which.max(info$mtime)]
  
  readRDS(newest)
}

# Helper function to style tables for pdf output
style_table <- function(df, tbl_head){
  grid::grid.newpage() 
  tg <- gridExtra::tableGrob(df, rows = NULL, 
                             theme = gridExtra::ttheme_default(colhead = list(
                               fg_params = list(fontface = "bold", col = "black"),
                               bg_params = list(fill = "white"))))
  
  gridExtra::grid.arrange(tg, top = grid::textGrob(tbl_head,
                                                   gp = grid::gpar(fontsize = 14, 
                                                                   fontface = "bold")), 
                          newpage = FALSE)
}

# ---- 2. Prepare GeneID Excel sheet ----
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

GO_combined_data <- import_newest_object("GO_combined_data")
GeneID_combined_data <- GO_combined_data[, c("Description", "geneID", "Region")]
write_xlsx(GeneID_combined_data, file.path("results", 
                                           "supplementary_tables", 
                                           paste0("GeneID_combined_data_", 
                                                  timestamp,
                                                  ".xlsx")))
save_object(GeneID_combined_data, "GeneID_combined_data")

# ---- 3. Significant QTL table ----
sigQTL <- read_csv(file.path("data", "processed", "sigQTL.csv"))
qtl_table <- sigQTL[, c("qtl", "CHROM", "start", "end", "length", "nSNPs", "avgDeltaSNP")]

# ---- 4. VCF statistics table ---- 
VCF_File_Statistics_Combined <- read_csv(file.path("data", "raw",
                                                   "VCF_File_Statistics_Combined.csv"))

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

# ---- 6. Write additional tables to pdf ----
out_pdf <- file.path("results", "supplementary_tables", 
                     paste0("03_Supplementary_tables_", timestamp, ".pdf"))
pdf(out_pdf, width = 14, height = 8)

style_table(qtl_table, "Table of significant QTL regions")

style_table(VCF_File_Statistics_Combined, 
            "Variant calling summary table")

style_table(QC_rawdata, "Sequencing quality control information")

dev.off()
###EOF