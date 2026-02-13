##############################################
# 01_QTLseq_Analysis.R
# QTL-seq analysis for Drosophila ananassae
##############################################

# ---- 1. Load packages ----
source(file.path("scripts", "_install_packages.R"))
source(file.path("scripts", "_save_objects.R"))

# ---- 2. Load input data ----
pool_Tolerant_Offspring <- read.delim(file.path("data", "raw", "fast_O.table"))
pool_Sensitive_Offspring <- read.delim(file.path("data", "raw", "slow_O.table"))
pool_Tolerant_Parent <- read.delim(file.path("data", "raw", "fast_P.table"))
pool_Sensitive_Parent <- read.delim(file.path("data", "raw", "slow_P.table"))

# ---- 3. Rename columns for clarity ----
rename_columns <- function(df, prefix) {
  names(df)[names(df) == "sample01.AD"] <- paste0(prefix, ".AD")
  names(df)[names(df) == "sample01.DP"] <- paste0(prefix, ".DP")
  names(df)[names(df) == "sample01.GT"] <- paste0(prefix, ".GT")
  return(df)
}

pool_Tolerant_Offspring <- rename_columns(pool_Tolerant_Offspring, "TolerantOffspring")
pool_Sensitive_Offspring <- rename_columns(pool_Sensitive_Offspring, "SensitiveOffspring")
pool_Tolerant_Parent <- rename_columns(pool_Tolerant_Parent, "TolerantParent")
pool_Sensitive_Parent <- rename_columns(pool_Sensitive_Parent, "SensitiveParent")

# ---- 4. Merge offspring pools ----
offspring_pools <- full_join(pool_Sensitive_Offspring, pool_Tolerant_Offspring)
write_tsv(offspring_pools, file.path("data", "processed", "offspring_pools.table"))

# ---- 5. Define analysis parameters ----
file <- file.path("data", "processed", "offspring_pools.table")
Chroms <- c("CM029943.2", "CM029944.2", "CM033063.1", "CM033064.1", "CM029945.2", "CM029946.2")
HighBulk <- "TolerantOffspring"
LowBulk <- "SensitiveOffspring"

# ---- 6. Import SNP data ----
df <- QTLseqr::importFromGATK(file = file, highBulk = HighBulk, 
                     lowBulk = LowBulk, chromList = Chroms)

# ---- 7. Visualize depth and allele frequency distributions ----
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_pdf <- file.path("results", "plots", 
                     paste0("01_qtlseq_plots_", timestamp, ".pdf"))
pdf(out_pdf)

ggplot(df) + geom_histogram(aes(x = DP.HIGH + DP.LOW)) + 
  theme_minimal() + 
  xlim(0, 1000) +
  ggtitle("Distribution of total read depth per SNP (DP.HIGH + DP.LOW)")
ggplot(df) + geom_histogram(aes(x = REF_FRQ)) + 
  theme_minimal()+
  ggtitle("Overall reference allele frequency across bulks (REF_FRQ)")

# ---- 8. Filter SNPs ----
df_filt <- filterSNPs(SNPset = df, 
                      refAlleleFreq = 0.20, 
                      depthDifference = 100, 
                      maxTotalDepth = 400, 
                      verbose = TRUE)
df_nona <- na.omit(df_filt)

# ---- 9. Visualize filtered data ----
ggplot(df_nona) + geom_histogram(aes(x = REF_FRQ)) + 
  theme_minimal()+
  ggtitle("Reference allele frequency after filtering")
ggplot(df_nona) + geom_histogram(aes(x = DP.HIGH + DP.LOW)) + 
  theme_minimal() + 
  xlim(0, 1000)+
  ggtitle("Total read depth after filtering")
ggplot(df_nona) + geom_histogram(aes(x = SNPindex.HIGH)) + 
  theme_minimal()+
  ggtitle("Alternative allele frequency in the HIGH bulk (SNP-index)")
ggplot(df_nona) + geom_histogram(aes(x = SNPindex.LOW)) + 
  theme_minimal()+
  ggtitle("Alternative allele frequency in the LOW bulk (SNP-index)")

# ---- 10. Run QTLseq and G' analyses ----
chromosomes <- c("Chr 2L", "Chr 2R", "Chr 3L", "Chr 3R", "Chr XL", "Chr XR")

qtl_results <- QTLseqr::runQTLseqAnalysis(df_nona, 
                                 windowSize = 1e6, 
                                 popStruc = "RIL", 
                                 bulkSize = 100, 
                                 replications = 1e6, 
                                 intervals = c(95, 99))
qtl_results$CHROM <- factor(qtl_results$CHROM,
                            levels = Chroms,
                            labels = chromosomes)

gprime_results <- QTLseqr::runGprimeAnalysis(df_nona, 
                                    windowSize = 1e6, 
                                    outlierFilter = "deltaSNP", 
                                    filterThreshold = 0.05)
gprime_results$CHROM <- factor(gprime_results$CHROM,
                               levels = Chroms,
                               labels = chromosomes)

# ---- 11. Plot QTL statistics ----
QTLseqr::plotQTLStats(qtl_results, var = "deltaSNP", plotIntervals = TRUE) + 
  theme_minimal() + 
  scale_color_manual(values = c("coral2", "blue"))+
  ggtitle("Tri-cube weighted delta SNP-index")
QTLseqr::plotQTLStats(qtl_results, var = "nSNPs") + 
  theme_minimal() +
  ggtitle("Distribution of SNPs used to calculate G'")

# ---- 12. Significant regions ----
sigRegions_qtl <- QTLseqr::getSigRegions(qtl_results, method = "QTLseq")
sigRegions_gprime <- QTLseqr::getSigRegions(gprime_results, method = "Gprime")

# ---- 13. Save key objects ----
save_object(qtl_results, "qtl_results")
save_object(gprime_results, "gprime_results")
save_object(sigRegions_qtl, "sigRegions_qtl")
save_object(sigRegions_gprime, "sigRegions_gprime")

dev.off()
###EOF