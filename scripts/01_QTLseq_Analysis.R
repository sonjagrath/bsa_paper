##############################################
# 01_QTLseq_Analysis.R
# QTL-seq analysis for Drosophila ananassae
##############################################

# ---- 1. Load packages define helper functions----
source(file.path("scripts", "_install_packages.R"))
source(file.path("scripts", "_save_objects.R"))

rename_columns <- function(df, prefix) {
  names(df)[names(df) == "sample01.AD"] <- paste0(prefix, ".AD")
  names(df)[names(df) == "sample01.DP"] <- paste0(prefix, ".DP")
  names(df)[names(df) == "sample01.GT"] <- paste0(prefix, ".GT")
  return(df)
}

# ---- 2. Create pool file if it does not exist ----
pool_file <- file.path("data", "processed", "offspring_pools.table")

if(!file.exists(pool_file)){
  # Import parent and offspring pools
  pool_Tolerant_Offspring <- read.delim(file.path("data", "raw", "fast_O.table"))
  pool_Sensitive_Offspring <- read.delim(file.path("data", "raw", "slow_O.table"))
  pool_Tolerant_Parent <- read.delim(file.path("data", "raw", "fast_P.table"))
  pool_Sensitive_Parent <- read.delim(file.path("data", "raw", "slow_P.table"))

  # Rename columns for clarity
  pool_Tolerant_Offspring <- rename_columns(pool_Tolerant_Offspring, "TolerantOffspring")
  pool_Sensitive_Offspring <- rename_columns(pool_Sensitive_Offspring, "SensitiveOffspring")
  pool_Tolerant_Parent <- rename_columns(pool_Tolerant_Parent, "TolerantParent")
  pool_Sensitive_Parent <- rename_columns(pool_Sensitive_Parent, "SensitiveParent")

  # Merge offspring pools
  a <- full_join(pool_Tolerant_Offspring, pool_Tolerant_Parent)
  b <- full_join(pool_Sensitive_Offspring, pool_Sensitive_Parent)
  offspring_pools <- full_join(a, b)
  message(paste0("writing offspring pools file to ", pool_file))
  write_tsv(offspring_pools, pool_file)
}

# ---- 3. Define analysis parameters ----
Chroms <- c("CM029943.2", "CM029944.2", "CM033063.1", "CM033064.1", "CM029945.2", "CM029946.2")
HighBulk <- "TolerantOffspring"
LowBulk <- "SensitiveOffspring"

# ---- 4. Import SNP data ----
SNPset <- QTLseqr::importFromGATK(file = pool_file, highBulk = HighBulk, 
                     lowBulk = LowBulk, chromList = Chroms)

# ---- 5. Visualize depth and allele frequency distributions ----
p1 <- ggplot(SNPset) + geom_histogram(aes(x = DP.HIGH + DP.LOW)) + 
  theme_minimal() + 
  xlim(0, 1000) +
  theme(panel.grid = element_blank())+
  ggtitle("Distribution of total read depth per SNP (DP.HIGH + DP.LOW)")
p2 <- ggplot(SNPset) + geom_histogram(aes(x = REF_FRQ)) + 
  theme_minimal()+
  theme(panel.grid = element_blank())+
  ggtitle("Overall reference allele frequency across bulks (REF_FRQ)")

# ---- 6. Filter SNPs ----
SNPset_filt <- QTLseqr::filterSNPs(SNPset = SNPset, 
                      refAlleleFreq = 0.30, 
                      depthDifference = 50, minTotalDepth = 40,
                      maxTotalDepth = 400, 
                      verbose = TRUE)%>%
  na.omit()

# ---- 7. Visualize filtered data ----
p3 <- ggplot(SNPset_filt) + geom_histogram(aes(x = REF_FRQ)) + 
  theme_minimal()+
  theme(panel.grid = element_blank())+
  ggtitle("Reference allele frequency after filtering")
p4 <- ggplot(SNPset_filt) + geom_histogram(aes(x = DP.HIGH + DP.LOW)) + 
  theme_minimal() + 
  xlim(0, 1000)+
  theme(panel.grid = element_blank())+
  ggtitle("Total read depth after filtering")
p5 <- ggplot(SNPset_filt) + geom_histogram(aes(x = SNPindex.HIGH)) + 
  theme_minimal()+
  theme(panel.grid = element_blank())+
  ggtitle("Alternative allele frequency in the HIGH bulk (SNP-index)")
p6 <- ggplot(SNPset_filt) + geom_histogram(aes(x = SNPindex.LOW)) + 
  theme_minimal()+
  theme(panel.grid = element_blank())+
  ggtitle("Alternative allele frequency in the LOW bulk (SNP-index)")

# ---- 8. Run QTLseq and G' analyses ----
chromosomes <- c("Chr 2L", "Chr 2R", "Chr 3L", "Chr 3R", "Chr XL", "Chr XR")

RNGkind("Mersenne-Twister", "Inversion")
set.seed(125)
qtl_results <- QTLseqr::runQTLseqAnalysis(SNPset_filt, 
                                 windowSize = 1e6, 
                                 popStruc = "RIL", 
                                 bulkSize = 100, 
                                 intervals = c(95, 99))

qtl_results$CHROM <- factor(qtl_results$CHROM,
                            levels = Chroms,
                            labels = chromosomes)

gprime_results <- QTLseqr::runGprimeAnalysis(SNPset_filt, 
                                    windowSize = 1e6, 
                                    outlierFilter = "deltaSNP", 
                                    filterThreshold = 0.05)
gprime_results$CHROM <- factor(gprime_results$CHROM,
                               levels = Chroms,
                               labels = chromosomes)

# ---- 9. Plot QTL statistics ----
p7 <- QTLseqr::plotQTLStats(qtl_results, var = "deltaSNP", plotIntervals = TRUE) + 
  theme_minimal() + 
  scale_color_manual(values = c("darkgrey", "red"))+
  ggtitle("deltaSNP-index plot across the genome")
p8 <- QTLseqr::plotQTLStats(qtl_results, var = "nSNPs") + 
  theme_minimal() +
  ggtitle("SNP density across genome within 1 Mb windows")
p9 <- QTLseqr::plotQTLStats(gprime_results, var = "Gprime", plotThreshold = TRUE) + 
  theme_minimal() +
  ggtitle("G` value across the genome within 1 Mb windows")

# ---- 10. Significant regions ----
sigRegions_qtl <- QTLseqr::getSigRegions(qtl_results, method = "QTLseq")
sigRegions_gprime <- QTLseqr::getSigRegions(gprime_results, method = "Gprime")

# ---- 11. Save key objects and plots  ----
save_object(qtl_results, "qtl_results")
save_object(gprime_results, "gprime_results")
save_object(sigRegions_qtl, "sigRegions_qtl")
save_object(sigRegions_gprime, "sigRegions_gprime")

qtl_path <- file.path("data", "processed", "sigQTL.csv")
getQTLTable(qtl_results, method = "QTLseq", export = TRUE, fileName = qtl_path)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_pdf <- file.path("results", "plots", 
                     paste0("01_qtlseq_plots_", timestamp, ".pdf"))
pdf(out_pdf, width = 10, height = 5)

print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
print(p7)
print(p8)
print(p9)

dev.off()
###EOF