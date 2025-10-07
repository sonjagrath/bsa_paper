##############################################
# 01_QTLseq_Analysis.R
# QTL-seq analysis for Drosophila ananassae
##############################################

# ---- 1. Load packages ----
source(file.path("scripts", "_install_packages.R"))
source(file.path("scripts", "_save_objects.R"))

library(vcfR)
library(QTLseqr)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggpubr)

# ---- 2. Load input data ----
pool_Tolerant_Offspring <- read_tsv(file.path("data", "fast_O.table"))
pool_Sensitive_Offspring <- read_tsv(file.path("data", "slow_O.table"))
pool_Tolerant_Parent <- read_tsv(file.path("data", "fast_P.table"))
pool_Sensitive_Parent <- read_tsv(file.path("data", "slow_P.table"))

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
df <- importFromGATK(file = file, highBulk = HighBulk, lowBulk = LowBulk, chromList = Chroms)

# ---- 7. Visualize depth and allele frequency distributions ----
ggplot(df) + geom_histogram(aes(x = DP.HIGH + DP.LOW)) + theme_minimal() + xlim(0, 1000)
ggplot(df) + geom_histogram(aes(x = REF_FRQ)) + theme_minimal()

# ---- 8. Filter SNPs ----
df_filt <- filterSNPs(SNPset = df, refAlleleFreq = 0.20, depthDifference = 100, maxTotalDepth = 400, verbose = TRUE)
df_nona <- na.omit(df_filt)

# ---- 9. Visualize filtered data ----
ggplot(df_nona) + geom_histogram(aes(x = REF_FRQ)) + theme_minimal()
ggplot(df_nona) + geom_histogram(aes(x = DP.HIGH + DP.LOW)) + theme_minimal() + xlim(0, 1000)
ggplot(df_nona) + geom_histogram(aes(x = SNPindex.HIGH)) + theme_minimal()
ggplot(df_nona) + geom_histogram(aes(x = SNPindex.LOW)) + theme_minimal()

# ---- 10. Run QTLseq and G' analyses ----
qtl_results <- runQTLseqAnalysis(df_nona, windowSize = 1e6, popStruc = "RIL", bulkSize = 100, replications = 1e6, intervals = c(95, 99))
qtl_results$CHROM <- factor(qtl_results$CHROM,
                            levels = Chroms,
                            labels = c("Chr 2L", "Chr 2R", "Chr 3L", "Chr 3R", "Chr XL", "Chr XR"))

gprime_results <- runGprimeAnalysis(df_nona, windowSize = 1e6, outlierFilter = "deltaSNP", filterThreshold = 0.05)
gprime_results$CHROM <- factor(gprime_results$CHROM,
                               levels = Chroms,
                               labels = c("Chr 2L", "Chr 2R", "Chr 3L", "Chr 3R", "Chr XL", "Chr XR"))

# ---- 11. Plot QTL statistics ----
plotQTLStats(qtl_results, var = "deltaSNP", plotIntervals = TRUE) + theme_minimal() + scale_color_manual(values = c("coral2", "blue"))
plotQTLStats(qtl_results, var = "nSNPs") + theme_minimal()

# ---- 12. Significant regions ----
sigRegions_qtl <- getSigRegions(qtl_results, method = "QTLseq")
sigRegions_gprime <- getSigRegions(gprime_results, method = "Gprime")

# ---- 13. Save key objects ----
save_object(qtl_results, "qtl_results")
save_object(gprime_results, "gprime_results")
save_object(sigRegions_qtl, "sigRegions_qtl")
save_object(sigRegions_gprime, "sigRegions_gprime")

