##############################################
# 02_QTL_Annotation_and_GO_Enrichment.R
# Annotation and Gene Ontology enrichment analysis
##############################################

# ---- 1. Load packages ----
source(file.path("scripts", "_install_packages.R"))
library(vcfR)
library(QTLseqr)
library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(ggplot2)
library(ggpubr)
library(GenomicRanges)
library(rtracklayer)
library(AnnotationDbi)
library(org.Dm.eg.db)
library(clusterProfiler)
library(enrichR)
library(openxlsx)
source(file.path("scripts", "_save_objects.R"))

# ---- 2. Load data ----
refseq_flybase <- read_tsv(file.path("data", "REFSEQ_FLYBASE_Dana.txt"))
dmel_dana_ortho <- read_excel(file.path("data", "dmel_dana_orthologs.xlsx"))
annotation <- import(file.path("data", "genomic.gtf"))
sigQTL <- read_csv(file.path("data", "sigQTL.csv"))

# ---- 3. Prepare ortholog mapping ----
orthologs <- AnnotationDbi::select(org.Dm.eg.db,
                                   keys = dmel_dana_ortho$ortholog,
                                   keytype = "FLYBASE",
                                   columns = c("FLYBASE", "ENTREZID")) %>%
  right_join(dmel_dana_ortho, by = c("FLYBASE" = "ortholog"))

# ---- 4. Map chromosomes ----
chromosome_mapping <- c(
  "Chr 2L" = "NC_057927.1",
  "Chr 2R" = "NC_057928.1",
  "Chr 3L" = "NC_057929.1",
  "Chr 3R" = "NC_057930.1",
  "Chr XL" = "NC_057931.1",
  "Chr XR" = "NC_057932.1"
)
sigQTL <- sigQTL %>%
  mutate(CHROM = recode(CHROM, !!!chromosome_mapping))

# ---- 5. Separate positive and negative QTL regions ----
pos_QTL <- sigQTL[which(sigQTL$avgDeltaSNP > 0), 1:4]
neg_QTL <- sigQTL[which(sigQTL$avgDeltaSNP < 0), 1:4]

# ---- 6. Function to process QTL regions and annotate genes ----
process_qtl <- function(sigQTL_subset, annotation, refseq_flybase, orthologs, output_file) {
  
  granges_sig_regions <- makeGRangesFromDataFrame(sigQTL_subset[,-2], keep.extra.columns = TRUE)
  overlaps <- findOverlaps(granges_sig_regions, annotation)
  overlapping_genes <- annotation[subjectHits(overlaps)]
  
  unique_transcript_ids <- as.data.frame(unique(overlapping_genes@elementMetadata@listData[["transcript_id"]]))
  colnames(unique_transcript_ids) <- "transcript_id"
  
  filtered_transcripts <- unique_transcript_ids %>%
    separate(transcript_id, c("REFSEQ_TYPE", "REFSEQ_ID"), "_") %>%
    filter(REFSEQ_TYPE == "XM") %>%
    unite("REFSEQ", REFSEQ_TYPE, REFSEQ_ID, sep = "_") %>%
    mutate(REFSEQ = sub("\\.\\d+$", "", REFSEQ))
  
  annotated_transcripts <- merge(x = filtered_transcripts, y = refseq_flybase, by = "REFSEQ")
  colnames(annotated_transcripts) <- c("REFSEQ", "Dana_ID")
  
  annotated_orthologs <- merge(x = annotated_transcripts, y = orthologs, by = "Dana_ID")
  
  write_tsv(annotated_orthologs, output_file)
  return(annotated_orthologs)
}

# ---- 7. Annotate QTLs ----
orthologs_positive <- process_qtl(pos_QTL, annotation, refseq_flybase, orthologs, file.path("results", "enrichment", "positive_annotated_orthologs.tsv"))
orthologs_negative <- process_qtl(neg_QTL, annotation, refseq_flybase, orthologs, file.path("results", "enrichment", "negative_annotated_orthologs.tsv"))

# ---- 8. GO enrichment analysis ----
perform_GO_analysis <- function(result_data, org_db, region_label) {
  gseGO_BP <- perform_enrichGO(result_data$ENTREZID, "BP", org_db)
  gseGO_MF <- perform_enrichGO(result_data$ENTREZID, "MF", org_db)
  gseGO_CC <- perform_enrichGO(result_data$ENTREZID, "CC", org_db)
  
  BP <- process_enrich_result(gseGO_BP, "Biological Process", region_label)
  MF <- process_enrich_result(gseGO_MF, "Molecular Function", region_label)
  CC <- process_enrich_result(gseGO_CC, "Cellular Component", region_label)
  
  combined_data <- rbind(BP, MF, CC)
  return(combined_data)
}

perform_enrichGO <- function(entrez_ids, ont, org_db) {
  enrichGO(gene = entrez_ids,
           OrgDb = org_db,
           keyType = "ENTREZID",
           ont = ont,
           readable = TRUE,
           pAdjustMethod = "fdr",
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)
}

process_enrich_result <- function(enrich_result, ontology, region_label, p_adjust_cutoff = 0.05) {
  result_df <- as.data.frame(enrich_result@result)
  result_df <- result_df[result_df$p.adjust <= p_adjust_cutoff, ]
  result_df$Ontology <- ontology
  result_df$Region <- region_label
  return(result_df)
}

GO_positive_data <- perform_GO_analysis(orthologs_positive, org.Dm.eg.db, "Positive")
GO_negative_data <- perform_GO_analysis(orthologs_negative, org.Dm.eg.db, "Negative")
GO_combined_data <- rbind(GO_positive_data, GO_negative_data)

# ---- 9. Save key objects ----
save_object(orthologs_positive, "orthologs_positive")
save_object(orthologs_negative, "orthologs_negative")
save_object(GO_combined_data, "GO_combined_data")

