##############################################
# 02_QTL_Annotation_and_GO_Enrichment.R
# Annotation and Gene Ontology enrichment analysis
##############################################

# ---- 1. Load packages and define helper functions ----
source(file.path("scripts", "_install_packages.R"))
source(file.path("scripts", "_save_objects.R"))

# Function to process QTL regions and annotate genes
process_qtl <- function(sigQTL_subset, 
                        annotation, 
                        refseq_flybase, 
                        orthologs, 
                        output_file) {
  
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

# Function to perform GO analysis
perform_GO_analysis <- function(result_data, 
                                org_db, 
                                region_label) {
  gseGO_BP <- perform_enrichGO(result_data$ENTREZID, "BP", org_db)
  gseGO_MF <- perform_enrichGO(result_data$ENTREZID, "MF", org_db)
  gseGO_CC <- perform_enrichGO(result_data$ENTREZID, "CC", org_db)
  
  BP <- process_enrich_result(gseGO_BP, "Biological Process", region_label)
  MF <- process_enrich_result(gseGO_MF, "Molecular Function", region_label)
  CC <- process_enrich_result(gseGO_CC, "Cellular Component", region_label)
  
  combined_data <- rbind(BP, MF, CC)
  return(combined_data)
}

# Wrapper function to call clusterProfiler::enrichGO
perform_enrichGO <- function(entrez_ids, 
                             ont, 
                             org_db) {
  clusterProfiler::enrichGO(gene = entrez_ids,
                            OrgDb = org_db,
                            keyType = "ENTREZID",
                            ont = ont,
                            readable = TRUE,
                            pAdjustMethod = "fdr",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05)
}

# Wrapper function to process GO results
process_enrich_result <- function(enrich_result, 
                                  ontology, 
                                  region_label, 
                                  p_adjust_cutoff = 0.05) {
  result_df <- as.data.frame(enrich_result@result)
  result_df <- result_df[result_df$p.adjust <= p_adjust_cutoff, ]
  result_df$Ontology <- ontology
  result_df$Region <- region_label
  return(result_df)
}

# Functions to plot GO enrich results
plot_GO_bar <- function(enrichResult, title = "GO plot", n_terms = NULL){
  df <- enrichResult%>%
    separate(GeneRatio, into = c("Gene", "Ratio"), sep = "/")%>%
    mutate(GeneRatio = as.numeric(Gene) / as.numeric(Ratio))
  
  if(!is.null(n_terms)){
    df <- df%>% arrange(df, p.adjust)%>% slice_head(n = n_terms)
    title <- paste0(title, ", showing ", n_terms, " terms")
  }
  
  df <- arrange(df, desc(GeneRatio))
  
  ontology_labels <- c(`Biological Process` = "Biological\nProcess",
                       `Molecular Function` = "Molecular\nFunction",
                       `Cellular Component` = "Cellular\nComponent")
  
  ggplot(df, aes(x = reorder(Description, GeneRatio), y = GeneRatio*100,
                 fill = Ontology))+
    geom_bar(stat = "identity") +
    theme_bw() + 
    facet_grid(. ~ Ontology, scales = "free_x", space = "free_x",
               labeller = labeller(Ontology = ontology_labels))+
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 35))+
    scale_fill_manual(values = c("#BC3C29FF", "#0072B5FF", "#E18727FF"))+
    theme_bw()+
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
          panel.grid = element_blank(),
          strip.background = element_rect(fill = "white", colour = "black"))+
    labs(y = "GeneRatio [%]", x = "", title = title)
}

plot_GO_point <- function(enrichResult, n_terms = NULL, title = "GO plot"){
  df <- enrichResult%>%
    separate(GeneRatio, into = c("Gene", "Ratio"), sep = "/")%>%
    mutate(GeneRatio = as.numeric(Gene) / as.numeric(Ratio))  
  
  if(!is.null(n_terms)){
      df <- df%>% arrange(df, p.adjust)%>% slice_head(n = n_terms)
      title <- paste0(title, ", showing ", n_terms, " terms")
    }
  
  df <- arrange(df, desc(GeneRatio))
  
  ontology_labels <- c(`Biological Process` = "Biological\nProcess",
                       `Molecular Function` = "Molecular\nFunction",
                       `Cellular Component` = "Cellular\nComponent")
  
  ggplot(df, aes(y = reorder(Description, GeneRatio), x = GeneRatio*100,
                 size = Count, color = p.adjust))+
    geom_point()+
    facet_grid(Ontology ~ ., scales = "free_y", space = "free_y",
               labeller = labeller(Ontology = ontology_labels))+
    scale_color_gradient(low = "red", high = "blue")+
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 35))+
    theme_bw()+
    theme(legend.position = "bottom", 
          axis.text.y = element_text(lineheight = 1.1),
          panel.grid = element_blank(),
          strip.background = element_rect(fill = "white", colour = "black"))+
    labs(x = "GeneRatio [%]", y = "", title = title)
    
}

# ---- 2. Load data ----
refseq_flybase <- read_tsv(file.path("data", "raw", "REFSEQ_FLYBASE_Dana.txt"))
dmel_dana_ortho <- read_excel(file.path("data", "raw", "dmel_dana_orthologs.xlsx"))
annotation <- rtracklayer::import(file.path("data", "raw", "genomic.gtf"))
sigQTL <- read_csv(file.path("data", "processed", "sigQTL.csv"))

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

# ---- 7. Annotate QTLs ----
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

orthologs_all <- process_qtl(sigQTL, 
                             annotation, 
                             refseq_flybase, 
                             orthologs, 
                             file.path("results", "supplementary_tables", 
                                       paste0("orthologs_dana_dmel_all_regions_", 
                                              timestamp, ".tsv")))

orthologs_positive <- process_qtl(pos_QTL, 
                                  annotation, 
                                  refseq_flybase, 
                                  orthologs, 
                                  file.path("results", "supplementary_tables", 
                                            paste0("orthologs_dana_dmel_positive_regions_", 
                                            timestamp, ".tsv")))
orthologs_negative <- process_qtl(neg_QTL, 
                                  annotation, 
                                  refseq_flybase, 
                                  orthologs, 
                                  file.path("results", "supplementary_tables", 
                                            paste0("orthologs_dana_dmel_negative_regions_", 
                                            timestamp, ".tsv")))

# ---- 8. GO enrichment analysis ----
GO_all_regions <- perform_GO_analysis(orthologs_all, org.Dm.eg.db, "All")
GO_positive_regions <- perform_GO_analysis(orthologs_positive, org.Dm.eg.db, "Positive")
GO_negative_regions <- perform_GO_analysis(orthologs_negative, org.Dm.eg.db, "Negative")
GO_combined_data <- rbind(GO_positive_regions, GO_negative_regions)

# ---- 9. Visualize enriched GO-terms ----
p1 <- plot_GO_bar(GO_all_regions,
                  title = "Figure 4: GO terms of genes within significant QTL regions",
                  n_terms = 25)

p2 <- plot_GO_bar(GO_positive_regions,
                 title = "Figure S4: GO terms of genes in regions with positive deltaSNP",
                 n_terms = 25)

p3 <- plot_GO_bar(GO_negative_regions,
                  title = "Figure S5: GO terms of genes in regions with negative deltaSNP",
                  n_terms = 25)

p4 <- plot_GO_point(GO_all_regions,
                    title = "GO terms of genes within significant QTL regions",
                    n_terms = 25)

p5 <- plot_GO_point(GO_positive_regions,
                    title = "GO terms of genes in regions with positive deltaSNP",
                    n_terms = 25)

p6 <- plot_GO_point(GO_negative_regions,
                    title = "GO terms of genes in regions with negative deltaSNP",
                    n_terms = 25)

# ---- 9. Save key objects and plots----
save_object(orthologs_all, "orthologs_all")
save_object(orthologs_positive, "orthologs_positive")
save_object(orthologs_negative, "orthologs_negative")
save_object(GO_combined_data, "GO_combined_data")

out_pdf <- file.path("results", "plots", 
                     paste0("02_GO_enrichment_plots_", timestamp, ".pdf"))
pdf(out_pdf, width = 14, height = 12)

print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)

dev.off()
##EOF