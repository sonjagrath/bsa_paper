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
GeneID_combined_data <- GO_combined_data[, c("ID", "Description", "geneID", 
                                             "p.adjust", "Ontology", "Region")]

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

# ---- 7. Table with CCRT values, main and supplementary CCRT plots ----
order <- c("BKK5","BKK6","BKK10","BKK12","BKK13","BKK16","BKK17","BKK18",
           "KATH14","KATH19","KATH23",
           "RIL57","RIL23","RIL22","RIL41","RIL50","RIL25","RIL20",
           "RIL30","RIL15","RIL47","RIL7","RIL81","RIL14","RIL58",
           "RIL93","RIL90", "RIL80")

seq_lines <- c("BKK12", "BKK13")

CCRT <- read_excel("data/raw/phenotype_RIL_IL.xlsx", sheet = "CCRT")%>%
  mutate(Population = ifelse(grepl("BKK", RIL), "BKK",
                             ifelse(grepl("RIL", RIL), "RIL", "KATH")))
CCRT$RIL <- factor(CCRT$RIL, levels = order)

fig1 <- ggplot(CCRT, aes(x = as.factor(RIL), y = Time, color = Population))+
  geom_boxplot()+
  ggh4x::facet_nested(.~Sex + Population, scales = "free_x", space = "free_x")+
  scale_x_discrete(labels = function(x){
    ifelse(x %in% seq_lines,
           paste0("<b>", x, "</b>"), x)})+
  labs(x = "Line",
       y = "Chill coma recovery time [min]",
       title = "Figure 1: Chill coma recovery time (in minutes) of iso-female and recombinant inbred lines.")+
  theme_bw()+
  scale_y_continuous(breaks = seq(5, 95, by = 10))+
  scale_color_manual(values = c("darkblue", "seagreen", "yellowgreen"))+
  theme(panel.grid = element_blank(),
        axis.text.x = ggtext::element_markdown(angle = 90, vjust = 0.5, hjust = 1),
        strip.background = element_rect(fill = "white"),
        legend.position = "bottom")

CCRT_table <- CCRT%>%
  group_by(RIL, Sex)%>%
  summarize(mean_T = mean(Time),
            sd_T = sd(Time))

figS1 <- ggplot(CCRT_table, aes(x = Sex, y = mean_T))+
  geom_boxplot()+
  geom_point(aes(colour = Sex))+
  geom_line(aes(group = RIL), color = "darkgrey")+
  geom_signif(comparisons = list(c("Male", "Female")),
              map_signif_level = FALSE, test = "t.test",
              test.args = list(paired = TRUE))+
  scale_y_continuous(breaks = seq(5, 95, by = 10))+
  labs(y = "Mean CCRT [min]",
       title = "Mean CCRT ")+
  scale_color_manual(values = c("purple", "orange"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = "none")

CCRT_table <- CCRT_table %>%
  pivot_wider(names_from = Sex, values_from = c(mean_T, sd_T), names_sep = "_" )
colnames(CCRT_table) <- c("Line", "Female mean CCRT [min]", "Female St Dev", 
                          "Male mean CCRT [min]", "Male St Dev")

# ---- 8. Table with Mortality values, main and supplementary Mortality plots ----
MORT <- read_excel("data/raw/phenotype_RIL_IL.xlsx", sheet = "CS")%>%
  mutate(mortality = Mortality/10)%>%
  mutate(Population = Root)
MORT$Line <- factor(MORT$Line, levels = order)

fig2 <- ggplot(MORT%>%mutate(mortality=mortality*100), aes(x = Line, y = mortality, fill = Population, group = Line))+
  geom_boxplot(aes(color = Population),fill="white")+
  geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.03,
               dotsize = 60)+
  ggh4x::facet_nested(.~Sex + Population, scales = "free_x", space = "free_x")+
  scale_x_discrete(labels = function(x){
    ifelse(x %in% seq_lines,
           paste0("<b>", x, "</b>"), x)})+
  scale_y_continuous(breaks = seq(0.0, 100, by = 10))+
  scale_fill_manual(values = c("darkblue", "seagreen", "yellowgreen"))+
  scale_color_manual(values = c("darkblue", "seagreen", "yellowgreen"))+
  theme_bw()+
  labs(x = "Line", y = "Mortality [%]",
       title = "Figure 2: Percent mortality of iso-female and recombinant inbred lines after 8 h cold shock")+
  theme(panel.grid = element_blank(),
        axis.text.x = ggtext::element_markdown(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom",
        strip.background = element_rect(fill = "white"))
MORT_table <- MORT%>%
  group_by(Line, Sex)%>%
  summarize(mean_M = mean(mortality),
            sd_M = sd(mortality))

figS2 <- ggplot(MORT_table%>%mutate(mean_M=mean_M*100), aes(x = Sex, y = mean_M))+
  geom_boxplot()+
  geom_point(aes(colour = Sex))+
  geom_line(aes(group = Line), color = "darkgrey")+
  geom_signif(comparisons = list(c("Male", "Female")),
              map_signif_level = FALSE, test = "t.test",
              test.args = list(paired = TRUE))+
  labs(y = "Mean mortality [%]",
       title = "Mean mortality upon 8-hour cold shock")+
  scale_color_manual(values = c("purple", "orange"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = "none")

MORT_table <- MORT_table %>%
  mutate(mean_M = round(mean_M * 10, 2))%>%
  mutate(sd_M = (sd_M * 10))%>%
  select(Line, Sex, mean_M, sd_M)%>%
  pivot_wider(names_from = Sex, values_from = c(mean_M, sd_M), names_sep = "_" )%>%
  select(Line, mean_M_Female, sd_M_Female, mean_M_Male, sd_M_Male)
colnames(MORT_table) <- c("Line", "Female mean Mortality", "Female St Dev", 
                          "Male mean Mortality", "Male St Dev")

# ---- 9. Calculation of LTi50 values, plot and table ----
lti_ril <- read_excel("BSA_manuscript/phenotype_RIL_IL.xlsx", sheet = "RIL_Mortality")
colnames(lti_ril)[1]<-"Line"
lti_il <- read_excel("BSA_manuscript/phenotype_RIL_IL.xlsx", sheet = "IL_Mortality")

LTI <- rbind(lti_ril, lti_il%>%filter(Time != 1 &  Time != 3))%>%
  mutate(Total = 10)%>%
  mutate(Line = factor(Line, levels = order))%>%
  arrange(Line)%>% na.omit()

LTI_table <- data.frame()
for (line in unique(LTI$Line)){
  for (sex in unique(LTI$Sex)){
    data <- LTI%>%filter(Line == line & Sex == sex)
    
    mod <- glm(cbind(Mortality, Total - Mortality) ~ Time, data, family="quasibinomial")
    anova(mod, test="LRT")
    
    LT50 <- -mod$coef[1]/mod$coef[2]
    grad <- c(-1/mod$coef[2], mod$coef[1]/(mod$coef[2]^2))
    CI <- sqrt(t(grad) %*% vcov(mod) %*% grad) * 1.96
    
    LTI_table <- rbind(LTI_table, data.frame(Line = line, Sex = sex, LT50 = LT50, CI = CI))
  }
}

LTI_table$Line <- factor(LTI_table$Line, levels = order)
LTI_table <- mutate(LTI_table, Population = ifelse(grepl("BKK", Line), "BKK",
                                     ifelse(grepl("RIL", Line), "RIL", "KATH")))

fig3 <- ggplot(LTI_table, aes(x = Line, y = LT50, color = Population))+
  ggh4x::facet_nested(.~Sex + Population, scales = "free_x", space = "free_x")+
  scale_color_manual(values = c("darkblue", "seagreen", "yellowgreen"))+
  geom_point()+
  geom_errorbar(aes(ymax = LT50 + CI, ymin = LT50 - CI), width = 0.4)+
  scale_x_discrete(labels = function(x){
    ifelse(x %in% seq_lines,
           paste0("<b>", x, "</b>"), x)})+
  labs(x = "Line", y = "LTi50 [h]")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        strip.background = element_rect(fill = "white"),
        axis.text.x = ggtext::element_markdown(angle = 90, vjust = 0.5, hjust = 1))
  
figS3 <- ggplot(LTI_table, aes(x = Sex, y = LT50))+
    geom_boxplot()+
    geom_point(aes(colour = Sex))+
    geom_line(aes(group = Line), color = "darkgrey")+
    geom_signif(comparisons = list(c("Male", "Female")),
                map_signif_level = FALSE, test = "t.test",
                test.args = list(paired = TRUE))+
    labs(y = "LTi50 values [h]",
         title = "LTi50 values")+
    scale_color_manual(values = c("purple", "orange"))+
    theme_bw()+
    theme(panel.grid = element_blank(),
          legend.position = "none",
          axis.text.x = ggtext::element_markdown())

LTI_table <- pivot_wider(LTI_table, names_from = Sex, values_from = c(LT50, CI))%>%
  select(Line, LT50_Female, CI_Female, LT50_Male, CI_Male)%>%
  arrange(Line)%>%
  mutate(LT50_Female = round(LT50_Female,2),
         LT50_Male = round(LT50_Male,2))
colnames(LTI_table) <- c("Line", "Female LTi50", "Female confidence interval",
                         "Male LTi50", "Male confidence interval")

# ---- 10. Write additional tables and plots to xlsx ----
out_xlsx <- file.path("results", "supplementary_tables", 
                     paste0("03_Supplementary_tables_", timestamp, ".xlsx"))

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Table S5")
openxlsx::writeData(wb, "Table S5", "Table S5: Sequencing quality control summary",
                    startRow = 1, startCol = 1)
openxlsx::writeData(wb, "Table S5", QC_rawdata,
                    startRow = 3, startCol = 1)

openxlsx::addWorksheet(wb, "Table S6")
openxlsx::writeData(wb, "Table S6", "Table S6: Summary of variant calling on the RIL sequencing data",
                    startRow = 1, startCol = 1)
openxlsx::writeData(wb, "Table S6", VCF_File_Statistics_Combined,
                    startRow = 3, startCol = 1)

openxlsx::addWorksheet(wb, "Table S7")
openxlsx::writeData(wb, "Table S7", "Table S7: Quantitative trait loci determined by QTLseqr",
                    startRow = 1, startCol = 1)
openxlsx::writeData(wb, "Table S7", qtl_table,
                    startRow = 3, startCol = 1)

openxlsx::addWorksheet(wb, "Table S8")
openxlsx::writeData(wb, "Table S8", "Table S8: Gene ontology term enrichment results for terms enriched in positive and negative deltaSNP regions determined by QTLseqr.",
                    startRow = 1, startCol = 1)
openxlsx::writeData(wb, "Table S8", GeneID_combined_data,
                    startRow = 3, startCol = 1)

openxlsx::addWorksheet(wb, "Table S9")
openxlsx::writeData(wb, "Table S9", "Table S5: Raw cold shock recovery times for all lines examined during the study, separated by sex.",
                    startRow = 1, startCol = 1)
openxlsx::writeData(wb, "Table S9", CCRT,
                    startRow = 3, startCol = 1)

#openxlsx::addWorksheet(wb, "Table S6")
#openxlsx::writeData(wb, "Table S6", "Table S6: Cold shock recovery times with standard devaition separated by sex for all strains examined during the study.",
#                    startRow = 1, startCol = 1)
#openxlsx::writeData(wb, "Table S6", CCRT_table,
#                    startRow = 3, startCol = 1)

openxlsx::addWorksheet(wb, "Table S10")
openxlsx::writeData(wb, "Table S10", "Table S10: Raw mortality data for 2, 4, 6, 8, 12 and 24 h cold shock separated by sex.",
                    startRow = 1, startCol = 1)
openxlsx::writeData(wb, "Table S10", LTI,
                    startRow = 3, startCol = 1)

#openxlsx::addWorksheet(wb, "Table S8")
#openxlsx::writeData(wb, "Table S8", "Table S8: Mortality upon 8h cold shock for all lines examined during the study.",
#                    startRow = 1, startCol = 1)
#openxlsx::writeData(wb, "Table S8", MORT_table,
#                    startRow = 3, startCol = 1)

#openxlsx::addWorksheet(wb, "Table S9")
#openxlsx::writeData(wb, "Table S9", "Table S9: LTi50 values and confidence intervals for males and females of all strains examined during the study.",
#                    startRow = 1, startCol = 1)
#openxlsx::writeData(wb, "Table S9", LTI_table,
#                    startRow = 3, startCol = 1)

saveWorkbook(wb, out_xlsx, overwrite = TRUE)

out_pdf <- file.path("results", "plots", 
                     paste0("03_Phenotyping_plots_", timestamp, ".pdf"))
pdf(out_pdf, width = 15, height = 8)

print(fig1)
print(fig2)
print(fig3)
print(figS1 + figS2 + figS3 +
        patchwork::plot_annotation(tag_levels = "A",
                                   title = "Boxplots comparing Male and Female mean CCRT, mean mortality and LTi50"))

dev.off()
###EOF

f1 <- figS1 + theme(legend.position = "none",
                   plot.title = element_blank(),
                   plot.margin = margin(0,0,0,0))+labs(x="",title="")
f2 <- figS2 + theme(legend.position = "none",
                   plot.title = element_blank(),
                   plot.margin = margin(0,0,0,0))+labs(x="",title="")
f3 <- figS3 + theme(legend.position = "none",
                   plot.title = element_blank(),
                   plot.margin = margin(0,0,0,0))+labs(x="",title="")

p <- f1 + f2 + f3+
  patchwork::plot_annotation(tag_levels = "A")
p

pdf("fig1.pdf", width = 7, height = 8)
p
dev.off()
