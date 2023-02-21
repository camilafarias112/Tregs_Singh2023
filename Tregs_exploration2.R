# In this script I investigate the expression of FOXP3 in human cutaneous leishmaniasis lesions,
# it's impact in gene expression, S. aureus burden and clinical outcome.
# This is for Tej's project, where he observed that upon colonization with S. aureus,
# and depletion of Tregs, Ifn-g production increases and it is associated with worst pathology.

# Here I explore this information from our LeishMicrobiome dataset (Amorim et al. 2023, 3rd lesion).

# Libraries ----
library(ggpubr)
library(ggthemes)
library(ggbreak)
library(ggrepel)
library(gplots)
library(patchwork)
library(RColorBrewer)
library(ggforce)
library(Seurat)
library(limma)
library(edgeR)
library(tidyverse)

#"%ni%" <- Negate("%in%")

# Ratio palette ----
myratiocol <- colorRampPalette(colors=c("red","white","blue"))(10)
myheatcol2 <- colorRampPalette(colors=c("gray","blue"))(100)

# Import the dataset and metadata ----
# Gene expression matrix:
geneExpression <- read_delim("~/Dropbox/CamilaAnalysis/Lesion3rd_dataset_Microbiome/Amorim_LeishMicrobiome/Public_submission/RNA-seq/GSE214397_Amorim_LeishMicrobiome_processed_norm_filt_log2cpm.txt.gz", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)
geneExpression <- as.matrix(column_to_rownames(geneExpression,"geneSymbol"))

# Study design:
studyDesign <- read_delim("~/Dropbox/CamilaAnalysis/Lesion3rd_dataset_Microbiome/Amorim_LeishMicrobiome/Public_submission/RNA-seq/GSE214397_Amorim_LeishMicrobiome_SupplemmentalTable1_LeishOmics_StudyDesign.txt.gz", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)
studyDesign <- studyDesign %>%
  mutate(LTCP_patient_ID = str_replace_all(LTCP_patient_ID, "^", "CL"))

studyDesign <- studyDesign %>% # Just biopsy samples
  filter(`Bacterial load (qPCR)` != "NA")

# FOXP3 expression ----
datagene <- t(geneExpression)
datagene <- rownames_to_column(as.data.frame(datagene),"LTCP_patient_ID")

datagene <- datagene %>%
  dplyr::select(LTCP_patient_ID, IFNG, FOXP3, IL1A, IL1B, LORICRIN, FLG, COL1A1, EPCAM,
         KRT1,
         KRT2,KRT10,KRT14,MKI67, DCN, CXCL5, IL1RL1, IL33,IL10,CTLA4) %>%
  left_join(studyDesign) %>%
  mutate(group = case_when(
    str_detect(LTCP_patient_ID, "HS") ~ "Healthy \nskin",
    TRUE ~ "CL \nlesions")) %>%
  mutate(FOXP3_groups = case_when(
    FOXP3 < 0.95 ~ "FOXP3_low",
    FOXP3 < 1.6 ~ "FOXP3_medium",
    TRUE ~ "FOXP3_high")) %>%
  dplyr::rename(cpm_saureus = "S. aureus abundance (CPM)") %>%
  dplyr::rename(cpm_saureus_cat = "S. aureus abundance - categorical") %>%
  mutate(cpm_saureus = case_when(cpm_saureus == 0 ~ (cpm_saureus + 0.0000001), TRUE ~ cpm_saureus)) %>%
  mutate(cpm_saureus_log = log10(cpm_saureus)) %>%
  mutate(cpm_saureus_cat = case_when(cpm_saureus_cat == "Saureus_detec" ~ "+",
                                     TRUE ~ "-"))

datagene %>%
  ggplot(., aes(x=fct_rev(group), y=FOXP3, color=group)) +
  geom_violin(size = 0.7, trim = T, adjust = 0.4) +
  geom_sina(size=2, shape=21, stroke = 1) +
  theme_stata() + scale_color_manual(values = c("black","red")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, vjust = 0.5), plot.background = element_rect(fill="white"),
        axis.text.y = element_text(size = 13, angle = 0), axis.title = element_text(size = 13),
        legend.position="none", legend.title = element_text(size = 17)) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.3, color="black") +
  #stat_compare_means(method = "wilcox.test", label.y = 0.9,
  #                   aes(label = ..p.signif..),
  #                   label.x = 1.4, size = 9, color="black") +
  xlab("")+ ylab(paste("FOXP3")) + 
  geom_hline(yintercept = 1.6) + geom_hline(yintercept = 0.95) + 
  NULL

# PCA ----
CL_datageneX <- datagene %>%
  filter(group == "CL \nlesions") %>%
  filter(FOXP3_groups != "FOXP3_medium")

pca.res <- prcomp(t(geneExpression[,CL_datageneX$LTCP_patient_ID]), scale.=F, retx=T)
pc.var<-pca.res$sdev^2
pc.per<-round(pc.var/sum(pc.var)*100, 1)

as.data.frame(pca.res$x) %>%
  rownames_to_column("LTCP_patient_ID") %>%
  select(LTCP_patient_ID, PC1, PC2, PC3) %>%
  left_join(datagene, by="LTCP_patient_ID") %>%
  ggplot(., aes(x=PC1, y=PC2, color=FOXP3_groups)) +
  theme_stata() + scale_color_manual(values = c("red","blue")) +
  geom_point(size=4) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 19, vjust = 0.5), plot.background = element_rect(fill="white"),
        axis.text.y = element_text(size = 19, angle = 0), axis.title = element_text(size = 19),
        legend.position="none", legend.title = element_text(size = 17)) +
  stat_ellipse(level = 0.95) +
  xlab(paste("PC1 -",pc.per[1],"%")) +
  ylab(paste("PC2 -",pc.per[2],"%")) +
  #annotate("text", x = -60, y = -60, label = "Permanova\nPr(>F)=0.005**", size=4, fontface=3) +
  coord_fixed()

library(vegan)
perm_vegan <- how(nperm = 199)

dist_vegan <- vegdist(t(2^geneExpression[,CL_datageneX$LTCP_patient_ID]), method = "bray")
vegan <- adonis2(dist_vegan ~ CL_datageneX$FOXP3_groups, permutations = perm_vegan, method="bray")
vegan

# Correlation IFNG and FOXP3 ----
datagene %>%
  mutate(group = case_when(
    PatientID %in% NA_character_ ~ "Healthy \nskin",
    TRUE ~ "CL \nlesions")) %>%
  filter(group == "CL \nlesions") %>%
  ggplot(., aes(x=IFNG, y=FOXP3)) +
  geom_smooth(method=lm, color = "red") +
  geom_point(size=2, shape=21, stroke=1) +
  theme_stata() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, vjust = 0.5), plot.background = element_rect(fill="white"),
        axis.text.y = element_text(size = 13, angle = 0), axis.title = element_text(size = 13),
        legend.position="none", legend.title = element_text(size = 17)) +
  #stat_cor(method = "pearson",
  #         label.sep = ", ", label.x.npc = 0.1, label.y.npc = 0.9,
  #         output.type = "expression",
  #         geom = "text", position = "identity", na.rm = FALSE,
  #         inherit.aes = TRUE) + 
  scale_x_break(c(-0.8, 2.8)) + 
  NULL

# Correlation between S. aureus and genes ----
datagene %>%
  mutate(group = case_when(
    PatientID %in% NA_character_ ~ "Healthy \nskin",
    TRUE ~ "CL \nlesions")) %>%
  filter(group == "CL \nlesions") %>%
  ggplot(., aes(x=cpm_saureus_cat, y=FOXP3)) +
  geom_jitter(position=position_jitter(0.1), size=2, shape = 21, stroke = 1) +
  theme_stata() + scale_fill_gradientn(colours = myratiocol) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, vjust = 0.5), plot.background = element_rect(fill="white"),
        axis.text.y = element_text(size = 13, angle = 0), axis.title = element_text(size = 13),
        legend.position="none", legend.title = element_text(size = 17)) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.3, color="black") +
  stat_compare_means(method = "wilcox.test", label.y = 2,
                     aes(label = ..p.signif..),
                     label.x = 1.3, size = 7, color="black") +
  xlab("S. aureus") + #ylab(paste("FOXP3/IFNG \nratio (log2)")) + 
  NULL

datagene %>%
  mutate(group = case_when(
    PatientID %in% NA_character_ ~ "Healthy \nskin",
    TRUE ~ "CL \nlesions")) %>%
  filter(group == "CL \nlesions") %>%
  filter(FOXP3_groups != "FOXP3_medium") %>%
  ggplot(., aes(x=FOXP3_groups, y=cpm_saureus_log)) +
  geom_jitter(position=position_jitter(0.1), size=2, shape = 21, stroke = 1) +
  theme_stata() + scale_fill_gradientn(colours = myratiocol) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, vjust = 0.5), plot.background = element_rect(fill="white"),
        axis.text.y = element_text(size = 13, angle = 0), axis.title = element_text(size = 13),
        legend.position="none", legend.title = element_text(size = 17)) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.3, color="black") +
  stat_compare_means(method = "wilcox.test", label.y = 2,
                     aes(label = ..p.format..),
                     label.x = 1.3, size = 7, color="black") +
  #xlab("S. aureus") + #ylab(paste("FOXP3/IFNG \nratio (log2)")) + 
  scale_y_break(c(-6.75, -1)) + 
  NULL

datagene %>%
  mutate(group = case_when(
    PatientID %in% NA_character_ ~ "Healthy \nskin",
    TRUE ~ "CL \nlesions")) %>%
  filter(group == "CL \nlesions") %>%
  filter(FOXP3_groups != "FOXP3_medium") %>%
  ggplot(., aes(x=FOXP3_groups, y=`Bacterial load (qPCR)`)) +
  geom_jitter(position=position_jitter(0.1), size=2, shape = 21, stroke = 1) +
  theme_stata() + scale_fill_gradientn(colours = myratiocol) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, vjust = 0.5), plot.background = element_rect(fill="white"),
        axis.text.y = element_text(size = 13, angle = 0), axis.title = element_text(size = 13),
        legend.position="none", legend.title = element_text(size = 17)) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.3, color="black") +
  stat_compare_means(method = "wilcox.test", label.y = 2,
                     aes(label = ..p.format..),
                     label.x = 1.3, size = 7, color="black") +
  #xlab("S. aureus") + #ylab(paste("FOXP3/IFNG \nratio (log2)")) + 
  #scale_y_break(c(-6.75, -1)) + 
  NULL

datagene %>%
  mutate(group = case_when(
    PatientID %in% NA_character_ ~ "Healthy \nskin",
    TRUE ~ "CL \nlesions")) %>%
  filter(group == "CL \nlesions") %>%
  #filter(FOXP3_groups != "FOXP3_medium") %>%
  ggplot(., aes(x=cpm_saureus_log, y=FOXP3)) +
  geom_smooth(method=lm, color = "purple") +
  geom_point(size=2, shape=21, stroke=1) +
  theme_stata() + 
  scale_fill_gradientn(colours = myratiocol) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, vjust = 0.5), plot.background = element_rect(fill="white"),
        axis.text.y = element_text(size = 13, angle = 0), axis.title = element_text(size = 13),
        legend.position="none", legend.title = element_text(size = 17)) +
  stat_cor(method = "spearman",
           label.sep = ", ", label.x.npc = 0.7, label.y.npc = 0.9,
           output.type = "expression",
           geom = "text", position = "identity", na.rm = FALSE,
           inherit.aes = TRUE) + 
  #ylim(-1,6.5) +
  scale_x_break(c(-6.75, -1)) + 
  ylab(paste("FOXP3 (log2CPM)")) + 
  #facet_wrap(. ~ FOXP3_groups) +
  #facet_wrap(. ~ cpm_saureus_cat) +
  NULL

datagene %>%
  mutate(group = case_when(
    PatientID %in% NA_character_ ~ "Healthy \nskin",
    TRUE ~ "CL \nlesions")) %>%
  filter(group == "CL \nlesions") %>%
  #filter(FOXP3_groups != "FOXP3_medium") %>%
  ggplot(., aes(x=`Bacterial load (qPCR)`, y=FOXP3)) +
  geom_smooth(method=lm, color = "purple") +
  geom_point(size=2, shape=21, stroke=1) +
  theme_stata() + 
  scale_fill_gradientn(colours = myratiocol) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, vjust = 0.5), plot.background = element_rect(fill="white"),
        axis.text.y = element_text(size = 13, angle = 0), axis.title = element_text(size = 13),
        legend.position="none", legend.title = element_text(size = 17)) +
  stat_cor(method = "spearman",
           label.sep = ", ", label.x.npc = 0.5, label.y.npc = 0.9,
           output.type = "expression",
           geom = "text", position = "identity", na.rm = FALSE,
           inherit.aes = TRUE) + 
  #ylim(-1,6.5) +
  #scale_x_break(c(-6.75, -1)) + 
  ylab(paste("FOXP3 (log2CPM)")) + 
  #facet_wrap(. ~ FOXP3_groups) +
  #facet_wrap(. ~ cpm_saureus_cat) +
  NULL

# Correlation of random genes and FOXP3 ----
datagene %>%
  mutate(group = case_when(
    PatientID %in% NA_character_ ~ "Healthy \nskin",
    TRUE ~ "CL \nlesions")) %>%
  filter(group == "CL \nlesions") %>%
  ggplot(., aes(x=IL1A, y=FOXP3, fill=FOXP3)) +
  geom_smooth(method=lm, color = "purple") +
  geom_point(size=2, shape=21, stroke=1) +
  theme_stata() + 
  scale_fill_gradientn(colours = myratiocol) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, vjust = 0.5), plot.background = element_rect(fill="white"),
        axis.text.y = element_text(size = 13, angle = 0), axis.title = element_text(size = 13),
        legend.position="none", legend.title = element_text(size = 17)) +
  stat_cor(method = "spearman",
           label.sep = ", ", label.x.npc = 0.2, label.y.npc = 0.8,
           output.type = "expression",
           geom = "text", position = "identity", na.rm = FALSE,
           inherit.aes = TRUE) + 
  #ylim(-1,6.5) +
  #scale_x_break(c(-6.75, -1)) + 
  ylab(paste("FOXP3/IFNG \nratio (log2)")) + 
  NULL

datagene %>%
  mutate(group = case_when(
    PatientID %in% NA_character_ ~ "Healthy \nskin",
    TRUE ~ "CL \nlesions")) %>%
  filter(group == "CL \nlesions") %>%
  ggplot(., aes(x=IL1B, y=FOXP3, fill=FOXP3)) +
  geom_smooth(method=lm, color = "purple") +
  geom_point(size=2, shape=21, stroke=1) +
  theme_stata() + 
  scale_fill_gradientn(colours = myratiocol) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, vjust = 0.5), plot.background = element_rect(fill="white"),
        axis.text.y = element_text(size = 13, angle = 0), axis.title = element_text(size = 13),
        legend.position="none", legend.title = element_text(size = 17)) +
  stat_cor(method = "spearman",
           label.sep = ", ", label.x.npc = 0.2, label.y.npc = 0.8,
           output.type = "expression",
           geom = "text", position = "identity", na.rm = FALSE,
           inherit.aes = TRUE) + 
  #ylim(-1,6.5) +
  #scale_x_break(c(-6.75, -1)) + 
  ylab(paste("FOXP3/IFNG \nratio (log2)")) + 
  NULL

datagene %>%
  mutate(group = case_when(
    PatientID %in% NA_character_ ~ "Healthy \nskin",
    TRUE ~ "CL \nlesions")) %>%
  filter(group == "CL \nlesions") %>%
  ggplot(., aes(x=LORICRIN, y=FOXP3, fill=FOXP3)) +
  geom_smooth(method=lm, color = "purple") +
  geom_point(size=2, shape=21, stroke=1) +
  theme_stata() + 
  scale_fill_gradientn(colours = myratiocol) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, vjust = 0.5), plot.background = element_rect(fill="white"),
        axis.text.y = element_text(size = 13, angle = 0), axis.title = element_text(size = 13),
        legend.position="none", legend.title = element_text(size = 17)) +
  stat_cor(method = "spearman",
           label.sep = ", ", label.x.npc = 0.2, label.y.npc = 0.8,
           output.type = "expression",
           geom = "text", position = "identity", na.rm = FALSE,
           inherit.aes = TRUE) + 
  #ylim(-1,6.5) +
  #scale_x_break(c(-6.75, -1)) + 
  ylab(paste("FOXP3/IFNG \nratio (log2)")) + 
  NULL +
  
  datagene %>%
  mutate(group = case_when(
    PatientID %in% NA_character_ ~ "Healthy \nskin",
    TRUE ~ "CL \nlesions")) %>%
  filter(group == "CL \nlesions") %>%
  ggplot(., aes(x=FLG, y=FOXP3, fill=FOXP3)) +
  geom_smooth(method=lm, color = "purple") +
  geom_point(size=2, shape=21, stroke=1) +
  theme_stata() + 
  scale_fill_gradientn(colours = myratiocol) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, vjust = 0.5), plot.background = element_rect(fill="white"),
        axis.text.y = element_text(size = 13, angle = 0), axis.title = element_text(size = 13),
        legend.position="none", legend.title = element_text(size = 17)) +
  stat_cor(method = "spearman",
           label.sep = ", ", label.x.npc = 0.2, label.y.npc = 0.8,
           output.type = "expression",
           geom = "text", position = "identity", na.rm = FALSE,
           inherit.aes = TRUE) + 
  #ylim(-1,6.5) +
  #scale_x_break(c(-6.75, -1)) + 
  ylab(paste("FOXP3/IFNG \nratio (log2)")) + 
  NULL

datagene %>%
  mutate(group = case_when(
    PatientID %in% NA_character_ ~ "Healthy \nskin",
    TRUE ~ "CL \nlesions")) %>%
  filter(group == "CL \nlesions") %>%
  ggplot(., aes(x=IFNG, y=FOXP3, fill=FOXP3)) +
  geom_smooth(method=lm, color = "purple") +
  geom_point(size=2, shape=21, stroke=1) +
  theme_stata() + 
  scale_fill_gradientn(colours = myratiocol) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, vjust = 0.5), plot.background = element_rect(fill="white"),
        axis.text.y = element_text(size = 13, angle = 0), axis.title = element_text(size = 13),
        legend.position="none", legend.title = element_text(size = 17)) +
  stat_cor(method = "spearman",
           label.sep = ", ", label.x.npc = 0.2, label.y.npc = 0.8,
           output.type = "expression",
           geom = "text", position = "identity", na.rm = FALSE,
           inherit.aes = TRUE) + 
  #ylim(-1,6.5) +
  #scale_x_break(c(-6.75, -1)) + 
  ylab(paste("FOXP3/IFNG \nratio (log2)")) + 
  NULL

datagene %>%
  mutate(group = case_when(
    PatientID %in% NA_character_ ~ "Healthy \nskin",
    TRUE ~ "CL \nlesions")) %>%
  filter(group == "CL \nlesions") %>%
  ggplot(., aes(x=COL1A1, y=FOXP3, fill=FOXP3)) +
  geom_smooth(method=lm, color = "purple") +
  geom_point(size=2, shape=21, stroke=1) +
  theme_stata() + 
  scale_fill_gradientn(colours = myratiocol) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, vjust = 0.5), plot.background = element_rect(fill="white"),
        axis.text.y = element_text(size = 13, angle = 0), axis.title = element_text(size = 13),
        legend.position="none", legend.title = element_text(size = 17)) +
  stat_cor(method = "spearman",
           label.sep = ", ", label.x.npc = 0.2, label.y.npc = 0.8,
           output.type = "expression",
           geom = "text", position = "identity", na.rm = FALSE,
           inherit.aes = TRUE) + 
  #ylim(-1,6.5) +
  #scale_x_break(c(-6.75, -1)) + 
  ylab(paste("FOXP3/IFNG \nratio (log2)")) + 
  NULL

datagene %>%
  mutate(group = case_when(
    PatientID %in% NA_character_ ~ "Healthy \nskin",
    TRUE ~ "CL \nlesions")) %>%
  filter(group == "CL \nlesions") %>%
  ggplot(., aes(x=LORICRIN, y=IFNG, fill=FOXP3)) +
  geom_smooth(method=lm, color = "purple") +
  geom_point(size=2, shape=21, stroke=1) +
  theme_stata() + 
  scale_fill_gradientn(colours = myratiocol) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, vjust = 0.5), plot.background = element_rect(fill="white"),
        axis.text.y = element_text(size = 13, angle = 0), axis.title = element_text(size = 13),
        legend.position="none", legend.title = element_text(size = 17)) +
  stat_cor(method = "spearman",
           label.sep = ", ", label.x.npc = 0.2, label.y.npc = 0.3,
           output.type = "expression",
           geom = "text", position = "identity", na.rm = FALSE,
           inherit.aes = TRUE) + 
  #ylim(-1,6.5) +
  #scale_x_break(c(-6.75, -1)) + 
  NULL +
  
  datagene %>%
  mutate(group = case_when(
    PatientID %in% NA_character_ ~ "Healthy \nskin",
    TRUE ~ "CL \nlesions")) %>%
  filter(group == "CL \nlesions") %>%
  ggplot(., aes(x=FLG, y=IFNG, fill=FOXP3)) +
  geom_smooth(method=lm, color = "purple") +
  geom_point(size=2, shape=21, stroke=1) +
  theme_stata() + 
  scale_fill_gradientn(colours = myratiocol) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, vjust = 0.5), plot.background = element_rect(fill="white"),
        axis.text.y = element_text(size = 13, angle = 0), axis.title = element_text(size = 13),
        legend.position="none", legend.title = element_text(size = 17)) +
  stat_cor(method = "spearman",
           label.sep = ", ", label.x.npc = 0.2, label.y.npc = 0.3,
           output.type = "expression",
           geom = "text", position = "identity", na.rm = FALSE,
           inherit.aes = TRUE) + 
  #ylim(-1,6.5) +
  #scale_x_break(c(-6.75, -1)) + 
  NULL

datagene %>%
  mutate(group = case_when(
    PatientID %in% NA_character_ ~ "Healthy \nskin",
    TRUE ~ "CL \nlesions")) %>%
  filter(group == "CL \nlesions") %>%
  ggplot(., aes(x=IL10, y=FOXP3, fill=FOXP3)) +
  geom_smooth(method=lm, color = "purple") +
  geom_point(size=2, shape=21, stroke=1) +
  theme_stata() + 
  scale_fill_gradientn(colours = myratiocol) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, vjust = 0.5), plot.background = element_rect(fill="white"),
        axis.text.y = element_text(size = 13, angle = 0), axis.title = element_text(size = 13),
        legend.position="none", legend.title = element_text(size = 17)) +
  stat_cor(method = "spearman",
           label.sep = ", ", label.x.npc = 0.2, label.y.npc = 0.8,
           output.type = "expression",
           geom = "text", position = "identity", na.rm = FALSE,
           inherit.aes = TRUE) + 
  #ylim(-1,6.5) +
  #scale_x_break(c(-6.75, -1)) + 
  ylab(paste("FOXP3/IFNG \nratio (log2)")) + 
  NULL +

datagene %>%
  mutate(group = case_when(
    PatientID %in% NA_character_ ~ "Healthy \nskin",
    TRUE ~ "CL \nlesions")) %>%
  filter(group == "CL \nlesions") %>%
  ggplot(., aes(x=CTLA4, y=FOXP3, fill=FOXP3)) +
  geom_smooth(method=lm, color = "purple") +
  geom_point(size=2, shape=21, stroke=1) +
  theme_stata() + 
  scale_fill_gradientn(colours = myratiocol) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, vjust = 0.5), plot.background = element_rect(fill="white"),
        axis.text.y = element_text(size = 13, angle = 0), axis.title = element_text(size = 13),
        legend.position="none", legend.title = element_text(size = 17)) +
  stat_cor(method = "spearman",
           label.sep = ", ", label.x.npc = 0.2, label.y.npc = 0.8,
           output.type = "expression",
           geom = "text", position = "identity", na.rm = FALSE,
           inherit.aes = TRUE) + 
  #ylim(-1,6.5) +
  #scale_x_break(c(-6.75, -1)) + 
  ylab(paste("FOXP3/IFNG \nratio (log2)")) + 
  NULL

# And with individual genes ----
datagene %>%
ggplot(., aes(x=cpm_saureus_log, y=IFNG)) +
  geom_smooth(method=lm, color = "red") +
  geom_point(size=2) +
  theme_classic() +
  #scale_color_manual(values = Dark24) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11, vjust = 0.5, hjust = 1), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="right") +
  stat_cor(method = "spearman",
           label.sep = ", ", label.x.npc = 0.7, label.y.npc = "bottom",
           output.type = "expression",
           geom = "text", position = "identity", na.rm = FALSE,
           inherit.aes = TRUE) + 
  #ylim(-1,6.5) +
  scale_x_break(c(-6.75, -1)) + 
  NULL

datagene %>%
  ggplot(., aes(x=cpm_saureus_log, y=FOXP3)) +
  geom_smooth(method=lm, color = "red") +
  geom_point(size=2) +
  theme_classic() +
  #scale_color_manual(values = Dark24) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11, vjust = 0.5, hjust = 1), 
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 11),
        legend.position="right") +
  stat_cor(method = "spearman",
           label.sep = ", ", label.x.npc = 0.7, label.y.npc = "bottom",
           output.type = "expression",
           geom = "text", position = "identity", na.rm = FALSE,
           inherit.aes = TRUE) + 
  #ylim(-1,6.5) +
  scale_x_break(c(-6.75, -1)) + 
  NULL




# Metadata and FOXP3 ----
datagene %>%
  filter(treatment_other_drug == "antimony") %>%
  ggplot(., aes(x=treatment_outcome, y=FOXP3, color=treatment_outcome)) +
  geom_jitter(position=position_jitter(0.1), size=3) +
  theme_classic() + scale_color_manual(values = c("dark gray","red")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 20, vjust = 1, hjust= 1,angle = 45), 
        axis.text.y = element_text(size = 20), axis.title = element_text(size = 20),
        legend.position="none", legend.text = element_text(size = 20), legend.title = element_blank()) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.3, color="black") +
  stat_compare_means(method = "wilcox.test", label.y = 0,
                     aes(label = ..p.format..),
                     label.x = 1.5, size = 7, color="black") +
  xlab("") + ylab(paste("FOXP3/IFNG \nratio (log2)")) + 
  NULL

datagene %>%
  filter(treatment_other_drug == "antimony") %>%
  ggplot(., aes(x=healing_time_days, y=FOXP3)) +
  geom_smooth(method=lm, color = "red") +
  geom_point(size=2) +
  theme_classic() +
  #scale_color_manual(values = Dark24) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 20, vjust = 0.5, hjust = 1), 
        axis.text.y = element_text(size = 20), axis.title = element_text(size = 20),
        legend.position="right") +
  stat_cor(method = "spearman",
           label.sep = ", ", label.x.npc = 0.2, label.y.npc = 0.8,
           output.type = "expression",
           geom = "text", position = "identity", na.rm = FALSE,
           inherit.aes = TRUE) + 
  ylab(paste("FOXP3/IFNG \nratio (log2)")) + 
  NULL


datagene %>%
  filter(treatment_other_drug == "antimony") %>%
  ggplot(., aes(x=size_lesion_mm2, y=FOXP3)) +
  geom_smooth(method=lm, color = "red") +
  geom_point(size=2) +
  geom_smooth(method=lm, color = "red") +
  geom_point(size=2, shape=21, stroke=1) +
  theme_stata() + scale_fill_gradientn(colours = myratiocol) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, vjust = 0.5), plot.background = element_rect(fill="white"),
        axis.text.y = element_text(size = 13, angle = 0), axis.title = element_text(size = 13),
        legend.position="none", legend.title = element_text(size = 17)) +
  stat_cor(method = "spearman",
           label.sep = ", ", label.x.npc = 0.3, label.y.npc = 0.8,
           output.type = "expression",
           geom = "text", position = "identity", na.rm = FALSE,
           inherit.aes = TRUE) + 
  ylab(paste("FOXP3")) +
  xlab(paste("Lesion size (mm2)")) + 
  NULL

# Parasite loads and Ratios ----
datagene %>%
  filter(treatment_other_drug == "antimony") %>%
  ggplot(., aes(x=log(`L.braziliensis load (qPCR)`,2), y=FOXP3)) +
  geom_smooth(method=lm, color = "purple") +
  geom_point(size=2, shape=21, stroke=1) +
  theme_stata() + scale_fill_gradientn(colours = myratiocol) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, vjust = 0.5), plot.background = element_rect(fill="white"),
        axis.text.y = element_text(size = 13, angle = 0), axis.title = element_text(size = 13),
        legend.position="right", legend.title = element_text(size = 17)) +
  stat_cor(method = "spearman",
           label.sep = ", ", label.x.npc = 0.2, label.y.npc = 0.8,
           output.type = "expression",
           geom = "text", position = "identity", na.rm = FALSE,
         inherit.aes = TRUE) + 
  ylab(paste("FOXP3")) + xlab("L.braziliensis burden")
  NULL

# DE with FOXP3 ----
CL_datagene <- datagene %>%
  filter(group == "CL \nlesions")
design <- model.matrix(~0 + factor(CL_datagene$FOXP3_groups))
colnames(design) <- levels(factor(CL_datagene$FOXP3_groups))
head(design)
dim(design)

v.DEGList.filtered.norm <- voom(2^geneExpression[,colnames(geneExpression) %in% CL_datagene$LTCP_patient_ID], design)
fit <- lmFit(v.DEGList.filtered.norm, design)
contrast.matrix <- makeContrasts(High_Low = FOXP3_high - FOXP3_low,
                                 levels=design)
fits <- contrasts.fit(fit, contrast.matrix)
ebFit <- eBayes(fits)

myTopHits0 <- topTable(ebFit, adjust ="BH", number=40000, sort.by="logFC")
genelist1 <- c("ADRA2A", "ZFP36L1", "SERPING1", "CD34", "COL1A1", "COL5A1", "CCN2",
               "F13A1", "FBLN1", "GLI3", "GPR4", "CCN1", "MAP1B", "MMP2", "PDGFB",
               "SAA1", "TGFB1", "ZFP36", "MYL9", "TXN2", "FKBP10", "PEAR1", "APLNR", "DMTN", "COX2", "PRKCD", "THBS1")

genelist2 <- c("IFNG", "KLRB1", "KLRC1", "KLRC3", "KLRD1",
               "SH2D1A", "PRF1", "KLRC4", "KLRK1", "CADM1", "IL21", "GZMH", "GZMA")
# Volcano:
myTopHits <- myTopHits0 %>%
rownames_to_column("geneSymbol") %>%
  mutate(siggenes3 = case_when(
    geneSymbol %in% "FOXP3" ~ "FOXP3",
    TRUE ~ NA_character_)) %>%
  mutate(siggenes1 = case_when(
    logFC > 0.59 & adj.P.Val < 0.05 ~ "FOXP3_high",
    logFC < -0.59 & adj.P.Val < 0.05 ~ "FOXP3_low",
    TRUE ~ "not_sig")) %>%
  mutate(siggenes2 = case_when(
    logFC > 0.59 & adj.P.Val < 0.05 | logFC < -0.59 & adj.P.Val < 0.05 ~ geneSymbol,
    TRUE ~ NA_character_)) %>%
  mutate(siggenes2 = case_when(
    geneSymbol %in% "FOXP3" ~ NA_character_,
    TRUE ~ siggenes2)) %>%
  
  # Genes asked by Tej:
  mutate(siggenes3 = case_when(
    geneSymbol %in% c(genelist1,genelist2) & geneSymbol != "IFNG" ~ geneSymbol,
    TRUE ~ NA_character_)) %>%
  mutate(siggenes4 = case_when(
    geneSymbol %in% "IFNG" ~ "IFNG",
    TRUE ~ NA_character_)) %>%
  mutate(siggenes5 = case_when(
    geneSymbol %in% genelist1 ~ "WH",
    geneSymbol %in% genelist2 ~ "CD",
    TRUE ~ "other"))
  
# ggplot:
myTopHits %>%
  arrange(siggenes5) %>%
  ggplot(., aes(y=-log10(adj.P.Val), x=logFC,
                color = siggenes5, size = siggenes5, alpha = siggenes5)) +
  geom_jitter() +
  theme_stata() +
  scale_color_manual(values = c("red","gray","blue")) +
  scale_alpha_manual(values = c(1,0.4,1)) +
  scale_size_manual(values = c(2,1,2)) +
  geom_text_repel(aes(label = siggenes4), size = 5, fontface=4, color="red",
                  xlim = c(-1.5, -0.8),
                  ylim = c(0.5, 0.8)) +
  geom_text_repel(aes(label = siggenes3), size = 4.5, fontface=3, color="black",
                  max.overlaps = 50,
                  NULL, force = 70) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 17, vjust = 0.5), plot.background = element_rect(fill="white"),
        axis.text.y = element_text(size = 17, angle = 0), axis.title = element_text(size = 17),
        legend.position="none", legend.title = element_text(size = 17), axis.line = element_line(size = 1)) +
  #geom_hline(yintercept = -log10(0.05)) + 
  #geom_vline(xintercept = -0.59) + 
  #xlim(0,0.85) +
  #ylim(1,3.35) +
  ylim(0,3.35) +
  xlab(paste("logFC FOXP3high vs. FOXP3low"))

# For exporting:

myTopHits_export <- myTopHits0 %>%
  rownames_to_column("geneSymbol") %>%
  mutate(siggenes1 = case_when(
    logFC > 0.59 & adj.P.Val < 0.05 ~ "FOXP3_high",
    logFC < -0.59 & adj.P.Val < 0.05 ~ "FOXP3_low",
    TRUE ~ "not_sig")) %>% filter(siggenes1 != "not_sig")

FOXP3_high <- myTopHits_export %>% filter(siggenes1 %in% "FOXP3_high")
FOXP3_high <- FOXP3_high$geneSymbol

FOXP3_low <- myTopHits_export %>% filter(siggenes1 %in% "FOXP3_low")
FOXP3_low <- FOXP3_low$geneSymbol

write_tsv(myTopHits, "outputs/myTopHits.txt")
write_tsv(as.data.frame(FOXP3_high), "outputs/FOXP3_high.txt")
write_tsv(as.data.frame(FOXP3_low), "outputs/FOXP3_low.txt")

# Cell estimations ----
#library(immunedeconv)
# MCP counter:
res_mcp_counter <- deconvolute(2^geneExpression, "mcp_counter")

res_mcp_counter_g <- res_mcp_counter %>%
  gather(LTCP_patient_ID, cell_score, -cell_type) %>%
  left_join(datagene, by="LTCP_patient_ID")


res_mcp_counter_g %>%
    mutate(cell_type = case_when(
      cell_type %in% "Cancer associated fibroblast" ~ "Fibroblasts",
      cell_type %in% "cytotoxicity score" ~ "CTLs",
      TRUE ~ cell_type)) %>%
  filter(FOXP3_groups != "FOXP3_medium") %>%
  filter(group != "Healthy \nskin") %>%
    #filter(cell_type %in% c("T cell", "CTLs","NK cell")) %>%
    ggplot(., aes(x=fct_rev(FOXP3_groups), y=log(cell_score,2), color=FOXP3_groups)) +
    geom_violin(size = 0.7, trim = T, adjust = 0.4) +
    geom_sina(size=2, shape=21, stroke = 1) +
    theme_stata() + scale_color_manual(values = c("red","blue")) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = 13, vjust = 0.5), plot.background = element_rect(fill="white"),
          axis.text.y = element_text(size = 13, angle = 0), axis.title = element_text(size = 13),
          legend.position="none", legend.title = element_text(size = 17), strip.background = element_rect(fill = NA)) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.3, color="black") +
  stat_compare_means(method = "wilcox.test", #label.y = 3,
                     aes(label = ..p.signif..),
                     label.x = 1.4, size = 9, color="black") +  
  ylab(paste("MCP counter \nscore")) +  xlab("FOXP3 groups") +
    facet_wrap(. ~ cell_type, scales = "free")

res_mcp_counter_g %>%
  mutate(cell_type = case_when(
    cell_type %in% "Cancer associated fibroblast" ~ "Fibroblasts",
    cell_type %in% "cytotoxicity score" ~ "CTLs",
    TRUE ~ cell_type)) %>%
  #filter(FOXP3_groups != "FOXP3_medium") %>%
  filter(group != "Healthy \nskin") %>%
  #filter(cell_type %in% c("T cell", "CTLs","NK cell")) %>%
  ggplot(., aes(x=FOXP3, y=log(cell_score,2))) +
  geom_smooth(method=lm, color = "purple") +
  geom_point(size=2, shape=21, stroke=1) +
  theme_stata() + 
  scale_fill_gradientn(colours = myratiocol) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, vjust = 0.5), plot.background = element_rect(fill="white"),
        axis.text.y = element_text(size = 13, angle = 0), axis.title = element_text(size = 13),
        legend.position="none", legend.title = element_text(size = 17)) +
  stat_cor(method = "spearman",
           label.sep = ", ", label.x.npc = 0.2, label.y.npc = 0.8,
           output.type = "expression",
           geom = "text", position = "identity", na.rm = FALSE,
           inherit.aes = TRUE) +  
  ylab(paste("MCP counter \nscore")) +  xlab("FOXP3 groups") +
  facet_wrap(. ~ cell_type, scales = "free")


#xCell:
res_xcell <- deconvolute(2^geneExpression, "xcell")
  
res_xcell_g <- res_xcell %>%
    gather(LTCP_patient_ID, cell_score, -cell_type) %>%
    left_join(datagene, by="LTCP_patient_ID")
  
res_xcell_g %>%
  filter(FOXP3_groups != "FOXP3_medium") %>%
  filter(group != "Healthy \nskin") %>%
  #filter(cell_type %in% c("T cell", "CTLs","NK cell")) %>%
  ggplot(., aes(x=fct_rev(FOXP3_groups), y=cell_score, color=FOXP3_groups)) +
  geom_violin(size = 0.7, trim = T, adjust = 0.4) +
  geom_sina(size=2, shape=21, stroke = 1) +
  theme_stata() + scale_color_manual(values = c("red","blue")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, vjust = 0.5), plot.background = element_rect(fill="white"),
        axis.text.y = element_text(size = 13, angle = 0), axis.title = element_text(size = 13),
        legend.position="none", legend.title = element_text(size = 17), strip.background = element_rect(fill = NA)) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.3, color="black") +
  stat_compare_means(method = "wilcox.test", #label.y = 3,
                     aes(label = ..p.signif..),
                     label.x = 1.4, size = 9, color="black") +  
  ylab(paste("MCP counter \nscore")) +  xlab("FOXP3 groups") +
  facet_wrap(. ~ cell_type, scales = "free")

res_xcell_g %>%
  filter(FOXP3_groups != "FOXP3_medium") %>%
  filter(group != "Healthy \nskin") %>%
  filter(str_detect(cell_type, "regulatory")) %>%
  ggplot(., aes(x=fct_rev(FOXP3_groups), y=cell_score, color=FOXP3_groups)) +
  geom_violin(size = 0.7, trim = T, adjust = 0.4) +
  geom_sina(size=2, shape=21, stroke = 1) +
  theme_stata() + scale_color_manual(values = c("red","blue")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, vjust = 0.5), plot.background = element_rect(fill="white"),
        axis.text.y = element_text(size = 13, angle = 0), axis.title = element_text(size = 13),
        legend.position="none", legend.title = element_text(size = 17), strip.background = element_rect(fill = NA)) +
  stat_summary(fun = median,
               geom = "crossbar", width = 0.3, color="black") +
  stat_compare_means(method = "wilcox.test", label.y = 0.02,
                     aes(label = ..p.signif..),
                     label.x = 1.4, size = 9, color="black") +  
  ylab(paste("MCP counter \nscore")) +  xlab("FOXP3 groups") +
  facet_wrap(. ~ cell_type, scales = "free")

# Quantiseq
res_quantiseq <- deconvolute(2^geneExpression, "quantiseq", tumor = FALSE)

res_quantiseq %>%
  gather(LTCP_patient_ID, cell_score, -cell_type) %>%
  #left_join(datagene, by="LTCP_patient_ID") %>%
  ggplot(aes(x = LTCP_patient_ID, y = cell_score, fill = cell_type)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_brewer(palette = "Paired") +
  scale_x_discrete(limits = rev(levels(res_quantiseq)))

res_quantiseq_g <- res_quantiseq %>%
  gather(LTCP_patient_ID, cell_score, -cell_type) %>%
  left_join(datagene, by="LTCP_patient_ID")

res_quantiseq_g %>%
  filter(FOXP3_groups != "FOXP3_medium") %>%
  filter(group != "Healthy \nskin") %>%
  #filter(cell_type %in% c("T cell", "CTLs","NK cell")) %>%
  ggplot(., aes(x=fct_rev(FOXP3_groups), y=cell_score, color=FOXP3_groups)) +
  geom_violin(size = 0.7, trim = T, adjust = 0.4) +
  geom_sina(size=2, shape=21, stroke = 1) +
  theme_stata() + scale_color_manual(values = c("red","blue")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, vjust = 0.5), plot.background = element_rect(fill="white"),
        axis.text.y = element_text(size = 13, angle = 0), axis.title = element_text(size = 13),
        legend.position="none", legend.title = element_text(size = 17), strip.background = element_rect(fill = NA)) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.3, color="black") +
  stat_compare_means(method = "wilcox.test", #label.y = 3,
                     aes(label = ..p.signif..),
                     label.x = 1.4, size = 9, color="black") +  
  ylab(paste("QuanTIseq \nscore")) +  xlab("FOXP3 groups") +
  facet_wrap(. ~ cell_type, scales = "free")

# CIBERSORT - COULDN'T MAKE IT WORK
set_cibersort_binary("/path/to/CIBERSORT.R")
set_cibersort_mat("/path/to/LM22.txt")
res_cibersort <- deconvolute(2^geneExpression, "cibersort")

res_cibersort %>%
  gather(LTCP_patient_ID, cell_score, -cell_type) %>%
  #left_join(datagene, by="LTCP_patient_ID") %>%
  ggplot(aes(x = LTCP_patient_ID, y = cell_score, fill = cell_type)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_brewer(palette = "Paired") +
  scale_x_discrete(limits = rev(levels(res_cibersort)))

res_cibersort_g <- res_cibersort %>%
  gather(LTCP_patient_ID, cell_score, -cell_type) %>%
  left_join(datagene, by="LTCP_patient_ID")

res_cibersort_g %>%
  filter(FOXP3_groups != "FOXP3_medium") %>%
  filter(group != "Healthy \nskin") %>%
  #filter(cell_type %in% c("T cell", "CTLs","NK cell")) %>%
  ggplot(., aes(x=fct_rev(FOXP3_groups), y=cell_score, color=FOXP3_groups)) +
  geom_violin(size = 0.7, trim = T, adjust = 0.4) +
  geom_sina(size=2, shape=21, stroke = 1) +
  theme_stata() + scale_color_manual(values = c("red","blue")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, vjust = 0.5), plot.background = element_rect(fill="white"),
        axis.text.y = element_text(size = 13, angle = 0), axis.title = element_text(size = 13),
        legend.position="none", legend.title = element_text(size = 17), strip.background = element_rect(fill = NA)) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.3, color="black") +
  stat_compare_means(method = "wilcox.test", #label.y = 3,
                     aes(label = ..p.signif..),
                     label.x = 1.4, size = 9, color="black") +  
  ylab(paste("cibersort \nscore")) +  xlab("FOXP3 groups") +
  facet_wrap(. ~ cell_type, scales = "free")

# DE Cell estimations:
design2 <- model.matrix(~0 + factor(CL_datageneX$FOXP3_groups))
colnames(design2) <- levels(factor(CL_datageneX$FOXP3_groups))

res_xcell_matrix <- as.matrix(column_to_rownames(res_xcell, "cell_type"))
v.cell <- voom(res_xcell_matrix[,colnames(res_xcell_matrix) %in% CL_datageneX$LTCP_patient_ID], design2)
fit.cell <- lmFit(v.cell, design2)
contrast.matrix.cell <- makeContrasts(High_Low = FOXP3_high - FOXP3_low,
                                 levels=design2)
fits.cell <- contrasts.fit(fit.cell, contrast.matrix.cell)
ebFit.cell <- eBayes(fits.cell)

myTopHits.cell0 <- topTable(ebFit.cell, adjust ="BH", number=40000, sort.by="logFC")
myTopHits.cell0 %>% filter(P.Value < 0.05)

# Something is wrong here!


# None of thee have Keratinocytes - Tej's focus. So I'm gonna get Keratinocyte gene markers from the 
# L.major scRNAseq dataset.
  
# GSVA for Cell death and wound healing scoring ----
library(GSVA)
library(GSEABase)

TejSig_gmt <- getGmt("Tej_geneSig.gmt", geneIdType=SymbolIdentifier()) #Importing gmt

GSVA.res <- gsva(geneExpression, TejSig_gmt, mx.diff=FALSE, method="ssgsea")

#my_comparisons <- list(c("+","-"))
as.data.frame(GSVA.res) %>%
  rownames_to_column("geneset") %>%
  gather(LTCP_patient_ID, ss_score, -geneset) %>%
  full_join(datagene) %>%
  mutate(group = case_when(
    PatientID %in% NA_character_ ~ "Healthy \nskin",
    TRUE ~ cpm_saureus_cat),
    group = factor(group, levels = c("Healthy \nskin","-","+"))) %>%
  #filter(group == "CL \nlesions") %>%
  ggplot(., aes(x=group, y=ss_score)) +
  geom_jitter(position=position_jitter(0.1), size=2, shape = 21, stroke = 1) +
  theme_stata() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, vjust = 0.5), plot.background = element_rect(fill="white"),
        axis.text.y = element_text(size = 13, angle = 0), axis.title = element_text(size = 13),
        legend.position="none", legend.title = element_text(size = 17)) +
  stat_summary(fun = median,
               geom = "crossbar", width = 0.6, color="red") +
  #stat_compare_means(method = "wilcox.test", label.y = 0.6,
  #                   aes(label = ..p.signif..), comparisons = my_comparisons,
  #                   label.x = 1.3, size = 7, color="black") +
  xlab("") + ylab("Cell death score \n(ssGSEA)") +
  facet_wrap(. ~ geneset)

as.data.frame(GSVA.res) %>%
  rownames_to_column("geneset") %>%
  gather(LTCP_patient_ID, ss_score, -geneset) %>%
  full_join(datagene) %>%
  mutate(group = case_when(
    PatientID %in% NA_character_ ~ "Healthy \nskin",
    TRUE ~ cpm_saureus_cat),
    group = factor(group, levels = c("Healthy \nskin","-","+"))) %>%
  filter(group != "Healthy \nskin") %>%
  ggplot(., aes(x=group, y=ss_score)) +
  geom_jitter(position=position_jitter(0.1), size=2, shape = 21, stroke = 1) +
  theme_stata() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, vjust = 0.5), plot.background = element_rect(fill="white"),
        axis.text.y = element_text(size = 13, angle = 0), axis.title = element_text(size = 13),
        legend.position="none", legend.title = element_text(size = 17)) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.6, color="red") +
  stat_compare_means(method = "wilcox.test", label.y = 0.20,
                     aes(label = ..p.signif..), comparisons = my_comparisons,
                     label.x = 1.3, size = 7, color="black") +
  #scale_y_break(c(-0.15, 0.25)) + 
  xlab("") + ylab("Cell death score \n(ssGSEA)") +
  facet_wrap(. ~ geneset, scales = "free")

# Export raw values to be plotted in GraphPad ----
library("writexl")
write_xlsx(as.data.frame(log(t(res_xcell %>%
                                 column_to_rownames("cell_type")),2)) %>%
             rownames_to_column("LTCP_patient_ID") %>%
             left_join(datagene) %>%
             rename(Lbraziliensis_burdenlog2 = "L.braziliensis load (qPCR)") %>%
             mutate(Lbraziliensis_burdenlog2 = log(Lbraziliensis_burdenlog2,2)) %>%
             select(-cpm_saureus) %>%
             left_join(as.data.frame(pca.res$x) %>%
                         rownames_to_column("LTCP_patient_ID") %>%
                         select(PC1,PC2, LTCP_patient_ID)),"outputs/RawValues_FOXP3_project_version6.xlsx")

write_tsv(as.data.frame(log(t(res_xcell %>%
                                column_to_rownames("cell_type")),2)) %>%
            rownames_to_column("LTCP_patient_ID") %>%
            left_join(datagene) %>%
            rename(Lbraziliensis_burdenlog2 = "L.braziliensis load (qPCR)") %>%
            mutate(Lbraziliensis_burdenlog2 = log(Lbraziliensis_burdenlog2,2)) %>%
            select(-cpm_saureus) %>%
            left_join(as.data.frame(pca.res$x) %>%
                        rownames_to_column("LTCP_patient_ID") %>%
                        select(PC1,PC2, LTCP_patient_ID)) %>%
            select(1,41:77) %>%
            mutate(group = case_when(
              str_detect(LTCP_patient_ID, "HS") ~ "HS",
              TRUE ~ "CL")),"outputs/RawValues_FOXP3_project_version6.txt")

write_tsv(as.data.frame(log(t(res_mcp_counter %>%
                                column_to_rownames("cell_type")),2)) %>%
            rownames_to_column("LTCP_patient_ID") %>%
            rename("Fibroblasts" = "Cancer associated fibroblast",
                   "CTLs" = "cytotoxicity score") %>%
            left_join(datagene) %>%
            rename(Lbraziliensis_burdenlog2 = "L.braziliensis load (qPCR)") %>%
            mutate(Lbraziliensis_burdenlog2 = log(Lbraziliensis_burdenlog2,2)) %>%
            select(-cpm_saureus) %>%
            left_join(as.data.frame(pca.res$x) %>%
                        rownames_to_column("LTCP_patient_ID") %>%
                        select(PC1,PC2, LTCP_patient_ID)) %>%
            mutate(group = case_when(
              str_detect(LTCP_patient_ID, "HS") ~ "HS",
              TRUE ~ "CL")),"outputs/RawValues_FOXP3_project_version7.txt")

write_tsv(as.data.frame(log(t(res_mcp_counter %>%
                                column_to_rownames("cell_type")),2)) %>%
            rownames_to_column("LTCP_patient_ID") %>%
            rename("Fibroblasts" = "Cancer associated fibroblast",
                   "CTLs" = "cytotoxicity score") %>%
            left_join(datagene) %>%
            dplyr::rename(Lbraziliensis_burdenlog2 = "L.braziliensis load (qPCR)") %>%
            mutate(Lbraziliensis_burdenlog2 = log(Lbraziliensis_burdenlog2,2)) %>%
            dplyr::select(-cpm_saureus) %>%
            
            # PCA
            left_join(as.data.frame(pca.res$x) %>%
                        rownames_to_column("LTCP_patient_ID") %>%
                        dplyr::select(PC1,PC2, LTCP_patient_ID)) %>%
            mutate(group = case_when(
              str_detect(LTCP_patient_ID, "HS") ~ "HS",
              TRUE ~ "CL")) %>%
            
            #ssGSEA
            left_join(as.data.frame(GSVA.res) %>%
                        rownames_to_column("geneset") %>%
                        gather(LTCP_patient_ID, ss_score, -geneset) %>%
                        pivot_wider(names_from = "geneset", values_from = ss_score)),
          "outputs/RawValues_FOXP3_project_version8.txt")