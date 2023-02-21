# In this script, I investigate 3 scRNA-seq datasets from inflammatory conditions.
# The goal is to extract gene expression associated with especific cells, for example, Tregs.

# I downloaded the Seurat object from the publication:
# Classification of human chronic inflammatory skin disease based on single-cell immune profiling. Liu et al. 2022.
# https://www.science.org/doi/10.1126/sciimmunol.abl9165.

# Load libraries ----
library(Seurat)
library(tidyverse)
library(ggrepel)
library(RColorBrewer)

# Seurat objects: https://zenodo.org/record/6470377#.Y1gAa-zMIUH ----
Three_together <- readRDS("Three_together.rds")

Three_together <- SetIdent(Three_together, value = "Ident2")

DimPlot(Three_together, reduction="umap", label=F, pt.size = .1, label.size=4.5, #group.by = "Ident1",
        repel = T)+
  scale_color_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(41))+
  theme(text=element_text(family="Arial",size=15)) +
  theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',family="Arial",size=18),axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(0.6,"lines"),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=15),
        axis.title.y=element_text(colour='black', size=15),
        axis.text=element_text(colour='black',size=15),
        #legend.title=element_blank(),
        #legend.text=element_text(family="Arial", size=15),
        #legend.key=element_blank(),
        legend.position = "none")+
  theme(plot.title = element_text(size=22,colour = "black",face = "bold"))+ 
  guides(colour = guide_legend(override.aes = list(size=5)))

# Tregs analysis ----
## Combine Treg IDs:
Combined_Tregs <- factor(str_replace_all(Three_together@meta.data[["Ident2"]],
                                  "eTreg1|Treg-c|cmTreg|eTreg2",  "Tregs"))

Three_together@meta.data[["Combined_Tregs"]] <- Combined_Tregs
Three_together <- SetIdent(Three_together, value = "Combined_Tregs")

# Find Treg markers:
myTopHits_Tregs <- FindMarkers(Three_together, ident.1 = "Tregs")

myTopHits_Tregs1 <- myTopHits_Tregs %>%
  rownames_to_column("geneSymbol") %>%
  mutate(siggenes1 = case_when(
    avg_log2FC >= 1 & p_val_adj == 0 ~ geneSymbol,
    TRUE ~ NA_character_)) %>%
  mutate(siggenes2 = case_when(
    avg_log2FC >= 1 & p_val_adj == 0 ~ "sig",
    TRUE ~ "not_sig")) %>%
  
  arrange(desc(siggenes2))
  
  # ggplot:
myTopHits_Tregs1 %>%
  ggplot(., aes(y=pct.1*100, x=2^avg_log2FC,
                color = siggenes2, size = pct.1,
                alpha = siggenes2,
                NULL)) +
  geom_jitter() +
  theme_classic() +
  scale_color_manual(values = c("dark gray","red")) +
  scale_alpha_manual(values = c(0.3,1)) +
  geom_text_repel(aes(label = siggenes1), size = 3.5, fontface="bold",color="black",
                  max.overlaps = 200,
                  NULL) +
  theme(legend.position="right",legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15), panel.border = element_rect(fill = NA, size = 1.2)) +
  #geom_vline(xintercept = 0) + 
  xlab("avg_log2FC Treg vs. other cells") +
  xlim(1.1,6) + ylim(20,100)

# Dotplot of markers:
markers.to.plot <- na.omit(myTopHits_Tregs1$siggenes1)

library(viridis)
library(scales)
DotPlot(Three_together, features = markers.to.plot, cols=c("#5F4B8BFF", "#ED2B33FF"), assay = "RNA", col.min = 0.3,
        col.max = 0.8, dot.min=0.12, dot.scale = 1,
        cluster.idents=T)+
  scale_size(range = c(0, 5))+ 
  scale_size_area(max_size = 5)+ 
  #scale_color_viridis_c(name = 'log2 (count + 1)') + 
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 9),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size = 9),
        legend.text = element_text(size=9),
        legend.title = element_text(size = 9)) +
  scale_color_gradientn(colours = viridis::magma(20), limits = c(0,1), oob = scales::squish, name = 'log2 (count + 1)') +
  coord_flip()

myheatcol <- colorRampPalette(colors=c("gray","red"))(100)

FeaturePlot(Three_together, 
            reduction = "umap", 
            features = c("AC133644.2"),
            pt.size = 0.6, 
            order = TRUE,
            #split.by = "treatment",
            min.cutoff = 'q10',
            label = FALSE) +
  scale_color_gradientn(colours = myheatcol)

FeaturePlot(Three_together, 
            reduction = "umap", 
            features = c("FOXP3"),
            pt.size = 0.6, 
            order = TRUE,
            #split.by = "treatment",
            min.cutoff = 'q10',
            label = FALSE) +
  scale_color_gradientn(colours = myheatcol)

FeaturePlot(Three_together, 
            reduction = "umap", 
            features = c("TTN"),
            pt.size = 0.6, 
            order = TRUE,
            #split.by = "treatment",
            min.cutoff = 'q10',
            label = FALSE) +
  scale_color_gradientn(colours = myheatcol)
