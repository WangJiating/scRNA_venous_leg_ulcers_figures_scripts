rm(list = ls())
library(Seurat)
library(harmony)
library(tidyverse)
library(patchwork)
library(dittoSeq)
library(cowplot)
library(paletteer) 

scobj <- readRDS(file = "output/scobj_combined.rds")

#### umap ----
Idents(scobj) <- "Condition"
plot1 <- DimPlot(scobj, label = F, pt.size = 0.1, 
                 cols = c("#5F9EA0","#DE2D26"))+labs(x = "UMAP1", y = "UMAP2") +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
plot1
Idents(scobj) <- "Celltype"
pal <- paletteer_d("ggsci::nrc_npg")
plot2 <- DimPlot(scobj, label = T, pt.size = 0.1, repel = F,shuffle = T, split.by = "Condition",
                 cols = pal)+labs(x = "UMAP1", y = "UMAP2") +
  NoLegend()+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
plot2

##### marker genes ----
table(scobj$Celltype)
Idents(scobj) <- "Celltype"
DefaultAssay(scobj) <- "SCT"
marker_genes <- c("KRT1","KRT5","KRT15",
                  "TYRP1","PMEL","DCT",
                  "COL1A1","COL3A1","DCN",
                  "RERGL","MYH11","TAGLN",
                  "RGS5","PDGFRB","CALD1",
                  "PECAM1","CLDN5","TM4SF1",
                  "TPSB2","TPSAB1","CPA3",
                  "LYZ","CD14","CD68",
                  "IL7R","CD3D","TRAC")
data.usage <- DotPlot(scobj, features = marker_genes,cols = c("white", "#CD1076"), col.min = 0
)+ RotatedAxis() +
  theme(panel.border = element_rect(color="black", size=0.5, linetype="solid"))
data.usage$data
data.anno <- data.frame(features.plot = unique(data.usage$data$features.plot),
                        label= c(rep("Keratinocyte",3), rep("Melanocyte",3),
                                 rep("Fibroblast",3),rep("Smooth muscle",3),
                                 rep("Pericyte",3), rep("Endothelium",3),
                                 rep("Mast",3),rep("Myeloid",3),
                                 rep("Lymphoid",3)))
df.plot <- plyr::join(data.usage$data,data.anno)
df.plot$id <- factor(df.plot$id,levels = c("Keratinocyte","Melanocyte","Fibroblast",
                                           "Smooth muscle","Pericyte","Endothelium",
                                           "Mast","Myeloid","Lymphoid"))
df.plot$id
df.plot$label <- factor(df.plot$label,levels = c("Keratinocyte","Melanocyte","Fibroblast",
                                                 "Smooth muscle","Pericyte","Endothelium",
                                                 "Mast","Myeloid","Lymphoid"))
p <- ggplot(df.plot,aes(x=features.plot,y =   as.numeric(id),size = pct.exp, fill = avg.exp.scaled) )+
  geom_point(shape= 21, stroke = 1,color= "black") + 
  scale_size("Percent Expressed", range = c(0,7)) + 
  coord_cartesian(#xlim = c(0, 25),
    ylim = c(0.5, 9.5),expand = TRUE)+
  scale_fill_gradientn(colours = c("grey99", "#CD1076"), 
                       guide = guide_colorbar(ticks.colour = "black",frame.colour = "black"),
                       name = "Average expression") +
  cowplot::theme_cowplot() + 
  ylab("") + xlab("") + 
  scale_y_continuous(breaks = 1:length(levels(df.plot$id)),labels = levels(df.plot$id),sec.axis = dup_axis())+ 
  facet_grid(~label, scales="free_x",space = "free")+theme_classic() +
  theme_bw() + 
  theme(
    panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
    axis.text.x = element_text(size=12, angle=45, hjust=1, color="black"),
    axis.text.y = element_text(size=12,colour = "black"),  
    axis.title.x = element_text(size=14,colour = 'black',vjust = -0.8,hjust = 0.5),
    
    axis.ticks.y.right = element_blank(),
    axis.text.y.right = element_blank(),
    axis.line = element_line(colour = 'grey30',size = 0.5), 
    panel.spacing=unit(0, "mm"), 
    strip.text.x = element_text(size=10, face="bold",color = "#FFFFFF",
                                vjust = 0.5,margin = margin(b = 3,t=3)),
    strip.background = element_rect(colour="grey30", fill="grey60",size = 0.5)
  )
p
g <- ggplot_gtable(ggplot_build(p))
strips <- which(grepl('strip-', g$layout$name))
cols <- paletteer_d("ggsci::nrc_npg")
for (i in seq_along(strips)) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- cols[i]
}
plot(g)

##### cell type nDEGs ----
scobj$celltype_condition <- paste(scobj$Celltype, scobj$Condition, sep = "_")
Idents(scobj) <- "celltype_condition"
metadata <- scobj@meta.data
scobj <- PrepSCTFindMarkers(object = scobj)

marker_condition=data.frame()
for ( ci in sort(as.character(unique(scobj@meta.data$Celltype))) ) {
  tmp.marker <- FindMarkers(
    scobj, logfc.threshold = 0.25, min.pct = 0.1,
    only.pos = F, test.use = "wilcox",
    ident.1=paste0(ci,"_VU"),ident.2=paste0(ci,"_NS")
  )
  
  tmp.marker$gene=rownames(tmp.marker)
  tmp.marker$condition=ifelse(tmp.marker$avg_log2FC > 0,paste0(ci,"_VU"),paste0(ci,"_NS"))
  tmp.marker$cluster=ci
  
  tmp.marker=tmp.marker%>%filter(p_val_adj < 0.01)
  tmp.marker=as.data.frame(tmp.marker)
  tmp.marker=tmp.marker%>%arrange(desc(avg_log2FC))
  
  marker_condition=marker_condition%>%rbind(tmp.marker)
}
FB_Condtion <- as.data.frame(table(marker_condition$condition))
write.csv(FB_Condtion,file = "FB_DIFF_condition.csv")

df1 = as.data.frame(table(marker_condition$cluster,marker_condition$condition))
colnames(df1) = c('celltype','celltype_condition','nDEGs')
df2 <- df1[df1$nDEGs != 0,]
colnames(df2) = c('Celltype','Celltype_condition','nDEGs')

df_up <- df2[df2$Celltype_condition %in% c("Keratinocyte_VU","Melanocyte_VU","Fibroblast_VU",
                                           "Smooth muscle_VU","Smooth muscle_VU","Pericyte_VU",
                                           "Endothelium_VU","Mast_VU",
                                           "Myeloid_VU","Lymphoid_VU"),]
df_down <- df2[! df2$Celltype_condition %in% c("Keratinocyte_VU","Melanocyte_VU","Fibroblast_VU",
                                               "Smooth muscle_VU","Smooth muscle_VU","Pericyte_VU",
                                               "Endothelium_VU","Mast_VU",
                                               "Myeloid_VU","Lymphoid_VU"),]

scobj$upregulated_nDEGs = df_up[match(scobj$celltype, df_up$Celltype),3]
scobj$downregulated_nDEGs = df_down[match(scobj$celltype, df_down$Celltype),3]
metadata <- scobj@meta.data

library(patchwork)
Idents(scobj) <- "celltype"
p1 <- FeaturePlot(scobj,'upregulated_nDEGs',cols = c("lightgrey", "#ef3b2c"))+labs(x = "UMAP1", y = "UMAP2") +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p2 <- FeaturePlot(scobj,'downregulated_nDEGs',cols = c("lightgrey", "#427bbd"))+labs(x = "UMAP1", y = "UMAP2") +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p1|p2



##### cell proportion ----
source("R/Singlecellratio_plotstat.R")
my_comparisons <- list(c("NS", "VU"))
Singlecellratio_plotstat(scobj, group_by = "Condition",
                         meta.include = c("Condition","orig.ident"),
                         comparisons = my_comparisons, color_by = 'Condition',
                         group_by.point = "orig.ident",label.x = 1, pt.size = 3,
                         label = 'p.format', ncol =3)
