rm(list = ls())
library(Seurat)
library(harmony)
library(tidyverse)
library(patchwork)
library(dittoSeq)
library(cowplot)
library(paletteer)
library(ComplexHeatmap)
library(circlize)
library(scRNAtoolVis)
library(RColorBrewer)
library(enrichplot)
library(clusterProfiler)

scobj <- readRDS(file = "output/scobj_subEC_SMC.rds")

#### umap ----
Idents(scobj) <- "Condition"
plot1 <- DimPlot(scobj, label = F, pt.size = 0.5, 
                 cols = c("#5F9EA0","#DE2D26"))+labs(x = "UMAP1", y = "UMAP2") +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
plot1

Idents(scobj) <- "celltype"
Idents(scobj) <- factor(Idents(scobj), 
                        levels = c('Arterial EC','Venular EC','Capillary EC','Lymphatic EC',
                                   'Smooth muscle','Pericyte'))
plot2 <- DimPlot(scobj, label = T, pt.size = 0.5, repel = F,shuffle = T,split.by = "Condition",
                 cols = c("#20B2AA","#FFA500","#9370DB","#98FB98",
                          "#3C5488FF",'#F39B7FFF'))+labs(x = "UMAP1", y = "UMAP2") +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
plot2

#### marker genes ----
library(Nebulosa)
DefaultAssay(scobj) <- "RNA"
marker_genes <- c("SEMA3G","SELE",'APLN','PROX1',
                  "MYH11","DCN")
plot_density(scobj, features = marker_genes, size = 0.5, legend_title = FALSE) +
  plot_layout(ncol = 3) +
  labs(x = "UMAP1", y = "UMAP2") +
  theme(
    panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"  
  )

#### GSEA ----
sce_vessel <- subset(scobj,idents = c('Venular EC','Arterial EC','Capillary EC'))
Idents(sce_vessel) <- "Condition"
diff <- FindMarkers(object = sce_vessel,
                    ident.1 = "VU",
                    ident.2 = "NS",
                    assay = "SCT",
                    recorrect_umi = FALSE,
                    logfc.threshold = 0, 
                    min.pct = 0.1,
                    test.use = "wilcox")
gene_df <- diff
geneList <- gene_df$avg_log2FC
names(geneList) =  rownames(gene_df)
geneList = sort(geneList, decreasing = TRUE)

head(geneList)

genesets  <- read.gmt("resource/c5.go.bp.v2023.2.Hs.symbols.gmt")
y <- GSEA(geneList,TERM2GENE = genesets)
yd <- as.data.frame(y)
library(GseaVis)
geneSetID = c('GOBP_RESPONSE_TO_CYTOKINE',
              "GOBP_CYTOKINE_PRODUCTION",
              'GOBP_LYMPHOCYTE_ACTIVATION',
              "GOBP_T_CELL_ACTIVATION",
              'GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_II',
              "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN"
)
# add NES and Pvalue
gseaNb(object = y,
       geneSetID = geneSetID,
       subPlot = 2,
       termWidth = 35,
       legend.position = c(0.8,0.8),
       #addGene = gene,
       addPval = T,
       curveCol =  jjAnno::useMyCol('paired',6),
       pvalX = 0.05,pvalY = 0.05)

#### EC DEGs ----
scobj_EC <- subset(scobj, idents = c('Arterial EC','Venular EC','Capillary EC'))
markers <- c("CCL2",'CXCL2','CXCL12','CXCL14','IL6','IL33','ICAM1','TNFAIP3','NFKBIA',
             'CD9','CD36','CD55','CD74','HLA-B','HLA-DMA','HLA-DPA1','HLA-DPB1','HLA-DQA1',
             'HLA-DQB1','HLA-DRA','HLA-DRB1','HLA-F','HSP90AB1','HSPA1A','HSPA1B','HSPD',
             'COL4A1','COL4A2','SPARC','ENG','POSTN','PLVAP',
             'FABP5','FLNB','APLN','ITGB1','LAMA4','LAMB1','VCL','VWA1','VWF',"VCAM1"
)
DefaultAssay(scobj_EC) <- "SCT"
gene_cell_exp <- AverageExpression(scobj_EC,
                                   features = markers,
                                   group.by = 'orig.ident',
                                   slot = 'data') 
gene_cell_exp <- as.data.frame(gene_cell_exp$SCT)
marker_exp <- t(scale(t(gene_cell_exp),scale = T,center = T))
exp <- t(marker_exp)
df <- data.frame(colnames(gene_cell_exp))
colnames(df) <- 'class'
df$condition <- c(rep('NS',5), rep('VU',4))
df <- as.data.frame(df[,-1])
colnames(df) <- 'condition'
annotation_row = data.frame(
  group = c(rep('NS',5), rep('VU',4)))

row.names(annotation_row) <- rownames(exp)
ann_color = list(group = c('NS'="#5F9EA0",
                           'VU'="#DE2D26"),
                 Change= c('up'='#F3B1A0',
                           'down'='#C1E6F3'),
                 Function = c("Inflammation"='#FEE066',
                              "Immune"='#FF5E00',
                              "Remodeling"="#5CB3BA"))
col = colorRamp2(c(-2,0,2),c("#96CCCB",'grey98',"#F6CAE5"))
annotation_col = data.frame(
  Change = c(rep("down",9), rep("down",16), rep("up", 16)),
  Function = c(rep("Inflammation",9), rep("Immune",16),rep("Remodeling",16))
)
annotation_col$Change <- factor(annotation_col$Change, levels = c('up','down'))
annotation_col$Function <- factor(annotation_col$Function,levels = c("Remodeling",
                                                                     "Immune",
                                                                     "Inflammation" ))
row.names(annotation_col) <- colnames(exp)
str(annotation_col)
pheatmap(exp, scale = "column",
         show_rownames = F,
         show_colnames = T,
         cluster_rows = F,
         cluster_cols = F,
         col = col,
         border= "grey",
         border_color = NA,
         annotation_col = annotation_col, 
         annotation_row = annotation_row,
         annotation_colors = ann_color,
         row_split = annotation_row$group,
         column_split = annotation_col$Function,
         annotation_names_row = F,
         annotation_names_col = F ,
         column_title = NULL,
         row_title = NULL)
