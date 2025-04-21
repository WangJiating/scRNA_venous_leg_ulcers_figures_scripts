rm(list = ls())
library(Seurat)
library(harmony)
library(tidyverse)
library(patchwork)
library(dittoSeq)
library(cowplot)
library(paletteer)
library(scRNAtoolVis)
library(Nebulosa)
library(clusterProfiler)
library(GseaVis)
library(pheatmap)
library(ComplexHeatmap)
gc()

scobj <- readRDS("output/scobj_subFB.rds")

#### subFB umap ----
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
plot2 <- DimPlot(scobj, label = T, pt.size = 0.5, repel = F,shuffle = T,
                 cols = c("#AEC7E8","#FFBB78","#9467BD",
                          "#FF7F0E","#C49C94"))+labs(x = "UMAP1", y = "UMAP2") +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
plot1+plot2
  
#### cell density ----
library(ggpointdensity) 
scRNA <- readRDS("output/scobj_subFB.rds")
data <- cbind(Embeddings(object=scRNA[['umap']]),
              FetchData(scRNA,c("Condition","orig.ident")))
#plot
ggplot(data = data, mapping = aes(x = UMAP_1,y = UMAP_2)) + 
  geom_pointdensity(adjust=6) + 
  geom_density_2d(bins = 5, colour="black")+
  facet_wrap(~Condition, scales = 'free_y')+
  theme_bw()+
  theme(
    axis.text=element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    legend.text = element_text(size =10),
    aspect.ratio=1,
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.x = element_text(size=12, color = "black",
                                vjust = 0.5,margin = margin(b = 3,t=3)), 
    strip.background = element_blank(),
    plot.margin=unit(c(1, 1, 1, 1),'cm'),
    legend.margin = margin(-0.2,-0.2,0,0,'cm'),
    legend.key.height = unit(1,'cm'),
    legend.key.width = unit(0.4,'cm'))+
  scale_color_gradientn(colours = c('#5749a0', '#0f7ab0', '#00bbb1',
                                    '#bef0b0', '#fdf4af', '#f9b64b',
                                    '#ec840e', '#ca443d', '#a51a49'), 
                        guide = guide_colorbar(frame.colour = "black", 
                                               ticks.colour = NA),
                        name = "Density",
                        labels = c('low',"high"),
                        breaks = c(100,700)) 

#### subFB proportion ----
source("R/Singlecellratio_plotstat.R")
my_comparisons <- list(c("NS", "VU"))
head(Idents(scobj))
Idents(scobj) <- "celltype"
scobj$celltype <- Idents(scobj)
Singlecellratio_plotstat(scobj, group_by = "Condition",
                         meta.include = c("Condition","orig.ident"),
                         comparisons = my_comparisons, color_by = 'Condition',
                         group_by.point = "orig.ident",label.x = 1, pt.size = 3,
                         label = 'p.format', ncol =5)

#### subFB marker genes ----
marker_genes <- c("ADAM12","ASPN","CCL19",'APOE',"APCDD1","COL18A1",
                  "JUN","IER5",
                  "MKI67","PCLAF")
scobj$celltype <- factor(scobj$celltype, levels = c("subFB5","subFB4","subFB3","subFB2","subFB1"))
scobj$Celltype_condition <- paste(scobj$sub_FB,scobj$Condition,sep = "_")
jjDotPlot(object =scobj,
          gene = marker_genes,
          rescale.min = -2,
          rescale.max = 2,
          midpoint = 0,
          dot.col = c("#1B7837","#F7F7F7", "#762A83"), #,"#91CAE8","#990000"
          id = 'Celltype_condition',
          anno = T,
          ytree = FALSE,
          cluster.order = cr,
          split.by.aesGroup = T)

marker_genes <- c("ADAM12","CCL19","APCDD1","JUN",
                  "MKI67","ASPN","APOE","COL18A1","IER5","PCLAF")
plot_density(scobj,features = marker_genes,size = 1,legend_title = FALSE,
             pal =  'magma')+plot_layout(ncol = 5)&NoLegend()&labs(x = "UMAP1", y = "UMAP2") &
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

#### GSEA ----
Idents(scobj) <- "Condition"
scobj <- PrepSCTFindMarkers(object = scobj)
gene_df <- FindMarkers(scobj, 
                        ident.1 = "VU",
                        ident.2 = "NS", 
                        logfc.threshold = 0,
                        test.use = "wilcox")
geneList <- gene_df$avg_log2FC
names(geneList) =  rownames(gene_df)
geneList = sort(geneList, decreasing = TRUE)
head(geneList)

genesets <- clusterProfiler::read.gmt("resource/h.all.v2022.1.Hs.symbols.gmt")
y <- GSEA(geneList,TERM2GENE = genesets)
yd <- as.data.frame(y)
geneSetID = c('HALLMARK_TNFA_SIGNALING_VIA_NFKB',
              "HALLMARK_IL2_STAT5_SIGNALING",
              "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
gseaNb(object = y,
       geneSetID = geneSetID,
       subPlot = 2,
       termWidth = 35,
       legend.position = c(0.8,0.8),
       #addGene = gene,
       addPval = T,
       curveCol =  c("#4370B4",'#549f9A','#C30078'),
       pvalX = 0.05,pvalY = 0.05)

#### DEGs heatmap ----
markers <- c("COL1A1","COL1A2","COL3A1","COL4A1","COL4A2","COL5A1","COL5A2","COL5A3","COL6A1","COL6A2",
             "COL7A1","COL8A1","COL10A1","COL11A1","COL12A1","COL27A1",   #  up
             "COL6A5","COL6A6","COL14A1","COL15A1","COL23A1",             #  down
             'IL11','LRRC15',"MMP1","MMP2","MMP3","MMP11","MMP14",
             'WNT2','WNT5A','GREM1','PTK7','CPE','PTEN','IGFBP2','SEMA5A',
             'MT1A','MT1E','MT1G','MT1M','MT1X','MT2A',
             'CCL3','CCL26',"CXCL1","CXCL3","CXCL5","CXCL6","CXCL8","CXCL13",'CHI3L1','CHI3L2', 'IL24', #  up
             'CCL2','CCL19','CXCL2',"CXCL12","CXCL14",
             'C3',"IL15",'IL33','APOD','APOE','TNC','NAMPT','NFKBIA','NFIA','NFIB','NFIX','PI16',
             'FAP','FABP5','FN1','POSTN','COMP','ASPN','ADAM12','LUM','CTHRC1','ACTA2','MFAP2'
)
DefaultAssay(scobj) <- "SCT"
gene_cell_exp <- AverageExpression(scobj,
                                   features = markers,
                                   group.by = 'orig.ident',
                                   slot = 'data') 
gene_cell_exp <- as.data.frame(gene_cell_exp$SCT)
marker_exp <- t(scale(t(gene_cell_exp),scale = T,center = T))
exp <- t(marker_exp)
annotation_row = data.frame(
  group = c(rep('NS',5), rep('VU',4)))
cols =RColorBrewer::brewer.pal(9, "Pastel1")
row.names(annotation_row) <- rownames(exp)
ann_color = list(group = c('NS'="#5F9EA0",
                           'VU'="#DE2D26"),
                 Change= c('up'='#F3B1A0',
                           'down'='#C1E6F3'),
                 Function = c("Collagen"="#B3CDE3",
                              'ECM' = "#CCEBC5",
                              'Metallothionein' = "#DECBE4",
                              'Wnt' = "#FED9A6",
                              "Immunosuppression"="#FFFFCC",
                              "Inflammation"="#E5D8BD"))

annotation_col = data.frame(
  Change = c(rep("up",16), rep("down", 5),rep("up",15), rep("down", 6),
             rep("up",11), rep("down", 17),
             rep("up",11)),
  Function = c(rep("Collagen",21), rep("Immunosuppression",7),rep("Wnt",8),
               rep("Metallothionein",6),rep("Inflammation",28),rep("ECM",11))
)
annotation_col$Change <- factor(annotation_col$Change, levels = c('up','down'))
annotation_col$Function<- factor(annotation_col$Function, 
                                 levels = c("Collagen","ECM",'Metallothionein',
                                            'Wnt','Immunosuppression',
                                            "Inflammation"
                                 ))
str(annotation_col)
col = colorRamp2(c(-2,0,2),c("#148f28","grey98","#EA71AE"))
pheatmap(exp, scale = "col",
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
         column_split = annotation_col$Function,
         row_split = annotation_row$group,
         annotation_names_row = F,
         annotation_names_col = F ,
         column_title = NULL,
         row_title = NULL)

#### GO function ----
scobj$celltype <- factor(scobj$celltype,c("subFB1", "subFB2", "subFB3", "subFB4", "subFB5"))
Idents(scobj) <- "celltype"
scobj <- PrepSCTFindMarkers(object = scobj)
scobj.markers <- FindAllMarkers(object = scobj,
                                logfc.threshold = 0.25,
                                min.pct = 0.25,
                                test.use = "wilcox",
                                only.pos = TRUE)
top50_subFB <-scobj.markers %>%
  group_by(cluster) %>%
  top_n(n=-50,wt=p_val_adj)
gid <- bitr(unique(top50_subFB$gene), 'SYMBOL', 'ENTREZID', OrgDb= 'org.Hs.eg.db')
colnames(gid)[1] <- "gene"
markers <- merge(top50_subFB, gid, by='gene')
library(org.Hs.eg.db)
y = compareCluster(ENTREZID ~ cluster, OrgDb="org.Hs.eg.db",data = markers, fun='enrichGO',
                   ont = "BP")
go <- as.data.frame(y@compareClusterResult)
go <-  go %>% 
  dplyr::filter(Count>5)
go_top5 <- go[c(1:2,4,7,10,39,46,48:50,74,80,81,97,112,134,137:140,152:156),]
go_top5$Description <- factor(go_top5$Description)
go_top5$group = c(25:1)
go_top5$Description <- reorder(go_top5$Description,go_top5$group)
go_top5$Description <- factor(go_top5$Description ,levels = 
                                rev(c("organelle fission", "nuclear division", "chromosome segregation", 
                                      "nuclear chromosome segregation", "spindle organization", "cellular response to oxidative stress", 
                                      "response to steroid hormone", "intrinsic apoptotic signaling pathway", 
                                      "RNA splicing", "response to heat", "canonical Wnt signaling pathway", 
                                      "cell-cell signaling by wnt", "Wnt signaling pathway", "response to hypoxia", 
                                      "leukocyte migration", "leukocyte chemotaxis", "regulation of leukocyte migration", 
                                      "mononuclear cell migration", "regulation of lipid metabolic process", 
                                      "ossification", "collagen metabolic process", 
                                      "collagen fibril organization", "extracellular structure organization",
                                      "extracellular matrix organization"
                                )))
go_top5$Cluster <- factor(go_top5$Cluster,levels = rev(levels(go_top5$Cluster)))
ggplot(go_top5, aes(Cluster, y = Description)) +
  geom_point(aes(color=qvalue, size=Count))+
  theme_bw()+
  ylab("GO term") + xlab("subFB") + 
  scale_color_continuous(low="#6EE2FFFF", high="#D0DFE6FF")+
  theme(
    axis.text.x = element_text(size=12,angle=45,hjust=1, color="black"),
    axis.text.y = element_text(size=12,colour = "black"), 
    axis.title.x = element_text(size=12,colour = 'black',vjust = -0.8,hjust = 0.5),
    axis.title.y = element_text(size=12,colour = 'black',vjust = -0.8,hjust = 0.5),
    legend.position="right", legend.box = ""
  )+ scale_size(range=c(4, 6))+
  coord_flip()
