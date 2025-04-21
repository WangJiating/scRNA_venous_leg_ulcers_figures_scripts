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
library(ClusterGVis)
library(org.Hs.eg.db)

scobj <- readRDS("output/scobj_KC.rds")

##### marker genes ----
Idents(sce) <- sce$celltype
levels(Idents(sce)) 
sce <- PrepSCTFindMarkers(object = sce)
sce.markers.all <- Seurat::FindAllMarkers(sce,
                                          only.pos = TRUE,
                                          min.pct = 0.1,
                                          logfc.threshold = 0.25)

# get top 40 genes
sce.markers <- sce.markers.all %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 40, wt = avg_log2FC)

# check
head(sce.markers)

# prepare data from seurat object
## ??ClusterGVis::prepareDataFromscRNA
st.data <- prepareDataFromscRNA(object = sce,
                                assays = "SCT",
                                diffData = sce.markers,
                                showAverage = TRUE)

# check
str(st.data)

# enrich for clusters
enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = 5,
                        seed = 5201314)

head(enrich)
markGenes <- c("COL17A1","LAMB4","COL7A1","DST","KRT15",
               "TOP2A","MKI67","NUSAP1","UBE2C","HMGB2",
               "KRT6C","KRT6A","LY6D","KRT6B","FGF7","GMNN",
               "IVL","KLK7","DSC1","ARL5A","KRT78")

head(sce.markers)
# line plot
visCluster(object = st.data,
           plot.type = "line")
visCluster(object = st.data,
           plot.type = "heatmap",
           #ncol = 20,
           column_names_rot = 45,
           markGenes = markGenes,
           cluster.order = c(1:8))
dev.off()

Idents(sce)
# add bar plot
??ClusterGVis::visCluster

order <-  c('Basal','Proliferating',
            'Spinous',"Granular")
library(ggsci)
library(scRNAtoolVis)
library(RColorBrewer)

cols = c( "#FF95A8FF" ,"#FF6348FF", "#84D7E1FF" ,"#ADE2D0FF" )
#pdf('sc3.pdf',height = 10,width = 14,onefile = F)
visCluster(object = st.data,
           plot.type = "both",
           ctAnno.col = cols,
           boxcol = cols,
           column_names_rot = 45,
           show_row_dend = F,
           markGenes = markGenes,
           markGenes.side = "left",
           annoTerm.data = enrich,
           line.side = "left",
           sample.order = order,
           sample.col = cols,
           cluster.order = c(1:4),
           go.col = rep(cols,each = 5),
           add.bar = T)
dev.off()

#### KC DEGs ----
Idents(sce) <- 'celltype'
DimPlot(sce)
marker <- c("KRT1","KRT2","KRT10","KRT15","KRT31",
            "KRT6A","KRT6B","KRT6C","KRT14","KRT16","KRT17","KRT18",
            "S100A2","S100A6","S100A7","S100A7A","S100A8",
            "S100A9","S100A11","S100A13","S100A16",
            'CCNB1','CENPF','H2AZ1','H1-1','H4C3','MKI67','PTTG1','UBE2S','CDC20',
            "HLA-A", "HLA-B","HLA-F","B2M",
            "CCL2","CCL27","CXCL14","IL18","IL34",
            "NFIB","NFIX","NFKBIA","NFKBIZ",
            "TNFRSF18","TNFRSF19","TNFAIP3"
)
DefaultAssay(sce) <- "SCT"
gene_cell_exp <- AverageExpression(sce,
                                   features = marker,
                                   group.by = 'orig.ident',
                                   slot = 'data') 
gene_cell_exp <- as.data.frame(gene_cell_exp$SCT)
marker_exp <- t(scale(t(gene_cell_exp),scale = T,center = T))
exp <- marker_exp
df <- data.frame(colnames(gene_cell_exp))
colnames(df) <- 'class'
df$condition <- c(rep('NS',5), rep('VU',4))
df <- as.data.frame(df[,-1])
colnames(df) <- 'condition'
top_anno = HeatmapAnnotation(df = df,
                             border = F,
                             show_annotation_name = F,
                             gp = gpar(col = NA),
                             col = list(condition = c('NS'="#5F9EA0",
                                                      'VU'="#DE2D26")))
range(marker_exp)
col = colorRamp2(c(-2,0,2),c("#91CAE8","grey98","#F48892"))
annotation_col = data.frame(
  group = c(rep('NS',5), rep('VU',4)))
rownames(annotation_col) <- colnames(marker_exp)
ann_color = list(group = c('NS'="#5F9EA0",
                           'VU'="#DE2D26"),
                 Change= c('up'='#F3B1A0',
                           'down'='#C1E6F3'),
                 Function = c("Keratin"=  "#BD9AAD",
                              "AMP"= "#E8D2B3" ,
                              "HLA"="#6D65A3",
                              'Proliferation' ="#9193B4",
                              "Inflammation"='#FEE066'))

annotation_row = data.frame(
  Change = c(rep("down",5), rep("up", 7),rep("up",18), rep("down", 16)),
  Function = c(rep("Keratin",12), rep("AMP",9),rep("Proliferation",9),rep("HLA",4),rep("Inflammation",12))
)
annotation_row$Function <- factor(annotation_row$Function,levels = c("Keratin","AMP",
                                                                     "Proliferation",
                                                                     "HLA","Inflammation"))
rownames(annotation_row) <- rownames(marker_exp)
str(annotation_row)
pheatmap(marker_exp, scale = "row",
         show_rownames = T,
         show_colnames = F,
         cluster_rows = F,
         cluster_cols = F,
         col = col,
         border= "grey",
         border_color = NA,
         annotation_col = annotation_col, 
         annotation_row = annotation_row,
         annotation_colors = ann_color,
         row_split = annotation_row$Function,
         column_split = annotation_col$group,
         annotation_names_row = F,
         annotation_names_col = F ,
         column_title = NULL,
         row_title = NULL)

