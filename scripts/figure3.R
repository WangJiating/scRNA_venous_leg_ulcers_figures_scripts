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
library(ggsci)
library(RColorBrewer)

sce <- readRDS("output/scobj_IM_outKC.rds")

#### immune cells UMAP ----
Idents(sce) <- "Condition"
plot1 <- DimPlot(sce, label = F, pt.size = 0.5, split.by = "Condition",
                 cols = c("#5F9EA0","#DE2D26"))+labs(x = "UMAP1", y = "UMAP2") +
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
plot1

Idents(sce) <- "Celltype"
sce$Celltype <- factor(sce$Celltype,levels = c("LC",'cDC1','cDC2','DC3','pDC',
                                               "Macro1","Macro2","Macro3","Mast","Neutro",
                                               "Th","Treg","Tc",'NK','B','Plasma'))
cols = c(paletteer_d("ggsci::nrc_npg")[7],pal_igv()(20))
plot2 <- DimPlot(sce, label = F, pt.size = 0.5, repel = F,shuffle = T,split.by = "Condition",
                 cols = c(paletteer_d("ggsci::nrc_npg")[7],pal_igv()(20)))+labs(x = "UMAP1", y = "UMAP2") +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
plot2

#### marker genes ----
marker_genes <- c("CD1A","CD207",'FCGBP',             
                  "CLEC9A","CLIC2","CPNE3",           
                  "FCGR2B","CD1C","CLEC10A",          
                  'LAMP3','CCL22','FSCN1',           
                  'IL3RA',"TCF4","SOX4",              
                  'C1QA','C1QB','C1QC',             
                  'MMP9','MMP19','CTSB',              
                  'IL1B','THBS1','EREG',             
                  "TPSB2","TPSAB1","CPA3",            
                  "FCGR3B","CSF3R","CXCL8",           
                  "CD40LG","IL7R","CD3D",             
                  "FOXP3","CTLA4","TIGIT",            
                  "CD8A","GZMK","GZMA",                
                  "GNLY","KLRD1",'NKG7',               
                  "CD79B","CD79A","MS4A1",           
                  "MZB1","IGKC","IGHG1"              
                  
)
data.usage <- DotPlot(sce, features = marker_genes,cols = c("white", "#CD1076"), col.min = 0
)+ RotatedAxis() +
  theme(panel.border = element_rect(color="black", size=0.5, linetype="solid"))
data.anno <- data.frame(features.plot = unique(data.usage$data$features.plot),
                        label= c(rep("LC",3),
                                 rep("cDC1",3), rep("cDC2",3),
                                 rep("DC3",3),
                                 rep("pDC",3),rep("Macro1",3),
                                 rep("Macro2",3), rep("Macro3",3),
                                 rep("Mast",3), 
                                 rep("Neutro",3),rep("Th",3),
                                 rep("Treg",3), rep("Tc",3),
                                 rep("NK",3),rep("B",3),
                                 rep("Plasma",3)))
df.plot <- plyr::join(data.usage$data,data.anno)
df.plot$id <- factor(df.plot$id,levels =   c("LC",'cDC1','cDC2','DC3','pDC',
                                             "Macro1","Macro2","Macro3","Mast","Neutro",
                                             "Th","Treg","Tc",'NK','B','Plasma'))
df.plot$id
df.plot$label <- factor(df.plot$label,levels =  c("LC",'cDC1','cDC2','DC3','pDC',
                                                  "Macro1","Macro2","Macro3","Mast","Neutro",
                                                  "Th","Treg","Tc",'NK','B','Plasma'))


p <- ggplot(df.plot,aes(x=features.plot,y =   as.numeric(id),size = pct.exp, fill = avg.exp.scaled) )+
  geom_point(shape= 21, stroke = 1,color= "black") + 
  scale_size("Percent Expressed", range = c(0,6)) + 
  coord_cartesian(#xlim = c(0, 25),
    ylim = c(1, 16),expand = TRUE)+
  scale_fill_gradientn(colours = c("grey99", "#CD1076"), 
                       guide = guide_colorbar(ticks.colour = "black",frame.colour = "black"),
                       name = "Average expression") +
  cowplot::theme_cowplot() + 
  ylab("") + xlab("") + 
  scale_y_continuous(breaks = 1:length(levels(df.plot$id)),labels = levels(df.plot$id),sec.axis = dup_axis())+ #复制 y轴 代替边框效果
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
    strip.background = element_rect(colour="grey30", fill="grey60",size = 0.5),
    legend.position="top", legend.box = ""
  )
p
g <- ggplot_gtable(ggplot_build(p))
strips <- which(grepl('strip-', g$layout$name))
cols <- cols
for (i in seq_along(strips)) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- cols[i]
}
plot(g)


#### cell proportion ----
source("R/Singlecellratio_plotstat.R")
my_comparisons <- list(c("NS", "VU"))
sce$celltype <- sce$Celltype
Singlecellratio_plotstat(sce, group_by = "Condition",
                         meta.include = c("Condition","orig.ident"),
                         comparisons = my_comparisons, color_by = 'Condition',
                         group_by.point = "orig.ident",label.x = 1, pt.size = 3,
                         label = 'p.format', ncol =8)


#### DC cells ----
sce_DC <- subset(sce, Celltype %in% c('LC','cDC1','cDC2','DC3','pDC'))
Idents(sce_DC) <- "Condition"
DC.markers <- FindMarkers(object = sce_DC,
                          ident.1 = "VU",
                          ident.2 = "NS",
                          logfc.threshold = 0.25, min.pct = 0.1,
                          test.use = "wilcox")
DC.markers$gene <- rownames(DC.markers)
features <- c("ITGAX","CD1A","CD1B","CD1C","CD1E","CD207",
              "HLA-DQA1", "HLA-DRB1", "HLA-DQB2","HLA-DQB1", 
              "HLA-DMB","IL18","IL23A","TNF")
whitePurple = c("#F7F7F7", '#00bbb1')
jjDotPlot(object =sce_DC,
          gene = features,
          dot.col = whitePurple,
          id = 'celltype_condition',
          split.by = 'Condition',
          ytree = FALSE,
          split.by.aesGroup = T)


genesets <- read.gmt("resource/c5.go.bp.v2023.2.Hs.symbols.gmt")
signatures <- split(genesets$gene,genesets$term)
features <- signatures[c("GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION",
                         "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_ENDOGENOUS_LIPID_ANTIGEN_VIA_MHC_CLASS_IB",
                         "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_II")]
Idents(sce_DC) <- "Condition"
sce_DC <- AddModuleScore(sce_DC, 
                         features = features, 
                         nbin = 15,
                         seed = 1234,
                         name = "seurat", 
                         assay = "RNA")
test <- sce_DC@meta.data
names(sce_DC@meta.data)[grep("seurat\\d",names(sce_DC@meta.data))] <- names(features)
library(ggpubr)
df<- FetchData(sce_DC,vars = c("Condition",
                               "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION"))
df$Condition <- factor(df$Condition,levels = c("NS","VU"))
colnames(df) = c('Condition','Pathway')
my_comparisons <- list(c("VU", "NS")) 
ggplot(df,aes(x=Condition ,y=Pathway, fill=Condition))+
  geom_violin(color='black',size=0.3,color=NA,alpha=0.5)+
  theme_classic() + 
  labs(title = "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION")+
  theme(text = element_text(size=10, colour = "black")) + 
  theme(plot.title = element_text(hjust = 0.5, size = 10)) +
  theme(legend.position="none") +  
  geom_jitter(width = 0.35,
              height =0, size = 1.5, shape = 21, 
              stroke = 0.3, 
              show.legend = FALSE)+
  geom_boxplot(width=0.2, outlier.shape = NA) +
  scale_fill_manual(values = c("VU" = "#F48892", "NS" = "#91CAE8"))+
  stat_compare_means(method="wilcox",hide.ns = F,
                     comparisons =c(my_comparisons),
                     label="p.signif",
                     bracket.size=0.8,
                     size=6)


#### Macrophage -----
sce_Macro <- subset(sce, Celltype %in% c('Macro1','Macro2','Macro3'))
Idents(sce_Macro) <- "Condition"
Macro.markers <- FindMarkers(object = sce_Macro,
                             ident.1 = "VU",
                             ident.2 = "NS",
                             logfc.threshold = 0.25, min.pct = 0.1,
                             test.use = "wilcox")
features <- c('IL10','TGFBI',
              'CD68','CD163',
              'C1QA','C1QB','C1QC',
              'CTSB','CTSH','CTSS',
              'MMP9','MMP19','FN1', 'SPP1',
              'IL1A','IL1B','IL6',
              'IL23A','CXCL8','THBS1','PTGS2','CD86','EREG','AREG','NLRP3')
jjDotPlot(object = sce_Macro,
          gene = features,
          id = 'celltype',
          split.by = 'Condition',
          split.by.aesGroup = T,
          dot.col = c( "#1B7837","#F7F7F7","#762A83"), 
          ytree = FALSE,
          point.geom = F,
          tile.geom = T,
          legend.position="top")

#### T cells ----
sce_Tcell <- subset(sce, Celltype %in% c('Th','Treg','Tc','NK'))
Idents(sce_Tcell) <- "Condition"
T.markers <- FindMarkers(object = sce_Tcell,
                         ident.1 = "VU",
                         ident.2 = "NS",
                         logfc.threshold = 0.25, min.pct = 0.1,
                         test.use = "wilcox")
features <- c('AREG','PRF1','NKG7','GNLY','KLRD1',
              'GZMA','GZMB','GZMK',"GZMH",
              "IFI16",
              'IFNG','CCL4','CCL5',
              'CCL20',
              'CTLA4','TIGIT',"LAG3",
              "CXCL13","TNFRSF9","EOMES",
              "TOX",
              "ITGAE",
              "HAVCR2","CCL27",'TNF',
              "SELL","LEF1","CCR7",'TCF7','IL7R',
              "GATA3","ANXA1","ANXA2","LTB")
whitePurple = c("#F7F7F7", '#00bbb1')
jjDotPlot(object =sce_Tcell,
          gene = features,
          dot.col = whitePurple,
          id = 'Condition',
          split.by = 'Condition',
          ytree = FALSE,
          split.by.aesGroup = T)

genesets <- read.gmt("resource/c5.go.bp.v2023.2.Hs.symbols.gmt")
signatures <- split(genesets$gene,genesets$term)
features <- signatures[c("GOBP_T_CELL_CYTOKINE_PRODUCTION",
                         'GOBP_T_CELL_ACTIVATION')]
features <- c(features,list("T cell Exhaustion" =c('PRF1',
                                                   'GZMA','GZMB',"GZMH",
                                                   "IFI16",
                                                   'IFNG',
                                                   'CTLA4','TIGIT',"LAG3",
                                                   "CXCL13","TNFRSF9","CD160","CD244","EOMES",
                                                   "CD38","TOX",
                                                   "ITGAE","PDCD1","ENTPD1",
                                                   "HAVCR2")))
sce_Tcell <- AddModuleScore(sce_Tcell, 
                            features = features, 
                            nbin = 15,
                            seed = 1234,
                            name = "seurat", 
                            assay = "RNA")
test <- sce_Tcell@meta.data
names(sce_Tcell@meta.data)[grep("seurat\\d",names(sce_Tcell@meta.data))] <- names(features)
df<- FetchData(sce_Tcell,vars = c("Condition",
                                  "GOBP_T_CELL_CYTOKINE_PRODUCTION"))
df$Condition <- factor(df$Condition,levels = c("NS","VU"))
colnames(df) = c('Condition','Pathway')
my_comparisons <- list(c("VU", "NS")) 
ggplot(df,aes(x=Condition ,y=Pathway, fill=Condition))+
  geom_violin(color='black',size=0.3,color=NA,alpha=0.5)+
  theme_classic() + 
  labs(title = "GOBP_T_CELL_CYTOKINE_PRODUCTION")+
  theme(text = element_text(size=10, colour = "black")) + 
  theme(plot.title = element_text(hjust = 0.5, size = 10)) +
  theme(legend.position="none") +  
  geom_jitter(width = 0.35,
              height =0, size = 1.5, shape = 21, 
              stroke = 0.3,
              show.legend = FALSE)+
  geom_boxplot(width=0.2, outlier.shape = NA) +
  scale_fill_manual(values = c("VU" = "#F48892", "NS" = "#91CAE8"))+
  stat_compare_means(method="wilcox",hide.ns = F,
                     comparisons =c(my_comparisons),
                     label="p.signif",
                     bracket.size=0.8,
                     size=6)
