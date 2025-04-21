rm(list = ls())
library(CellChat)
library(RColorBrewer)
library(ggplot2)
gc()

VU_cellchat<-readRDS("./CellChat/VUcellchat_2.rds")
VU_cellchat <- aggregateNet(VU_cellchat)
NS_cellchat<-readRDS("./CellChat/NScellchat_2.rds")
NS_cellchat <- aggregateNet(NS_cellchat)


#### macrophage to fibroblast ----
NS.net <- subsetCommunication(NS_cellchat)
unique(NS.net$source)
Macro_net <- subset(NS.net, source %in% c("Macro1","Macro2","Macro3"))
Macro_net <- subset(Macro_net, target %in% c("subFB1","subFB2","subFB3",
                                             "subFB4","subFB5"))
Macro_net <- subset(Macro_net, annotation %in% c("Secreted Signaling"))
Macro_net$cell_inter <- paste0(Macro_net$source,"->",Macro_net$target)
Macro_net$group <- "NS"


VU.net <- subsetCommunication(VU_cellchat)
Macro_net2 <- subset(VU.net, source %in% c("Macro1","Macro2","Macro3"))
Macro_net2 <- subset(Macro_net2, target %in% c("subFB1","subFB2","subFB3",
                                               "subFB4","subFB5"))
Macro_net2 <- subset(Macro_net2, annotation %in% c("Secreted Signaling"))
Macro_net2$cell_inter <- paste0(Macro_net2$source,"->",Macro_net2$target)
Macro_net2$group <- "VU"

# merge
df.net <- rbind(Macro_net, Macro_net2)
ggplot(df.net,aes(x=cell_inter,y=interaction_name)) +
  geom_point(aes(size=prob,color=prob)) +
  geom_point(shape=21,aes(size=prob))+
  scale_size(range=c(2, 6))+ 
  facet_wrap(~group)+
  scale_color_gradientn('Communication\nProbability', 
                        colors=colorRampPalette(rev(brewer.pal(9, "Spectral")))(100)) +
  theme_bw() +
  theme(axis.text=element_text(size=10, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=10),
        axis.text.y = element_text(size=8, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
