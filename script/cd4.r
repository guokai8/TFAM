setwd("ram/resn/")
library(Seurat)
load("cd4.rdata")
#####
DefaultAssay(cd4)<-"integrated"
cd4 <- RunPCA(cd4, verbose = T)
#####
cd4 <- RunUMAP(cd4, dims = 1:20)
cd4 <- FindNeighbors(cd4, dims = 1:20)
cd4 <- FindClusters(
  object = cd4, 
  reduction.type = "pca", 
  dims.use = 1:20, 
  resolution = seq(0.1,1.2,0.1), 
  print.output = FALSE, 
  save.SNN = TRUE
)
sapply(grep("^integrated_snn_res",colnames(cd4@meta.data),value = TRUE),
       function(x) length(unique(cd4@meta.data[,x])))
cd4<-RunTSNE(cd4,dims = 1:20)
###use 0.7
cd4 <- FindClusters(cd4, resolution = 0.7)
DotPlot(cd4,features = c("Sell","Lef1","Tcf7","Ccr7","S100a4","Lgals1","S100a10","Tnfsf8","Lag3","Tbc1d4","Ikzf2","Cd74","Il2ra","Icos","Rgs1","Ccl5","Gzma","Gzmk","Eomes"))+coord_flip()+scale_color_viridis_c()
#####
library(tidyverse)
library(future)
plan("multiprocess", workers = 40)
DefaultAssay(cd4) <- "RNA"
cd4markers<-FindAllMarkers(cd4)
cd4markers %>% group_by(cluster_id) %>% top_n(n = 2, wt = WT_avg_log2FC)
### assign the celltype
clusterid<-c("Naive","Naive","Naive","Cytotoxic","TEM","aTregs","Exhausted","nTregs","Naive","Naive","Exhausted")
names(clusterid) <- levels(cd4)
cd4 <- RenameIdents(cd4, clusterid)
cd4$celltype<-Idents(cd4)
mycol<-c("#6CB2E2","#88918A","#BB66A8","#F6BF39","#B55622","#006D4A")
names(mycol)<-levels(cd4$celltype)
###generate figures Figure 1C
DimPlot(cd4,label=T,cols=mycol,reduction = "tsne",split.by = "Group")
dev.print(pdf,file="cd4_tsne_group.pdf")
###generate circle figures Figure S2A
meta<-cd4@meta.data
meta%>%group_by(Group,celltype)%>%summarise(count=n())%>%mutate(pro=count/sum(count))%>%filter(Group=="WT")%>%
  ggplot(aes(x=2,pro,fill=celltype))+geom_bar(width = 1, size = 0.1, color = "white", stat = "identity")+scale_fill_manual(values = mycol)+
  geom_text_repel(aes(label = paste0(round(100*pro,1), "%")), position = position_stack(vjust = 0.5)) + 
  xlim(0.5,2.5)+ coord_polar("y")+ guides(fill=guide_legend(ncol=2)) +  theme_classic()+
  labs(x = NULL, y = NULL, fill = NULL,
       title = "WT") +theme(axis.line = element_blank(),
                            axis.text = element_blank(),
                            axis.ticks = element_blank(),
                            plot.title = element_text(hjust = 0.5, color = "#666666"))

dev.print(pdf,file="WT_CD4_circle.pdf")
meta%>%group_by(Group,celltype)%>%summarise(count=n())%>%mutate(pro=count/sum(count))%>%filter(Group=="KO")%>%
  ggplot(aes(x=2,pro,fill=celltype))+geom_bar(width = 1, size = 0.1, color = "white", stat = "identity")+scale_fill_manual(values = mycol)+
  geom_text_repel(aes(label = paste0(round(100*pro,1), "%")), position = position_stack(vjust = 0.5)) + 
  xlim(0.5,2.5)+ coord_polar("y")+ guides(fill=guide_legend(ncol=2)) +  theme_classic()+
  labs(x = NULL, y = NULL, fill = NULL,
       title = "KO") +theme(axis.line = element_blank(),
                            axis.text = element_blank(),
                            axis.ticks = element_blank(),
                            plot.title = element_text(hjust = 0.5, color = "#666666"))

dev.print(pdf,file="KO_CD4_circle.pdf")
###########################
##generate figures for proportion changes
###Figure 1E
meta%>%group_by(Group,celltype)%>%summarise(count=n())%>%spread(Group,count)%>%mutate(Change=(KO-WT)/WT)
meta%>%group_by(Group,celltype)%>%summarise(count=n())%>%mutate(pro=count/sum(count))%>%select(celltype,Group,pro)%>%spread(Group,pro)%>%
  mutate(Change=(KO-WT)/WT)%>%ggplot(aes(celltype,Change*100,fill=celltype))+
  geom_bar(stat="identity")+geom_text(aes(label=paste0(round(Change*100,1),"(%)")))+scale_fill_manual(values=mycol)+ylim(-80,300)+
    theme_classic(base_size = 15)+coord_flip()+ylab("Relative Change (%)")+xlab("")
dev.print(pdf,file="CD4_relative_change.pdf")
#######################
################################################
####calculate the DEGs
cd4$condition<-paste0(cd4$Group,cd4$celltype)
Idents(cd4)<-"condition"
deg<-lapply(as.character(unique(cd4$celltype)), function(x)FindMarkers(cd4,ident.1 = paste("KO",x,sep=""),ident.2 = paste("WT",x,sep=""),test.use = "MAST",logfc.threshold = 0))
names(deg)<-as.character(unique(cd4$celltype))
###write out
sapply(names(deg), function(x)write.csv(deg[[x]],file=paste(x,"_KOvsWT_DEG.csv",sep="")))
##########
##GSEA
filename<-list.files(pattern="KOvsWT")
deg<-lapply(filename, function(x)read.csv(x,row.names = 1))
library(richR)
mmko<-buildAnnot(species = "mouse",keytype = "SYMBOL",anntype = "KEGG",builtin = F)
gseak<-function(x){
  fc<-x$avg_log2FC
  names(fc)<-rownames(x)
  res<-richGSEA(fc,mmko,minSize = 5)
  return(res)
}
kegg<-lapply(deg, function(x)gseak(x))
sapply(names(kegg), function(x)write.csv(result(kegg[[x]]),file=paste(x,"KOVS_WT_GSEA_KEGG.csv",sep="_")))
dfk<-do.call(rbind,lapply(kegg,result))
dfk$Group<-sub('\\..*','',rownames(dfk))
dfk$Group<-factor(dfk$Group,levels = c("Naive","TEM","Exhausted","nTregs","aTregs","Cytotoxic"))
###make figures
##  Figure 1F
selp<-c("Antigen processing and presentation","Calcium signaling pathway","cAMP signaling pathway",
        "Chemokine signaling pathway","Cytokine-cytokine receptor interaction",
        "Hematopoietic cell lineage","IL-17 signaling pathway","JAK-STAT signaling pathway",
        "Natural killer cell mediated cytotoxicity","Oxidative phosphorylation",
        "PD-L1 expression and PD-1 checkpoint pathway in cancer",
        "Ribosome","T cell receptor signaling pathway",
        "TNF signaling pathway")
ggplot(subset(dfk,pval<0.01&pathway%in%selp),aes(Group,pathway,color=NES,size=-log10(padj)))+
  geom_point(alpha=0.75)+scale_color_gradient2(high='red',low='cyan4',mid='white',midpoint = 0)+theme_minimal(base_size = 14)+
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust = 0.5))+xlab("")+ylab("")
dev.print(pdf,file="Cd4_GSEA_KOvsWT.pdf")
### for all pathways 
ggplot(subset(dfk,pval<0.01),aes(Group,pathway,color=NES,size=-log10(padj)))+
  geom_point(alpha=0.75)+scale_color_gradient2(high='red',low='cyan4',mid='white',midpoint = 0)+theme_minimal(base_size = 14)+
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust = 0.5))+xlab("")+ylab("")

#####
####Figure S2B
gene=c("S100a8","S100a9","Cdk6","Rora","Cdkn1b","Gzmb")
VlnPlot(cd4,features = gene,split.by = "Group")
dev.print(pdf,file="cd4_select_gene_vln.pdf")

