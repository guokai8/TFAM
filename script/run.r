setwd("ram/resn/")
library(Seurat)
library(cowplot)
options(future.globals.maxSize = 500000 * 1024^2)
WT1.data<-Read10X("/home/guokai8/ram/res/75_1/outs/filtered_feature_bc_matrix/")
WT2.data<-Read10X("/home/guokai8/ram/res/75_2/outs/filtered_feature_bc_matrix/")
KO1.data<-Read10X("/home/guokai8/ram/res/75_3/outs/filtered_feature_bc_matrix/")
KO2.data<-Read10X("/home/guokai8/ram/res/75_4/outs/filtered_feature_bc_matrix/")
############
wt1<-CreateSeuratObject(counts = WT1.data, project = "WT1",min.cells = 3, min.features = 200)
wt2<-CreateSeuratObject(counts = WT2.data, project = "WT2",min.cells = 3, min.features = 200)
ko1<-CreateSeuratObject(counts = KO1.data, project = "KO1",min.cells = 3, min.features = 200)
ko2<-CreateSeuratObject(counts = KO2.data, project = "KO2",min.cells = 3, min.features = 200)

#################
#assign group
#######
wt1$group<-"WT1"
wt2$group<-"WT2"
ko1$group<-"KO1"
ko2$group<-"KO2"

#############check mt
wt1[["percent.mt"]]<-PercentageFeatureSet(wt1,pattern = "^mt-")
pp1<-VlnPlot(wt1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.1)
wt2[["percent.mt"]]<-PercentageFeatureSet(wt2,pattern = "^mt-")
pp2<-VlnPlot(wt2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.1)
ko1[["percent.mt"]]<-PercentageFeatureSet(ko1,pattern = "^mt-")
pp3<-VlnPlot(ko1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.1)
ko2[["percent.mt"]]<-PercentageFeatureSet(ko2,pattern = "^mt-")
pp4<-VlnPlot(ko2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.1)
##

########################################################################
###subset
wt1 <- subset(wt1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 25)
wt2 <- subset(wt2, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 25)
ko1 <- subset(ko1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 25)
ko2 <- subset(ko2, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 25)
##########
##Remove doublet with Scrublet
library(reticulate)
use_python("/usr/bin/python3")
reticulate::py_config()
scr<-import("scrublet")
wt1d<-scr$Scrublet(counts_matrix = t(as.matrix(wt1@assays$RNA@counts)),expected_doublet_rate=0.1)
wt1d<-wt1d$scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=50L)
###
wt1s<-as.vector(wt1d[[1]])
wt1p<-as.vector(wt1d[[2]])
wt2d<-scr$Scrublet(counts_matrix = t(as.matrix(wt2@assays$RNA@counts)),expected_doublet_rate=0.1)
wt2d<-wt2d$scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=50L)
### 
wt2s<-as.vector(wt2d[[1]])
wt2p<-as.vector(wt2d[[2]])
ko1d<-scr$Scrublet(counts_matrix = t(as.matrix(ko1@assays$RNA@counts)),expected_doublet_rate=0.1)
ko1d<-ko1d$scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=50L)
###
ko1s<-as.vector(ko1d[[1]])
ko1p<-as.vector(ko1d[[2]])
ko2d<-scr$Scrublet(counts_matrix = t(as.matrix(ko2@assays$RNA@counts)),expected_doublet_rate=0.1)
ko2d<-ko2d$scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=50L)
### 
ko2s<-as.vector(ko2d[[1]])
ko2p<-as.vector(ko2d[[2]])
#### assign back
wt1$doublet_score<-wt1s
wt1$doublet<-ifelse(wt1p,"Doublet","Singlet")
wt2$doublet_score<-wt2s
wt2$doublet<-ifelse(wt2p,"Doublet","Singlet")
ko1$doublet_score<-ko1s
ko1$doublet<-ifelse(ko1p,"Doublet","Singlet")
ko2$doublet_score<-ko2s
ko2$doublet<-ifelse(ko2p,"Doublet","Singlet")
####only keep Singlet
####
#https://github.com/satijalab/seurat/issues/1679
# NormalizeData prior to CellCycleScoring. 
# This learns cell cycle scores that can be added to the vars.to.regress parameter in 
# SCTransform. For all downstream analyses, you can use the SCT assay.
sam<-list(WT1=subset(wt1,doublet=="Singlet"),WT2=subset(wt2,doublet=="Singlet"),KO1=subset(ko1,doublet=="Singlet"),KO2=subset(ko2,doublet=="Singlet"))
### cell cycle scoring
cell_cycle_genes<-read.csv("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Mus_musculus.csv")
library(GenomicFeatures)
library(AnnotationHub)
library(tidyverse)
ah <- AnnotationHub()
# Access the Ensembl database for organism
ahDb <- query(ah, 
              pattern = c("Mus musculus", "EnsDb"), 
              ignore.case = TRUE)

# Acquire the latest annotation files
id <- ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)

# Download the appropriate Ensembldb database
edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb, 
                     return.type = "data.frame")

# Select annotations of interest
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)

cell_cycle_markers <- dplyr::left_join(cell_cycle_genes, annotations, by = c("geneID" = "gene_id"))

# Acquire the S phase genes
s_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "S") %>%
  pull("gene_name")

# Acquire the G2M phase genes        
g2m_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "G2/M") %>%
  pull("gene_name")
  ###Normalized data
for (i in 1:length(sam)) {
  sam[[i]] <- NormalizeData(sam[[i]], verbose = T)
}

#######
for (i in 1:length(sam)) {
  sam[[i]] <- CellCycleScoring(sam[[i]], s.features = s_genes, g2m.features = g2m_genes)
}

########SCTransformed
for (i in 1:length(sam)) {
  sam[[i]] <- SCTransform(sam[[i]])
}
##########
sam.features <- SelectIntegrationFeatures(object.list = sam, nfeatures = 3000)
sam <- PrepSCTIntegration(object.list = sam, anchor.features = sam.features, 
                                   verbose = FALSE)
###
sam.anchors <- FindIntegrationAnchors(object.list = sam, normalization.method = "SCT",
                                           anchor.features = sam.features, verbose = T)
sam.integrated <- IntegrateData(anchorset = sam.anchors,  normalization.method = "SCT",
                                     verbose = T)
###
sam.integrated <- RunPCA(sam.integrated, verbose = T)
####
####PC numbers
ElbowPlot(sam.integrated,ndims = 30)
dev.print(pdf,file="sample_all_elbow.pdf")
###
###############################################
sam.integrated <- RunUMAP(sam.integrated, dims = 1:30)
sam.integrated <- FindNeighbors(sam.integrated, dims = 1:30)
sam.integrated <- FindClusters(
  object = sam.integrated, 
  reduction.type = "pca", 
  dims.use = 1:30, 
  resolution = seq(0.1,1.2,0.1), 
  print.output = FALSE, 
  save.SNN = TRUE
)
sapply(grep("^integrated_snn_res",colnames(sam.integrated@meta.data),value = TRUE),
       function(x) length(unique(sam.integrated@meta.data[,x])))

sam.integrated <- FindClusters(sam.integrated, resolution = 0.7)
sam.integrated$Group<-sub('\\d+','',sam.integrated$group)
#########################################################
###
DimPlot(sam.integrated,split.by = "Group",cols=c(distcolor[1:14],lightcolor[1:12]),label = T,repel = T)+NoLegend()
dev.print(pdf,file="sample_umap_new.pdf")
##########################################################
DimPlot(sam.integrated, reduction = "umap", label=T,cols=c(distcolor[1:14],lightcolor[1:12]),repel = T)+NoLegend()
dev.print(pdf,file="sample_umap_label_all_new.pdf")
##############################################################
#############################################################
library(future)
plan("multiprocess", workers = 40)
DefaultAssay(sam.integrated) <- "RNA"
get_conserved <- function(cluster){
  FindConservedMarkers(sam.integrated,
                       ident.1 = cluster,
                       grouping.var = "Group") %>%
    rownames_to_column(var = "gene") %>%
    cbind(cluster_id = cluster, .)
}
conserved_markers <- map_dfr(0:26, get_conserved)
conserved_markers %>% group_by(cluster_id) %>% top_n(n = 2, wt = WT_avg_log2FC)
#sample.combined<-sample.combined.bk
sam.integrated.bk<-sam.integrated
DefaultAssay(sam.integrated) <- "RNA"
###################################################################
clusterid<-c("BC1","CD8+ Naive","CD4+ Naive","Erythroid","CD8+","CD4+","NK",
             "PC1","BC2","Erythroid","Plasma","BC3",
             "BC4","Neu","Plasma","PC2","Mac","Mono","Erythroid","Plasma",
             "DC","TAC","Neu","iDC","PC3","pDC")
names(clusterid) <- levels(sam.integrated)
sam.integrated <- RenameIdents(sam.integrated, clusterid)
sam.integrated$celltype<-Idents(sam.integrated)
mycol<-c("#A6761D" ,"#D95F02","deepskyblue" ,"#E7298A" ,"#7570B3","#E6AB02", "#A6CEE3" ,"#1F78B4","#B2DF8A","#33A02C","slateblue1","darkgreen" ,"darkred","plum1","darkmagenta" ,
         "hotpink2","magenta2", "#57C3F3","#E59CC4" ,"#23452F" ,"grey60")
names(mycol)<-c("BC1","CD4+ Naive", "BC2","CD8+ Naive","Erythroid","CD8+","PC1","CD4+","NK","BC3","Neu","TAC","BC4","Plasma","Mac","Mono","DC","PC2","iDC","PC3","pDC")
#######
gene1<-c("Hbb-bt","Hbb-bs")
samf<-subset(sam.integrated,celltype!="Erythroid")
fcol<-mycol[levels(unique(samf$celltype))]
DimPlot(samf,cols=fcol,split.by = "Group",reduction = "tsne",label=T)
dev.print(pdf,file="tsne_label_group.pdf")
#######################
genes<-c("Cd79a","Ms4a1","Cd19","Ebf1","Il7r",
         "Cd3d","Cd3g","Cd3e",
         "Tcf7","Cd79b","Cd52","Cd74","Trbc2",
         "Cd8a","Cd8b1","Sell","Mki67",
         "H2afz","Hist1h2ae","S100a4","Nkg7","Ccl5","Klrd1",
         "Cd38","Cd22","Pax5","S100a8","S100a9","Ran","Nme1","Eif5a",
         "Jchain","Igha","C1qa","Mafb","Mrc1",
         "Ly6c2","Ccr2","Lyz2","Ighg1","H2-Aa","Cst3","H2-Ab1",
         "H2-Eb1","Bst2","Irf8","Tcf4")
DotPlot(sam.integrated,features = genes,cluster.idents = T,cols=c("lightgrey","red"))+coord_flip()+theme(axis.text.x=element_text(angle=90,hjust = 1,vjust = 0.5))
dev.print(pdf,file="marker.pdf")
########################
meta<-samf@meta.data
meta%>%group_by(Group,celltype)%>%summarise(count=n())%>%mutate(pro=count/sum(count))%>%
  ggplot(aes(Group,pro,fill=celltype))+geom_bar(stat="identity")+scale_fill_manual(values=fcol)+
  theme_light(base_size = 15)+coord_flip()+theme(legend.position = "top")+labs(fill="")
dev.print(pdf,file="cell_prop.pdf")
meta%>%group_by(Group,celltype)%>%summarise(count=n())%>%mutate(pro=count/sum(count))%>%dplyr::filter(Group=="WT")%>%
ggplot(aes(x=2,pro,fill=celltype))+geom_bar(width = 1, size = 0.1, color = "white", stat = "identity")+scale_fill_manual(values = fcol)+
  geom_text_repel(aes(label = paste0(round(100*pro,1), "%")), position = position_stack(vjust = 0.5)) + 
  xlim(0.5,2.5)+ coord_polar("y")+ guides(fill=guide_legend(ncol=2)) +  theme_classic()+
  labs(x = NULL, y = NULL, fill = NULL,
       title = "WT") +theme(axis.line = element_blank(),
                           axis.text = element_blank(),
                           axis.ticks = element_blank(),
                           plot.title = element_text(hjust = 0.5, color = "#666666"))

dev.print(pdf,file="cell_prop_WT.pdf")
##################################
meta%>%group_by(Group,celltype)%>%summarise(count=n())%>%mutate(pro=count/sum(count))%>%dplyr::filter(Group=="KO")%>%
  ggplot(aes(x=2,pro,fill=celltype))+geom_bar(width = 1, size = 0.1, color = "white", stat = "identity")+scale_fill_manual(values = fcol)+
  geom_text_repel(aes(label = paste0(round(100*pro,1), "%"),max.overlaps=100), position = position_stack(vjust = 0.5)) + 
  xlim(0.5,2.5)+ coord_polar("y")+ guides(fill=guide_legend(ncol=2)) +  theme_classic()+
  labs(x = NULL, y = NULL, fill = NULL,
       title = "KO") +theme(axis.line = element_blank(),
                            axis.text = element_blank(),
                            axis.ticks = element_blank(),
                            plot.title = element_text(hjust = 0.5, color = "#666666"))

dev.print(pdf,file="cell_prop_KO.pdf")
#######################
meta%>%group_by(Group,celltype)%>%summarise(count=n())%>%
  ggplot(aes(Group,count,fill=celltype))+geom_bar(stat="identity")+scale_fill_manual(values=fcol)+
  theme_light(base_size = 15)+coord_flip()+theme(legend.position = "top")+labs(fill="")+xlab("")+ylab("")
dev.print(pdf,file="cell_number.pdf")
#######################################
##
library(SeuratWrappers)
wt1loom<-ReadVelocity("/home/guokai8/ram/res/75_1/velocyto/75_1.loom")
wt2loom<-ReadVelocity("/home/guokai8/ram/res/75_2/velocyto/75_2.loom")
ko1loom<-ReadVelocity("/home/guokai8/ram/res/75_3/velocyto/75_3.loom")
ko2loom<-ReadVelocity("/home/guokai8/ram/res/75_4/velocyto/75_4.loom")
#####
saml<-list()
for(i in names(wt1loom)){
  colnames(wt1loom[[i]])<-sub('x','-1_1',sub('.*:','',colnames(wt1loom[[i]])))
  colnames(wt2loom[[i]])<-sub('x','-1_2',sub('.*:','',colnames(wt2loom[[i]])))
  colnames(ko1loom[[i]])<-sub('x','-1_3',sub('.*:','',colnames(ko1loom[[i]])))
  colnames(ko2loom[[i]])<-sub('x','-1_4',sub('.*:','',colnames(ko2loom[[i]])))
  saml[[i]]<-cbind(wt1loom[[i]],wt2loom[[i]],ko1loom[[i]],ko2loom[[i]])
}

##############################
for (i in names(x = saml)) {
  ### Store assay in a new variable
  assay <- saml[[i]]
  ### Subset to filtered cells in Seurat object
  assay <- assay[,colnames(samf)]
  ### Add assay to Seurat object
  samf[[i]] <- CreateAssayObject(counts = assay)
}
#############################################
geness<-c("Ms4a1","Cd19","Ebf1","Cd79a","Il7r","Cd3d","Cd3g","Tcf7","Cd74","Trbc2","Cd8a","Cd8b1","Mki67",
         "H2afz","Hist1h2ae","Nkg7","Ccl5","Klrd1","S100a8","S100a9","Ran","Nme1","Eif5a",
         "Jchain","Igha","Xbp1","C1qa","Mafb","Mrc1","Ly6c2","Ccr2","Lyz2","H2-Aa",
         "H2-Ab1","H2-Eb1","Cst3","Bst2","Irf8","Tcf4")

samf$celltype<-factor(as.vector(samf$celltype),levels=c("BC1","BC2","BC3","BC4","CD4+","CD4+ Naive","CD8+","CD8+ Naive",
                                                        "PC1","PC2","PC3","NK","Neu","TAC","Plasma","Mac","Mono","DC","iDC","pDC"))
Idents(samf)<-"celltype"
DotPlot(samf,features = geness,cluster.idents = F)+theme(axis.text.x=element_text(angle=90,hjust = 1,vjust = 0.5))+scale_color_viridis_c()
dev.print(pdf,file="marker_label.pdf")
###
cd4<-subset(samf,celltype%in%c("CD4+","CD4+ Naive"))
#####
###DEGs and GSEA
samf$condition<-paste(samf$Group,samf$celltype,sep="_")
Idents(samf)<-"condition"
library(future)
plan("multiprocess", workers = 40)
deg<-lapply(as.character(unique(samf$celltype)), function(x)FindMarkers(samf,ident.1 = paste("KO",x,sep="_"),ident.2 = paste("WT",x,sep="_"),test.use = "MAST",logfc.threshold = 0))
names(deg)<-as.character(unique(samf$celltype))
library(richR)
mmko<-buildAnnot(species = "mouse",keytype = "SYMBOL",anntype = "KEGG",builtin = F)
gseak<-function(x){
  fc<-x$avg_log2FC
  names(fc)<-x$gene
  res<-richGSEA(fc,mmko,minSize = 5)
  return(res)
}
kegg<-lapply(deg, function(x)gseak(x))
sapply(names(kegg), function(x)write.csv(result(kegg[[x]]),file=paste(x,"GSEA_KEGG.csv",sep="_")))
selp<-c("Antigen processing and presentation","B cell receptor signaling pathway",
        "Chemokine signaling pathway","Cytokine-cytokine receptor interaction",
        "Hematopoietic cell lineage","IL-17 signaling pathway",
        "Natural killer cell mediated cytotoxicity","Oxidative phosphorylation",
        "Parkinson disease","PD-L1 expression and PD-1 checkpoint pathway in cancer",
        "Ribosome","TNF signaling pathway")
#######
res<-do.call(rbind,lapply(kegg, function(x)result(x)))
res$Group<-sub('\\..*','',rownames(res))
res<-subset(res,padj<0.05)

res$Group<-factor(res$Group,levels=c("B","NB1","NB2","Follicular B","CD4+ Memory","CD4+ Naive","CD8+","CD8+ Naive",
                                     "PC1","PC2","PC3","NK","Neu","TAC","Plasma","Mac","Mono","DC","iDC","pDC"))

ggplot(subset(res,pathway%in%c(selp)),aes(Group,pathway,color=NES,size=-log10(padj)))+geom_point()+
  scale_color_gradient2(low="cyan4",high="red",mid="white",midpoint = 0)+theme_minimal(base_size = 14)+
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust = 0.5))+xlab("")

