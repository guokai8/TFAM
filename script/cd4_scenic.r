library(Seurat)
library(SCENIC)
library(RcisTarget)
load("cd4.rdata")
exprMat <- cd4@assays$RNA@counts
exprMat <- as.matrix(exprMat)
##meta information
cellInfo <- data.frame(cd4@meta.data)
colnames(cellInfo)[which(colnames(cellInfo)=="celltype")] <- "celltype"
colnames(cellInfo)[which(colnames(cellInfo)=="Group")] <- "Group"
cellInfo <- cellInfo[,c("celltype","Group")]
colnames(cellInfo)[1]<-"CellType"
#####
dbs <- c('mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather', 'mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather')
dbDir='/home/guokai8/cisTarget_databases/' # RcisTarget databases location
scenicOptions <- initializeScenic(org="mgi", dbDir=dbDir,dbs=dbs, nCores=40)
genesKept <- geneFiltering(exprMat, scenicOptions, minCountsPerGene = 3 * 0.01 * ncol(exprMat), minSamples = ncol(exprMat) * 0.01)
exprMat_filtered <- exprMat[genesKept,]
####
write.csv(exprMat_filtered,file="cd4_expr.csv")
###run run_scenic.py
library(reticulate)
### run run_scenic.py in RStudio
####
binary<-py$binary_mtx
aucmtx<-py$auc_mtx
thres<-py$auc_thresholds
###################
#############generate figures Figure S4
tsne<-Embeddings(cd4,reduction = "tsne")
tsne<-cbind(tsne,meta[rownames(tsne),c("Group","celltype")])
tsne<-cbind(tsne,aucmtx[rownames(tsne),])
#
ggplot(tsne,aes(tSNE_1,tSNE_2))+
  geom_point(aes(color=`Nfkb1(+)`),size=0.5,alpha=0.75)+
  scale_color_gradient2(high="red",low = "white",mid="lightgrey",midpoint = thres["Nfkb1(+)"])+facet_wrap(~Group)+
  theme_classic()
dev.print(pdf,file="Nfkb1_aub_tsne.pdf")
#
ggplot(tsne,aes(tSNE_1,tSNE_2))+
  geom_point(aes(color=`Irf8(+)`),size=0.5,alpha=0.75)+
  scale_color_gradient2(high="red",low = "white",mid="lightgrey",midpoint = thres["Irf8(+)"])+facet_wrap(~Group)+
  theme_classic()
dev.print(pdf,file="Irf8_aub_tsne.pdf")
#
ggplot(tsne,aes(tSNE_1,tSNE_2))+
  geom_point(aes(color=`Runx2(+)`),size=0.5,alpha=0.75)+
  scale_color_gradient2(high="red",low = "white",mid="lightgrey",midpoint = thres["Runx2(+)"])+facet_wrap(~Group)+
  theme_classic()
dev.print(pdf,file="Runx2_aub_tsne.pdf")
#
ggplot(tsne,aes(tSNE_1,tSNE_2))+
  geom_point(aes(color=`Bcl11a(+)`),size=0.5,alpha=0.75)+
  scale_color_gradient2(high="red",low = "white",mid="lightgrey",midpoint = thres["Bcl11a(+)"])+facet_wrap(~Group)+
  theme_classic()
dev.print(pdf,file="Bcl11a_aub_tsne.pdf")
#########
ann_colors<-list(celltype=mycol,Group = c(WT = "darkgreen", KO = "#D95F02"))
anncol<-tsne[,c("Group","celltype")]
#### read regulon dataframe
regulons<-read.csv("df_to_regulons.csv",skip=2)
reglist<-split(regulons,regulons$TF)
reg<-lapply(reglist, function(x)unique(gsub("\\'",'',unlist(lapply(lapply(lapply(lapply(strsplit(x$X.6,"\\)"),function(x)sub('\\, \\(','',x)),function(x)sub('\\[\\(','',x)),function(x)sub('\\]','',x)),function(x)sub('\\,.*','',x))))))
regn<-cbind(unlist(lapply(reg, function(x)length(x))))
regn<-as.data.frame(regn)
regn$name<-paste0(rownames(regn),"(",regn$V1,"g)")
colnames(binary)<-regn$name
######### heatmap Figure 3A
anncol<-anncol[order(anncol$Group,anncol$celltype),]
pheatmap::pheatmap(t(binary[rownames(anncol),]),  fontsize_row=3,cluster_cols = F,
                   color = colorRampPalette(c("white","black"))(100), breaks=seq(0, 1, length.out = 100),
                   border_color=NA,annotation_col = anncol,annotation_colors=ann_colors,show_colnames = F)
dev.print(pdf,file="cd4_binary_heatmap_celltype.pdf")
############ Figure 3B
DotPlot(subset(cd4,Group=="KO"),features = c("Prdm1","Nfkb1","Nfkb2","Ahr"))+scale_color_viridis_c()
dev.print(pdf,file="regulon_dot_KO.pdf")
DotPlot(subset(cd4,Group=="WT"),features = c("Prdm1","Nfkb1","Nfkb2","Ahr"))+scale_color_viridis_c()
dev.print(pdf,file="regulon_dot_WT.pdf")




