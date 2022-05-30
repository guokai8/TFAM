library(monocle)
library(Seurat)
load("cd4.rdata")
cd4c <- as(as.matrix(cd4@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = cd4@meta.data)
fData <- data.frame(gene_short_name = row.names(cd4c), row.names = row.names(cd4c))
fd <- new('AnnotatedDataFrame', data = fData)
###
cds <- newCellDataSet(cd4c,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
#
cds <- detectGenes(cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))
diff_test_res <- differentialGeneTest(cds[ expressed_genes, ], fullModelFormulaStr="~ celltype",cores=40)
#only use the genes with differential expressed
ordering_genes <- row.names(subset(diff_test_res, qval < 0.001))
length(ordering_genes) 
##reduce dimension
cds <- setOrderingFilter(cds, ordering_genes)  
cds <- reduceDimension(cds, max_components = 2, method='DDRTree')
#order cells
cds <- orderCells(cds)
cell_colors <- mycol
###generate trajectory figures Figure 1G
plot_cell_trajectory(cds, color_by="celltype",cell_size = 0.5) + scale_color_manual(values=cell_colors) +facet_wrap(~Group)
###### Figure S2C
cds <- orderCells(cds, root_state = 6)
plot_cell_trajectory(cds, color_by = "Pseudotime",show_branch_points = F)
########## State color Figure S2D
scolor<-c("#B2DF8A", "#e6194b", "#F4D000", "#E58308", "cyan4","darkgreen","darkred")
plot_cell_trajectory(cds, color_by="State") + scale_color_manual(values=scolor) 
#### State figures
cs<-data.frame(Celltype=cds$celltype,State=cds$State,Group=cds$Group)
### State WT and cKO Figure S2E
subset(cs,Group=="KO")%>%group_by(State,Celltype)%>%summarise(count=n())%>%mutate(pro=count/sum(count))%>%ggplot(aes(State,pro,fill=Celltype))+geom_bar(stat="identity")+scale_fill_manual(values=mycol)+theme_classic(base_size=14)
subset(cs,Group=="WT")%>%group_by(State,Celltype)%>%summarise(count=n())%>%mutate(pro=count/sum(count))%>%ggplot(aes(State,pro,fill=Celltype))+geom_bar(stat="identity")+scale_fill_manual(values=mycol)+theme_classic(base_size=14)
###State group Figure S2D
cs%>%group_by(Group,State)%>%summarise(count=n())%>%mutate(pro=count/sum(count))%>%ggplot(aes(Group,pro,fill=State))+geom_bar(stat="identity")+scale_fill_manual(values=scolor)+theme_classic(base_size = 14)+ylab("Cell proportion (%)")+xlab("")
dev.print(pdf,file="State_group.pdf")
##########
###heatmap and enrichment analysis
###Figure 1I
pesudo<- differentialGeneTest(cds[ ordering_genes, ], fullModelFormulaStr = "~sm.ns(Pseudotime)",cores = 40)
pesudo<-pesudo[order(pesudo$pval),]
sig_gene_names <- row.names(subset(pesudo, qval < 0.01))
ph<-plot_pseudotime_heatmap(cds[ sig_gene_names, ], num_clusters = 4, cores=40, show_rownames=TRUE,hmcols = viridis(128),return_heatmap = T)
dev.print(pdf,file="traj_heatmap.pdf")
gs<-ph$tree_row$labels[ph$tree_row$order]
###save genes in each cluster
write.csv(cbind(gs),file="cell_traj_cluster.csv")
### Do enrichment analysis
clus<-read.csv("cell_traj_cluster.csv",row.names = 1)
clus<-split(clus$gene,paste0("C",clus$X.1))
names(clus)<-paste0("C",names(clus))
library(richR)
mmko<-buildAnnot(species = "mouse",keytype = "SYMBOL",anntype = "KEGG",builtin = F)
clusg<-lapply(clus, function(x)richKEGG(x,mmko,builtin = F))
sapply(names(clusg), function(x)write.csv(clusg[[x]],file=paste0(x,"KEGG.csv")))
####State KEGG
cd4$State<-pData(cds)[rownames(cd4@meta.data),"State"]
cd4$State<-paste0("S",cd4$State)
Idents(cd4)<-"State"
ds<-FindAllMarkers(cd4,test.use = "MAST",logfc.threshold = 0)
### only use up-regulated DEGs
###Figure S2F
dss<-subset(ds,p_val_adj<0.05 & avg_log2FC>0.25)
dsg<-split(dss$gene,dss$cluster)
dsgk<-lapply(dsg, function(x)richKEGG(x,mmko,builtin = F))
res<-do.call(rbind,lapply(dsgk, function(x)result(x)))
res$Group<-sub('\\..*','',rownames(res))
d<-res%>%filter(Pvalue<0.01)%>%group_by(Group)%>%slice(1:20)
ggplot(d,aes(Group,Term,color=-log10(Pvalue)))+geom_point(size=3)+scale_color_gradient(low="pink",high="red")+theme_minimal(base_size = 12)+xlab("")+ylab("")
dev.print(pdf,file="State_KEGG_top20.pdf")
###write out 
sapply(names(dsgk), function(x)write.csv(result(dsgk[[x]]),file=paste0(x,"_KEGG.csv")))
###############################################################
###########################


