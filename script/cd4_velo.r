library(Seurat)
#########
library(SeuratDisk)
load("cd4.rdata")
### convert data
SaveH5Seurat(cd4, filename = "cd4.h5Seurat")
Convert("cd4.h5Seurat", dest = "h5ad")
### write out splice and unsplice data
write.csv(t(as.matrix(cd4@assays$spliced@counts)),file="cd4_spliced.csv")
write.csv(t(as.matrix(cd4@assays$unspliced@counts)),file="cd4_unspliced.csv")
write.csv(t(as.matrix(cd4@assays$RNA@counts)),file='cd4_count.csv')
######
####Figure 1H
### use python to run following code
### import scvelo as scv
### import scanpy as sp
### cd4 = sp.read_h5ad("cd4.h5ad")
### s=scv.read('cd4_spliced.csv')
### u=scv.read('cd4_unspliced.csv')
### adata=s
### adata.layers['spliced']=s.X
### adata.layers['unspliced']=u.X
### adata.obsm['X_umap']=mg.obsm['X_umap']
### adata.obs['seurat_clusters']=cd4.obs['celltype']
### adata.obsm['X_tsne']=cd4.obsm['X_tsne']
### adata.obsm['X_pca']=cd4.obsm['X_pca']
### adata.obs['Group']=cd4.obs['Group']
### adata.write('adata.h5ad',compression="gzip")
#### use R to load the h5ad files
library(reticulate)
use_python("/usr/bin/python3")
reticulate::py_config()
scv<-import("scvelo")
sp<-import('scanpy')
adata<-sp$read_h5ad('adata.h5ad')
wt=sp$read_h5ad('wt.h5ad')
ko=sp$read_h5ad('ko.h5ad')
###wt
scv$pp$filter_genes(wt) ## filter
scv$pp$moments(wt) ## normalize and compute moments
scv$tl$recover_dynamics(wt) ## model
##############################
scv$tl$velocity(wt, mode='dynamical')
scv$tl$velocity_graph(wt)
scv$pl$velocity_embedding_stream(wt, basis='tsne',save = "wt_velo_tsne.svg",color="celltype",palette = mycol)
scv$pl$velocity_embedding_grid(wt, basis='tsne', color='celltype', save='wt_embedding_grid_tsne.pdf', title='', scale=0.25,palette = mycol)
#####ko
scv$pp$filter_genes(ko) ## filter
scv$pp$moments(ko) ## normalize and compute moments
scv$tl$recover_dynamics(ko) ## model
##############################
scv$tl$velocity(ko, mode='dynamical')
scv$tl$velocity_graph(ko)
scv$pl$velocity_embedding_stream(ko, basis='tsne',save = "ko_velo_tsne.svg",color="celltype",palette = mycol)
scv$pl$velocity_embedding_grid(ko, basis='tsne', color='celltype', save='ko_embedding_grid_tsne.pdf', title='', scale=0.25,palette = mycol)
####all cells
scv$pp$filter_genes(adata) ## filter
scv$pp$moments(adata) ## normalize and compute moments
scv$tl$recover_dynamics(adata) ## model


