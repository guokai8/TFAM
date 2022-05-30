# TFAM
Mitochondrial Transcription Factor A Regulates FOXP3+ T Cell Mediated Maintenance of CD4+ T Cell Landscapes and Immunological Aging
### Description
Data were collected after mapping all fastq files to the mm10 genome with the cellranger v6.0 software with default paramaters by using "cellranger count". 
### Required packages:
Please check the sessionInfo.txt

#### Read data, integration, clustering
run the __run.r__ for clustering, cell annotation, Differential expressed analysis and GSEA for all cell types. The metadata.csv and Embedding.csv in the data folder can be used to skip the integration, clustering analysis. The metadata.csv includes all cell information and the Embedding.csv includes UMAP and TSNE information.

#### CD4 
run the __cd4.r__ for the cd4 analysis. The cd4metadata.csv and cd4Embedding.csv in the data folder can be used to skip the integration, clustering analysis. The cd4metadata.csv includes all cell information and the cd4Embedding.csv includes UMAP and TSNE information.

#### Cell trajectory, RNA velocity, SCENIC
run the __cd4_traj.r__, __cd4_scenic.r__ and __cd4_velo.r__ for all other analysis. The cd4_traj.rdata and cd4_scenic.rdata in the data folder are the results from Monocle and pySCEINIC packages.

### 
