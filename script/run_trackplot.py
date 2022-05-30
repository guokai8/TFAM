##Figure 1D
import scanpy as sc
cd4=sc.read_h5ad('cd4.h5ad')
col={'Naive':"#6CB2E2",'TEM':"#88918A",'Exhausted':"#BB66A8",'rTregs':"#F6BF39",'aTregs':"#B55622",'Cytotoxic':"#006D4A"}
marker_genes_dict = {'Naive': ['Sell', 'Lef1','Tcf7'], 'Cytotoxic': ['Ccl5', 'Gzmb','Gzma'],  'TEM': ['S100a4','S100a6'],'aTregs': ['Icos', 'Rgs1'],'Exhausted': ['Tnfsf8','Tbc1d4'], 'nTregs': ['Ikzf2','Il2ra','Cd74']}
sc.pl.tracksplot(cd4,marker_genes_dict,groupby="celltype",palette=col,save='marker')