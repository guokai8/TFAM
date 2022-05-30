### load the packages
import os
import glob
import pickle
import pandas as pd
import numpy as np
from dask.diagnostics import ProgressBar
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2
from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell
### the path of data and database
db_fnames=glob.glob('/home/guokai8/cisTarget_databases/mm10*.feather')
dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]
MM_TFS_FNAME='mm_mgi_tfs.txt'
SC_EXP_FNAME='exprMat.csv'
ADJACENCIES_FNAME = os.path.join("./", "adjacencies.tsv")
MODULES_FNAME = os.path.join("./", "modules.p")
MOTIFS_FNAME = os.path.join("./", "motifs.csv")
REGULONS_FNAME = os.path.join("./", "regulons.p")
AUC_NAME=os.path.join("./","aucmtx.p")
#########
ex_matrix = pd.read_csv(SC_EXP_FNAME, sep=',', header=0, index_col=0).T
#### read the data
tf_names = load_tf_names(MM_TFS_FNAME)
### load data base
def name(fname):
    return os.path.splitext(os.path.basename(fname))[0]
dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]
####
MOTIF_ANNOTATIONS_FNAME="/home/guokai8/cisTarget_databases/motifs-v9-nr.mgi-m0.001-o0.0.tbl"
adjancencies = grnboost2(expression_data=ex_matrix, tf_names=tf_names, verbose=True,seed=777)
####
adjancencies.to_csv(ADJACENCIES_FNAME,index=False,sep='\t')
modules = list(modules_from_adjacencies(adjancencies, ex_matrix))
with open(MODULES_FNAME, 'wb') as f:
    pickle.dump(modules, f)
   
####
with ProgressBar():
    df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME)
# Create regulons from this table of enriched motifs.
regulons = df2regulons(df)
df.to_csv('df_to_regulons.csv')
with open(REGULONS_FNAME, "wb") as f:
    pickle.dump(regulons, f)   
####
auc_mtx = aucell(ex_matrix, regulons, num_workers=50,seed=777)
with open(AUC_NAME, "wb") as f:
    pickle.dump(auc_mtx, f)
##### binary
from pyscenic.binarization import binarize
binary_mtx, auc_thresholds = binarize( auc_mtx, num_workers=50 )
### We just need binary_mtx, auc_mtx and auc_thresholds and df_to_regulons.csv
####################################################
