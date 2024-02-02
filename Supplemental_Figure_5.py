# Import libraries
import numpy as np
import pandas as pd
import scanpy as sc

# Read in dataset from current work
df_AK = sc.read_h5ad("/panfs/accrepfs.vampire/home/kjara/10191-AK-analysis/adata_concat.h5ad")

# Read in Human Neural Organoid Atlas
df_HNOCA = sc.read_h5ad("/gpfs52/data/lippmann_lab/HNOCA/hnoca_pre-release_public_subset.h5ad")

# Prepare dataframes for ingest
var_names = df_HNOCA.var_names.intersection(df_AK.var_names)
df_HNOCA_sub = df_HNOCA[:, var_names]
df_AK_sub = df_AK[:, var_names]
sc.pp.neighbors(df_HNOCA_sub)
sc.pp.neighbors(df_AK_sub)

b = np.array(list(map(len, df_HNOCA_sub.obsp['distances'].tolil().rows))) # number of neighbors of each cell
df_HNOCA_sub = df_HNOCA_sub[np.where(b == 15-1)[0]] # select only those with the right number
sc.pp.neighbors(df_HNOCA_sub, 15) # rebuild the neighbor graph

# Transfer annotations
sc.tl.ingest(df_AK_sub, df_HNOCA_sub, obs='annot_level_4')
sc.tl.ingest(df_AK_sub, df_HNOCA_sub, obs='annot_level_2')
