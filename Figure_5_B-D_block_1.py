# Import packages
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scanpy.external as sce
from pandas import *

# Set working directory file
results_file = '2023.11.08_adata_concat.h5ad' 

# Import datasets
adata1 = sc.read_10x_mtx(
    "/data/lippmann_lab/10490-AK/10490-AK-0001_S1_L005_results/filtered_matrix/sensitivity_5",  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                              # write a cache file for faster subsequent reading
adata1.obs["condition"] = 'D90 GelMA_HC'
adata1.obs["day"] = 'D90'
adata1.obs["crosslinker"] = 'HC'
adata1.obs["matrix"] = 'GelMA'

adata2 = sc.read_10x_mtx(
    "/data/lippmann_lab/10490-AK/10490-AK-0002_S1_L005_results/filtered_matrix/sensitivity_5",  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                              # write a cache file for faster subsequent reading
adata2.obs["condition"] = 'D90 GelMA_LC'
adata2.obs["day"] = 'D90'
adata2.obs["crosslinker"] = 'LC'
adata2.obs["matrix"] = 'GelMA'

adata3 = sc.read_10x_mtx(
    "/data/lippmann_lab/10191-AK/10191-AK-0003_S1_L005_results/filtered_matrix/sensitivity_5",  # the directory with the `.mtx` file
    var_names='gene_symbols',  # use gene symbols for the variable names (variables-axis index)
    cache=True)  # write a cache file for faster subsequent reading
adata3.obs["condition"] = 'D90 GelMA-Cad_LC'  # Replace 'adata2' with 'adata3'
adata3.obs["day"] = 'D90'  # Replace 'adata2' with 'adata3'
adata3.obs["crosslinker"] = 'LC'  # Replace 'adata2' with 'adata3'
adata3.obs["matrix"] = 'GelMA-Cad'  # Replace 'adata2' with 'adata3'

adata4 = sc.read_10x_mtx(
    "/data/lippmann_lab/10191-AK/10191-AK-0004_S1_L005_results/filtered_matrix/sensitivity_5",  # the directory with the `.mtx` file
    var_names='gene_symbols',  # use gene symbols for the variable names (variables-axis index)
    cache=True)  # write a cache file for faster subsequent reading
adata4.obs["condition"] = 'D90 GelMA-Cad_HC'  # Replace 'adata2' with 'adata4'
adata4.obs["day"] = 'D90'  # Replace 'adata2' with 'adata4'
adata4.obs["crosslinker"] = 'HC'  # Replace 'adata2' with 'adata4'
adata4.obs["matrix"] = 'GelMA-Cad'  # Replace 'adata2' with 'adata4'
adata4.write(results_file)  # Replace 'adata2' with 'adata4'

adata5 = sc.read_10x_mtx(
    "/data/lippmann_lab/10191-AK/10191-AK-0005_S1_L005_results/filtered_matrix/sensitivity_5",  # the directory with the `.mtx` file
    var_names='gene_symbols',  # use gene symbols for the variable names (variables-axis index)
    cache=True)  # write a cache file for faster subsequent reading
adata5.obs["condition"] = 'D90 Matrigel'  # Replace 'adata2' with 'adata5'
adata5.obs["day"] = 'D90'  # Replace 'adata2' with 'adata5'
adata5.obs["crosslinker"] = 'NA'  # Replace 'adata2' with 'adata5'
adata5.obs["matrix"] = 'Matrigel'  # Replace 'adata2' with 'adata5'

adata6 = sc.read_10x_mtx(
    "/data/lippmann_lab/10579-AK/10579-AK-0001_S1_L005_results/filtered_matrix/sensitivity_5",  # the directory with the `.mtx` file
    var_names='gene_symbols',  # use gene symbols for the variable names (variables-axis index)
    cache=True)  # write a cache file for faster subsequent reading
adata6.obs["condition"] = 'D120 GelMA_HC'  # Replace 'adata2' with 'adata6'
adata6.obs["day"] = 'D120'  # Replace 'adata2' with 'adata6'
adata6.obs["crosslinker"] = 'HC'  # Replace 'adata2' with 'adata6'
adata6.obs["matrix"] = 'GelMA'  # Replace 'adata2' with 'adata6'

adata7 = sc.read_10x_mtx(
    "/data/lippmann_lab/10579-AK/10579-AK-0002_S1_L005_results/filtered_matrix/sensitivity_5",  # the directory with the `.mtx` file
    var_names='gene_symbols',  # use gene symbols for the variable names (variables-axis index)
    cache=True)  # write a cache file for faster subsequent reading
adata7.obs["condition"] = 'D120 GelMA_LC'  # Replace 'adata2' with 'adata7'
adata7.obs["day"] = 'D120'  # Replace 'adata2' with 'adata7'
adata7.obs["crosslinker"] = 'LC'  # Replace 'adata2'
adata7.obs["matrix"] = "GelMA"

adata8 = sc.read_10x_mtx(
    "/data/lippmann_lab/10579-AK/10579-AK-0003_S1_L005_results/filtered_matrix/sensitivity_5",  # the directory with the `.mtx` file
    var_names='gene_symbols',  # use gene symbols for the variable names (variables-axis index)
    cache=True)  # write a cache file for faster subsequent reading
adata8.obs["condition"] = 'D120 Matrigel'  # Replace 'adata2' with 'adata8'
adata8.obs["day"] = 'D120'  # Replace 'adata2' with 'adata8'
adata8.obs["crosslinker"] = 'NA'  # Replace 'adata2' with 'adata8'
adata8.obs["matrix"] = 'Matrigel'  # Replace 'adata2' with 'adata8'

adata9 = sc.read_10x_mtx(
    "/data/lippmann_lab/10191-AK/10191-AK-0002_S1_L005_results/filtered_matrix/sensitivity_5",  # the directory with the `.mtx` file
    var_names='gene_symbols',  # use gene symbols for the variable names (variables-axis index)
    cache=True)  # write a cache file for faster subsequent reading
adata9.obs["condition"] = 'D120 GelMA-Cad_HC'  # Replace 'adata2' with 'adata9'
adata9.obs["day"] = 'D120'  # Replace 'adata2' with 'adata9'
adata9.obs["crosslinker"] = 'HC'  # Replace 'adata2' with 'adata9'
adata9.obs["matrix"] = 'GelMA-Cad'  # Replace 'adata2' with 'adata9'

adata10 = sc.read_10x_mtx(
    "/data/lippmann_lab/10191-AK/10191-AK-0001_S1_L005_results/filtered_matrix/sensitivity_5",  # the directory with the `.mtx` file
    var_names='gene_symbols',  # use gene symbols for the variable names (variables-axis index)
    cache=True)  # write a cache file for faster subsequent reading
adata10.obs["condition"] = 'D120 GelMA-Cad_LC'  # Replace 'adata2' with 'adata10'
adata10.obs["day"] = 'D120'  # Replace 'adata2' with 'adata10'
adata10.obs["crosslinker"] = 'LC'  # Replace 'adata2' with 'adata10'
adata10.obs["matrix"] = 'GelMA-Cad'  # Replace 'adata2' with 'adata10'

adata_concat = ad.concat([adata1, adata2, adata3, adata4, adata5, adata6, adata7, adata8, adata9, adata10])
adata_concat.obs_names_make_unique()

# Filter cells and quality control measures
sc.pp.filter_cells(adata_concat, min_genes=2000)
sc.pp.filter_genes(adata_concat, min_cells=3)
adata_concat.var['mt'] = adata_concat.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata_concat, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata_concat = adata_concat[adata_concat.obs.n_genes_by_counts < 8000, :]
adata_concat = adata_concat[adata_concat.obs.pct_counts_mt < 5, :]
sc.pp.normalize_total(adata_concat, target_sum=1e4)
sc.pp.log1p(adata_concat)
sc.pp.highly_variable_genes(adata_concat, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata_concat.raw = adata_concat
adata_concat = adata_concat[:, adata_concat.var.highly_variable]

sc.pp.regress_out(adata_concat, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata_concat, max_value=10)
sc.tl.pca(adata_concat, svd_solver='arpack')
sce.pp.harmony_integrate(adata_concat, 'condition')

sc.pp.neighbors(adata_concat, n_neighbors=15, n_pcs=50, use_rep = 'X_pca_harmony')
sc.tl.tsne(adata_concat, use_rep = 'X_pca_harmony')

sc.tl.leiden(adata_concat, resolution = 1.7)

# Annotate cells
adata_concat.obs["cell_type"] = adata_concat.obs["leiden"]
cell_names = ["MC-1","ChP-1","RG-1","GB","ChP-2","N-1","RG-2","DLExN","ULExN","RG-3","N-2",
             "RG-4","RG-5","MC-2","RG-6","RG-7","RG-8","IPC-1","DLX5InH","RG-9","N-3",
             "IPC-2","RELNInH","IPC-3","MC-3","RG-10","MC-4","N-4","RG-11","N-5","MC-5",
             "Unknown-1","ChP-3","ChP-4","IPC-4","MC-6","MC-7","MC-8","IPC-5","Unknown-2","Unknown-3"]

# Add scores for modules
data = read_csv("hsa04724.csv")
sc.tl.score_genes(adata_concat, gene_list = data['x'].tolist(), score_name = 'Glutamatergic_synapse_module_score')

data = read_csv("hsa04727.csv")
sc.tl.score_genes(adata_concat, gene_list = data['x'].tolist(), score_name = 'GABAergic_synapse_module_score')

data = read_csv("hsa04360.csv")
sc.tl.score_genes(adata_concat, gene_list = data['x'].tolist(), score_name = 'Axon_guidance_module_score')
