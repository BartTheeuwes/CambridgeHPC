import sys
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import cellrank as cr
#import bbknn
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import random
random.seed(123)

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')
scv.set_figure_params()

Data = "/rds/project/bg200/rds-bg200-hphi-gottgens/users/bt392/mouse/Mixl1_KO/data/mixl_chim_atlas.h5"
Results = "/rds/project/bg200/rds-bg200-hphi-gottgens/users/bt392/mouse/Mixl1_KO/data/mixl_chim_atlas_analysis.h5"
out_folder = "/rds/project/bg200/rds-bg200-hphi-gottgens/users/bt392/mouse/Mixl1_KO/data/"


# distance chim from atlas instead of other way around
adata = sc.read(Results)
adata.obs['origin'] = np.where(adata.obs['tdTom'] != 'nan', 'chimaera', 'atlas')
X_pca = adata.obsm['X_pca']
pca_atlas = sc.AnnData(X_pca[adata.obs['origin'] == 'atlas']).X
pca_chimaera = sc.AnnData(X_pca[adata.obs['origin'] == 'chimaera']).X

k=15

print('calculating distances')
cell_distances_ct = cdist(pca_atlas,pca_chimaera, metric='euclidean')
print('finished distances')

Adj = np.zeros(cell_distances_ct.shape, dtype=float)
indices = np.zeros((cell_distances_ct.shape[0], k), dtype=np.int_)

for irow, row in enumerate(cell_distances_ct):
    idcs = np.argpartition(row, k)[:k]
    indices[irow] = idcs

for irow, row in enumerate(indices):
    Adj[irow, row] = 1

adjacency_score = np.sum(Adj, axis=0)

print('finished rest')

chim = adata[adata.obs['origin'] == 'chimaera']
chim.obs['adjacency_score_chim'] = adjacency_score


chim.obs.to_csv(out_folder + 'chimaera_distances_complete.csv')