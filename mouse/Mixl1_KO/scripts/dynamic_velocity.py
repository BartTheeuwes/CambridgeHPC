import numpy as np
import pandas as pd
import loompy
import scanpy as sc
import scvelo as scv

adata_in = '/rds/project/bg200/rds-bg200-hphi-gottgens/users/bt392/mouse/Mixl1_KO/data/chim_velocity_final.h5ad'
adata_out = '/rds/project/bg200/rds-bg200-hphi-gottgens/users/bt392/mouse/Mixl1_KO/data/chim_velocity_final_3000.h5ad'

adata = sc.read(adata_in)

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=3000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

scv.tl.recover_dynamics(adata)
print('Saving recovered dynamics')
adata.write(adata_out)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)
print('saving velocity graph')
adata.write(adata_out)