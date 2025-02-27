# Turn atlas from anndata to SCE object

import scanpy as sc
import pandas as pd

main = '/rds/project/rds-SDzz0CATGms/users/bt392/atlasses/extended/'
atlas_in = main + 'embryo_mnn_bbknn_shiny.h5ad'

atlas = sc.read_h5ad(atlas_in)

blood_meta = pd.read_csv(main + 'blood/blood_meta.csv')
atlas = atlas[atlas.obs['cell'].isin(blood_meta['cell'].to_numpy())]

import anndata2ri
from rpy2.robjects import r
from rpy2.robjects.conversion import localconverter
#anndata2ri.activate()

with localconverter(anndata2ri.converter):
    r.globalenv['adata'] = atlas

#sce = r('as(adata, "SingleCellExperiment")')
r("saveRDS(adata, paste0(main, 'atlas_sce_script.rds')")