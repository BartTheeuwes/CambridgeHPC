import scvi
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
sc.set_figure_params(figsize=(8, 8))

adata = sc.read_h5ad('/rds/project/rds-SDzz0CATGms/users/bt392/09_Eomes_invitro_blood/code/rna/mapping/scVI_test/anndata.h5ad')

# preprocessing
sc.pp.filter_genes(adata, min_counts=3)
adata.layers["counts"] = adata.X.copy() # preserve counts
adata.raw = adata # freeze the state in `.raw`

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2500,
    subset=True,
    layer="counts",
    flavor="cell_ranger",
    batch_key="stage"
)

#################
### Try 1
#################


scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    batch_key="stage",
    categorical_covariate_keys=["sample"]
)

model = scvi.model.SCVI(adata, n_latent=40)

model.train()

model.save("/rds/project/rds-SDzz0CATGms/users/bt392/09_Eomes_invitro_bloodcode/rna/mapping/scVI_test/scvi_model_try1/")

latent = model.get_latent_representation()
adata.obsm["X_scVI"] = latent
adata.layers["scvi_normalized"] = model.get_normalized_expression(
    library_size=10e4
)

adata.obs.tail()

# use scVI latent space for UMAP generation
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=35)
sc.tl.umap(adata, min_dist=0.4)

sc.pl.umap(
    adata,
    color=["celltype", "origin", "stage", "sample"],
    frameon=False,
    ncols=2, legend_loc='none',
    save = 'try1.png'
)

#################
### Try 2
#################


scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    batch_key="origin",
    categorical_covariate_keys=["stage", "sample"]
)

model = scvi.model.SCVI(adata, n_latent=40)

model.train()

model.save("/rds/project/rds-SDzz0CATGms/users/bt392/09_Eomes_invitro_bloodcode/rna/mapping/scVI_test/scvi_model_try2/")

latent = model.get_latent_representation()
adata.obsm["X_scVI"] = latent
adata.layers["scvi_normalized"] = model.get_normalized_expression(
    library_size=10e4
)

adata.obs.tail()

# use scVI latent space for UMAP generation
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=35)
sc.tl.umap(adata, min_dist=0.4)

sc.pl.umap(
    adata,
    color=["celltype", "origin", "stage", "sample"],
    frameon=False,
    ncols=2, legend_loc='none',
    save = 'try2.png'
)

#################
### Try 3
#################


scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    batch_key="sample",
    categorical_covariate_keys=["stage", "origin"]
)

model = scvi.model.SCVI(adata, n_latent=40)

model.train()

model.save("/rds/project/rds-SDzz0CATGms/users/bt392/09_Eomes_invitro_bloodcode/rna/mapping/scVI_test/scvi_model_try3/")

latent = model.get_latent_representation()
adata.obsm["X_scVI"] = latent
adata.layers["scvi_normalized"] = model.get_normalized_expression(
    library_size=10e4
)

adata.obs.tail()

# use scVI latent space for UMAP generation
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=35)
sc.tl.umap(adata, min_dist=0.4)

sc.pl.umap(
    adata,
    color=["celltype", "origin", "stage", "sample"],
    frameon=False,
    ncols=2, legend_loc='none',
    save = 'try3.png'
)




























