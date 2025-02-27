######################
## Import libraries ##
######################

import os
from re import search
from dfply import *
import seaborn as sns
from collections import Counter
import argparse

from metacells.core import Metacells
# from metacells import plot


################################
## Initialise argument parser ##
################################

p = argparse.ArgumentParser( description='' )
p.add_argument( '--text_outfile',               type=str,                required=True,           help='Output file to store the cell2metacell dataframe (.txt)' )
p.add_argument( '--anndata_outfile',               type=str,                required=True,           help='Output file to store the metacell anndata (.h5)' )
p.add_argument( '--samples',            type=str, nargs="+",             default="all",             help='samples to use ' )
p.add_argument( '--seed',                  type=int,                default=42,               help='Random seed' )
p.add_argument( '--number_metacells',            type=int,              default=1000,             help='Number of metacells' )
p.add_argument( '--n_iter',       type=int,              default=50,              help='Number of iterations')
args = p.parse_args()


## START TEST ##
# args = {}
# args["anndata_outfile"] = "/hps/nobackup2/research/stegle/users/ricard/gastrulation_multiome_10x/results/rna/metacells/rna_atac/anndata_metacell_rna_atac.h5ad"
# args["text_outfile"] = "/hps/nobackup2/research/stegle/users/ricard/gastrulation_multiome_10x/results/rna/metacells/rna_atac/cell2metacell_1000metacells.txt"
# args["samples"] = ["E8.5_rep2"]
# args["seed"] = 42
# args["number_metacells"] = 1000
# args["n_iter"] = 5
## END TEST ##

#####################
## Define settings ##
#####################

# Load default settings
if search("ricard", os.uname()[1]):
    exec(open('/Users/ricard/gastrulation_multiome_10x/settings.py').read())
    exec(open('/Users/ricard/gastrulation_multiome_10x/utils.py').read())
elif search("ebi", os.uname()[1]):
    exec(open('/homes/ricard/gastrulation_multiome_10x/settings.py').read())
    exec(open('/homes/ricard/gastrulation_multiome_10x/utils.py').read())
else:
    exit("Computer not recognised")

# I/O
io["pca_rna"] = io["basedir"] + "/results/rna/dimensionality_reduction/all_cells/E7.5_rep1-E7.5_rep2-E8.0_rep1-E8.0_rep2-E8.5_rep1-E8.5_rep2_pca_features2500_pcs30_batchcorrectionbysample.txt.gz"
io["pca_atac"] = io["basedir"] + "/results/atac/archR/dimensionality_reduction/PeakMatrix/all_cells/E7.5_rep1-E7.5_rep2-E8.0_rep1-E8.0_rep2-E8.5_rep1-E8.5_rep2_umap_nfeatures50000_ndims50_neigh45_dist0.45.txt.gz"

###################
## Load metadata ##
###################

metadata = (pd.read_table(io["metadata"]) >>
    # mask(X.pass_rnaQC==True, X.doublet_call==False, X["celltype.mapped"].isin(opts["celltypes"])) >>
    mask(X.pass_rnaQC==True, X.pass_atacQC==True, X.doublet_call==False, X["celltype.mapped"].isin(opts["celltypes"])) >>
    mask(X["sample"].isin(args.samples))
).set_index("cell", drop=False)


##################
## Load AnnData ##
##################

adata = load_adata(adata_file = io["anndata"], metadata_file = io["metadata"], normalise = True, cells = metadata.index.values)

##############################
## Dimensionality reduction ##
##############################

# Load precomputed PCA coordinates
# pca_mtx = pd.read_csv(io["pca_rna"]).set_index("cell", drop=True).loc[adata.obs.index].to_numpy()
pca_mtx = pd.read_csv(io["pca_atac"]).set_index("cell", drop=True).loc[adata.obs.index].to_numpy()
adata.obsm["X_pca"] = pca_mtx

# run PCA
# sc.tl.pca(adata, svd_solver='arpack')

# Plot PCA
# sc.pl.pca(adata, components=[1,2], color=["celltype.mapped","stage"], size=25, legend_loc=None)


###############
## Metacells ##
###############

# Fast kernel archetypal analysis. 
# - Finds archetypes and weights given annotated data matrix. 
# - Modifies annotated data matrix in place to include Metacell assignments in ad.obs['Metacell']
model = Metacells(adata, build_kernel_on="X_pca", n_metacells=args.number_metacells, verbose=True)

model.fit(n_iter=args.n_iter)

# Save cell2metacell assignment
to_save = adata.obs.loc[:,['Metacell',"cell"]]
to_save.to_csv(args.text_outfile, sep="\t", header=True, index=False)


# Plot (TO-DO)
# graph_model.plot_convergence()
# plot.plot_metacell_sizes(ad, bins=15)

##################################
## Create AnnData for metacells ##
##################################

# DOESNT WORK IN THE CLUSTER BECAUSE OF NA's in HDF5>=3: TypeError: Can't implicitly convert non-string objects to strings
# adata.raw = adata

# metacell_adata = model.summarize_by_metacell(aggregate_by='sum')

# # Normalise and log-transform
# # sc.pp.normalize_total(metacell_adata)
# # sc.pp.log1p(metacell_adata)

# # TO-DO
# # metacell_adata.uns["celltype.mapped_colors"] = adata.uns["celltype.mapped_colors"]
# # metacell_adata.uns["stage_colors"] = adata.uns["stage_colors"]

# Save
# metacell_adata.write_h5ad(args.anndata_outfile)
