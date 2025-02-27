
import os
from re import search
from dfply import *
import scvelo as scv


if search("ricard", os.uname()[1]):
    exec(open('/Users/ricard/gastrulation_multiome_10x/settings.py').read())
    exec(open('/Users/ricard/gastrulation_multiome_10x/utils.py').read())
elif search("ebi", os.uname()[1]):
    exec(open('/homes/ricard/gastrulation_multiome_10x/settings.py').read())
    exec(open('/homes/ricard/gastrulation_multiome_10x/utils.py').read())
else:
    exit("Computer not recognised")


io["outfile"] = io["basedir"] + "/processed/rna/velocyto/anndata_scvelo.h5ad"


opts["samples"] = [
	"E7.5_rep1",
	"E7.5_rep2",
	"E8.0_rep1",
	"E8.0_rep2",
	"E8.5_rep1",
	"E8.5_rep2"
]


# # Load metadata

metadata = (pd.read_table(io["metadata"]) >>
    mask(X.pass_rnaQC==True, X.doublet_call==False) >>
    mask(X["sample"].isin(opts["samples"]))
).set_index("cell", drop=False)
metadata.shape

metadata.index.values


# # Load anndata object

adata = load_adata(adata_file = io["anndata"], metadata_file = io["metadata"], normalise = False, cells = metadata.index.values)


# # Load spliced and unspliced counts from loom files

rename_dict = {
    "E8.5_rep1" : "multiome1",
    "E8.5_rep2" : "multiome2",
    "E8.0_rep1" : "E8_0_rep1_multiome.loom",
    "E8.0_rep2" : "E8_0_rep2_multiome.loom",
    "E7.5_rep1" : "rep1_L001_multiome",
    "E7.5_rep2" : "rep2_L002_multiome"
}

looms = [None for i in range(len(opts["samples"]))]
for i in range(len(opts["samples"])):
    io["loom_velocyto"] = io["basedir"] + "/processed/rna/velocyto/" + opts["samples"][i] + ".loom"
    looms[i] = sc.read_loom(io["loom_velocyto"], sparse=True, X_name='spliced', obs_names='CellID', obsm_names=None, var_names='Gene')
    looms[i].var_names_make_unique()
    looms[i].obs.index = looms[i].obs.index.str.replace(rename_dict[opts["samples"][i]]+":",opts["samples"][i]+"_").str.replace("x","-1")
    print(looms[i].shape)
    print(looms[i].obs.head())

adata_loom = anndata.AnnData.concatenate(*looms, join='inner', batch_key=None, index_unique=None)
del looms

# Remove non-used layers to save memory
del adata_loom.layers["ambiguous"]
del adata_loom.layers["matrix"]


# Merge anndata objects

adata_final = scv.utils.merge(adata, adata_loom)
del adata_loom
del adata
adata_final

adata_final.obs.index.name = None


# # Save anndata object

adata_final.write_h5ad(io["outfile"])
