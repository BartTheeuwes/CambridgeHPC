######################
## Import libraries ##
######################

import os
from re import search

###########################
## Load default settings ##
###########################

if search("BI2404M", os.uname()[1]):
    exec(open('/Users/argelagr/gastrulation_multiome_10x/settings.py').read())
    exec(open('/Users/argelagr/gastrulation_multiome_10x/utils.py').read())
elif search("pebble|headstone", os.uname()[1]):
    exec(open('/bi/group/reik/ricard/scripts/gastrulation_multiome_10x/settings.py').read())
    exec(open('/bi/group/reik/ricard/scripts/gastrulation_multiome_10x/utils.py').read())
else:
    exit("Computer not recognised")

################################
## Initialise argument parser ##
################################

p = argparse.ArgumentParser( description='' )
p.add_argument( '--anndata',               type=str,                required=True,           help='Anndata file')
# p.add_argument( '--metadata',               type=str,                required=True,           help='Cell metadata file')
p.add_argument( '--outdir',               type=str,                required=True,           help='Output directory')
p.add_argument( '--trajectory_name',               type=str,                required=True,           help='Trajectory')
args = p.parse_args()

# convert args to dictionary
args = vars(args)

#####################
## Parse arguments ##
#####################

# I/O
if not os.path.isdir(args["outdir"]): os.makedirs(args["outdir"])
args["outdir"] = Path(args["outdir"])

sc.settings.figdir = args["outdir"] / "pdf"

print("Infering trajectory %s using metacells..." % args["trajectory_name"])

print(args)

##################
## Load anndata ##
##################

# adata = read...(adata, SEACells_label='SEACell', summarize_layer='raw')
adata = load_adata(
    adata_file = args["anndata"], 
    # metadata_file = io["metadata"],
    # cells = metadata.cell.values, 
    normalise = True, 
    filter_lowly_expressed_genes = True
)
adata

print(adata.obs["celltype"].value_counts())

#########################################
## Normalisation and feature selection ##
#########################################

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=1500)

#########
## PCA ##
#########

sc.tl.pca(adata, n_comps=15)

#########
## kNN ##
#########

sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_pca')

########################
## Force-atlas layout ##
########################

sc.tl.draw_graph(adata, layout="fa", init_pos=None)
sc.pl.draw_graph(adata, color=["celltype"], size=150, save="_metacell_trajectory.pdf")

##########
## Save ##
##########

to_save = pd.DataFrame(adata.obsm["X_draw_graph_fa"], index=adata.obs_names, columns=["FA1","FA2"])
to_save.to_csv(io["outdir"] + "/metacell_trajectory.txt.gz", sep='\t')
