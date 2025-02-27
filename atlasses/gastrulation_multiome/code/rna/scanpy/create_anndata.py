# NOTE: WE CREATED THE ANNDATA FILE USING (...)/other/conversions/rna/convert_SingleCellExperiment_to_anndata.R

exec(open('../../settings.py').read())
exec(open('../../utils.py').read())

######################
## Define arguments ##
######################

p = argparse.ArgumentParser( description='' )
p.add_argument( '--input_directory',               type=str,                       help='Input directory for the 10x output files' )
p.add_argument( '--outfile',               type=str,                       help='Output file (anndata)' )
p.add_argument( '--metadata',               type=str,                       help='Metadata file' )
p.add_argument( '--samples',               type=str, nargs="+",                       help='Samples' )
args = p.parse_args()

#####################
## Define settings ##
#####################

## START TEST ##
args.outfile = io["basedir"]+"/processed/rna/anndata.h5ad"
args.input_directory = io["basedir"]+"/original"
args.metadata = io["basedir"]+"/results_new2/rna/mapping/sample_metadata_after_mapping.txt.gz"
args.samples = ["E7.5_rep1", "E7.5_rep2"]
## END TEST ##


###########################
## Create anndata object ##
###########################

print("Loading 10x output files...")

anndatas = [None for i in range(len(args.samples))]
for i in range(len(args.samples)):
    sample = args.samples[i]
    anndatas[i] = sc.read_10x_mtx(path=args.input_directory+"/"+sample+"/filtered_feature_bc_matrix", gex_only=True)
    # anndatas[i].obs["sample"] = sample
    anndatas[i].obs["cell"] = sample + "#" + anndatas[i].obs.index
    anndatas[i].obs.set_index("cell", inplace=True, drop=True)
    print(sample)
    print(anndatas[i].shape)


# Concatenate anndata experiments
# AnnData.concatenate(*adatas, join='inner', batch_key='batch', batch_categories=None, uns_merge=None, index_unique='-', fill_value=None)
adata = anndata.AnnData.concatenate(*anndatas, join='inner', batch_key=None, index_unique=None)
print(adata.shape)

# Sanity checks
assert adata.obs.index.duplicated().sum()==0
assert adata.var.index.duplicated().sum()==0

#################################################
## Add existing metadata to the AnnData object ##
#################################################

print("Loading metadata...")

# Load existing metadata
metadata = pd.read_csv(args.metadata, sep="\t")
metadata.set_index("cell", inplace=True, drop=True)


# foo = pd.merge(left=adata.obs, right=metadata, left_on=["cell"], right_on=['cell']).set_index("cell")
foo = pd.merge(left=adata.obs, right=metadata, left_index=True, right_index=True)
assert adata.shape[0] == foo.shape[0]
adata.obs = foo

##########
## Save ##
##########

print("Saving anndata object...")

adata.write_h5ad(args.outfile)