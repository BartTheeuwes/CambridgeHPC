#!/g/stegle/ricard/anaconda/envs/mofa2/bin/python
#SBATCH -N 1                        # number of nodes
#SBATCH -n 1                        # number of cores
#SBATCH --mem 10G                   # memory pool for all cores
#SBATCH -t 0-10:00                   # runtime limit (D-HH:MM:SS)
#SBATCH -o slurm.%N.%j.out          # STDOUT
#SBATCH -e slurm.%N.%j.err          # STDERR
#SBATCH --mail-type=END,FAIL        # notifications for job done & fail

from mofapy2.run.entry_point import entry_point
import pandas as pd
import numpy as np
import argparse
from scipy.io import mmread, mmwrite
from pathlib import Path

################################
## Initialise argument parser ##
################################

p = argparse.ArgumentParser( description='' )

# I/O options
p.add_argument( '--input_folder',          type=str,              required=True,          help='Input data file (matrix format)' )
p.add_argument( '--outfile',               type=str,              required=True,          help='Output file to store the model (.hdf5)' )

# Model options
p.add_argument( '--factors',               type=int,              default=25,             help='Number of factors' )

# Training options
p.add_argument( '--seed',                  type=int,              default=42,             help='Random seed' )
p.add_argument( '--cores',                  type=int,              default=42,             help='Number of cores' )
p.add_argument( '--test',               action="store_true",                           help='Do stochastic inference?' )
p.add_argument( '--convergence_mode',      type=str,              default="fast",       help='Convergence mode')
args = p.parse_args()

## START TEST ##
# args.input_folder = "/bi/group/reik/ricard/data/gastrulation_multiome_10x/results_new/rna_atac/mofa/all_cells"
# args.outfile = "/bi/group/reik/ricard/data/gastrulation_multiome_10x/results_new/rna_atac/mofa/all_cells/test.hdf5"
# args.factors = 25
# args.seed = 42
# args.convergence_mode = "fast"
## END TEST ##

###############
## Load data ##
###############

input_folder = Path(args.input_folder)

rna_mtx = mmread(input_folder/"rna.mtx").todense().T
atac_mtx = mmread(input_folder/"atac_tfidf.mtx").todense().T

rna_features = pd.read_csv(input_folder/"rna_features.txt", header=None)[0].tolist()
atac_features = pd.read_csv(input_folder/"atac_features.txt", header=None)[0].tolist()
cells = pd.read_csv(input_folder/"cells.txt", header=None)[0].tolist()

# sample_metadata = pd.read_csv(input_folder/"sample_metadata.txt.gz")
# assert sample_metadata.cell.tolist() == cells

########################
## Create MOFA object ##
########################

# initialise entry point    
ent = entry_point()

# Set data
ent.set_data_matrix(
	data = [[rna_mtx], [atac_mtx]], 
	views_names = ["RNA","ATAC"], 
	samples_names = [ cells ], 
	features_names = [ rna_features, atac_features ]
)

# Set data options
ent.set_data_options(use_float32 = True)

# Set model options
ent.set_model_options(factors=args.factors, spikeslab_factors=False, spikeslab_weights=False)

# Set training options
if args.test:
	ent.set_train_options(iter=5)
else:
	ent.set_train_options(convergence_mode=args.convergence_mode, seed=args.seed)

###############################
## Build and train the model ##
###############################


# Build the model
ent.build()

# Train the model
ent.run()

####################
## Save the model ##
####################

ent.save(args.outfile, save_data=True)

##########
## TEST ##
##########

# Check multithreading
# import os
# os.environ["OMP_NUM_THREADS"] = "4"
# os.environ["OPENBLAS_NUM_THREADS"] = "4"
# os.environ["MKL_NUM_THREADS"] = "6"
# os.environ["VECLIB_MAXIMUM_THREADS"] = "4"
# os.environ["NUMEXPR_NUM_THREADS"] = "6"