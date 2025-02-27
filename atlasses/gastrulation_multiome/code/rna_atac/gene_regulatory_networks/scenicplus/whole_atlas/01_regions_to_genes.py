# All things needed to mimic scenic wrapper
import dill
import scanpy as sc
import os
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import pyranges
# Set stderr to null to avoid strange messages from ray
import sys
_stderr = sys.stderr
null = open(os.devnull,'wb')
work_dir = '/rds/project/rds-SDzz0CATGms/users/bt392/atlasses/gastrulation_multiome/results/rna_atac/gene_regulatory_networks/scenicplus/whole_atlas/'
tmp_dir = '/home/bt392/ry/'

# load scenic object
scplus_obj = dill.load(open(os.path.join(work_dir, 'scenicplus/scplus_obj.pkl'), 'rb'))

# biomart host
biomart_host = "http://sep2019.archive.ensembl.org/"

# Wrapper imports
from scenicplus.scenicplus_class import SCENICPLUS, create_SCENICPLUS_object
from scenicplus.preprocessing.filtering import *
from scenicplus.cistromes import *
from scenicplus.enhancer_to_gene import get_search_space, calculate_regions_to_genes_relationships, GBM_KWARGS
from scenicplus.enhancer_to_gene import export_to_UCSC_interact 
from scenicplus.utils import format_egrns, export_eRegulons
from scenicplus.eregulon_enrichment import *
from scenicplus.TF_to_gene import *
from scenicplus.grn_builder.gsea_approach import build_grn
from scenicplus.dimensionality_reduction import *
from scenicplus.RSS import *
from scenicplus.diff_features import *
from scenicplus.loom import *
from typing import Dict, List, Mapping, Optional, Sequence
import os
import dill
import time

# wrapper inputs
scplus_obj = scplus_obj
variable = ['GEX_celltype']
species = 'mmusculus'
assembly = 'mm10'
tf_file = "/rds/project/rds-SDzz0CATGms/users/bt392/software/github/scenicplus/resources/allTFs_mm.txt"
save_path = os.path.join(work_dir, 'scenicplus')
biomart_host = biomart_host
upstream = [1000, 150000]
downstream = [1000, 150000]
calculate_TF_eGRN_correlation = True
calculate_DEGs_DARs = True
export_to_loom_file = True
export_to_UCSC_file = True
path_bedToBigBed = 'pbmc_tutorial'
n_cpu = 30
_temp_dir = os.path.join(tmp_dir, 'ray_spill')

# Create logger
level = logging.INFO
log_format = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'
handlers = [logging.StreamHandler(stream=sys.stdout)]
logging.basicConfig(level=level, format=log_format, handlers=handlers)
log = logging.getLogger('SCENIC+_wrapper')

start_time = time.time()
check_folder = os.path.isdir(save_path)

# untangled scenic wrapper:

#only keep the first two columns of the PCA embedding in order to be able to visualize this in SCope
scplus_obj.dr_cell['GEX_X_pca'] = scplus_obj.dr_cell['GEX_X_pca'].iloc[:, 0:2]

if not check_folder:
    os.makedirs(save_path)
    log.info("Created folder : "+ save_path)

else:
    log.info(save_path + " folder already exists.")

if 'Cistromes' not in scplus_obj.uns.keys():
    log.info('Merging cistromes')
    merge_cistromes(scplus_obj)


if 'search_space' not in scplus_obj.uns.keys():
    log.info('Getting search space')
    get_search_space(scplus_obj,
                 biomart_host = biomart_host,
                 species = species,
                 assembly = assembly, 
                 upstream = upstream,
                 downstream = downstream)

if 'region_to_gene' not in scplus_obj.uns.keys():
    log.info('Inferring region to gene relationships')
    calculate_regions_to_genes_relationships(scplus_obj, 
                    ray_n_cpu = n_cpu, 
                    _temp_dir = _temp_dir,
                    importance_scoring_method = 'GBM',
                    importance_scoring_kwargs = GBM_KWARGS)


# Save object
dill.dump(scplus_obj, open(os.path.join(work_dir, 'scenicplus/scplus_obj.pkl'), 'wb'), protocol=-1)
