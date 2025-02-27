#!/usr/bin/env python
# coding: utf-8

# # Running SCENIC

# In[3]:


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
work_dir = '/rds/project/rds-SDzz0CATGms/users/bt392/09_Eomes_invitro_blood/results/multiome_atlas/SCENICplus_sc/'
projDir = work_dir
tmp_dir = '/home/bt392/ry/'


# #### RNA 

# In[21]:
adata = sc.read_mtx(projDir + 'rna.mtx') # need to check if this is not transposed the wrong way
adata = adata.T
cluster_metadata = pd.read_csv(projDir + 'meta_subsampled.txt')
genes = pd.read_csv(projDir + 'genes.txt', sep = ' ')

cluster_metadata = cluster_metadata.set_index('cell')
genes = genes.set_index('x')

adata.obs = cluster_metadata
adata.var = genes

# process
adata.raw = adata
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)

sc.pp.pca(adata) #calculate PCA embedding


# #### ATAC

# In[6]:


# Load ATAC
cistopic_obj = dill.load(open(os.path.join(work_dir, 'pycistopic/cisTopicObject.pkl'), 'rb'))
menr = dill.load(open(os.path.join(work_dir, 'pycistopic/motifs/menr.pkl'), 'rb'))


# # Start SCENIC+ analysis

# In[14]:


#adata.obs


# In[15]:


#cistopic_obj.cell_data


# In[22]:


from scenicplus.scenicplus_class import create_SCENICPLUS_object
import numpy as np
scplus_obj = create_SCENICPLUS_object(
    GEX_anndata = adata.raw.to_adata(),
    cisTopic_obj = cistopic_obj,
    menr = menr,
    bc_transform_func = lambda x: f'{x}___cisTopic' #function to convert scATAC-seq barcodes to scRNA-seq ones
)
scplus_obj.X_EXP = np.array(scplus_obj.X_EXP.todense())
scplus_obj


# In[17]:


ensembl_version_dict = {'105': 'http://www.ensembl.org',
                        '104': 'http://may2021.archive.ensembl.org/'}

import pybiomart as pbm
def test_ensembl_host(scplus_obj, host, species):
    dataset = pbm.Dataset(name=species+'_gene_ensembl',  host=host)
    annot = dataset.query(attributes=['chromosome_name', 'transcription_start_site', 'strand', 'external_gene_name', 'transcript_biotype'])
    annot.columns = ['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
    annot['Chromosome'] = annot['Chromosome'].astype('str')
    filter = annot['Chromosome'].str.contains('CHR|GL|JH|MT')
    annot = annot[~filter]
    annot.columns=['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
    gene_names_release = set(annot['Gene'].tolist())
    ov=len([x for x in scplus_obj.gene_names if x in gene_names_release])
    print('Genes recovered: ' + str(ov) + ' out of ' + str(len(scplus_obj.gene_names)))
    return ov

n_overlap = {}
for version in ensembl_version_dict.keys():
    print(f'host: {version}')
    try:
        n_overlap[version] =  test_ensembl_host(scplus_obj, ensembl_version_dict[version], 'mmusculus')
    except:
        print('Host not reachable')
v = sorted(n_overlap.items(), key=lambda item: item[1], reverse=True)[0][0]
print(f"version: {v} has the largest overlap, use {ensembl_version_dict[v]} as biomart host")


# In[18]:


biomart_host = "http://may2021.archive.ensembl.org/"

# In[2]:


#get_ipython().system('wget -O /rds/project/rds-SDzz0CATGms/users/bt392/software/bedToBigBed http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed')
#get_ipython().system('chmod +x /rds/project/rds-SDzz0CATGms/users/bt392/software/bedToBigBed')


# In[23]:


#only keep the first two columns of the PCA embedding in order to be able to visualize this in SCope
scplus_obj.dr_cell['GEX_X_pca'] = scplus_obj.dr_cell['GEX_X_pca'].iloc[:, 0:2]
#scplus_obj.dr_cell['GEX_rep'] = scplus_obj.dr_cell['GEX_rep'].iloc[:, 0:2]


# In[ ]:


#get_ipython().system('wget -O pbmc_tutorial/data/utoronto_human_tfs_v_1.01.txt  http://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt')
# mouse variant?


# In[ ]:


from scenicplus.wrappers.run_scenicplus import run_scenicplus
try:
    run_scenicplus(
        scplus_obj = scplus_obj,
        variable = ['GEX_celltype'],
        species = 'mmusculus',
        assembly = 'mm10',
        tf_file = "/rds/project/rds-SDzz0CATGms/users/bt392/software/github/scenicplus/resources/allTFs_mm.txt",
        save_path = os.path.join(work_dir, 'scenicplus'),
        biomart_host = biomart_host,
        upstream = [1000, 150000],
        downstream = [1000, 150000],
        calculate_TF_eGRN_correlation = True,
        calculate_DEGs_DARs = True,
        export_to_loom_file = True,
        export_to_UCSC_file = True,
        path_bedToBigBed = '/rds/project/rds-SDzz0CATGms/users/bt392/software/',
        n_cpu = 30,
        _temp_dir = os.path.join(tmp_dir, 'ray_spill'))
except Exception as e:
    #in case of failure, still save the object
    dill.dump(scplus_obj, open(os.path.join(work_dir, 'scenicplus/scplus_obj.pkl'), 'wb'), protocol=-1)
    raise(e)


# In[ ]:


dill.dump(scplus_obj, open(os.path.join(work_dir, 'scenicplus/scplus_obj.pkl'), 'wb'), protocol=-1)

