#!/usr/bin/env python
# coding: utf-8

# # Import libraries

# In[1]:


import os
from re import search
from dfply import *


# # Load settings

# In[2]:


if search("ricard", os.uname()[1]):
    exec(open('/Users/ricard/gastrulation_multiome_10x/settings.py').read())
    exec(open('/Users/ricard/gastrulation_multiome_10x/utils.py').read())
elif search("ebi", os.uname()[1]):
    exec(open('/homes/ricard/gastrulation_multiome_10x/settings.py').read())
    exec(open('/homes/ricard/gastrulation_multiome_10x/utils.py').read())
else:
    exit("Computer not recognised")


# ## Define I/O

# In[3]:


io["outdir"] = io["basedir"] + "/..."


# ## Define options 

# scanpy options

# In[4]:


# %%capture
# sc.settings.verbosity = 3
# sc.logging.print_versions()
sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(8, 7), facecolor='white')


# In[5]:


opts["samples"] = [
	"E7.5_rep1",
	"E7.5_rep2",
	"E8.0_rep1",
	"E8.0_rep2",
	"E8.5_rep1",
	"E8.5_rep2"
]

opts["celltypes"] = [
   # "Nascent_mesoderm",
   "Mixed_mesoderm",
   # "Allantois",
   # "ExE_mesoderm",
   # "Mesenchyme",
   "Haematoendothelial_progenitors",
   "Endothelium",
   "Blood_progenitors_1",
   "Blood_progenitors_2",
   "Erythroid1"
]


# ## Load cell metadata

# In[6]:


metadata = (pd.read_table(io["metadata"]) >>
    mask(X["pass_rnaQC"]==True, X["doublet_call"]==False) >>
    mask(X["sample"].isin(opts["samples"]), X["celltype.mapped"].isin(opts["celltypes"]))
)
metadata.shape


# In[7]:


metadata.head()


# # Load anndata object

# In[8]:


adata = load_adata(
    adata_file = io["anndata"], 
    cells = metadata.cell.values, 
    normalise = True, 
    filter_lowly_expressed_genes = True
)
adata


# In[9]:


colPalette_celltypes = [opts["celltype_colors"][i.replace(" ","_")] for i in sorted(np.unique(adata.obs['celltype.mapped']))]
adata.uns['celltype.mapped_colors'] = colPalette_celltypes
colPalette_stages = [opts["stages_colors"][i.replace(" ","_")] for i in sorted(np.unique(adata.obs['stage']))]
adata.uns['stage_colors'] = colPalette_stages


# # Dimensionality reduction

# ## PCA

# Run PCA

# In[10]:


sc.tl.pca(adata, svd_solver='arpack')


# Plot PCA

# In[11]:


sc.pl.pca(adata, components=[1,2], color=["celltype.mapped"], size=25, legend_loc=None)


# ## k-NN graph

# Build kNN graph

# In[12]:


sc.pp.neighbors(adata, n_neighbors=15, n_pcs=15)


# ## UMAP

# Run UMAP

# In[13]:


sc.tl.umap(adata, min_dist=0.5, n_components=2)


# Plot UMAP

# In[18]:


sc.pl.umap(adata, color=["celltype.mapped","stage"], size=25, legend_loc="on data")


# ## Force-directed layout

# In[15]:


sc.tl.draw_graph(adata, layout="fa")


# In[17]:


sc.pl.draw_graph(adata, color=["celltype.mapped","stage"], size=25, legend_loc="on data")


# In[ ]:




