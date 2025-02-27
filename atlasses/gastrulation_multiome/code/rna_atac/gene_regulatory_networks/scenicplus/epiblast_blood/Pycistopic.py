#!/usr/bin/env python
# coding: utf-8

# # Pycistopic

# In[4]:


import warnings
warnings.simplefilter(action='ignore')
import pycisTopic
pycisTopic.__version__
import pickle

# In[5]:


# Project directory
projDir = '/rds/project/rds-SDzz0CATGms/users/bt392/atlasses/gastrulation_multiome/results/rna_atac/gene_regulatory_networks/scenicplus/blood/'
# Output directory
outDir = projDir + 'pycistopic/'
import os
if not os.path.exists(outDir):
    os.makedirs(outDir)

work_dir = outDir
# In[6]:


# Temp dir
tmpDir = '/home/bt392/ry/'


# ## 1. Creating a cisTopic object

# In[7]:


# Create cisTopic object
from pycisTopic.cistopic_class import *
#count_matrix=pd.read_csv(projDir+'epiblast_blood_atac_mtx.tsv', sep='\t')


# In[8]:


#count_matrix = count_matrix.set_index('cell')


# In[9]:


##path_to_blacklist='/staging/leuven/stg_00002/lcb/cbravo/Multiomics_pipeline/pycisTopic/blacklist/hg19-blacklist.v2.bed'
#cistopic_obj = create_cistopic_object(fragment_matrix=count_matrix) #, path_to_blacklist=path_to_blacklist)
## Adding cell information
#cistopic_obj = pickle.load(open(outDir+'cisTopicObject.pkl', 'rb')) # hash this out, was just to fix metadata!
#cell_data =  pd.read_csv(projDir+'epiblast_blood_sample_metadata.tsv', sep='\t')
#cell_data = cell_data.set_index('cell')
#cistopic_obj.add_cell_data(cell_data)


# In[10]:


#print(cistopic_obj)


# ## 2. Run models

# In[ ]:


import os
path_to_mallet_binary='mallet'
# Run models
#models= run_cgs_models(cistopic_obj,  # run_cgs_models_mallet(path_to_mallet_binary,  #  tmp_path
#                    n_topics=[2,4,10,15,25,35],
#                    n_cpu=30,
#                    n_iter=500,
#                    random_state=555,
#                    alpha=50,
#                    alpha_by_topic=True,
#                    eta=0.1,
#                    eta_by_topic=False,
#                    _temp_dir  =tmpDir, #Use SCRATCH if many models or big data set
#                    save_path=None)


# In[ ]:


# Save
#import pickle
#with open(outDir+'models_500.pkl', 'wb') as f:
#  pickle.dump(models, f)


# ## 3. Model selection

#cistopic_obj = pickle.load(open(outDir+'cisTopicObject.pkl', 'rb'))
#models = pickle.load(open(outDir+'models_500.pkl', 'rb'))


# In[ ]:

from pycisTopic.lda_models import *
if not os.path.exists(outDir+'models/'):
    os.mkdir(outDir+'models/')
    
#model=evaluate_models(models,
#                     select_model=None,
#                     return_model=True,
#                     metrics=['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
#                     plot_metrics=False,
#                     save= outDir + 'models/model_selection.pdf')


# In[ ]:


#from pycisTopic.lda_models import *
#model = evaluate_models(models,
#                       select_model=16,
#                       return_model=True,
#                       metrics=['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
#                       plot_metrics=False)


# In[ ]:


#cistopic_obj.add_LDA_model(model)
#pickle.dump(cistopic_obj,
#            open(os.path.join(outDir, 'epiblast_blood_cisTopicObject.pkl'), 'wb'))


# In[ ]:


# # Save
# with open(outDir + 'epiblast_blood_cisTopicObject.pkl', 'wb') as f:
#   pickle.dump(cistopic_obj, f)


# ## 4. Clustering & Visualisation

# In[ ]:


# Load cisTopic object
#import pickle
#infile = open(outDir + 'epiblast_blood_cisTopicObject.pkl', 'rb')
#cistopic_obj = pickle.load(infile)
#infile.close()


# In[ ]:


from pycisTopic.clust_vis import *
#find_clusters(cistopic_obj,
#                 target  = 'cell',
#                 k = 10,
#                 res = [0.6],
#                 prefix = 'pycisTopic_',
#                 scale = True,
#                 split_pattern = '-')


# In[ ]:


#run_umap(cistopic_obj,
#                 target  = 'cell', scale=True)


# In[ ]:


#run_tsne(cistopic_obj,
#                 target  = 'cell', scale=True)


# In[ ]:

#if not os.path.exists(outDir+'visualization/'):
#    os.mkdir(outDir+'visualization/')
#plot_metadata(cistopic_obj,
#                 reduction_name='UMAP',
#                 variables=['celltype', 'stage', 'pycisTopic_leiden_10_0.6'], # Labels from RNA and new clusters
#                 target='cell', num_columns=3,
#                text_size=10,
#                 dot_size=5,
#                 figsize=(15,5),
#                 save= outDir + 'visualization/dimensionality_reduction_label.pdf')


# In[ ]:


#annot_dict={}
#annot_dict['pycisTopic_leiden_10_0.6'] = {'0':'MM029 (0)', '1':'MM001 (1)', '2': 'MM034 (2)', '3': 'MM047 (3)', '4': 'MM011 (4)'}
#cistopic_obj.cell_data['pycisTopic_leiden_10_0.6'] = [annot_dict['pycisTopic_leiden_10_0.6'][x] for x in cistopic_obj.cell_data['pycisTopic_leiden_10_0.6'].tolist()]


# In[ ]:


#plot_metadata(cistopic_obj,
#                 reduction_name='UMAP',
#                 variables=['celltype', 'stage', 'pycisTopic_leiden_10_0.6'], # Labels from RNA and new clusters
#                 target='cell', num_columns=3,
#                 text_size=10,
#                 dot_size=5,
#                 figsize=(15,5),
#                 save= outDir + 'visualization/dimensionality_reduction_label.pdf')


# In[ ]:


#plot_topic(cistopic_obj,
#            reduction_name = 'UMAP',
#            target = 'cell',
#            num_columns=5,
#            save= outDir + 'visualization/dimensionality_reduction_topic_contr.pdf')


# In[ ]:


#cell_topic_heatmap(cistopic_obj,
#                     variables = ['celltype', 'stage'],
#                     scale = False,
#                     legend_loc_x = 1.05,
#                     legend_loc_y = -1.2,
#                     legend_dist_y = -1,
#                     figsize=(10,10),
#                     save = outDir + 'visualization/heatmap_topic_contr.pdf')


# In[ ]:


# Save
#with open(outDir + 'epiblast_blood_cisTopicObject.pkl', 'wb') as f:
#  pickle.dump(cistopic_obj, f)


# ## 5. Topic binarization & qc

# In[ ]:


# Load cisTopic object
#import pickle
#infile = open(outDir + 'epiblast_blood_cisTopicObject.pkl', 'rb')
#cistopic_obj = pickle.load(infile)
#infile.close()


# In[ ]:

if not os.path.exists(outDir+'topic_binarization/'):
    os.mkdir(outDir+'topic_binarization/')
    
#from pycisTopic.topic_binarization import *
#region_bin_topics = binarize_topics(cistopic_obj, method='otsu', ntop=3000, plot=True, num_columns=5, save= outDir + 'topic_binarization/otsu.pdf')


# In[ ]:


#binarized_cell_topic = binarize_topics(cistopic_obj, target='cell', method='li', plot=True, num_columns=5, nbins=60)


# In[ ]:


#from pycisTopic.topic_qc import *
#topic_qc_metrics = compute_topic_metrics(cistopic_obj)


# In[ ]:


#fig_dict={}
#fig_dict['CoherenceVSAssignments']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Log10_Assignments', var_color='Gini_index', plot=False, return_fig=True)
#fig_dict['AssignmentsVSCells_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Log10_Assignments', var_y='Cells_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
#fig_dict['CoherenceVSCells_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Cells_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
#fig_dict['CoherenceVSRegions_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Regions_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
#fig_dict['CoherenceVSMarginal_dist']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Marginal_topic_dist', var_color='Gini_index', plot=False, return_fig=True)
#fig_dict['CoherenceVSGini_index']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Gini_index', var_color='Gini_index', plot=False, return_fig=True)


# In[ ]:


# Plot topic stats in one figure
#fig=plt.figure(figsize=(40, 43))
#i = 1
#for fig_ in fig_dict.keys():
#    plt.subplot(2, 3, i)
#    img = fig2img(fig_dict[fig_]) #To convert figures to png to plot together, see .utils.py. This converts the figure to png.
#    plt.imshow(img)
#    plt.axis('off')
#    i += 1
#plt.subplots_adjust(wspace=0, hspace=-0.70)
#fig.savefig(outDir + 'topic_binarization/Topic_qc.pdf', bbox_inches='tight')
#plt.show()


# In[ ]:


#topic_annot = topic_annotation(cistopic_obj, annot_var='celltype', binarized_cell_topic=binarized_cell_topic, general_topic_thr = 0.2)


# In[ ]:


#topic_annot


# In[ ]:


#topic_qc_metrics = pd.concat([topic_annot[['celltype', 'Ratio_cells_in_topic', 'Ratio_group_in_population']], topic_qc_metrics], axis=1)


# In[ ]:


#topic_qc_metrics


# In[ ]:


# Save
#with open(outDir + 'topic_binarization/Topic_qc_metrics_annot.pkl', 'wb') as f:
#  pickle.dump(topic_qc_metrics, f)
#with open(outDir + 'topic_binarization/binarized_cell_topic.pkl', 'wb') as f:
#  pickle.dump(binarized_cell_topic, f)
#with open(outDir + 'topic_binarization/binarized_topic_region.pkl', 'wb') as f:
#  pickle.dump(region_bin_topics, f)


# ## 6. Differentially Accessible Regions (DARs)

# In[ ]:


# Load cisTopic object
import pickle
#infile = open(outDir + 'epiblast_blood_cisTopicObject.pkl', 'rb')
#cistopic_obj = pickle.load(infile)
#infile.close()


# In[ ]:


#from pycisTopic.diff_features import *
#imputed_acc_obj = impute_accessibility(cistopic_obj, selected_cells=None, selected_regions=None, scale_factor=10**6)


# In[ ]:


#normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)


# In[ ]:


if not os.path.exists(os.path.join(work_dir, 'DARs/')):
    os.makedirs(os.path.join(work_dir, 'DARs/'))
#variable_regions = find_highly_variable_features(normalized_imputed_acc_obj,
#                                           min_disp = 0.05,
#                                           min_mean = 0.0125,
#                                           max_mean = 3,
#                                           max_disp = np.inf,
#                                           n_bins=20,
#                                           n_top_features=None,
#                                           plot=True,
#                                           save= outDir + 'DARs/HVR_plot.pdf')


# In[ ]:


#len(variable_regions)


# In[ ]:


#markers_dict= find_diff_features(cistopic_obj,
#                      imputed_acc_obj,
#                      variable='celltype',
#                      var_features=variable_regions,
#                      contrasts=None,
#                      adjpval_thr=0.05,
#                      log2fc_thr=np.log2(1.5),
#                      n_cpu=30,
#                      _temp_dir=tmpDir,
#                      split_pattern = '-')


# In[ ]:


#from pycisTopic.clust_vis import *
#plot_imputed_features(cistopic_obj,
#                    reduction_name='UMAP',
#                    imputed_data=imputed_acc_obj,
#                    features=[markers_dict[x].index.tolist()[0] for x in ['Epiblast', 'Blood_progenitors_1', 'Endothelium', 'Erythroid3']],
#                    scale=False,
#                    num_columns=4,
#                    save= outDir + 'DARs/example_best_DARs.pdf')


# In[ ]:


#x = [print(x + ': '+ str(len(markers_dict[x]))) for x in markers_dict.keys()]


# In[ ]:





# In[ ]:


if not os.path.exists(os.path.join(work_dir, 'candidate_enhancers/')):
    os.makedirs(os.path.join(work_dir, 'candidate_enhancers/'))
#import pickle
#pickle.dump(markers_dict, open(os.path.join(work_dir, 'candidate_enhancers/markers_dict.pkl'), 'wb'))
#pickle.dump(imputed_acc_obj, open(outDir + 'DARs/Imputed_accessibility.pkl', 'wb'))


# # Motif enrichment analysis using pycistarget

# In[ ]:


import pickle
region_bin_topics_otsu = pickle.load(open(os.path.join(work_dir, 'topic_binarization/binarized_topic_region.pkl'), 'rb'))
#region_bin_topics_top3k = pickle.load(open(os.path.join(work_dir, 'scATAC/candidate_enhancers/region_bin_topics_top3k.pkl'), 'rb'))
markers_dict = pickle.load(open(os.path.join(work_dir, 'candidate_enhancers/markers_dict.pkl'), 'rb'))


# In[ ]:


import pyranges as pr
from pycistarget.utils import region_names_to_coordinates
region_sets = {}
region_sets['topics_otsu'] = {}
# region_sets['topics_top_3'] = {}
region_sets['DARs'] = {}
for topic in region_bin_topics_otsu.keys():
    regions = region_bin_topics_otsu[topic].index[region_bin_topics_otsu[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['topics_otsu'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
# for topic in region_bin_topics_top3k.keys():
#     regions = region_bin_topics_top3k[topic].index[region_bin_topics_top3k[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
#     region_sets['topics_top_3'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
for DAR in markers_dict.keys():
    regions = markers_dict[DAR].index[markers_dict[DAR].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['DARs'][DAR] = pr.PyRanges(region_names_to_coordinates(regions))


# In[ ]:


for key in region_sets.keys():
    print(f'{key}: {region_sets[key].keys()}')


# In[ ]:


db_fpath = "/rds/project/rds-SDzz0CATGms/users/bt392/software/cistarget_database/mm10"
motif_annot_fpath = "/rds/project/rds-SDzz0CATGms/users/bt392/software/cistarget_database/mm10"


# In[ ]:


rankings_db = os.path.join(db_fpath, 'mm10_screen_v10_clust.regions_vs_motifs.rankings.feather')
scores_db =  os.path.join(db_fpath, 'mm10_screen_v10_clust.regions_vs_motifs.scores.feather')
motif_annotation = os.path.join(motif_annot_fpath, 'motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl')


# In[ ]:


if not os.path.exists(os.path.join(work_dir, 'motifs')):
    os.makedirs(os.path.join(work_dir, 'motifs'))


# In[ ]:


from scenicplus.wrappers.run_pycistarget import run_pycistarget
run_pycistarget(
    region_sets = region_sets,
    species = 'mus_musculus',
    save_path = os.path.join(work_dir, 'motifs/'),
    ctx_db_path = rankings_db,
    dem_db_path = scores_db,
    path_to_motif_annotations = motif_annotation,
    run_without_promoters = True,
    n_cpu = 30,
    _temp_dir = tmpDir,
    annotation_version = 'v10nr_clust',
    )


# In[ ]:


#import dill
#menr = dill.load(open(os.path.join(work_dir, 'motifs/menr.pkl'), 'rb'))


# In[ ]:


#menr['DEM_topics_otsu_All'].DEM_results('Topic8')

