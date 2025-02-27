from celloracle.trajectory.oracle_utility import _adata_to_df
from celloracle.trajectory.oracle_utility import intersect
from tqdm.notebook import tqdm
from sklearn.linear_model import Ridge
from celloracle.visualizations.development_module_visualization import plot_cluster_cells_use
from celloracle.trajectory.oracle_utility import (_adata_to_matrix, _adata_to_df,
                             _adata_to_color_dict, _get_clustercolor_from_anndata,
                             _numba_random_seed, _linklist2dict,
                             _decompose_TFdict, _is_perturb_condition_valid)
from celloracle.trajectory.oracle_utility import (_adata_to_matrix, _adata_to_df )
from celloracle.visualizations.config import CONFIG
from celloracle.visualizations.development_module_visualization import *

def fit_GRN_for_simulation(self, GRN_unit="cluster", alpha=1, use_cluster_specific_TFdict=False):
        
        """
        Do GRN inference.
        Please see the paper of CellOracle paper for details.
        GRN can be constructed for the entire population or each clusters.
        If you want to infer cluster-specific GRN, please set [GRN_unit="cluster"].
        You can select cluster information when you import data.
        If you set [GRN_unit="whole"], GRN will be made using all cells.
        Args:
            GRN_unit (str): Select "cluster" or "whole"
            alpha (float or int): The strength of regularization.
                If you set a lower value, the sensitivity increases, and you can detect weaker network connections. However, there may be more noise.
                If you select a higher value, it will reduce the chance of overfitting.
        """
        # prepare data for GRN calculation
        def _getCoefMatrix(gem, TFdict, alpha=1):

            genes = gem.columns

            all_genes_in_dict = intersect(gem.columns, list(TFdict.keys()))
            zero_ = pd.Series(np.zeros(len(genes)+1), index=genes.append(pd.Index(['Intercept'])))

            def get_coef(target_gene):
                tmp = zero_.copy()

                # define regGenes
                reggenes = TFdict[target_gene]
                reggenes = intersect(reggenes, genes)

                if target_gene in reggenes:
                    reggenes.remove(target_gene)

                if len(reggenes) == 0 :

                    tmp[target_gene] = 1
                    return(tmp)


                # prepare learning data
                Data = gem[reggenes]
                Label = gem[target_gene]


                # model fitting
                model = Ridge(alpha=alpha, random_state=123)
                model.fit(Data, Label)

                tmp[reggenes] = model.coef_
                tmp['Intercept'] = model.intercept_
                return tmp




            li = []
            li_calculated = []

            for i in tqdm(genes):
                if not i in all_genes_in_dict:
                    tmp = zero_.copy()
                    tmp[i] = 1

                else:
                    tmp = get_coef(i)
                    li_calculated.append(i)

                li.append(tmp)


            coef_matrix = pd.concat(li, axis=1)

            coef_matrix.columns = genes

            print(f"genes_in_gem: {gem.shape[1]}")
            print(f"models made for {len(li_calculated)} genes")


            return coef_matrix #, li_calculated
        
        gem_imputed = _adata_to_df(self.adata, "imputed_count")

        self.adata.layers["simulation_input"] = self.adata.layers["imputed_count"].copy()
        self.alpha_for_trajectory_GRN = alpha
        self.GRN_unit = GRN_unit

        if use_cluster_specific_TFdict & (self.cluster_specific_TFdict is not None):
            self.coef_matrix_per_cluster = {}
            cluster_info = self.adata.obs[self.cluster_column_name]

            print(f"fitting GRN again...")
            for cluster in np.unique(cluster_info):
                print(f"calculating GRN in {cluster}")
                cells_in_the_cluster_bool = (cluster_info == cluster)
                gem_ = gem_imputed[cells_in_the_cluster_bool]
                self.coef_matrix_per_cluster[cluster] = _getCoefMatrix(gem=gem_,
                                                                       TFdict=self.cluster_specific_TFdict[cluster],
                                                                       alpha=alpha)


        else:
            if GRN_unit == "whole":
                self.coef_matrix = _getCoefMatrix(gem=gem_imputed, TFdict=self.TFdict, alpha=alpha)
            if GRN_unit == "cluster":
                self.coef_matrix_per_cluster = {}
                cluster_info = self.adata.obs[self.cluster_column_name]
                for cluster in np.unique(cluster_info):
                    print(f"calculating GRN in {cluster}")
                    cells_in_the_cluster_bool = (cluster_info == cluster)
                    gem_ = gem_imputed[cells_in_the_cluster_bool]
                    self.coef_matrix_per_cluster[cluster] = _getCoefMatrix(gem=gem_,
                                                                           TFdict=self.TFdict,
                                                                           alpha=alpha)

        self.extract_active_gene_lists(verbose=False)
        return self
    
    
def do_simulation_(coef_matrix,  gem, n_propagation,method ,n_diff_genes):
    """
    Silulate signal propagation in GRNs.
    Args:
        coef_matrix (pandas.DataFrame): 2d matrix that store GRN weights
        simulation_input (pandas.DataFrame): input for simulation
        gem (pandas.DataFrame): input for simulation
        n_propagation (int): number of propagation.
        method (str):'wilcoxon' 't-test_overestim_var' or 't-test'
        n_diff_genes (int): differential genes
    Returns:
        pandas.DataFrame: simulated gene expression matrix after signal progagation
    """
    X_simulated = gem.copy()
    n_genes = gem.shape[1]
    X_simulated = X_simulated.dot(coef_matrix.iloc[0:n_genes,])
    for i in coef_matrix.columns:
        X_simulated[i] =  X_simulated[i] + coef_matrix.loc['Intercept',i]
    # gene expression cannot be negative. adjust delta values to make sure that gene expression are not netavive values.
    gem_simulated = X_simulated

#     # find differential genes
#     X_diff = X_simulated.append(gem)
#     X_diff.index = pd.Index(X_simulated.index + '_sim').append(X_simulated.index)
#     adata_diff = sc.AnnData(X_diff)
#     adata_diff.obs['group'] = ['simulated'] * gem.shape[0] + ['original'] * gem.shape[0] 
#     sc.tl.rank_genes_groups(adata_diff, 'group', method=method, key_added = method)
#     genes_up = sc.get.rank_genes_groups_df(adata_diff, group='simulated', key=method)['names'][:n_diff_genes]
#     genes_down = sc.get.rank_genes_groups_df(adata_diff, group='original', key=method)['names'][:n_diff_genes]

#     # set delta_input
#     delta_input = X_simulated.copy()
#     for col in delta_input.columns:
#         delta_input[col].values[:] = 0
#     for col in genes_up:
#         delta_input[col].values[:] = 2
#     for col in genes_down:
#         delta_input[col].values[:] = -2
    
#     # propagate delta_input
#     delta_simulated = delta_input.copy()
#     for i in range(n_propagation):
#         delta_simulated = delta_simulated.dot(coef_matrix.iloc[0:n_genes,])
#         delta_simulated[delta_input != 0] = delta_input

#         # gene expression cannot be negative. adjust delta values to make sure that gene expression are not netavive values.
#         gem_tmp = gem + delta_simulated
#         gem_tmp[gem_tmp<0] = 0
#         delta_simulated = gem_tmp - gem

#     gem_simulated = gem + delta_simulated
    return gem_simulated


def do_simulation(coef_matrix, gem, n_propagation):
    """
    Simulate signal propagation in GRNs.
    Args:
        coef_matrix (pandas.DataFrame): 2d matrix that store GRN weights
        gem (pandas.DataFrame): input for simulation
        n_propagation (int): number of propagation.
    Returns:
        pandas.DataFrame: simulated gene expression matrix after signal progagation
    """
    n_genes = gem.shape[1]
    X_simulated = gem
    for i in range(n_propagation):
        X_simulated = X_simulated.dot(coef_matrix.iloc[0:n_genes,])
        for i in coef_matrix.columns:
            X_simulated[i] =  X_simulated[i] + coef_matrix.loc['Intercept',i]
        # gene expression cannot be negative. adjust delta values to make sure that gene expression are not netavive values.
        X_simulated[X_simulated<0] = 0
    return X_simulated


def simulate_future(self,GRN_unit=None,n_propagation=2, ignore_warning=False):       
    # ,method = 'wilcoxon',n_diff_genes=20
    # 1. load gene expression matrix (initiation information for the simulation)
    gem_imputed = _adata_to_df(self.adata, "imputed_count")

    # 2. do simulation for signal propagation within GRNs
    if GRN_unit == "whole":
        gem_simulated = do_simulation(self.coef_matrix,
                                       simulation_input,
                                       gem_imputed,
                                       n_propagation)

    elif GRN_unit == "cluster":
        simulated = []
        cluster_info = self.adata.obs[self.cluster_column_name]
        for cluster in np.unique(cluster_info):
            cells_in_the_cluster_bool = (cluster_info == cluster)
            gem_ = gem_imputed[cells_in_the_cluster_bool]
            simulated_in_the_cluster = do_simulation(
                                         self.coef_matrix_per_cluster[cluster],
                                         gem_,
                                         n_propagation,
                                         )
            #method,n_diff_genes
            simulated.append(simulated_in_the_cluster)
        gem_simulated = pd.concat(simulated, axis=0)
        gem_simulated = gem_simulated.reindex(gem_imputed.index)
        
    else:
        raise ValueError("GRN_unit shold be either of 'whole' or 'cluster'")
    
    # 4. store simulation results
    #  simulated future gene expression matrix
    self.adata.layers["simulated_count"] = gem_simulated.values

    #  difference between simulated values and original values
    self.adata.layers["delta_X"] = self.adata.layers["simulated_count"] - self.adata.layers["imputed_count"]
    return self

def plot_cluster(self, ax=None, s=CONFIG["s_scatter"], args=CONFIG["default_args"],color=None):

    if ax is None:
        ax = plt
    col_dict = _get_clustercolor_from_anndata(adata=self.adata,
                                                  cluster_name=color,
                                                  return_as="dict")
    col_array = np.array([col_dict[i] for i in self.adata.obs[color]])
    ax.scatter(self.embedding[:, 0], self.embedding[:, 1], c=col_array, s=s, **args)
    ax.axis("off")
    
def set_color(self,color=None):
    col_dict = _get_clustercolor_from_anndata(adata=self.adata,
                                                  cluster_name=color,
                                                  return_as="dict")
    self.colorandum = np.array([col_dict[i] for i in self.adata.obs[color]])
    return(self)

def visualize_development_module(self, scale_for_pseudotime=CONFIG["scale_dev"],
    scale_for_simulation=CONFIG["scale_simulation"], s=CONFIG["s_scatter"], s_grid=CONFIG["s_grid"], vm=1, show_background=True):


    fig, ax = plt.subplots(2, 2, figsize=[10, 10])
    ax_ = ax[0, 0]
    plot_reference_flow_on_grid(self, ax=ax_, scale=scale_for_pseudotime, show_background=show_background, s=s)
    ax_.set_title("Development flow")

    ##
    ax_ = ax[0,1]
    plot_simulation_flow_on_grid(self, ax=ax_, scale=scale_for_simulation, show_background=show_background, s=s)
    ax_.set_title("GRN based velocity")

    ####
    ax_ = ax[1, 0]
    plot_inner_product_on_grid(self, ax=ax_, vm=vm,s=s_grid, show_background=show_background)
    ax_.set_title("Inner product of \n GRN based velocity * Development flow")


    ax_ = ax[1, 1]
    plot_inner_product_on_pseudotime(self, ax=ax_, vm=vm, s=s_grid)
    
def visualize_simulation_module(self, scale_for_pseudotime=CONFIG["scale_dev"],
    scale_for_simulation=CONFIG["scale_simulation"], s=CONFIG["s_scatter"], s_grid=CONFIG["s_grid"], vm=1, show_background=True):


    fig, ax = plt.subplots(2, 2, figsize=[10, 10])

    ax_ = ax[0, 0]
    plot_reference_flow_on_grid(self, ax=ax_, scale=scale_for_pseudotime, show_background=show_background, s=s)
    ax_.set_title("Development flow")

    ##
    ax_ = ax[0,1]
    plot_simulation_flow_on_grid(self, ax=ax_, scale=scale_for_simulation, show_background=show_background, s=s)
    ax_.set_title("Perturbation simulation")

    ####
    ax_ = ax[1, 0]
    plot_inner_product_on_grid(self, ax=ax_, vm=vm,s=s_grid, show_background=show_background)
    ax_.set_title("Inner product of \n Perturbation simulation * Development flow")


    ax_ = ax[1, 1]
    plot_inner_product_on_pseudotime(self, ax=ax_, vm=vm, s=s_grid)


