def add_correlation(adjacencies: pd.DataFrame,
                    ex_mtx: pd.DataFrame,
                    rho_threshold=0.03,
                    mask_dropouts=False) -> pd.DataFrame:
    print("Calculating Pearson correlations.")
    COLUMN_NAME_TARGET = "target"
    COLUMN_NAME_WEIGHT = "importance"
    COLUMN_NAME_REGULATION = "regulation"
    COLUMN_NAME_CORRELATION = "rho"
    COLUMN_NAME_TF = "TF"
    ex_mtx = ex_mtx.T[~ex_mtx.columns.duplicated(keep='first')].T.astype(float)
    assert rho_threshold > 0, "rho_threshold should be greater than 0."
    if mask_dropouts:
        ex_mtx = ex_mtx.sort_index(axis=1)
        col_idx_pairs = _create_idx_pairs(adjacencies, ex_mtx)
        rhos = masked_rho4pairs(ex_mtx.values, col_idx_pairs, 0.0)
    else:
        genes = list(
            set(adjacencies[COLUMN_NAME_TF]).union(
                set(adjacencies[COLUMN_NAME_TARGET])))
        ex_mtx = ex_mtx[ex_mtx.columns[ex_mtx.columns.isin(genes)]]
        corr_mtx = pd.DataFrame(index=ex_mtx.columns,
                                columns=ex_mtx.columns,
                                data=np.corrcoef(ex_mtx.values.T))
        rhos = np.array([
            corr_mtx[s2][s1]
            for s1, s2 in zip(adjacencies.TF, adjacencies.target)
        ])

    regulations = (rhos > rho_threshold).astype(int) - (
        rhos < -rho_threshold).astype(int)
    return pd.DataFrame(
        data={
            COLUMN_NAME_TF: adjacencies[COLUMN_NAME_TF].values,
            COLUMN_NAME_TARGET: adjacencies[COLUMN_NAME_TARGET].values,
            COLUMN_NAME_WEIGHT: adjacencies[COLUMN_NAME_WEIGHT].values,
            COLUMN_NAME_REGULATION: regulations,
            COLUMN_NAME_CORRELATION: rhos
        })


def derive_regulons(
    motifs,
    db_names=('mm10__refseq-r80__10kb_up_and_down_tss.mc9nr',
              'mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr')):
    def contains(*elems):
        def f(context):
            return any(elem in context for elem in elems)

        return f

    # For the creation of regulons we only keep the 10-species databases and the activating modules. We also remove the
    # enriched motifs for the modules that were created using the method 'weight>50.0%' (because these modules are not part
    # of the default settings of modules_from_adjacencies anymore.
    motifs = motifs[
        np.fromiter(map(compose(op.not_, contains('weight>50.0%')), motifs.iloc[:,5]), dtype=np.bool) & \
        np.fromiter(map(contains(*db_names),motifs.iloc[:,5]), dtype=np.bool) & \
        np.fromiter(map(contains('activating'), motifs.iloc[:,5]), dtype=np.bool)]

    # We build regulons only using enriched motifs with a NES of 3.0 or higher; we take only directly annotated TFs or TF annotated
    # for an orthologous gene into account; and we only keep regulons with at least 10 genes.
    regulons = list(
        filter(
            lambda r: len(r) >= 10,
            df2regulons(motifs[
                (motifs.iloc[:, 1] >= 3.0)
                & ((motifs.iloc[:, 4] == 'gene is directly annotated')
                   |
                   (motifs.iloc[:, 4].str.startswith('gene is orthologous to')
                    & motifs.iloc[:, 4].str.endswith(
                        'which is directly annotated for motif')))])))

    # Rename regulons, i.e. remove suffix.
    return list(map(lambda r: r.rename(r.transcription_factor), regulons))



def add_scenic_metadata(adata: 'sc.AnnData',
                        auc_mtx: pd.DataFrame,
                        regulons,
                        bin_rep: bool = False,
                        copy: bool = False) -> 'sc.AnnData':
    """
    Add AUCell values and regulon metadata to AnnData object.
    :param adata: The AnnData object.
    :param auc_mtx: The dataframe containing the AUCell values (#observations x #regulons).
    :param bin_rep: Also add binarized version of AUCell values as separate representation. This representation
    is stored as `adata.obsm['X_aucell_bin']`.
    :param copy: Return a copy instead of writing to adata.
    :
    """
    # To avoid dependency with scanpy package the type hinting intentionally uses string literals.
    # In addition, the assert statement to assess runtime type is also commented out.
    #assert isinstance(adata, sc.AnnData)
    assert isinstance(auc_mtx, pd.DataFrame)
    assert len(auc_mtx) == adata.n_obs

    REGULON_SUFFIX_PATTERN = 'Regulon({})'

    result = adata.copy() if copy else adata

    # Add AUCell values as new representation (similar to a PCA). This facilitates the usage of
    # AUCell as initial dimensional reduction.
    result.obsm['X_aucell'] = auc_mtx.values.copy()
    if bin_rep:
        bin_mtx, _ = binarize(auc_mtx)
        result.obsm['X_aucell_bin'] = bin_mtx.values

    # Encode genes in regulons as "binary" membership matrix.
    if regulons is not None:
        genes = np.array(adata.var_names)
        data = np.zeros(shape=(adata.n_vars, len(regulons)), dtype=bool)
        for idx, regulon in enumerate(regulons):
            data[:, idx] = np.isin(genes, regulon.genes).astype(bool)
        regulon_assignment = pd.DataFrame(
            data=data,
            index=genes,
            columns=list(
                map(lambda r: REGULON_SUFFIX_PATTERN.format(r.name),
                    regulons)))
        result.var = pd.merge(result.var,
                              regulon_assignment,
                              left_index=True,
                              right_index=True,
                              how='left')

    # Add additional meta-data/information on the regulons.
    def fetch_logo(context):
        for elem in context:
            if elem.endswith('.png'):
                return elem
        return ""

    result.uns['aucell'] = {
        'regulon_names':
        auc_mtx.columns.map(lambda s: REGULON_SUFFIX_PATTERN.format(s)).values,
        'regulon_motifs':
        np.array([fetch_logo(reg.context)
                  for reg in regulons] if regulons is not None else [])
    }

    # Add the AUCell values also as annotations of observations. This way regulon activity can be
    # depicted on cellular scatterplots.
    mtx = auc_mtx.copy()
    mtx.columns = result.uns['aucell']['regulon_names']
    result.obs = pd.merge(result.obs,
                          mtx,
                          left_index=True,
                          right_index=True,
                          how='left')

    return result

