def umap_redu(pca):
    import umap
    import pandas as pd
    import numpy as np

    pca_py = pd.DataFrame(pca)   
    reducer = umap.UMAP(
        n_neighbors=30,
        min_dist=0.5,
        random_state=42)
    
    embedding = reducer.fit_transform(pca_py)
    return embedding
