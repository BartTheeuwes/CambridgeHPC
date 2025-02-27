# Create environment
#conda create --name mo python=3.8
#conda activate mo

# Install mamba
#conda install -c conda-forge mamba

# Install packages

mamba install -c bioconda r-seurat -y
mamba install -c bioconda bioconductor-genomicranges -y 
mamba install -c bioconda bioconductor-bsgenome -y
mamba install -c bioconda bioconductor-edger -y 
mamba install -c bioconda bioconductor-dropletutils -y 
mamba install -c bioconda bioconductor-biocparallel -y 
mamba install -c bioconda bioconductor-biomart -y
mamba install -c bioconda bioconductor-chromvar -y
mamba install -c bioconda bioconductor-rsamtools -y
mamba install -c bioconda bioconductor-biostrings -y
mamba install -c bioconda bioconductor-complexheatmap -y
mamba install -c bioconda bioconductor-genomicranges -y
mamba install -c bioconda bioconductor-batchelor -y
mamba install -c bioconda bioconductor-biocparallel -y
mamba install -c bioconda bioconductor-scater -y
mamba install -c bioconda bioconductor-genomicfeatures -y
mamba install -c bioconda bioconductor-annotationhub -y
mamba install -c conda-forge r-matrix -y
mamba install -c r r-irlba -y
mamba install -c conda-forge r-devtools -y
mamba install -c conda-forge r-uwot -y
mamba install -c bioconda r-harmony -y 
mamba install -c conda-forge r-ggplot2 -y
mamba install -c conda-forge r-matrixstats -y
mamba install -c conda-forge r-igraph -y
mamba install -c conda-forge r-reshape2 -y
mamba install -c conda-forge r-data.table -y
mamba install -c conda-forge r-knitr
mamba install -c bioconda bioconductor-scran -y
mamba install -c conda-forge r-reticulate -y
mamba install -c r r-gridextra -y
mamba install -c conda-forge r-rtsne -y
mamba install -c conda-forge r-viridis -y
mamba install -c conda-forge jupyterlab -y
mamba install -c conda-forge r-irkernel -y
mamba install -c conda-forge r-base=4.1.1