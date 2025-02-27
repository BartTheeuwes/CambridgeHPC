
umap_atac.df <- fread(paste0(basedir,"/umap_atac.txt.gz")) %>% setnames(c("cell","X","Y")) %>% tibble::column_to_rownames("cell") 
umap_rna.df <- fread(paste0(basedir,"/umap_rna.txt.gz")) %>% setnames(c("X","Y","cell")) %>% tibble::column_to_rownames("cell") 
TFs <- fread(paste0(basedir,"/TFs.txt"), header=F)[[1]]
celltypes <- fread(paste0(basedir,"/celltypes.txt"), header=F)[[1]]
genes <- fread(paste0(basedir,"/genes.txt"), header=F)[[1]]
paga <- readRDS(paste0(basedir,"/paga_network.rds"))
motifs <- readRDS(paste0(basedir,"/PWMatrixList.rds"))

sample_metadata <- fread(paste0(basedir,"/sample_metadata.txt.gz")) %>% 
  setnames("celltype.predicted","celltype") %>% 
  .[,celltype := factor(celltype, levels = names(celltype_colours), ordered = TRUE)] %>%
  .[,stage := factor(stage, levels = names(stage_colours), ordered = TRUE)] %>%
  # .[,sample = factor(sample, levels = names(samples), ordered = TRUE)
  setkey(cell)

genes <- fread(paste0(basedir,"/genes.txt"), header=F)[[1]]
cells_atac <- fread(paste0(basedir,"/cells_atac.txt"), header=F)[[1]]
cells_rna <- fread(paste0(basedir,"/cells_rna.txt"), header=F)[[1]]

# Load ATAC GeneScoreMatrix
link_gene_acc = HDF5Array(file = paste0(basedir,"/atac_GeneScoreMatrix.hdf5"), name = "atac_GeneScoreMatrix")
colnames(link_gene_acc) <- cells_atac
rownames(link_gene_acc) <- genes

# Load RNA expression matrix
link_rna_expr = HDF5Array(file = paste0(basedir,"/rna_expr.hdf5"), name = "rna_expr_logcounts")
colnames(link_rna_expr) <- cells_rna
rownames(link_rna_expr) <- genes

# Load RNA vs ATAC estimates
rna_vs_chromvar_pseudobulk.dt <- fread(paste0(basedir,"/rna_vs_chromvar_pseudobulk.txt.gz")) %>%
  .[,celltype := factor(celltype, levels = celltypes, ordered = TRUE)]
stopifnot(unique(rna_vs_chromvar_pseudobulk.dt$celltype)%in%celltypes)

gene_rna_vs_atac_pseudobulk.dt <- fread(paste0(basedir,"/gene_rna_vs_acc_pseudobulk.txt.gz")) %>%
  .[,celltype := factor(celltype, levels = celltypes, ordered = TRUE)]
stopifnot(unique(gene_rna_vs_atac_pseudobulk.dt$celltype)%in%celltypes)


cor_rna_vs_chromvar_per_gene.dt <- fread(paste0(basedir,"/cor_rna_vs_chromvar_per_gene.txt.gz")) %>%
  .[,log_pval:=-log10(padj_fdr+1e-100)] %>%
  .[,cor_sign:=as.factor(c("Repressor","Activator")[(r>0)+1])]

# in silico chip-seq results
insilico_chip_stats <- fread(paste0(basedir,"/insilico_chipseq/insilico_chip_stats.txt.gz"))



