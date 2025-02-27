########################
## Load ArchR Project ##
########################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
} else {
  stop("Computer not recognised")
}

#####################
## Define settings ##
#####################

# I/O
io$pca.rna <- paste0(io$basedir,"/results/rna/dimensionality_reduction/all_cells/E7.5_rep1-E7.5_rep2-E8.0_rep1-E8.0_rep2-E8.5_rep1-E8.5_rep2_pca_features2500_pcs30_batchcorrectionbysample.txt.gz")
io$pca.atac <- paste0(io$basedir,"/results/atac/archR/dimensionality_reduction/PeakMatrix/all_cells/E7.5_rep1-E7.5_rep2-E8.0_rep1-E8.0_rep2-E8.5_rep1-E8.5_rep2_lsi_features50000_ndims30.txt.gz")
io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/TF_activites")


# Options
opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E8.0_rep1",
  "E8.0_rep2",
  "E8.5_rep1",
  "E8.5_rep2"
)

opts$celltypes = c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  "PGC",
  "Anterior_Primitive_Streak",
  "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  "Mixed_mesoderm",
  "Intermediate_mesoderm",
  "Caudal_Mesoderm",
  "Paraxial_mesoderm",
  "Somitic_mesoderm",
  "Pharyngeal_mesoderm",
  "Cardiomyocytes",
  "Allantois",
  "ExE_mesoderm",
  "Mesenchyme",
  "Haematoendothelial_progenitors",
  "Endothelium",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
  "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm",
  "Visceral_endoderm",
  "ExE_endoderm",
  "ExE_ectoderm",
  "Parietal_endoderm"
)

# opts$motif_annotation <- "Motif_JASPAR2020_human"
opts$motif_annotation <- "Motif_cisbp"

opts$denoise <- TRUE
opts$knn <- 25

###################
## Load metadata ##
###################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & pass_rnaQC==TRUE & doublet_call==FALSE] %>%
  .[sample%in%opts$samples & celltype.predicted%in%opts$celltypes] %>%
  .[,celltype.predicted:=factor(celltype.predicted,levels=opts$celltypes)] 

##################
## Subset ArchR ##
##################

stopifnot(sample_metadata$cell %in% rownames(ArchRProject))
ArchRProject.filt <- ArchRProject[sample_metadata$cell,]

################################
## Load RNA and chromVAR data ##
################################

io$rna_chromvar.file <- sprintf("%s/chromvar_rna.txt.gz",io$outdir)

if (file.exists(io$rna_chromvar.file)) {
  
  sprintf("Found precomputed values, loading %s...",io$rna_chromvar.file)
  chromvar_rna_dt <- fread(io$rna_chromvar.file)
  
} else {
  
  if (grepl("ricard",Sys.info()['nodename'])) {
    source("/Users/ricard/gastrulation_multiome_10x/rna_atac/rna_vs_chromvar/load_rna_chromvar_single_cells.R")
  } else if (grepl("ebi",Sys.info()['nodename'])) {
    source("/homes/ricard/gastrulation_multiome_10x/rna_atac/rna_vs_chromvar/load_rna_chromvar_single_cells.R")
  } else {
    stop("Computer not recognised")
  }
  
  # Filter genes with low variability
  genes.to.keep.rna <- rna_dt[,.(var(expr)),by="gene"] %>% .[V1>0.001,gene] %>% as.character
  genes.to.keep.chromvar <- chromvar_dt[,.(var(chromvar_zscore)),by="gene"] %>% .[V1>0.10,gene] %>% as.character
  genes.to.keep <- intersect(genes.to.keep.rna,genes.to.keep.chromvar)
  rna_dt <- rna_dt[gene%in%genes.to.keep]
  chromvar_dt <- chromvar_dt[gene%in%genes.to.keep]
  
  # Merge
  chromvar_rna_dt <- merge(
    rna_dt,
    chromvar_dt,
    by = c("cell","gene")
  ) %>% .[,c("expr","chromvar_zscore"):=list(round(expr,2),round(chromvar_zscore,2))]
  
  # Calculate TF activities
  chromvar_rna_dt %>%
    .[,activity:=minmax.normalisation(expr)*minmax.normalisation(chromvar_zscore)] %>%
    .[,activity:=round(minmax.normalisation(activity),2)]
  
  length(unique(chromvar_rna_dt$gene))
  length(unique(chromvar_rna_dt$cell))
  
  fwrite(chromvar_rna_dt, io$rna_chromvar.file, quote=F, sep="\t", na="NA")
}

# Add sample metadata
chromvar_rna_dt <- chromvar_rna_dt %>%
  merge(sample_metadata[,c("cell","celltype.predicted")], by="cell")

##############################
## Dimensionality reduction ##
##############################

# Run PCA
activity.mtx <- chromvar_rna_dt %>% 
  dcast(cell~gene,value.var="activity") %>%
  matrix.please
pca.activity <- irlba::prcomp_irlba(activity.mtx, n=50)

rna.mtx <- chromvar_rna_dt %>% 
  dcast(cell~gene,value.var="expr") %>%
  matrix.please
pca.rna <- irlba::prcomp_irlba(rna.mtx, n=50)

chromvar.mtx <- chromvar_rna_dt %>% 
  dcast(cell~gene,value.var="chromvar_zscore") %>%
  matrix.please
pca.chromvar <- irlba::prcomp_irlba(chromvar.mtx, n=50)

# UMAP
umap.activity <- uwot::umap(pca.activity$x, n_neighbors=25, min_dist=0.3, metric="cosine")
umap.rna <- uwot::umap(pca.rna$x, n_neighbors=25, min_dist=0.3, metric="cosine")
umap.chromvar <- uwot::umap(pca.chromvar$x, n_neighbors=25, min_dist=0.3, metric="cosine")

umap_list <- list(
  "activity" = umap.activity, 
  "rna" = umap.rna, 
  "chromvar" = umap.chromvar
)

# Plot PCA
to.plot <- names(umap_list) %>% 
  map(function(i) umap_list[[i]] %>% as.data.table %>% .[,cell:=rownames(activity.mtx)] %>% .[,class:=i]) %>%
  rbindlist %>%
    merge(sample_metadata[,c("cell","celltype.predicted")])
  
ggplot(to.plot, aes_string(x="V1", y="V2", fill="celltype.predicted")) +
  facet_wrap(~class,nrow=1) +
  geom_point(size=1.5, shape=21, stroke=0.05) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme_classic() +
  ggplot_theme_NoAxes() +
  theme(
    legend.position = "none"
  )

# Overlay TF activities onto the UMAP
io$umap <- paste0(io$basedir,"/results/rna/dimensionality_reduction/all_cells_noExE/E7.5_rep1-E7.5_rep2-E8.0_rep1-E8.0_rep2-E8.5_rep1-E8.5_rep2_umap_features2500_pcs30_neigh25_dist0.3.txt.gz")
umap.dt <- fread(io$umap)

genes.to.plot <- unique(chromvar_rna_dt$gene)
# genes.to.plot <- cor.dt[sig==T & abs(r)>0.25,gene]

for (i in genes.to.plot) {
  
  outfile <- sprintf("%s/per_gene/%s_activity_allcells.png",io$outdir,i)
  
  if (file.exists(outfile)) {
    print(sprintf("file or %s already exists...",i))
  } else {
    
    # Plot UMAP
    to.plot <- chromvar_rna_dt[gene==i] %>%
      setnames(c("cell","gene","RNA expression", "Motif accessibility (chromVAR)", "TF activity", "celltype.predicted")) %>%
      melt(id.vars=c("cell","gene","celltype.predicted")) %>% 
      merge(umap.dt,by="cell")
    
    to.plot.umap <- to.plot %>% copy %>% .[,value:=minmax.normalisation(value),by="variable"]
    p1 <- ggplot(to.plot.umap, aes_string(x="UMAP1", y="UMAP2", fill="value")) +
      facet_wrap(~variable, nrow=1) +
      geom_point(size=1.5, shape=21, stroke=0.05) +
      scale_fill_gradient2(low = "gray50", mid="gray90", high = "red") +
      theme_classic() +
      ggplot_theme_NoAxes() +
      theme(
        legend.position = "none"
      )
    
    
    # Plot boxplot
    # to.plot.boxplot <- to.plot
    p2 <- ggboxplot(to.plot, x="celltype.predicted", y="value", fill="celltype.predicted", outlier.shape=NA) +
      facet_wrap(~variable, nrow=3, scales="free_y") +
      scale_fill_manual(values=opts$celltype.colors) +
      geom_hline(yintercept=0, linetype="dashed") +
      labs(x="", y="") +
      theme_classic() +
      guides(x = guide_axis(angle = 90)) +
      theme(
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()
      )
    
    p <- cowplot::plot_grid(plotlist=list(p1,p2), nrow = 1, rel_widths = c(1/2,1/2))
    
    # pdf(sprintf("%s/individual_genes/%s_rna_chromvar_vs_pseudotime_blood.pdf",args$outdir,i), width=8, height=6)
    png(outfile, width = 1200, height = 350)
    print(p)
    dev.off()
    
  }
}
