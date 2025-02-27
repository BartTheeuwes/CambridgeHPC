#####################
## Define settings ##
#####################

# load default setings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
} else {
  stop("Computer not recognised")
}


## Options

opts$number_metacells <- 1000

# Dot size
opts$size.mapped <- 0.18
opts$size.nomapped <- 0.1

# Dot transparency
opts$alpha.mapped <- 0.65
opts$alpha.nomapped <- 0.35

## I/O
io$cell2metacell <- sprintf("%s/results/rna_atac/metacells/cell2metacell_%smetacells.txt.gz",io$basedir,opts$number_metacells)
# io$umap <- paste0(io$basedir,"/results/rna/dimensionality_reduction/all_cells/E7.5_rep1-E7.5_rep2-E8.0_rep1-E8.0_rep2-E8.5_rep1-E8.5_rep2_umap_features2500_pcs30_neigh25_dist0.3.txt.gz")
io$umap <- paste0(io$basedir,"/results//atac/archR/dimensionality_reduction/PeakMatrix/all_cells/E7.5_rep1-E7.5_rep2-E8.0_rep1-E8.0_rep2-E8.5_rep1-E8.5_rep2_umap_nfeatures50000_ndims50_neigh30_dist0.45.txt.gz")
io$outdir <- paste0(io$basedir,"/results/rna/metacells/pdf"); dir.create(io$outdir, showWarnings = F)



############################
## User-defined functions ##
############################


plot.dimred <- function(plot_df, query.label, atlas.label = "Atlas") {
  
  # Define dot size  
  size.values <- c(opts$size.mapped, opts$size.nomapped)
  names(size.values) <- c(query.label, atlas.label)
  
  # Define dot alpha  
  alpha.values <- c(opts$alpha.mapped, opts$alpha.nomapped)
  names(alpha.values) <- c(query.label, atlas.label)
  
  # Define dot colours  
  colour.values <- c("red", "lightgrey")
  names(colour.values) <- c(query.label, atlas.label)
  
  # Plot
  ggplot(plot_df, aes(x=V1, y=V2)) +
    ggrastr::geom_point_rast(aes(size=mapped, alpha=mapped, colour=mapped)) +
    scale_size_manual(values = size.values) +
    scale_alpha_manual(values = alpha.values) +
    scale_colour_manual(values = colour.values) +
    # labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
    guides(colour = guide_legend(override.aes = list(size=6))) +
    theme_classic() +
    theme(
      legend.position = "top", 
      legend.title = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank()
    )
}


###################
## Load metadata ##
###################

sample_metadata <- fread(io$metadata) %>%
  # .[pass_rnaQC==TRUE & doublet_call==FALSE] %>%
  .[pass_rnaQC==TRUE & pass_atacQC==TRUE & doublet_call==FALSE] %>%
  .[sample%in%opts$samples & celltype.predicted%in%opts$celltypes]

########################
## Load cell2metacell ##
########################

cell2metacell <- fread(io$cell2metacell) %>% 
  merge(sample_metadata[,c("cell","sample","stage")] %>% copy %>% setnames("cell","Metacell"))

length(unique(cell2metacell$Metacell))

tmp <-cell2metacell[,c("Metacell","sample","stage")] %>% unique
table(tmp$stage)
table(tmp$sample)

###########################
## Load precomputed UMAP ##
###########################

umap.dt <- fread(io$umap) %>% setnames(c("cell","V1","V2"))

#########################################################
## Plot dimensionality reduction: one sample at a time ##
#########################################################

to.plot <- umap.dt %>% copy %>%
  .[,index:=match(cell, cell2metacell$Metacell )] %>% 
  .[,mapped:=as.factor(!is.na(index))] %>% 
  .[,mapped:=plyr::mapvalues(mapped, from = c("FALSE","TRUE"), to = c("Atlas","Metacell"))] %>%
  setorder(mapped) 

p <- plot.dimred(to.plot, query.label = "Metacell", atlas.label = "Atlas")

pdf(sprintf("%s/umap_metacell.pdf",io$outdir), width=8, height=6.5)
print(p)
dev.off()
