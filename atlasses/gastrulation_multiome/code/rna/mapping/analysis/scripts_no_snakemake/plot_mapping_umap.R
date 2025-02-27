#####################
## Define settings ##
#####################

# load default setings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/rna/mapping/analysis/plot_utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/rna/mapping/analysis/plot_utils.R")
} else {
  stop("Computer not recognised")
}

# I/O
io$outdir <- paste0(io$basedir,"/results/rna/mapping/pdf")
dir.create(io$outdir, showWarnings = F)
dir.create(paste0(io$outdir,"/per_sample"), showWarnings = F)
dir.create(paste0(io$outdir,"/per_stage"), showWarnings = F)

# Options

# Dot size
opts$size.mapped <- 0.18
opts$size.nomapped <- 0.1

# Transparency
opts$alpha.mapped <- 0.65
opts$alpha.nomapped <- 0.35

###################
## Load metadata ##
###################

sample_metadata <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE] %>%
  .[sample%in%opts$samples & celltype.predicted%in%opts$celltypes]

####################
## Load 10x atlas ##
####################

# Load atlas cell metadata
meta_atlas <- fread(io$rna.atlas.metadata) %>%
  .[celltype%in%opts$celltypes] %>%
  .[stripped==F & doublet==F]

# Extract precomputed dimensionality reduction coordinates
umap.dt <- meta_atlas[,c("cell","umapX","umapY","celltype")] %>%
  setnames(c("umapX","umapY"),c("V1","V2"))

####################
## Plot all cells ##
####################

to.plot <- umap.dt %>% copy %>%
  .[,index:=match(cell, sample_metadata[,closest.cell] )] %>% 
  .[,mapped:=as.factor(!is.na(index))] %>% 
  .[,mapped:=plyr::mapvalues(mapped, from = c("FALSE","TRUE"), to = c("scRNA-seq atlas","Multiome"))] %>%
  setorder(mapped) 

p <- plot.dimred(to.plot, query.label = "Multiome", atlas.label = "scRNA-seq atlas") +
  theme(
    legend.position = "none",
    axis.line = element_blank()
  )

pdf(sprintf("%s/umap_mapped_allcells.pdf",io$outdir), width=8, height=6.5)
print(p)
dev.off()

###############################
## Plot one sample at a time ##
###############################

for (i in opts$samples) {
  
  to.plot <- umap.dt %>% copy %>%
    .[,index:=match(cell, sample_metadata[sample==i,closest.cell] )] %>% 
    .[,mapped:=as.factor(!is.na(index))] %>% 
    .[,mapped:=plyr::mapvalues(mapped, from = c("FALSE","TRUE"), to = c("Atlas",i))] %>%
    setorder(mapped) 
  
  p <- plot.dimred(to.plot, query.label = i, atlas.label = "Atlas")
  
  pdf(sprintf("%s/per_sample/umap_mapped_%s.pdf",io$outdir,i), width=8, height=6.5)
  print(p)
  dev.off()
}

##############################
## Plot one stage at a time ##
##############################

for (i in opts$stages) {
  
  to.plot <- umap.dt %>% copy %>%
    .[,index:=match(cell, sample_metadata[stage==i,closest.cell] )] %>% 
    .[,mapped:=as.factor(!is.na(index))] %>% 
    .[,mapped:=plyr::mapvalues(mapped, from = c("FALSE","TRUE"), to = c("Atlas",i))] %>%
    setorder(mapped) 
  
  p <- plot.dimred(to.plot, query.label = i, atlas.label = "Atlas")
  
  pdf(sprintf("%s/per_stage/umap_mapped_%s.pdf",io$outdir,i), width=8, height=6.5)
  print(p)
  dev.off()
}
