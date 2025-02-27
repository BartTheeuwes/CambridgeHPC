#####################
## Define settings ##
#####################

# Load default settings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/utils.R")
} else {
  stop("Computer not recognised")
}

# I/O
io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/pseudobulk/TF_activites")

# Options
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

opts$motif_annotation <- "Motif_cisbp" # Motif_JASPAR2020_human

###################
## Load metadata ##
###################

# sample_metadata <- fread(io$metadata) %>%
#   .[pass_atacQC==TRUE & pass_rnaQC==TRUE & doublet_call==FALSE] %>%
#   .[celltype.mapped%in%opts$celltypes] %>%
#   .[,celltype.mapped:=factor(celltype.mapped,levels=opts$celltypes)] 

# subset celltypes with sufficient number of cells
# opts$min.cells <- 50
# sample_metadata <- sample_metadata %>%
#   .[,N:=.N,by=c("celltype.mapped")] %>% .[N>opts$min.cells] %>% .[,N:=NULL] %>% droplevels
# table(sample_metadata$celltype.mapped)

# opts$celltypes <- unique(sample_metadata$celltype.mapped) %>% as.character

################################################
## Load pseudobulk RNA and chromVAR estimates ##
################################################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/rna_atac/rna_vs_chromvar/pseudobulk/load_rna_chromvar_pseudobulk.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/rna_atac/rna_vs_chromvar/pseudobulk/load_rna_chromvar_pseudobulk.R")
} else {
  stop("Computer not recognised")
}

################
## Parse data ##
################

# Filter genes with low variability
genes.to.keep.rna <- rna_dt[,.(var(expr)),by="gene"] %>% .[V1>0.05,gene] %>% as.character
genes.to.keep.chromvar <- chromvar_dt[,.(var(chromvar_zscore)),by="gene"] %>% .[V1>0.05,gene] %>% as.character
genes.to.keep <- intersect(genes.to.keep.rna,genes.to.keep.chromvar)
rna_dt <- rna_dt[gene%in%genes.to.keep]
chromvar_dt <- chromvar_dt[gene%in%genes.to.keep]

# minmax scaling globally
# rna_dt[,expr:=minmax.normalisation(expr)]
# chromvar_dt[,chromvar_zscore:=minmax.normalisation(chromvar_zscore)]

# Minmax scaling by gene
rna_dt[,expr:=minmax.normalisation(expr),by="gene"]
chromvar_dt[,chromvar_zscore:=minmax.normalisation(chromvar_zscore),by="gene"]

###########
## Merge ##
###########

rna_chromvar.dt <- merge(
  rna_dt,
  chromvar_dt,
  by = c("celltype","gene")
)

#############################
## Calculate TF activities ##
#############################

rna_chromvar.dt %>%
  .[,activity:=expr*chromvar_zscore] %>%
  .[,activity:=minmax.normalisation(activity)]

# Save
fwrite(rna_chromvar.dt, paste0(io$outdir,"/tf_activities.txt.gz"), sep="\t", quote=F, na="NA")

#######################
## Barplots per gene ##
#######################

genes.to.plot <- unique(rna_chromvar.dt$gene)
# genes.to.plot <- cor.dt[sig==T & abs(r)>0.25,gene]

for (i in genes.to.plot) {
  
  outfile <- sprintf("%s/per_gene/%s_activity_pseudobulk.pdf",io$outdir,i)
  
  if (file.exists(outfile)) {
    print(sprintf("file or %s already exists...",i))
  } else {
    
    to.plot <- rna_chromvar.dt[gene==i] %>%
      melt(id.vars=c("celltype","gene")) %>%
      .[,celltype:=factor(celltype,levels=levels(sample_metadata$celltype.mapped))]
    
    p <- ggbarplot(to.plot, x="celltype", y="value", fill="celltype") +
      facet_wrap(~variable, nrow=3)  +
      scale_fill_manual(values=opts$celltype.colors) +
      guides(x = guide_axis(angle = 90)) +
      labs(x="", y="TF activity") +
      guides(fill=F) +
      theme(
        axis.text = element_text(size=rel(0.7))
      )
    
    pdf(outfile, width = 6, height = 8)
    print(p)
    dev.off()
  }
}

#############################
## Line plot per cell type ##
#############################

celltypes.to.plot <- unique(rna_chromvar.dt$celltype)
# genes.to.plot <- cor.dt[sig==T & abs(r)>0.25,gene]

for (i in celltypes.to.plot) {
  
  outfile <- sprintf("%s/per_celltype/%s_activity_pseudobulk.pdf",io$outdir,i)
  
  if (file.exists(outfile)) {
    print(sprintf("file or %s already exists...",i))
  } else {
    
    to.plot <- rna_chromvar.dt[celltype==i] %>%
      setorder(-activity) %>%
      head(n=50) %>%
      .[,gene:=factor(gene,levels=rev(gene))]
    
    p <- ggplot(to.plot, aes_string(x="gene", y="activity"), fill="gray70") +
      geom_point(size=2) +
      geom_segment(aes_string(xend="gene"), size=0.75, yend=0) +
      coord_flip(ylim = c(0.5,1.0)) +
      labs(x="", "TF activity") +
      theme_classic() +
      theme(
        axis.ticks.y = element_blank(),
        axis.text = element_text(size=rel(0.75), color="black")
      )
    
    pdf(outfile, width = 5, height = 8)
    print(p)
    dev.off()
  }
}

