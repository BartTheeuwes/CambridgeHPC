#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/utils.R")
} else if (grepl("Workstation",Sys.info()['nodename'])){
  source("/home/lijingyu/gastrulation/gastrulation_multiome_10x/settings.R")
  source("/home/lijingyu/gastrulation/gastrulation_multiome_10x/utils.R")
} else {
  stop("Computer not recognised")
}
  

# I/O
io$outdir <- paste0(io$basedir,"/results/rna_atac/GRN/coexpression")
io$rna.metacell.sce <-paste0(io$basedir,'/processed/rna/metacell/SingleCellExperiment_1000metacells.rds')
  
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

opts$motif_annotation <- "Motif_cisbp"
opts$remove.ExE.celltypes=TRUE


###################
## Load metadata ##
###################

sample_metadata <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE] %>%
  .[celltype.mapped%in%opts$celltypes] 
if (opts$remove.ExE.celltypes) {
  sample_metadata <- sample_metadata %>%
    .[!celltype.mapped%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
  opts$celltypes=opts$celltypes[!opts$celltypes%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}

#########################
## Load pseudobulk RNA ##
#########################

# Load SingleCellExperiment
rna.sce <- readRDS(io$rna.metacell.sce)
rna.sce <-rna.sce%>%.[,colData(.)$celltype.mapped %in% opts$celltypes]

# Rename genes
# rownames(rna.sce) <- toupper(rownames(rna.sce))

###########################
## Load motif annotation ##
###########################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/atac/archR/load_motif_annotation.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/atac/archR/load_motif_annotation.R")
} else if (grepl('Workstation',Sys.info()['nodename'])){
  source("/home/lijingyu/gastrulation/gastrulation_multiome_10x/atac/archR/load_motif_annotation.R")
} else{
  stop("Computer not recognised")
}

rna.sce.tf <- rna.sce[toupper(rownames(rna.sce))%in%motif2gene.dt$gene,]
rna.sce.target <- rna.sce
rownames(rna.sce.tf) <- toupper(rownames(rna.sce.tf))

##########################
## Correlation analysis ##
##########################

cor.mtx <- cor(t(logcounts(rna.sce.tf)),t(logcounts(rna.sce.target))) %>% round(2)

# Save
saveRDS(cor.mtx, paste0(io$outdir,"/correlation_matrix_tf2gene.rds"))

##########
## Plot ##
##########

cor.mtx <- readRDS(paste0(io$outdir,"/correlation_matrix_tf2gene.rds"))

i <- "HNF4A"
j <- "Afp"

to.plot <- data.table(
  TF = logcounts(rna.sce.tf[i,])[1,],
  target_gene = logcounts(rna.sce.target[j,])[1,],
  celltype=rna.sce.tf@colData$celltype.mapped
)


ggscatter(to.plot, x="TF", y="target_gene", fill="celltype", size=4, shape=21, 
               add="reg.line", add.params = list(color="black", fill="lightgray"), conf.int=TRUE) +
  stat_cor(method = "pearson") +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(x=sprintf("%s expression",i), y=sprintf("%s expression",j)) +
  guides(fill=F) +
  theme(
    plot.title = element_text(hjust = 0.5, size=rel(0.85)),
    axis.text = element_text(size=rel(0.7))
  )
