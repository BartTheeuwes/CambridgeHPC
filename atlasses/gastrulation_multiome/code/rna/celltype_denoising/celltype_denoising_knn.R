#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/rna/celltype_denoising/utils.R")
} else {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/rna/celltype_denoising/utils.R")
}

io$metadata <- paste0(io$basedir,"/results/rna/doublets/sample_metadata_after_doublets.txt.gz")
io$metadata.out <- paste0(io$basedir,"/results/rna/celltype_denoising/sample_metadata_after_celltype_denoising.txt.gz")
io$outdir <- paste0(io$basedir,"/results/rna/celltype_denoising")
opts$k <- 10

opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
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

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_QC==TRUE & doublet_call==FALSE] %>%
  .[celltype.mapped%in%opts$celltypes & sample%in%opts$samples]
table(sample_metadata$celltype.mapped)

################################################
## Load precomputed dimensionality reduction ##
################################################

# dimred.dt <- fread("/Users/ricard/data/gastrulation_multiome_10x/results/rna/dimensionality_reduction/umap.txt.gz") %>%
#   merge(sample_metadata[,c("cell","sample","celltype.mapped")],by="cell")

##############################################################
## Construct a similarity matrix using the UMAP coordinates ##
##############################################################

# TO-DO: TAKE MORE THAN 2 UMAP COORDINATES

# Euclidean similarity 
S <- 1 / as.matrix(dist(dimred.dt[,c("UMAP1","UMAP2")]))
rownames(S) <- colnames(S) <- dimred.dt$cell
diag(S) <- 0

################
## Define kNN ##
################

# get closest k cells
knn.matrix <- apply(S, 1, function(x) which.maxn(x,n=opts$k))
k.mapped  <- apply(knn.matrix, 2, function(x) dimred.dt$cell[x])

# get celltypes for the k nearest cells
celltypes <- apply(k.mapped, 2, function(x) dimred.dt$celltype.mapped[match(x, dimred.dt$cell)])

# For each cell calculate the most dominant cell type
celltype.denoised <- apply(celltypes, 2, function(x) getmode(x, 1:length(x)))

dimred.dt$celltype.denoised <- celltype.denoised

###########################
## Plot before denoising ##
###########################

p <- ggplot(dimred.dt, aes(x=UMAP1, y=UMAP2, fill=celltype.mapped)) +
  geom_point(size=1.5, shape=21, stroke=0.1) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position="none",
    legend.title=element_blank()
  )

pdf(paste0(io$outdir,"/umap_before_denoising.pdf"), width=6, height=4)
print(p)
dev.off()


###########################
## Plot after denoising ##
###########################

p <- ggplot(dimred.dt, aes(x=UMAP1, y=UMAP2, fill=celltype.denoised)) +
  geom_point(size=1.5, shape=21, stroke=0.1) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position="none",
    legend.title=element_blank()
  )

pdf(sprintf("%s/umap_after_denoising_k%s.pdf",io$outdir,opts$k), width=6, height=4)
print(p)
dev.off()

############################
## Update sample metadata ##
############################

sample_metadata.updated <- fread(io$metadata) %>%
  merge(dimred.dt[,c("cell","celltype.denoised")],by=c("cell"), all.x=T) %>%
  .[is.na(celltype.denoised),celltype.denoised:=celltype.mapped]
fwrite(sample_metadata.updated, io$metadata.out, sep="\t", na="NA", quote=F)
