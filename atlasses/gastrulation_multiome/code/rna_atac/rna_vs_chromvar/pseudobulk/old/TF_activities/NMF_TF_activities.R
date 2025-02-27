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
  # "Caudal_epiblast",
  # "PGC",
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
  # "Mesenchyme",
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
  # "Surface_ectoderm",
  "Visceral_endoderm",
  "ExE_endoderm",
  "ExE_ectoderm",
  "Parietal_endoderm"
)


###################################
## Load pseudobulk TF activities ##
###################################

TF_activites.dt <- fread("/Users/ricard/data/gastrulation_multiome_10x/results/rna_atac/rna_vs_chromvar/pseudobulk/TF_activities/tf_activities.txt.gz")

TF_activites.mtx <- TF_activites.dt %>% 
  .[celltype%in%opts$celltypes] %>%
  dcast(celltype~gene, value.var="activity") %>% 
  matrix.please %>% t

##########
## NMF ##
##########

library(NMF)

# nmfAlgorithm('brunet')
# nmfSeed('ica')

foo <- nmf(TF_activites.mtx, rank=5, method="brunet", seed="ica")

# Extract weights and factors
weights <- foo@fit@W
factors <- foo@fit@H %>% t

# Sort factors by variance explained
var.explained <- sapply(1:nrow(foo@fit@H), function(i) {
  sum((TF_activites.mtx - foo@fit@W[,i,drop=F] %*% foo@fit@H[i,,drop=F])**2)
})

factors <- factors[,order(var.explained, decreasing = T)]
weights <- weights[,order(var.explained, decreasing = T)]

##########
## Plot ##
##########

to.plot <- factors %>% 
  as.data.table(keep.rownames = T) %>%
  setnames("rn","sample")

ggscatter(to.plot, x="V1", y="V2", fill="sample", shape=21, size=5) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme(
    legend.position="none"
  )


opts$celltypes.to.rename <- c(
  "Erythroid1" = "Erythroid",
  "Erythroid2" = "Erythroid",
  "Erythroid3" = "Erythroid",
  "Blood_progenitors_1" = "Blood_progenitors_2",
  # "Parietal_endoderm" = "ExE_endoderm",
  # "Visceral_endoderm" = "ExE_endoderm",
  "Allantois" = "ExE_mesoderm",
  "Nascent_mesoderm" = "Mixed_mesoderm",
  "Paraxial_mesoderm" = "Somitic_mesoderm"
)
# opts$celltype.colors <- opts$celltype.colors[names(opts$celltype.colors)!="Mid_Hindbrain"]
# names(opts$celltype.colors) <- stringr::str_replace_all(names(opts$celltype.colors),opts$celltypes.to.rename)


to.plot.barplot <- to.plot %>% 
  melt(id.vars=c("sample"), variable.name="factor") %>%
  .[,sample:=stringr::str_replace_all(sample,opts$celltypes.to.rename)] %>%
  .[,.(value=sum(value)),by=c("factor","sample")]
  

ggbarplot(to.plot.barplot, x="factor", y="value", fill="sample") +
  scale_fill_manual(values=opts$celltype.colors) +
  theme(
    legend.position="none"
  )




#########
## PCA ##
#########

# Run PCA
pca.activity <- irlba::prcomp_irlba(activity.mtx, n=10)


# Lineplot of Var % values
to.plot <- data.table(
  chromvar_var = pca.chromvar$sdev**2 / pca.chromvar$totalvar,
  rna_var = pca.rna$sdev**2 / pca.rna$totalvar,
  activity_var = pca.activity$sdev**2 / pca.activity$totalvar
) %>% .[,pc:=paste0("PC",1:.N)] %>%
  melt(id.vars=c("pc"), variable.name="class")

ggline(to.plot, x="pc", y="value", color="class")

# Plot PCA

# to.plot <- pca.activity$x %>% as.data.table %>% 
# to.plot <- pca.rna$x %>% as.data.table %>% 
to.plot <- pca.chromvar$x %>% as.data.table %>% 
  .[,sample:=rownames(activity.mtx)] 

ggscatter(to.plot, x="PC1", y="PC2", fill="sample", shape=21, size=5) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme(
    legend.position="none"
  )
