here::i_am("rna/trajectories/infer_trajectory.R")

suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(destiny))
suppressPackageStartupMessages(library(reticulate))


################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--metadata',  type="character",              help='') 
p$add_argument('--sce',  type="character",              help='') 
p$add_argument('--trajectory_name',  type="character",              help='') 
p$add_argument('--celltype_label',  type="character",              help='') 
# p$add_argument('--genes_to_plot',  type="character", nargs="+", help='') 
p$add_argument('--outdir',          type="character",                help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

## START TEST ##
args <- list()
args$sce <-file.path(io$basedir,"processed/rna/SingleCellExperiment.rds")
# args$metadata <- file.path(io$basedir,"results/atac/archR/celltype_assignment/sample_metadata_after_celltype_assignment.txt.gz")
args$metadata <- file.path(io$basedir,"results/rna/mapping/sample_metadata_after_mapping.txt.gz")
args$trajectory_name <- "blood"
args$celltype_label <- "celltype.mapped"
args$outdir <- file.path(io$basedir,"results/rna/trajectories/blood")
## END TEST ##

# I/O
dir.create(args$outdir, showWarnings=F)

# Options
opts$celltype_trajectory_dic <- list(
  "blood" = c("Haematoendothelial_progenitors", "Blood_progenitors_1", "Blood_progenitors_2", "Erythroid1", "Erythroid2", "Erythroid3"),
  "ectoderm" = c("Epiblast", "Rostral_neurectoderm", "Forebrain_Midbrain_Hindbrain"),
  "endoderm" = c("Epiblast", "Anterior_Primitive_Streak", "Def._endoderm", "Gut"),
  "mesoderm" = c("Epiblast", "Primitive_Streak", "Nascent_mesoderm"),
  "nmp" = c("NMP","Caudal_Mesoderm", "Somitic_mesoderm", "Spinal_cord")
  )

stopifnot(args$trajectory_name%in%names(opts$celltype_trajectory_dic))
opts$celltypes <- opts$celltype_trajectory_dic[[args$trajectory_name]]

opts$genes2plot <- list(
  "blood" = c("Hbb-y","Etv2"),
  "ectoderm" = c("Utf1","Crabp1"),
  "endoderm" = c("Pou5f1","Krt8"),
  "mesoderm" = c("Dnmt3b","Mesp1"),
  "nmp" = c("Sox2","T")
)

stopifnot(args$trajectory_name%in%names(opts$genes2plot))
opts$genes_to_plot <- opts$genes2plot[[args$trajectory_name]]

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE]

stopifnot(args$celltype_label%in%colnames(sample_metadata))
sample_metadata <- sample_metadata   %>%
  .[,celltype:=eval(as.name(args$celltype_label))] %>%
  .[celltype%in%opts$celltypes] %>%
  .[,celltype:=factor(celltype,levels=opts$celltypes)]

table(sample_metadata$celltype)

# Save
fwrite(sample_metadata, file.path(args$outdir,sprintf("%s_sample_metadata.txt.gz",args$trajectory_name)))

#########################
## Load RNA expression ##
#########################
  
sce <- load_SingleCellExperiment(
  file = args$sce, 
  normalise = TRUE, 
  cells = sample_metadata$cell, 
  remove_non_expressed_genes = TRUE
)

colData(sce) <- sample_metadata %>% as.data.frame %>% tibble::column_to_rownames("cell") %>%
  .[colnames(sce),] %>% DataFrame()

#######################
## Feature selection ##
#######################

decomp <- modelGeneVar(sce)
decomp <- decomp[decomp$mean > 0.01,]
hvgs <- rownames(decomp)[decomp$p.value <= 0.01]

# Subset HVGs
sce_filt <- sce[hvgs,]
dim(sce_filt)

#########
## PCA ##
#########

sce_filt <- runPCA(sce_filt, ncomponents = 2, ntop=nrow(sce_filt))

# Define the starting point
foo <- mean(reducedDim(sce_filt,"PCA")[sce_filt$celltype==opts$celltypes[1],"PC1"])
if (foo>0) reducedDim(sce_filt,"PCA")[,"PC1"] <- -reducedDim(sce_filt,"PCA")[,"PC1"]


# Plot PCA
pdf(file.path(args$outdir,sprintf("%s_pca_celltype.pdf",args$trajectory_name)), width=8, height=5)
plotPCA(sce_filt, colour_by="celltype", ncomponents = c(1,2)) +
  scale_colour_manual(values=opts$celltype.colors) + theme(legend.position="none")
dev.off()

pdf(file.path(args$outdir,sprintf("%s_pca_stage.pdf",args$trajectory_name)), width=8, height=5)
plotPCA(sce_filt, colour_by="stage", ncomponents = c(1,2)) +
  scale_colour_manual(values=opts$stage.colors) + theme(legend.position="top", legend.title = element_blank())
dev.off()

###################
## Diffusion map ##
###################

set.seed(42)
dm <- DiffusionMap(sce_filt, n_pcs=2)

# Add to the SingleCellExperiment object
reducedDim(sce_filt, "DiffusionMap") <- dm@eigenvectors[,c(1,2)]

# Define the starting point
foo <- mean(reducedDim(sce_filt,"DiffusionMap")[sce_filt$celltype==opts$celltypes[1],"DC1"])
if (foo>0) reducedDim(sce_filt,"DiffusionMap")[,"DC1"] <- -reducedDim(sce_filt,"DiffusionMap")[,"DC1"]

# Plot
pdf(file.path(args$outdir,sprintf("%s_diffmap_celltype.pdf",args$trajectory_name)), width=8, height=5)
plotReducedDim(sce_filt, dimred = "DiffusionMap", colour_by="celltype", ncomponents = c(1,2)) +
  scale_colour_manual(values=opts$celltype.colors) + theme(legend.position="none")
dev.off()

pdf(file.path(args$outdir,sprintf("%s_diffmap_stage.pdf",args$trajectory_name)), width=8, height=5)
plotReducedDim(sce_filt, dimred = "DiffusionMap", colour_by="stage", ncomponents = c(1,2)) +
  scale_colour_manual(values=opts$stage.colors) + theme(legend.position="top", legend.title = element_blank())
dev.off()

####################
## Prepare output ##
####################

pseudotime.dt <- data.table(
  cell = colnames(sce_filt),
  PC1 = reducedDim(sce_filt,"PCA")[,"PC1"],
  DC1 = reducedDim(sce_filt,"DiffusionMap")[,"DC1"]
) %>% setorder(PC1)

# Save
fwrite(pseudotime.dt, file.path(args$outdir,sprintf("/%s_trajectory.txt.gz",args$trajectory_name)))

# Load precomputed pseudotime
# pseudotime.dt <- fread(paste0(args$outdir,"/blood_trajectory.txt.gz"))

################################################
## Boxplots of pseudotime values per celltype ##
################################################

to.plot <- pseudotime.dt %>%
  merge(sample_metadata[,c("cell","celltype","stage")], by="cell") %>%
  .[,pca:=minmax.normalisation(PC1)] %>% .[,diffmap:=minmax.normalisation(DC1)]

for (i in c("pca","diffmap")) {
  p <- ggboxplot(to.plot, x="celltype", y=i, fill="celltype", outlier.shape=NA) +
    scale_fill_manual(values=opts$celltype.colors[opts$celltypes]) +
    stat_summary(fun.data = give.n, geom = "text", size=4) +
    labs(x="", y=sprintf("Pseudotime (%s)",i)) +
    theme_classic() +
    theme(
      legend.position = "none",
      # legend.position = c(.8,.3),
      legend.title = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      # axis.text.x = element_text(color="black", size=rel(0.75)),
      axis.text.y = element_text(color="black", size=rel(0.9)),
      axis.title.y = element_text(color="black", size=rel(1)),
      axis.title.x = element_blank()
    )
  
  pdf(file.path(args$outdir,sprintf("%s_%s_boxplot_celltype.pdf",args$trajectory_name,i)), width=6, height=3.5)
  print(p)
  dev.off()
}

#########################################################################
## Scatterplots of pseudotime values versus expression of marker genes ##
#########################################################################

# Denoise
# pca.rna <- fread(io$pca.rna) %>% matrix.please %>% .[colnames(sce),]
# assay(sce_filt,"logcounts_denoised") <- smoother_aggregate_nearest_nb(mat=as.matrix(logcounts(sce_filt)), D=pdist(pca.rna), k=25)

rna.dt <- data.table(as.matrix(logcounts(sce)[opts$genes_to_plot,]), keep.rownames = T) %>% 
  setnames("rn","gene") %>%
  melt(id.vars="gene", variable.name="cell", value.name="expr") %>% 
  merge(sample_metadata[,c("cell","celltype.mapped_mnn","stage")])

to.plot2 <- to.plot %>% 
  merge(rna.dt[gene%in%opts$genes_to_plot],by="cell") %>%
  .[,gene:=factor(gene,levels=opts$genes_to_plot)] %>%
  .[,expr:=minmax.normalisation(expr),by="gene"]

for (i in c("pca","diffmap")) {
  p <- ggplot(to.plot2, aes_string(x=i, y="expr")) +
    geom_point(aes(fill=celltype), size=1.75, shape=21, stroke=0.1, alpha=0.5) +
    # ggrastr::geom_point_rast(aes(fill=celltype), size=1.5, shape=21, stroke=0.1, alpha=0.75) +
    stat_smooth(aes(fill=expr), method="loess", color="black", alpha=0.85, span=0.5) +
    geom_rug(aes(color=celltype), sides="b") +
    # viridis::scale_fill_viridis() +
    scale_color_manual(values=opts$celltype.colors) +
    scale_fill_manual(values=opts$celltype.colors) +
    facet_wrap(~gene, scales="free_y") +
    labs(x=sprintf("Pseudotime (%s)",i), y="Gene expression") +
    theme_classic() +
    guides(fill="none", color="none") +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      # axis.text.y = element_text(color="black"),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      strip.text = element_text(size=rel(1.25), color="black"),
      legend.position="right",
      legend.title = element_blank()
    )
  
  pdf(file.path(args$outdir,sprintf("%s_%s_vs_expression.pdf",args$trajectory_name,i)), width=7.5, height=4)
  print(p)
  dev.off()
}


###############################
## Save SingleCellExperiment ##
###############################

logcounts(sce) <- NULL
saveRDS(sce, file.path(args$outdir,sprintf("%s_SingleCellExperiment.rds",args$trajectory_name)))

##################
## Save Anndata ##
##################

# # sce <- readRDS(file.path(args$outdir,sprintf("%s_SingleCellExperiment.rds",args$trajectory_name)))
# 
# # reticulate::use_condaenv("celloracle_env", required = T)
# reticulate::use_condaenv("scanpy", required = T)
# sc <- import("scanpy")
# 
# # Create anndata
# adata_sce <- sc$AnnData(
#     X   = t(counts(sce)),
#     obs = as.data.frame(colData(sce)),
#     var = data.frame(gene=rownames(sce), row.names=rownames(sce))
# )
# 
# # Add cell type colors
# adata_sce$uns$update("celltype_colors" = opts$celltype.colors[sort(unique(as.character(adata_sce$obs[[args$celltype_label]])))])
# 
# # Save
# adata_sce$write_h5ad(file.path(args$outdir,sprintf("%s_adata.h5ad",args$trajectory_name)))

