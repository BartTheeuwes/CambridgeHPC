suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallel))

#########
## I/O ##
#########

io <- list()
io$basedir <- "/rds/project/rds-SDzz0CATGms/users/bt392/06_Runx1_RNA"
io$atlas.basedir <- "/rds/project/rds-SDzz0CATGms/users/bt392/atlasses/gastrulation/pijuansala2019_gastrulation10x"
io$gene_metadata <- "/rds/project/rds-SDzz0CATGms/users/bt392/atlasses/Mmusculus_genes_BioMart.87.txt"


io$metadata <- file.path(io$basedir,"processed/metadata.txt.gz")

# TFs
io$TFs <- file.path(io$basedir,"results_new/TFs.txt")

# RNA
io$rna.anndata <- file.path(io$basedir,"processed/anndata.h5ad")
io$rna.seurat <- file.path(io$basedir,"processed/seurat.rds")
io$rna.sce <- file.path(io$basedir,"processed/SingleCellExperiment.rds")
io$rna.differential <- file.path(io$basedir,"results/rna/differential")
io$rna.pseudobulk.sce <- file.path(io$basedir,"results/rna/pseudobulk/SingleCellExperiment.rds")

# RNA atlas (Pijuan-Sala2019)
io$rna.atlas.metadata <- file.path(io$atlas.basedir,"sample_metadata.txt.gz")
io$rna.atlas.marker_genes.up <- file.path(io$atlas.basedir,"results/marker_genes/all_stages/marker_genes.txt.gz")
# io$rna.atlas.marker_genes.all <- file.path(io$atlas.basedir,"results/marker_genes/all_stages/marker_genes_all.txt.gz")
io$rna.atlas.marker_TFs.up <- file.path(io$atlas.basedir,"results/differential/celltypes/TFs/TF_markers/marker_TFs_up.txt.gz")
# io$rna.atlas.marker_TFs.all <- file.path(io$atlas.basedir,"results/differential/celltypes/TFs/TF_markers/marker_TFs_all.txt.gz")
io$rna.atlas.differential <- file.path(io$atlas.basedir,"results/differential")
# io$rna.atlas.average_expression_per_celltype <- file.path(io$atlas.basedir,"results/marker_genes/all_stages/avg_expr_per_celltype_and_gene.txt.gz")
io$rna.atlas.sce.pseudobulk <- file.path(io$atlas.basedir,"results/pseudobulk/SingleCellExperiment_pseudobulk.rds")
io$rna.atlas.sce <- file.path(io$atlas.basedir,"processed/SingleCellExperiment.rds")
io$rna.atlas.celltype_proportions <- file.path(io$atlas.basedir,"results/celltype_proportions/celltype_proportions.txt.gz")

#############
## Options ##
#############

opts <- list()

opts$stages <- c(
  "E8.5"
)

opts$samples <- c(
    'SLX-21184_SITTA7_H5NKNDMXY',
    'SLX-21184_SITTB7_H5NKNDMXY',
    'SLX-21184_SITTC7_H5NKNDMXY',
    'SLX-21184_SITTD7_H5NKNDMXY',
    'SLX-21184_SITTE12_H5NKNDMXY',
    'SLX-21184_SITTH10_H5NKNDMXY'
)

opts$rename.samples <- c(
  "E8.75_rep1" = "e8_75_1_L002",
  "E8.75_rep2" = "e8_75_2_L002",
  "E8.5_rep1" = "multiome1",
  "E8.5_rep2" = "multiome2",
  "E8.0_rep1" = "E8_0_rep1_multiome",
  "E8.0_rep2" = "E8_0_rep2_multiome",
  "E7.5_rep1" = "rep1_L001_multiome",
  "E7.5_rep2" = "rep2_L002_multiome"
)

opts$sample2stage <- c(
  "21184_SITTA7_H5NKNDMXY" = "E8.5",
  "21184_SITTB7_H5NKNDMXY" = "E8.5",
  "21184_SITTC7_H5NKNDMXY" = "E8.5",
  "21184_SITTD7_H5NKNDMXY" = "E8.5",
  "21184_SITTE12_H5NKNDMXY" = "E8.5",
  "21184_SITTH10_H5NKNDMXY" = "E8.5"
)


opts$celltypes <- c(
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

opts$stage.colors = c(
   "E8.5" = "#440154FF",
   "E8.25" = "#472D7BFF",
   "E8.0" = "#3B528BFF",
   "E7.75" = "#2C728EFF",
   "E7.5" = "#21908CFF",
   "E7.25" = "#27AD81FF",
   "E7.0" = "#5DC863FF",
   "E6.75" = "#AADC32FF",
   "E6.5" = "#FDE725FF"
 )
opts$stage.colors <- viridis::viridis(n=length(opts$stages))
names(opts$stage.colors) <- rev(opts$stages)

opts$celltype.colors = c(
  "Epiblast" = "#635547",
  "Primitive_Streak" = "#DABE99",
  "Caudal_epiblast" = "#9e6762",
  "PGC" = "#FACB12",
  "Anterior_Primitive_Streak" = "#c19f70",
  "Notochord" = "#0F4A9C",
  "Def._endoderm" = "#F397C0",
  "Gut" = "#EF5A9D",
  "Nascent_mesoderm" = "#C594BF",
  "Mixed_mesoderm" = "#DFCDE4",
  "Intermediate_mesoderm" = "#139992",
  "Caudal_Mesoderm" = "#3F84AA",
  "Paraxial_mesoderm" = "#8DB5CE",
  "Somitic_mesoderm" = "#005579",
  "Pharyngeal_mesoderm" = "#C9EBFB",
  "Cardiomyocytes" = "#B51D8D",
  "Allantois" = "#532C8A",
  "ExE_mesoderm" = "#8870ad",
  "Mesenchyme" = "#cc7818",
  "Haematoendothelial_progenitors" = "#FBBE92",
  "Endothelium" = "#ff891c",
  "Blood_progenitors" = "#c9a997",
  "Blood_progenitors_1" = "#f9decf",
  "Blood_progenitors_2" = "#c9a997",
  "Erythroid" = "#EF4E22",
  "Erythroid1" = "#C72228",
  "Erythroid2" = "#f79083",
  "Erythroid3" = "#EF4E22",
  "NMP" = "#8EC792",
  "Neurectoderm" = "#65A83E",
  "Rostral_neurectoderm" = "#65A83E",
  "Caudal_neurectoderm" = "#354E23",
  "Neural_crest" = "#C3C388",
  "Forebrain_Midbrain_Hindbrain" = "#647a4f",
  "Spinal_cord" = "#CDE088",
  "Surface_ectoderm" = "#f7f79e",
  "Visceral_endoderm" = "#F6BFCB",
  "ExE_endoderm" = "#7F6874",
  "ExE_ectoderm" = "#989898",
  "Parietal_endoderm" = "#1A1A1A"
)

opts$chr <- paste0("chr",c(1:19,"X","Y"))
