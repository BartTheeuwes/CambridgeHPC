suppressMessages(library(SingleCellExperiment))
suppressMessages(library(data.table))
suppressMessages(library(purrr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(stringr))
suppressMessages(library(argparse))
suppressMessages(library(parallel))


#########
## I/O ##
#########

io <- list()
io$basedir <- "/rds/project/rds-SDzz0CATGms/users/bt392/10_Eomes_invitro_gut"
io$atlas.basedir <- "/rds/project/rds-SDzz0CATGms/users/bt392/atlasses/gastrulation/pijuansala2019_gastrulation10x"
io$gene_metadata <- "/rds/project/rds-SDzz0CATGms/users/bt392/atlasses/Mmusculus_genes_BioMart.87.txt"
io$archR.directory <- file.path(io$basedir,"processed/atac/archR")

io$metadata <- file.path(io$basedir,"sample_metadata.txt.gz")

# TFs
# io$TFs <- file.path(io$basedir,"results/TFs.txt")

# RNA
io$rna.anndata <- file.path(io$basedir,"processed/rna/anndata.h5ad")
io$rna.seurat <- file.path(io$basedir,"processed/rna/seurat.rds")
io$rna.sce <- file.path(io$basedir,"processed/rna/SingleCellExperiment.rds")
io$rna.differential <- file.path(io$basedir,"results/rna/differential")
io$rna.pseudobulk.sce <- file.path(io$basedir,"results/rna/pseudobulk/SingleCellExperiment.rds")

# RNA atlas (Pijuan-Sala2019)
io$rna.atlas.metadata <- file.path(io$atlas.basedir,"sample_metadata.txt.gz")
io$rna.atlas.marker_genes <- file.path(io$atlas.basedir,"results/marker_genes/all_stages/marker_genes.txt.gz")
io$rna.atlas.marker_TFs <- file.path(io$atlas.basedir,"results/differential/celltypes/TFs/TF_markers/marker_TFs_up.txt.gz")
io$rna.atlas.differential <- file.path(io$atlas.basedir,"results/differential")
io$rna.atlas.sce.pseudobulk <- file.path(io$atlas.basedir,"results/pseudobulk/SingleCellExperiment_pseudobulk.rds")
io$rna.atlas.sce <- file.path(io$atlas.basedir,"processed/SingleCellExperiment.rds")
io$rna.atlas.celltype_proportions <- file.path(io$atlas.basedir,"results/celltype_proportions/celltype_proportions.txt.gz")

# motifmatchr
# io$motifmatcher.se <- sprintf("%s/Annotations/Motif_cisbp-Matches-In-Peaks.rds",io$archR.directory)
io$motifmatcher.se <- sprintf("%s/Annotations/Motif_cisbp_lenient-Scores.rds",io$archR.directory)
io$motifmatcher_positions.se <- sprintf("%s/Annotations/Motif_cisbp-Positions-In-Peaks.rds",io$archR.directory)

# ATAC: archR
# io$atac.peak.annotation <- file.path(io$basedir,"/original/atac_peak_annotation.tsv")
io$archR.projectMetadata <- file.path(io$archR.directory,"projectMetadata.rds")
io$archR.peakSet.granges <- file.path(io$archR.directory,"PeakSet.rds")
io$archR.bgdPeaks <- file.path(io$archR.directory,"Background-Peaks.rds")
io$archR.peakSet.bed <- file.path(io$archR.directory,"PeakCalls/bed/peaks_archR_macs2.bed.gz")
io$archR.GeneScoreMatrix.se <- file.path(io$archR.directory,"/GeneScoreMatrix_no_distal_summarized_experiment.rds")
io$archR.peakMatrix.se <- file.path(io$archR.directory,"PeakCalls/PeakMatrix_summarized_experiment.rds")
io$archR.peak.variability <- file.path(io$basedir,"results/atac/archR/variability/peak_variability.txt.gz")
io$archR.peak.differential.dir <- file.path(io$basedir,"results/atac/archR/differential/PeakMatrix")
io$archR.peak.metadata <- file.path(io$archR.directory,"PeakCalls/peak_metadata.tsv.gz")
io$archR.peak.stats <- file.path(io$basedir,"results/atac/archR/peak_calling/peak_stats.txt.gz")
io$archR.peak2gene.all <- file.path(io$basedir,"results/atac/archR/peak_calling/peaks2genes/peaks2genes_all.txt.gz")
io$archR.peak2gene.nearest <- file.path(io$basedir,"results/atac/archR/peak_calling/peaks2genes/peaks2genes_nearest.txt.gz")
io$archr.chromvar.dir <- file.path(io$basedir,"results/atac/archR/chromvar")
io$archR.pseudobulk.GeneScoreMatrix.se <- file.path(io$archR.directory,"pseudobulk/pseudobulk_GeneScoreMatrix_summarized_experiment.rds")
io$archR.pseudobulk.peakMatrix.se <- file.path(io$archR.directory,"pseudobulk/pseudobulk_PeakMatrix_summarized_experiment.rds")
io$archR.pseudobulk.deviations.se <- file.path(io$basedir,"results/atac/archR/chromvar/pseudobulk/chromVAR_deviations_Motif_cisbp_lenient_archr_chip.rds")

# paga
io$paga.connectivity <- file.path(io$atlas.basedir,"results/paga/paga_connectivity.csv")
io$paga.coordinates <- file.path(io$atlas.basedir,"results/paga/paga_coordinates.csv")

# PCA
# io$pca.rna <- file.path(io$basedir,"results/rna/dimensionality_reduction/all_cells/rv_eo_deg_day3_5_control-rv_eo_deg_day3_5_dtag-E8.0_rep1-E8.0_rep2-E8.5_rep1-E8.5_rep2_pca_features2500_pcs30_batchcorrectionbysample.txt.gz")
# io$pca.atac <- file.path(io$basedir,"results/atac/archR/dimensionality_reduction/PeakMatrix/all_cells/rv_eo_deg_day3_5_control-rv_eo_deg_day3_5_dtag-E8.0_rep1-E8.0_rep2-E8.5_rep1-E8.5_rep2_lsi_features50000_ndims50.txt.gz")

# UMAP
# io$umap.rna <- file.path(io$basedir,"results/rna/dimensionality_reduction/all_cells/rv_eo_deg_day3_5_control-rv_eo_deg_day3_5_dtag-E8.0_rep1-E8.0_rep2-E8.5_rep1-E8.5_rep2_umap_features2500_pcs30_neigh25_dist0.3.txt.gz")
# io$umap.atac <- file.path(io$basedir,"results/atac/archR/dimensionality_reduction/PeakMatrix/all_cells/rv_eo_deg_day3_5_control-rv_eo_deg_day3_5_dtag-E8.0_rep1-E8.0_rep2-E8.5_rep1-E8.5_rep2_umap_nfeatures50000_ndims50_neigh30_dist0.45.txt.gz")

#############
## Options ##
#############

opts <- list()

opts$days <- c(
  "D2.5",
  "D3",
  "D4.5"
)


opts$samples <- c(
  # first batch
  "RBG43471",
  "RBG43472",
  "RBG43473",
  "RBG43474",
  "RBG43475",
  "RBG43476",
  "RBG43477",
  "RBG43478"
)

opts$sample2day <- c(
  "RBG43471" = "D2.5",
  "RBG43472" = "D2.5",
  "RBG43473" = "D3",
  "RBG43474" = "D3",
  "RBG43475" = "D3",
  "RBG43476" ="D3",
  "RBG43477" = "D4.5",
  "RBG43478" = "D4.5"
)

opts$sample2genotype <- c(
  "RBG43471" = "WT",
  "RBG43472" = "WT",
  "RBG43473" = "Eomes_KO",
  "RBG43474" = "Eomes_KO",
  "RBG43475" = "WT",
  "RBG43476" ="WT",
  "RBG43477" = "Eomes_KO",
  "RBG43478" = "WT"
)

opts$sample2alias <- c(
  "RBG43471" = "D2.5_WT_rep1",
  "RBG43472" = "D2.5_WT_rep2",
  "RBG43473" = "D3_KO_rep1",
  "RBG43474" = "D3_KO_rep2",
  "RBG43475" = "D3_WT_rep1",
  "RBG43476" = "D3_WT_rep2",
  "RBG43477" = "D4.5_KO_rep1",
  "RBG43478" = "D4.5_WT_rep1"
)

opts$genotypes <- c("WT","Eomes_KO")

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

opts$days.colors = c(
  "D5" = "#fde725",
  "D4.5" = "#5dc863",
  "D4" = "#21908c",
  "D3.5" = "#3b528b",
  "D3" = "#3d0154",
  "D2.5" = "#540134"
)

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

opts$celltype_ext.colors = c(
"Epiblast" = "#635547",
"Primitive Streak" = "#DABE99",
"Caudal epiblast" = "#9e6762",

"PGC" = "#FACB12",

"Anterior Primitive Streak" = "#c19f70",
"Node"="#153b3d",
"Notochord" = "#0F4A9C",



"Gut tube" = "#EF5A9D",
"Hindgut" = "#F397C0",
"Midgut" = "#ff00b2",
"Foregut" = "#ffb7ff",
"Pharyngeal endoderm"="#95e1ff",
"Thyroid primordium"="#97bad3",

"Nascent mesoderm" = "#C594BF",
"Intermediate mesoderm" = "#139992",
"Caudal mesoderm" = "#3F84AA",
"Lateral plate mesoderm" = "#F9DFE6",
"Limb mesoderm" = "#e35f82",
"Forelimb" = "#d02d75",
"Kidney primordium" = "#e85639",
"Presomitic mesoderm"="#5581ca",#"#0000ff",#blue
"Somitic mesoderm" = "#005579",
"Posterior somitic tissues" = "#5adbe4",#"#40e0d0",#turquoise
"Mesenchyme" = "#cc7818",
"Mixed mesoderm" = "#DFCDE4",



"Paraxial mesoderm" = "#8DB5CE",
"Cranial mesoderm" = "#456722",#"#006400",#darkgreen
"Anterior somitic tissues"= "#d5e839",
"Sclerotome" = "#e3cb3a",#"#ffff00",#yellow
"Dermomyotome" = "#00BFC4",#"#a52a2a",#brown



"Pharyngeal mesoderm" = "#C9EBFB",
"Cardiopharyngeal progenitors" = "#556789",
"Anterior cardiopharyngeal progenitors"="#683ed8",



"Allantois" = "#532C8A",
"Mesenchyme" = "#cc7818",
"YS mesothelium" = "#ff7f9c",
"Epicardium"="#f79083",
"Embryo proper mesothelium" = "#ff487d",



"Cardiopharyngeal progenitors FHF"="#d780b0",
"Cardiomyocytes FHF 1"="#a64d7e",
"Cardiomyocytes FHF 2"="#B51D8D",



"Cardiopharyngeal progenitors SHF"="#4b7193",
"Cardiomyocytes SHF 1"="#5d70dc",
"Cardiomyocytes SHF 2"="#332c6c",


                     
"Haematoendothelial progenitors" = "#FBBE92",
"Blood progenitors" = "#6c4b4c",
"Blood progenitors 1" = "#f9decf",
"Blood progenitors 2" = "#c9a997",
"Erythroid" = "#C72228",
"Erythroid1" = "#C72228",
"Erythroid2" = "#f79083",
"Erythroid3" = "#EF4E22",
"Chorioallantoic-derived erythroid progenitors"="#E50000",
"Megakaryocyte progenitors"="#e3cb3a",
"MEP"="#EF4E22",
"EMP"="#7c2a47",


"Endothelium" = "#ff891c",
"YS endothelium"="#ff891c",
"YS mesothelium-derived endothelial progenitors"="#AE3F3F",
"Allantois endothelium"="#2f4a60",
"Embryo proper endothelium"="#90e3bf",
"Venous endothelium"="#bd3400",
"Endocardium"="#9d0049",

"NMPs/Mesoderm-biased" = "#89c1f5",
"NMPs" = "#8EC792",

"Ectoderm" = "#ff675c",



"Optic vesicle" = "#bd7300",

"Ventral forebrain progenitors"="#a0b689",
"Early dorsal forebrain progenitors"="#0f8073",
"Late dorsal forebrain progenitors"="#7a9941",
"Midbrain/Hindbrain boundary"="#8ab3b5",
"Midbrain progenitors"="#9bf981",
"Dorsal midbrain neurons"="#12ed4c",
"Ventral hindbrain progenitors"="#7e907a",
"Dorsal hindbrain progenitors"="#2c6521",
"Hindbrain floor plate"="#bf9da8",
"Hindbrain neural progenitors"="#59b545",



"Neural tube"="#233629",



"Migratory neural crest"="#4a6798",
"Branchial arch neural crest"="#bd84b0",
"Frontonasal mesenchyme"="#d3b1b1",



"Spinal cord progenitors"="#6b2035",
"Dorsal spinal cord progenitors"="#e273d6",

"Non-neural ectoderm" = "#f7f79e",
"Surface ectoderm" = "#fcff00",
"Epidermis" = "#fff335",
"Limb ectoderm" = "#ffd731",
"Amniotic ectoderm" = "#dbb400",



"Placodal ectoderm" = "#ff5c00",



"Otic placode"="#f1a262",
"Otic neural progenitors"="#00b000",

"Visceral endoderm" = "#F6BFCB",
"ExE endoderm" = "#7F6874",
"ExE ectoderm" = "#989898",
"Parietal endoderm" = "#1A1A1A"
)

opts$genotype.colors = c(
    'Eomes_KO' = "red", 
    'WT' = "black"
)