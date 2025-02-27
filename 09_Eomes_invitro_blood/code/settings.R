suppressMessages(library(SingleCellExperiment))
suppressMessages(library(data.table))
suppressMessages(library(purrr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(stringr))
suppressMessages(library(argparse))
suppressMessages(library(parallel))
suppressPackageStartupMessages(library(BSgenome.Mmusculus.UCSC.mm10))


#########
## I/O ##
#########

io <- list()
io$basedir <- "/rds/project/rds-SDzz0CATGms/users/bt392/09_Eomes_invitro_blood"
io$atlas.basedir <- "/rds/project/rds-SDzz0CATGms/users/bt392/atlasses/gastrulation/pijuansala2019_gastrulation10x"
io$atlas.extended.sce = "/rds/project/rds-SDzz0CATGms/users/bt392/atlasses/extended/embryo_sce.rds"
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
  "D3",
  "D3.5",
  "D4",
  "D4.5",
  "D5"
)


opts$samples <- c(
  # first batch
  "rv_eo_deg_day3_5_control",
  "rv_eo_deg_day3_5_dtag",
  "rv_eo_deg_day4_5_control",
  "rv_eo_deg_day4_5_dtag",
  "rv_eo_deg_day4_control",
  "rv_eo_deg_day4_dtag",

  # second batch
  "1A_Eo_DEG_G9_day3",
  "1B_Eo_DEG_G9_day3",
  "2_Eo_DEG_G9_day3_5_dTAG",
  "2_Eo_DEG_G9_day3_5_VC",
  "2_Eo_DEG_G9_day4_dTAG",
  "2_Eo_DEG_G9_day4_VC",
  "2_Eo_DEG_G9_day5_dTAG",
  "2_Eo_DEG_G9_day5_VC"
)

opts$sample2day <- c(
  "rv_eo_deg_day3_5_control" = "D3.5",
  "rv_eo_deg_day3_5_dtag" = "D3.5",
  "rv_eo_deg_day4_5_control" = "D4.5",
  "rv_eo_deg_day4_5_dtag" = "D4.5",
  "rv_eo_deg_day4_control" = "D4",
  "rv_eo_deg_day4_dtag" ="D4",
  "1A_Eo_DEG_G9_day3" = "D3",
  "1B_Eo_DEG_G9_day3" = "D3",
  "2_Eo_DEG_G9_day3_5_dTAG" = "D3.5",
  "2_Eo_DEG_G9_day3_5_VC" = "D3.5",
  "2_Eo_DEG_G9_day4_dTAG" = "D4",
  "2_Eo_DEG_G9_day4_VC" = "D4",
  "2_Eo_DEG_G9_day5_dTAG" = "D5",
  "2_Eo_DEG_G9_day5_VC" = "D5"
)

opts$sample2genotype <- c(
  "rv_eo_deg_day3_5_control" = "WT",
  "rv_eo_deg_day3_5_dtag" = "KO",
  "rv_eo_deg_day4_5_control" = "WT",
  "rv_eo_deg_day4_5_dtag" = "KO",
  "rv_eo_deg_day4_control" = "WT",
  "rv_eo_deg_day4_dtag" ="KO",
  "1A_Eo_DEG_G9_day3" = "WT",
  "1B_Eo_DEG_G9_day3" = "WT",
  "2_Eo_DEG_G9_day3_5_dTAG" = "KO",
  "2_Eo_DEG_G9_day3_5_VC" = "WT",
  "2_Eo_DEG_G9_day4_dTAG" = "KO",
  "2_Eo_DEG_G9_day4_VC" = "WT",
  "2_Eo_DEG_G9_day5_dTAG" = "KO",
  "2_Eo_DEG_G9_day5_VC" = "WT"
)

opts$sample2exp <- c(
  "rv_eo_deg_day3_5_control" = "exp1",
  "rv_eo_deg_day3_5_dtag" = "exp1",
  "rv_eo_deg_day4_5_control" = "exp1",
  "rv_eo_deg_day4_5_dtag" = "exp1",
  "rv_eo_deg_day4_control" = "exp1",
  "rv_eo_deg_day4_dtag" ="exp1",
  "1A_Eo_DEG_G9_day3" = "exp2",
  "1B_Eo_DEG_G9_day3" = "exp2",
  "2_Eo_DEG_G9_day3_5_dTAG" = "exp2",
  "2_Eo_DEG_G9_day3_5_VC" = "exp2",
  "2_Eo_DEG_G9_day4_dTAG" = "exp2",
  "2_Eo_DEG_G9_day4_VC" = "exp2",
  "2_Eo_DEG_G9_day5_dTAG" = "exp2",
  "2_Eo_DEG_G9_day5_VC" = "exp2"
)

opts$sample2alias <- c(
  "rv_eo_deg_day3_5_control" = "day3.5_WT_rep1",
  "rv_eo_deg_day3_5_dtag" = "day3.5_KO_rep1",
  "rv_eo_deg_day4_5_control" = "day4.5_WT_rep1",
  "rv_eo_deg_day4_5_dtag" = "day4.5_KO_rep1",
  "rv_eo_deg_day4_control" = "day4_WT_rep1",
  "rv_eo_deg_day4_dtag" ="day4_KO_rep1",
  "1A_Eo_DEG_G9_day3" = "day3_WT_rep1",
  "1B_Eo_DEG_G9_day3" = "day3_WT_rep2",
  "2_Eo_DEG_G9_day3_5_VC" = "day3.5_WT_rep2",
  "2_Eo_DEG_G9_day3_5_dTAG" = "day3.5_KO_rep2",
  "2_Eo_DEG_G9_day4_VC" = "day4_WT_rep2",
  "2_Eo_DEG_G9_day4_dTAG" = "day4_KO_rep2",
  "2_Eo_DEG_G9_day5_VC" = "day5_WT_rep1",
  "2_Eo_DEG_G9_day5_dTAG" = "day5_KO_rep1"
)

opts$sample2rep <- c(
  "rv_eo_deg_day3_5_control" = "rep1",
  "rv_eo_deg_day3_5_dtag" = "rep1",
  "rv_eo_deg_day4_5_control" = "rep1",
  "rv_eo_deg_day4_5_dtag" = "rep1",
  "rv_eo_deg_day4_control" = "rep1",
  "rv_eo_deg_day4_dtag" ="rep1",
  "1A_Eo_DEG_G9_day3" = "rep1",
  "1B_Eo_DEG_G9_day3" = "rep2",
  "2_Eo_DEG_G9_day3_5_VC" = "rep2",
  "2_Eo_DEG_G9_day3_5_dTAG" = "rep2",
  "2_Eo_DEG_G9_day4_VC" = "rep2",
  "2_Eo_DEG_G9_day4_dTAG" = "rep2",
  "2_Eo_DEG_G9_day5_VC" = "rep2",
  "2_Eo_DEG_G9_day5_dTAG" = "rep2"
)

opts$genotypes <- c("WT","KO")

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

opts$stage.colors.original = c(
                  "mixed_gastrulation" = "#A9A9A9",
                  "E6.5" = "#D53E4F",
                  "E6.75" = "#F46D43",
                  "E7.0" = "#FDAE61",
                  "E7.25" = "#FEE08B",
                  "E7.5" = "#FFFFBF",
                  "E7.75" = "#E6F598",
                  "E8.0" = "#ABDDA4",
                  "E8.25" = "#66C2A5",
                  "E8.5" = "#3288BD"
 )

# opts$days.colors = c(
#   "D5" = "#21908CFF",
#   "D4.5" = "#27AD81FF",
#   "D4" = "#5DC863FF",
#   "D3.5" = "#AADC32FF",
#   "D3" = "#FDE725FF"
# )

opts$days.colors = c(
  "D5" = "#fde725",
  "D4.5" = "#5dc863",
  "D4" = "#21908c",
  "D3.5" = "#3b528b",
  "D3" = "#3d0154"
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

opts$celltype_extended.colors = c(
"Epiblast" = "#635547",
"Primitive_Streak" = "#DABE99",
"Caudal_epiblast" = "#9e6762",

"PGC" = "#FACB12",

"Anterior_Primitive Streak" = "#c19f70",
"Node"="#153b3d",
"Notochord" = "#0F4A9C",



"Gut_tube" = "#EF5A9D",
"Hindgut" = "#F397C0",
"Midgut" = "#ff00b2",
"Foregut" = "#ffb7ff",
"Pharyngeal_endoderm"="#95e1ff",
"Thyroid_primordium"="#97bad3",

"Nascent_mesoderm" = "#C594BF",
"Intermediate_mesoderm" = "#139992",
"Caudal_mesoderm" = "#3F84AA",
"Lateral_plate_mesoderm" = "#F9DFE6",
"Limb_mesoderm" = "#e35f82",
"Forelimb" = "#d02d75",
"Kidney_primordium" = "#e85639",
"Presomitic_mesoderm"="#5581ca",#"#0000ff",#blue
"Somitic_mesoderm" = "#005579",
"Posterior_somitic tissues" = "#5adbe4",#"#40e0d0",#turquoise



"Paraxial_mesoderm" = "#8DB5CE",
"Cranial_mesoderm" = "#456722",#"#006400",#darkgreen
"Anterior_somitic_tissues"= "#d5e839",
"Sclerotome" = "#e3cb3a",#"#ffff00",#yellow
"Dermomyotome" = "#00BFC4",#"#a52a2a",#brown



"Pharyngeal_mesoderm" = "#C9EBFB",
"Cardiopharyngeal_progenitors" = "#556789",
"Anterior_cardiopharyngeal_progenitors"="#683ed8",



"Allantois" = "#532C8A",
"Mesenchyme" = "#cc7818",
"YS_mesothelium" = "#ff7f9c",
"Epicardium"="#f79083",
"Embryo_proper_mesothelium" = "#ff487d",



"Cardiopharyngeal_progenitors_FHF"="#d780b0",
"Cardiomyocytes_FHF_1"="#a64d7e",
"Cardiomyocytes_FHF_2"="#B51D8D",



"Cardiopharyngeal_progenitors_SHF"="#4b7193",
"Cardiomyocytes_SHF_1"="#5d70dc",
"Cardiomyocytes_SHF_2"="#332c6c",



"Haematoendothelial_progenitors" = "#FBBE92",
"Blood_progenitors" = "#6c4b4c",
"Erythroid" = "#C72228",
"Chorioallantoic-derived_erythroid_progenitors"="#E50000",
"Megakaryocyte_progenitors"="#e3cb3a",
"MEP"="#EF4E22",
"EMP"="#7c2a47",



"YS_endothelium"="#ff891c",
"YS_mesothelium-derived_endothelial_progenitors"="#AE3F3F",
"Allantois_endothelium"="#2f4a60",
"Embryo_proper_endothelium"="#90e3bf",
"Venous_endothelium"="#bd3400",
"Endocardium"="#9d0049",

"NMPs_Mesoderm-biased" = "#89c1f5",
"NMPs" = "#8EC792",

"Ectoderm" = "#ff675c",



"Optic_vesicle" = "#bd7300",

"Ventral_forebrain_progenitors"="#a0b689",
"Early_dorsal_forebrain_progenitors"="#0f8073",
"Late_dorsal_forebrain_progenitors"="#7a9941",
"Midbrain_Hindbrain_boundary"="#8ab3b5",
"Midbrain_progenitors"="#9bf981",
"Dorsal_midbrain_neurons"="#12ed4c",
"Ventral_hindbrain_progenitors"="#7e907a",
"Dorsal_hindbrain_progenitors"="#2c6521",
"Hindbrain_floor_plate"="#bf9da8",
"Hindbrain_neural_progenitors"="#59b545",



"Neural_tube"="#233629",



"Migratory_neural_crest"="#4a6798",
"Branchial_arch_neural_crest"="#bd84b0",
"Frontonasal_mesenchyme"="#d3b1b1",



"Spinal_cord_progenitors"="#6b2035",
"Dorsal_spinal_cord_progenitors"="#e273d6",

"Non-neural_ectoderm" = "#f7f79e",
"Surface_ectoderm" = "#fcff00",
"Epidermis" = "#fff335",
"Limb_ectoderm" = "#ffd731",
"Amniotic_ectoderm" = "#dbb400",



"Placodal_ectoderm" = "#ff5c00",



"Otic_placode"="#f1a262",
"Otic_neural progenitors"="#00b000",

"Visceral_endoderm" = "#F6BFCB",
"ExE_endoderm" = "#7F6874",
"ExE_ectoderm" = "#989898",
"Parietal_endoderm" = "#1A1A1A"
)

opts$genotype.colors = c(
    'KO' = "red", 
    'WT' = "black"
)

# opts$clusters2celltypes = c(
#     '10' = 'Primitive_Streak', 
#     '11' = 'Primitive_Streak', 
#     '2' = 'Primitive_Streak',   
#     '0' = 'Mesoderm_Eomes+_PdgfRa+',
#     '1' = 'Mesoderm_Eomes+_PdgfRa+',
#     '8' = 'Mesoderm_Eomes+_PdgfRa+',
#     '17' = 'Mesoderm_Eomes+_PdgfRa+',
#     '5' = 'Mesoderm_Eomes-_T+',
#     '6' = 'Mesoderm_Eomes-_T+',
#     '15' = 'Mesoderm_Eomes-_T+',
#     '18' = 'Mesoderm_Meis2+',    #_Hand1+_Greb1l+, but sort of mutually exclusive?
#     '4' = 'Posterior_Mesoderm_Cdx4+_Hoxa9+',
#     '13' = 'Posterior_Mesoderm_Cdx4+_Hoxa9+',
#     '21' = 'Posterior_Mesoderm_Cdx4+_Hoxa9+',
#     '9' = 'Mesenchyme_Podxl+_Foxf1+',
#     '14' = 'Allantois_Precursor_Tbx4+_Pitx1-',
#     '19' = 'Allantois',
#     '3' = 'Allantois',
#     '12' = 'HE_precursor',
#     '7' = 'HE',
#     '16' = 'Endothelium', 
#     '20' = 'blood'
# )

opts$clusters2celltypes = c(
    '0' = 'Early_Mes_EOd',
    '1' = 'Early_Mes_EOd',
    '2' = 'Primitive_Streak',
    '3' = 'HE',
    '4' = 'Posterior_Mes',
    '5' = 'Mesenchyme',
    '6' = 'Early_Mes_EOi',
    '7' = 'Primitive_Streak',
    '8' = 'Early_Mes_EOi',
    '9' = 'Early_Mes_EOi',
    '10' = 'Allantois',
    '11' = 'Primitive_Streak',
    '12' = 'Allantois_Precursor',
    '13' = 'Posterior_Mes',
    '14' = 'HE_Precursor',
    '15' = 'Endothelium',
    '16' = 'Primitive_Streak',
    '17' = 'Blood_Progenitor'    
)

opts$celltype_levels = c(
    'Primitive_Streak',
    'Early_Mes_EOi',
    'Early_Mes_EOd',
    'HE_Precursor',
    'HE',
    'Blood_Progenitor',
    'Endothelium',
    'Allantois',
    'Allantois_Precursor',
    'Mesenchyme',
    'Posterior_Mes'
)
opts$celltype_levels2 = c(
    'Allantois',
    'Allantois_Precursor',
    'Mesenchyme',
    'Posterior_Mes',
    'Early_Mes_EOi',
    'Primitive_Streak',
    'Early_Mes_EOd',
    'HE_Precursor',
    'HE',
    'Endothelium', 
    'Blood_Progenitor' 
    )

opts$celltype_v1.colors = c(ArchR::ArchRPalettes$stallion, ArchR::ArchRPalettes$stallion2)[1:length(unique(opts$clusters2celltypes))]
names(opts$celltype_v1.colors) = unique(opts$clusters2celltypes)

opts$celltype_v2.colors = c(
    'PGC' = "#FACB12",
    'Allantois' = "#532C8A",
    'Allantois_Precursor' = "#8870ad", 
    'Mesenchyme' = "#cc7818", 
    'Posterior_Mes' = "#ffbfd0",  # "#F9DFE6"
    'Early_Mes_EOi' = "#a0d0fa", 
    'Primitive_Streak' = "#DABE99", 
    'Early_Mes_EOd' = "#fc4c5c",  # "#ffffba",  
    'HE_Precursor' = "#917057", 
    'HE' = "#FBBE92", 
    'Endothelium' = "#ff891c",
    'Blood_Progenitor' = "#f9decf"
)

opts$celltype_v2_figure = c(
'PGC' = 'PGC',
'Allantois' = 'Allantois',
'Allantois_Precursor' =  'Allantois Precursor',
'Mesenchyme' = 'Mesenchyme',
'Posterior_Mes' = 'Posterior Mes',
'Early_Mes_EOi' = 'Early-Mes #1',
'Primitive_Streak' = 'Primitive Streak',
'Early_Mes_EOd' = 'Early-Mes #2',
'HE_Precursor' = 'HE Precursor',
'HE' = 'HE',
'Endothelium' = 'Endothelium',
'Blood_Progenitor' = 'Blood Progenitor' 
)

opts$celltype_v2.colors2 = c(
    'PGC' = "#FACB12",
    'Allantois' = "#532C8A",
    'Allantois Precursor' = "#8870ad", 
    'Mesenchyme' = "#cc7818", 
    'Posterior Mes' = "#ffbfd0",  # "#F9DFE6"
    'Early-Mes #1' = "#a0d0fa", 
    'Primitive Streak' = "#DABE99", 
    'Early-Mes #2' = "#fc4c5c",  # "#ffffba",  
    'HE Precursor' = "#917057", 
    'HE' = "#FBBE92", 
    'Endothelium' = "#ff891c",
    'Blood Progenitor' = "#f9decf"
)