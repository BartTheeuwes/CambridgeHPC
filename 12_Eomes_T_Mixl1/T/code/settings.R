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
io$basedir <- "/rds/project/rds-SDzz0CATGms/users/bt392/12_Eomes_T_Mixl1/T/"
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
io$rna.atlas.metadata <- file.path(io$atlas.basedir,"sample_metadata_extended.txt.gz")
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
    'sample_1',
    'sample_2',
    'sample_5',
    'sample_6',
    'sample_7',
    'sample_8',
    'sample_9',
    'sample_10'
)

opts$rename.samples <- c(
'sample_1' = 'sample1KO',
'sample_2' = 'sample2WT',
'sample_5' = 'sample5KO',
'sample_6' = 'sample6WT',
'sample_7' = 'sample7KO',
'sample_8' = 'sample8WT',
'sample_9' = 'sample9KO',
'sample_10' = 'sample10WT'
)

opts$sample2stage <- c(
    'sample_1' = 'E8.5',
    'sample_2' = 'E8.5',
    'sample_5' = 'E8.5',
    'sample_6' = 'E8.5',
    'sample_7' = 'E8.5',
    'sample_8' = 'E8.5',
    'sample_9' = 'E8.5',
    'sample_10' = 'E8.5'
)


opts$sample2tomato = c(
    'sample_1' = TRUE,
    'sample_2' = FALSE,
    'sample_5' = TRUE,
    'sample_6' = FALSE,
    'sample_7' = TRUE,
    'sample_8' = FALSE,
    'sample_9' = TRUE,
    'sample_10' = FALSE
)


opts$sample2pool <- c(
    'sample_1' = 'pool1',
    'sample_2' = 'pool1',
    'sample_5' = 'pool2',
    'sample_6' = 'pool2',
    'sample_7' = 'pool3',
    'sample_8' = 'pool3',
    'sample_9' = 'pool4',
    'sample_10' = 'pool4'
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

opts$celltypes_extended <- c(
'Epiblast', 
'Primitive_Streak', 
'ExE_ectoderm', 
'Visceral_endoderm', 
'ExE_endoderm', 
'Non-neural_ectoderm', 
'Nascent_mesoderm', 
'Parietal_endoderm', 
'Ectoderm', 
'Anterior_Primitive_Streak', 
'Haematoendothelial_progenitors', 
'Caudal_epiblast', 
'Blood_progenitors', 
'Intermediate_mesoderm', 
'Paraxial_mesoderm', 
'Lateral_plate_mesoderm', 
'Mesenchyme', 
'PGC', 
'Node', 
'Gut_tube', 
'Embryo_proper_endothelium', 
'Cardiopharyngeal_progenitors_SHF', 
'Notochord', 
'Amniotic_ectoderm', 
'Venous_endothelium', 
'Presomitic_mesoderm', 
'Cardiomyocytes_FHF_1', 
'Allantois', 
'Cranial_mesoderm', 
'EMP', 
'Limb_mesoderm', 
'Anterior_somitic_tissues', 
'Pharyngeal_mesoderm', 
'Allantois_endothelium', 
'Thyroid_primordium', 
'Erythroid', 
'Hindbrain_neural_progenitors', 
'Cardiomyocytes_SHF 1', 
'NMPs', 
'Pharyngeal_endoderm', 
'Dorsal_spinal_cord_progenitors', 
'Anterior_cardiopharyngeal_progenitors', 
'Placodal_ectoderm', 
'Optic_vesicle', 
'Ventral_forebrain_progenitors', 
'Spinal_cord_progenitors', 
'Hindgut', 
'Caudal_mesoderm', 
'Embryo_proper_mesothelium', 
'Neural_tube', 
'Midbrain_Hindbrain_boundary', 
'Posterior_somitic_tissues', 
'Midgut', 
'Migratory_neural_crest', 
'Ventral_hindbrain_progenitors', 
'Surface_ectoderm', 
'YS_mesothelium', 
'Limb_ectoderm', 
'Somitic_mesoderm', 
'NMPs_Mesoderm-biased', 
'Cardiomyocytes_FHF_2', 
'Foregut', 
'Dermomyotome', 
'Kidney_primordium', 
'Otic_placode', 
'Cardiomyocytes_SHF_2', 
'Midbrain_progenitors', 
'Cardiopharyngeal_progenitors_FHF', 
'Epicardium', 
'Hindbrain_floor_plate', 
'Late_dorsal_forebrain_progenitors', 
'Dorsal_hindbrain_progenitors', 
'Sclerotome', 
'YS_endothelium', 
'Endocardium', 
'MEP', 
'Epidermis', 
'Megakaryocyte_progenitors', 
'Early_dorsal_forebrain_progenitors', 
'Cardiopharyngeal_progenitors', 
'Chorioallantoic-derived_erythroid_progenitors', 
'YS mesothelium-derived_endothelial_progenitors', 
'Dorsal_midbrain_neurons', 
'Branchial_arch_neural_crest', 
'Forelimb', 
'Frontonasal_mesenchyme', 
'Otic neural_progenitors'
)

opts$stage.colors = c(
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

opts$tdTom.color = c(
    'TRUE' = "red", 
    'FALSE' = "black"
)

opts$chr <- paste0("chr",c(1:19,"X","Y"))
