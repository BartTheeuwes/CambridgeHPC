# suppressMessages(library(SingleCellExperiment))
suppressMessages(library(data.table))
suppressMessages(library(purrr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(stringr))
# suppressMessages(library(Seurat))

#########
## I/O ##
#########

io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020"
  io$gene_metadata <- "/Users/ricard/data/ensembl/mouse/v87/BioMart/all_genes/Mmusculus_genes_BioMart.87.txt"
  io$encode_blacklist <- "/Users/ricard/data/mm10_regulation/encode_blacklisted/mm10.blacklist.bed.gz"
  io$atlas.basedir <- "/Users/ricard/data/gastrulation10x"
} else if (grepl("ebi",Sys.info()['nodename'])) {
  io$basedir <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020"
  io$gene_metadata <- "/hps/nobackup2/research/stegle/users/ricard/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
  io$encode_blacklist <- ""
  io$atlas.basedir <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x"
} else if (grepl("pebble|headstone", Sys.info()['nodename'])) {
  # io$basedir <- "/bi/scratch/Stephen_Clark/gastrulation_multiome_10x/"
} else if (grepl("bi2228m",Sys.info()['nodename'])) {
  # io$basedir <- "/Users/clarks/data/10X_multiome/"
} else {
  stop("Computer not recognised")
}

io$metadata <- paste0(io$basedir,"/cell_metadata.txt.gz")
io$matrix <- paste0(io$basedir,"/data/matrix.mtx.gz")
io$features <- paste0(io$basedir,"/data/features.tsv.gz")
io$fragments <- paste0(io$basedir,"/data/fragments.tsv.gz")
io$peak.metadata <- paste0(io$basedir,"/peak_metadata.csv.gz")
io$barcodes <- paste0(io$basedir,"/data/barcodes.tsv.gz")
# io$seurat <- paste0(io$basedir,"/data/processed/seurat.rds")
io$cistopic <- paste0(io$basedir,"/results/cisTopic/cistopic.rds")

# archR stuff
io$archR.directory <- paste0(io$basedir,"/data/processed/archR")
# io$archR.projectMetadata <- paste0(io$archR.directory,"/projectMetadata.rds")
# io$archR.peakSet.granges <- paste0(io$archR.directory,"/PeakSet.rds")
io$archR.bgdPeaks <- paste0(io$archR.directory,"/Background-Peaks.rds")
io$archR.peakSet.bed <- paste0(io$archR.directory,"/PeakCalls/bed/peaks_archR_macs2.bed.gz")
# io$archR.peak.variability <- paste0(io$basedir,"/results/atac/archR/variability/peak_variability.txt.gz")
io$archR.peak.differential.dir <- paste0(io$basedir,"/results/differential/PeakMatrix")
# io$archR.peak.metadata <- paste0(io$archR.directory,"/PeakCalls/peak_metadata.tsv.gz")
# io$archR.peak.stats <- paste0(io$basedir,"/results/atac/archR/peak_calling/peak_stats.txt.gz")
# io$archR.peak2gene.all <- paste0(io$basedir,"/results/atac/archR/peak_calling/peaks2genes/peaks2genes_all.txt.gz")
# io$archR.peak2gene.nearest <- paste0(io$basedir,"/results/atac/archR/peak_calling/peaks2genes/peaks2genes_nearest.txt.gz")
# io$archr.chromvar.dir <- paste0(io$basedir,"/results/atac/archR/chromvar")
io$archR.pseudobulk.peakMatrix.se <- paste0(io$archR.directory,"/pseudobulk/pseudobulk_PeakMatrix_summarized_experiment.rds")
io$archR.pseudobulk.TileMatrix.se <- paste0(io$archR.directory,"/pseudobulk/pseudobulk_TileMatrix_summarized_experiment.rds")
io$archR.pseudobulk.GeneScoreMatrix.se <- paste0(io$archR.directory,"/pseudobulk/pseudobulk_GeneScoreMatrix_summarized_experiment.rds")
io$archR.pseudobulk.deviations.se <- sprintf("%s/pseudobulk/pseudobulk_DeviationMatrix_Motif_cisbp_summarized_experiment.rds",io$archR.directory)


# RNA atlas (PijuanSala2019)
io$rna.atlas.metadata <- paste0(io$atlas.basedir,"/sample_metadata.txt.gz")
io$rna.atlas.marker_genes <- paste0(io$atlas.basedir,"/results/marker_genes/all_stages/marker_genes.txt.gz")
io$rna.atlas.differential <- paste0(io$atlas.basedir,"/results/differential")
io$rna.atlas.average_expression_per_celltype <- paste0(io$atlas.basedir,"/results/marker_genes/all_stages/avg_expr_per_celltype_and_gene.txt.gz")
io$rna.atlas.sce <- paste0(io$atlas.basedir,"/processed/SingleCellExperiment.rds")

#############
## Options ##
#############

opts <- list()

opts$stages <- c("E8.25")

# opts$batches <- c("","")

opts$nuclei_type <- c("all", "2n", "4n")

opts$celltypes = c(
  "Surface_ectoderm",
  "Notochord",
  "Gut",
  "Cardiomyocytes",
  "Mid_Hindbrain",
  "Endothelium",
  "Paraxial_mesoderm",
  "Spinal_cord",
  "Somitic_mesoderm",
  "Erythroid",
  "Neural_crest",
  "Mixed_mesoderm",
  "NMP",
  "Forebrain",
  "ExE_endoderm",
  "Allantois",
  "Mesenchyme",
  "Pharyngeal_mesoderm"
)

opts$celltype.colors = c(
  "Notochord" = "#0F4A9C",
  "Gut" = "#EF5A9D",
  "Mixed_mesoderm" = "#DFCDE4",
  "Paraxial_mesoderm" = "#8DB5CE",
  "Somitic_mesoderm" = "#005579",
  "Pharyngeal_mesoderm" = "#C9EBFB",
  "Cardiomyocytes" = "#B51D8D",
  "Allantois" = "#532C8A",
  "Mesenchyme" = "#cc7818",
  "Endothelium" = "#ff891c",
  "Erythroid" = "#EF4E22",
  "NMP" = "#8EC792",
  "Neural_crest" = "#C3C388",
  "Forebrain" = "#647a4f",
  "Mid_Hindbrain" = "black",
  "Spinal_cord" = "#CDE088",
  "Surface_ectoderm" = "#f7f79e",
  "ExE_endoderm" = "#7F6874"
)

opts$rna.celltype.colors = c(
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
	"Blood_progenitors_1" = "#f9decf",
	"Blood_progenitors_2" = "#c9a997",
	"Erythroid1" = "#C72228",
	"Erythroid2" = "#f79083",
	"Erythroid3" = "#EF4E22",
	"NMP" = "#8EC792",
	"Rostral_neurectoderm" = "#65A83E",
	"Caudal_neurectoderm" = "#354E23",
	"Neural_crest" = "#C3C388",
	"Forebrain_Midbrain_Hindbrain" = "#647a4f",
	"Spinal_cord" = "#CDE088",
	"Surface_ectoderm" = "#f7f79e",
	"Visceral_endoderm" = "#F6BFCB",
	"ExE_endoderm" = "#7F6874",
	"ExE_ectoderm" = "#989898",
	"Parietal_endoderm" = "#1A1A1A",
	
	# Additional
	"Erythroid" = "#EF4E22",
	"Blood_progenitors" = "#c9a997",
	"Neurectoderm" = "#65A83E"
)

###################
## Load metadata ##
###################

sample_metadata <- fread(io$metadata) %>%
  .[,cell:=paste("E8.5_gastrulation",cell,sep="#")] %>%
  .[,celltype:=stringr::str_replace_all(celltype," ","_")] %>%
  .[,celltype:=stringr::str_replace_all(celltype,"/","_")] 