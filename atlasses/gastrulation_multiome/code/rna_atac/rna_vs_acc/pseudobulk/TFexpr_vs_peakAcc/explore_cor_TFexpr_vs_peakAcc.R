
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
io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc")

###############
## Load data ##
###############

io$file <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc/cor_TFexpr_vs_peakAcc_SummarizedExperiment.rds")
tf2peak_cor.se <- readRDS(io$file)

io$file <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc/cor_TFexpr_vs_peakAcc.txt.gz")
cor_dt <- fread(io$file) %>%
  .[!is.na(cor)] %>%
  .[,cor_sign:=c("-","+")[(cor>0)+1]]

#####################
## Load PeakMatrix ##
#####################

# peakMatrix.se <- readRDS(io$archR.pseudobulk.peakMatrix.se)

# Load peak metadata
peak_metadata.dt <- fread(io$archR.peak.metadata) %>% 
  .[,peak:=sprintf("%s:%s-%s",chr,start,end)]

# Define peak names
# peak_names <- rowData(peakMatrix.se) %>% as.data.table %>% .[,idx:=sprintf("%s:%s-%s",seqnames,start,end)] %>% .$id
# rownames(peakMatrix.se) <- peak_names


###############################
## Load motifmatcher results ##
###############################

motifmatcher.se <- readRDS(sprintf("%s/Annotations/Motif_cisbp-Matches-In-Peaks.rds",io$archR.directory))
colnames(motifmatcher.se) <- colnames(motifmatcher.se) %>% toupper %>% stringr::str_split(.,"_") %>% map_chr(1)
motifmatcher.se <- motifmatcher.se[,!duplicated(colnames(motifmatcher.se))]
tmp <- rowRanges(motifmatcher.se)
rownames(motifmatcher.se) <- sprintf("%s:%s-%s",seqnames(tmp), start(tmp), end(tmp))
motifmatcher.se <- motifmatcher.se[unique(cor_dt$peak),]

#############################################
## Shared TFs between blood and ExE tissue ##
#############################################

io$archR.pseudobulk.deviations.se <- "/Users/ricard/data/gastrulation_multiome_10x/results/atac/archR/chromvar/pseudobulk/chromVAR_deviations_summarized_experiment_Motif_cisbp_pseudobulk_correlated_peaks.rds"
source("/Users/ricard/gastrulation_multiome_10x/rna_atac/rna_vs_chromvar/pseudobulk/load_rna_chromvar_pseudobulk.R")


rna_chromvar.dt <- merge(
  rna_dt,
  chromvar_dt,
  by = c("celltype","gene")
)

length(unique(rna_chromvar.dt$gene))

ggscatter(celltype=="Blood_")