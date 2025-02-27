here::i_am("rna_atac/virtual_chipseq_library/trajectories/create_virtual_chipseq_library_trajectory.R")

suppressPackageStartupMessages(library(furrr))

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
# p$add_argument('--trajectory',   type="character",    help='File with the trajectory')
p$add_argument('--trajectory_name',   type="character",    help='Name of the trajectory')
p$add_argument('--tf2peak_cor',  type="character",              help='Correlations between TF RNA expression and peak accessibility') 
# p$add_argument('--atac_peak_stats',  type="character",              help='ATAC Peak stats') 
p$add_argument('--atac_peak_matrix_pseudobulk',  type="character",              help='ATAC Peak matrix (pseudobulk)') 
p$add_argument('--peak_metadata',             type="character",        help='Peak metadata file')
# p$add_argument('--motif2gene',  type="character",              help='Motif annotation')
p$add_argument('--motifmatcher',  type="character",              help='Motif annotation') 
p$add_argument('--motif_annotation',  type="character",              help='Motif annotation') 
p$add_argument('--max_acc_score',     type="double",    default=1.5,    help='Maximum peak accessibility')
p$add_argument('--min_number_peaks',     type="integer",    default=50,    help='Minimum number of peaks per TF')
p$add_argument('--min_peak_score',     type="integer",    default=15,    help='Minimum peak score')
p$add_argument('--outdir',          type="character",                help='Output directory')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')
p$add_argument('--test_mode',  action="store_true",  help='Test mode?')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
args <- list()
args$trajectory_name <- "NMP"
args$motif_annotation <- "CISBP"
args$tf2peak_cor <- file.path(io$basedir,sprintf("results/rna_atac/rna_vs_acc/trajectories/%s/TFexpr_vs_peakAcc_%s_%s.rds",args$trajectory_name,args$motif_annotation,args$trajectory_name))
# args$atac_peak_stats <- file.path(io$basedir,"results/atac/archR/feature_stats/PeakMatrix_stats.txt.gz")
args$atac_peak_matrix_pseudobulk <- file.path(io$basedir,"results/atac/archR/pseudobulk/celltype.mapped_mnn/pseudobulk_PeakMatrix_summarized_experiment.rds")
args$peak_metadata <- file.path(io$archR.directory,"PeakCalls/peak_metadata.tsv.gz")
args$max_acc_score <- 1.5
args$min_number_peaks <- 50
args$min_peak_score <- 15
args$motifmatcher <- sprintf("%s/Annotations/%s-Scores.rds",io$archR.directory,args$motif_annotation)
# args$motif2gene <- sprintf("%s/Annotations/%s_TFs.txt.gz",io$archR.directory,args$motif_annotation)
args$outdir <- file.path(io$basedir,sprintf("results/rna_atac/virtual_chipseq/%s",args$motif_annotation))
args$threads <- 1
args$test_mode <- TRUE
## END TEST ##


# I/O
dir.create(args$outdir, showWarnings = F)

# Options
opts$celltype_trajectory_dic <- list(
  "blood" = c("Haematoendothelial_progenitors", "Blood_progenitors_1", "Blood_progenitors_2", "Erythroid1", "Erythroid2", "Erythroid3"),
  "ectoderm" = c("Epiblast", "Rostral_neurectoderm", "Forebrain_Midbrain_Hindbrain"),
  "endoderm" = c("Epiblast", "Anterior_Primitive_Streak", "Def._endoderm", "Gut"),
  "mesoderm" = c("Epiblast", "Primitive_Streak", "Nascent_mesoderm"),
  "NMP" = c("Spinal_cord", "NMP", "Somitic_mesoderm", "Caudal_mesoderm")
)

# Parallel processing
# plan(sequential)
plan(multisession, workers=args$threads)

###############################
## Load TF2peak correlations ##
###############################

tf2peak_cor.se <- readRDS(args$tf2peak_cor)

###############
## Load ATAC ##
###############

# Load pseudobulk ATAC peak matrix
atac_PeakMatrix_pseudobulk.se <- readRDS(args$atac_peak_matrix)

# Subset cell types
atac_PeakMatrix_pseudobulk.se <- atac_PeakMatrix_pseudobulk.se[,opts$celltype_trajectory_dic[[args$trajectory_name]]]

# Subset peaks
stopifnot(rownames(tf2peak_cor.se)%in%rownames(atac_PeakMatrix_pseudobulk.se))
atac_PeakMatrix_pseudobulk.se <- atac_PeakMatrix_pseudobulk.se[rownames(tf2peak_cor.se),]

# Load peak metadata
atac_peak_metadata.dt <- fread(args$peak_metadata) %>%
  .[score>=args$min_peak_score] %>%
  .[,idx:=sprintf("%s:%s-%s",chr,start,end)] %>%
  .[idx%in%rownames(tf2peak_cor.se)]

######################
## Load motifmatchR ##
######################

# source(here::here("load_motifmatchR.R"))
motifmatcher.se <- readRDS(args$motifmatcher)[atac_peak_metadata.dt$idx,]

# Subset peaks
# stopifnot(sort(rownames(motifmatcher.se))==sort(rownames(atac_pseudobulk_PeakMatrix.se)))
# motifmatcher.se <- motifmatcher.se[rownames(atac_pseudobulk_PeakMatrix.se),]

###########################
## Load motif annotation ##
###########################

opts$motif_annotation <- args$motif_annotation
source(here::here("atac/archR/load_motif_annotation.R"))

motif2gene.dt <- motif2gene.dt[gene%in%colnames(tf2peak_cor.se)]
# (TO FIX) Make sure that there is a one-to-one match between TFs and genes
# tmp <- motif2gene.dt[,.N,by=c("gene")] %>% .[N>1]
# motif2gene.dt[gene=="TFAP2B"]
# motif2gene.dt[gene=="TFAP2A"]
foo <- motif2gene.dt[,N:=.N,by=c("gene")] %>% .[N==1] %>% .[,N:=NULL] 
bar <- motif2gene.dt[,N:=.N,by=c("gene")] %>% .[N>1] %>% .[,N:=NULL] %>% .[,head(.SD,n=1), by="gene"]
motif2gene.dt <- rbind(foo,bar)

##################
## Prepare data ##
##################

peaks <- unique(atac_peak_metadata.dt$idx)
motifs <- motif2gene.dt$motif
TFs <- motif2gene.dt$gene

tf2peak_cor.mtx <- assay(tf2peak_cor.se,"cor")[peaks,TFs]
motifmatcher.mtx <- assay(motifmatcher.se,"motifScores")[peaks,motifs]
atac_PeakMatrix_pseudobulk.se <- assay(atac_PeakMatrix_pseudobulk.se[peaks,]) %>% round(3)

######################################
## Create virtual chip-seq library ##
######################################

print("Predicting TF binding sites...")

if (args$test_mode) {
  print("Testing mode activated...")
  TFs <- TFs %>% head(n=3)
} 

# sanity checks
stopifnot(motif2gene.dt$motif %in% colnames(motifmatcher.mtx))

# i <- "GATA1"
# virtual_chip.dt <- TFs %>% future_map(function(i) {
virtual_chip.dt <- TFs %>% map(function(i) {
# for (i in TFs) {
  print(i)
  motif <- motif2gene.dt[gene==i,motif]
  stopifnot(length(motif)==1)

  # peaks <- names(which(abs(assay(tf2peak_cor.se[peaks,i],"cor")[,1])>=0))
  peaks <- names(which(abs(tf2peak_cor.mtx[,i])>0)) # we only consider chromatin activators
  
  if (length(peaks)>=args$min_number_peaks) {
    
    # calculate accessibility score
    max_accessibility_score <- apply(atac_PeakMatrix_pseudobulk.se[peaks,],1,max) %>% round(2)
    max_accessibility_score[max_accessibility_score>=args$max_acc_score] <- args$max_acc_score
    
    # calculate correlation score
    # correlation_score <- assay(tf2peak_cor.se[,i],"cor")[peaks,1] %>% round(2)
    correlation_score <- tf2peak_cor.mtx[peaks,i] %>% round(2)
    # correlation_score[correlation_score==0] <- NA
    
    # calculate motif score
    # motif_score <- assay(motifmatcher.se[peaks,motif2gene.dt[gene==i,motif]],"motifScores")[,1]
    motif_score <- motifmatcher.mtx[peaks,motif]
    motif_score <- round(motif_score/max(motif_score),2)
    
    # calculate motif counts
    # motif_counts <- assay(motifmatcher.se[peaks,i],"motifCounts")[,1] %>% round(2)
    
    # predicted_score <- minmax.normalisation(max_accessibility_score * correlation_score * motif_score * motif_counts) 
    predicted_score <- correlation_score * minmax.normalisation(max_accessibility_score * motif_score)
    # predicted_score <- predicted_score/max(predicted_score,na.rm=T) %>% round(2)
    
    tmp <- data.table(
      peak = peaks, 
      correlation_score = correlation_score,
      max_accessibility_score = max_accessibility_score,
      motif_score = motif_score,
      # motif_counts = motif_counts,
      score = round(predicted_score,2)
    ) %>% sort.abs("score") %>% 
      # .[motif_score==0,motif_score:=NA] %>%
      .[,c("peak","score","correlation_score","max_accessibility_score","motif_score")]
    
    bed.dt <- tmp %>%
      .[,peak:=str_replace(peak,":","-")] %>%
      .[,chr:=strsplit(peak,"-") %>% map_chr(1)] %>%
      .[,start:=strsplit(peak,"-") %>% map_chr(2)] %>%
      .[,end:=strsplit(peak,"-") %>% map_chr(3)] %>%
      .[,c("chr","start","end","score")]
      
    # Save 
    fwrite(tmp, sprintf("%s/%s.txt.gz",args$outdir,i), sep="\t", quote=F, col.names = T)
    fwrite(bed.dt, sprintf("%s/%s.bed.gz",args$outdir,i), sep="\t", quote=F, col.names = F)
    
    to_return.dt <- tmp[!is.na(score),c("peak","score")] %>% .[,tf:=i]
    return(to_return.dt)
  }
}) %>% rbindlist


########################
## Save sparse matrix ##
########################

virtual_chip.mtx <- virtual_chip.dt %>% 
  data.table::dcast(peak~tf, value.var="score", fill=0) %>%
  matrix.please %>% Matrix::Matrix(.)

saveRDS(virtual_chip.mtx, file.path(args$outdir,"virtual_chip.mtx"))

###################################################################
## Update motifmatchr results using the virtual ChIP-seq library ##
###################################################################

print("Updating motifmatchr results using the virtual ChIP-seq library...")

# virtual_chip.mtx <- readRDS(file.path(args$outdir,"virtual_chip_mtx.rds"))

opts$min_chip_score <- 0.15

# i <- "CLOCK"
for (i in colnames(virtual_chip.mtx)) {
  
  motif <- motif2gene.dt[gene==i,motif]
  
  # peaks <- which(virtual_chip.mtx[,i]>0) %>% names
  peaks <- which(abs(virtual_chip.mtx[,i])>=opts$min_chip_score) %>% names
  # print(sprintf("%s: %d/%d",i,length(peaks),sum(assay(motifmatcher.se,"motifMatches")[,i])))
  
  assay(motifmatcher.se,"motifMatches")[,motif][!rownames(assay(motifmatcher.se))%in%peaks] <- FALSE
  assay(motifmatcher.se,"motifScores")[,motif][!rownames(assay(motifmatcher.se))%in%peaks] <- 0
  assay(motifmatcher.se,"motifCounts")[,motif][!rownames(assay(motifmatcher.se))%in%peaks] <- 0
}

saveRDS(motifmatcher.se, file.path(args$outdir,"motifmatchr_virtual_chip.rds"))
