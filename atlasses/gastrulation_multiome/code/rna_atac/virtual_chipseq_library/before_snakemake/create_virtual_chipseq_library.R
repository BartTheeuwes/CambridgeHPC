suppressPackageStartupMessages(library(furrr))

#####################
## Define settings ##
#####################

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

# I/O
io$outdir <- file.path(io$basedir,"results/rna_atac/virtual_chipseq")#; dir.create(io$outdir, showWarnings = F)

# opts$motif_annotation <- "Motif_cisbp"
opts$motif_annotation <- "Motif_cisbp_lenient"

# Parallel processing
# plan(sequential)
plan(multisession, workers=4)

###############################
## Load TF2peak correlations ##
###############################

io$tf2peak_cor.se <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc/cor_TFexpr_vs_peakAcc_SummarizedExperiment_lenient_v2.rds")
tf2peak_cor.se <- readRDS(io$tf2peak_cor.se)
# assay(tf2peak_cor.se,"cor") <- dropNA2matrix(assay(tf2peak_cor.se,"cor"))

# i <- "chr1:3035578-3036178"
# assay(tf2peak_cor.se[i,],"cor")[1,]["SOX21"]

#######################################
## Load pseudobulk RNA and ATAC data ##
#######################################

# io$archR.pseudobulk.deviations.se <- sprintf("%s/pseudobulk/pseudobulk_DeviationMatrix_%s_summarized_experiment.rds",io$archR.directory,opts$motif_annotation)
source(here::here("rna_atac/load_rna_atac_pseudobulk.R"))

######################
## Load motifmatchR ##
######################

source(here::here("load_motifmatchR.R"))

#################
## Filter data ##
#################

opts$min.peak_score <- 50

peaks <- peak_metadata.dt[score>=opts$min.peak_score,idx]

tf2peak_cor.mtx <- assay(tf2peak_cor.se,"cor")[peaks,]
motifmatcher.mtx <- assay(motifmatcher.se,"motifScores")[peaks,]
atac_pseudobulk_peakMatrix.se <- atac_pseudobulk_peakMatrix.se[peaks,]

######################################
## Create virtual chip-seq library ##
######################################

opts$max.acc.score <- 1.5
opts$min_number_peaks <- 50

opts$TFs <- colnames(tf2peak_cor.se)
# opts$TFs <- c("T", "ZIC2", "SOX2", "TAL1", "GATA1", "FOXA2", "GATA4", "NKX2-5")

# i <- "CLOCK"
opts$TFs %>% future_walk(function(i) {
# foreach(i=opts$TFs[1:16]) %dopar% {       # stopCluster(cl)
# for (i in opts$TFs) {
  # print(i)
  # peaks <- names(which(abs(assay(tf2peak_cor.se[peaks,i],"cor")[,1])>=0))
  peaks <- names(which(abs(tf2peak_cor.mtx[,i])>0))
  
  if (length(peaks)>=opts$min_number_peaks) {
    
    # calculate accessibility score
    max_accessibility_score <- apply(assay(atac_pseudobulk_peakMatrix.se[peaks,]),1,max) %>% round(2)
    max_accessibility_score[max_accessibility_score>=opts$max.acc.score] <- opts$max.acc.score
    
    # calculate correlation score
    # correlation_score <- assay(tf2peak_cor.se[,i],"cor")[peaks,1] %>% round(2)
    correlation_score <- tf2peak_cor.mtx[peaks,i] %>% round(2)
    correlation_score[correlation_score==0] <- NA
    
    # calculate motif score
    # motif_score <- assay(motifmatcher.se[peaks,i],"motifScores")[,1]
    motif_score <- motifmatcher.mtx[peaks,i]
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
    # ) %>% setorder(-score,na.last = T) %>% 
    ) %>% sort.abs("score") %>% 
      .[motif_score==0,motif_score:=NA] %>%
      .[,peak:=str_replace(peak,":","-")] %>%
      .[,chr:=strsplit(peak,"-") %>% map_chr(1)] %>%
      .[,start:=strsplit(peak,"-") %>% map_chr(2)] %>%
      .[,end:=strsplit(peak,"-") %>% map_chr(3)] %>%
      .[,peak:=NULL] %>%
      .[,c("chr","start","end","score","correlation_score","max_accessibility_score","motif_score")]
    
    # Save 
    # tmp[chr=="chr17" & start==64600366]
    fwrite(tmp, sprintf("%s/%s.txt.gz",io$outdir,i), sep="\t", quote=F, col.names = T)
    # fwrite(tmp[score>=0.10,c("chr","start","end","score")], sprintf("%s/%s_filt.bed.gz",io$outdir,i), sep="\t", quote=F, col.names = F)
    fwrite(tmp[!is.na(score),c("chr","start","end","score")], sprintf("%s/%s.bed.gz",io$outdir,i), sep="\t", quote=F, col.names = F)
  }
})


########################
## Save sparse matrix ##
########################

# opts$TFs <- c("FOXA2","MIXL1","TAL1")
opts$TFs <- list.files(io$virtual_chip.dir, ".bed.gz") %>% str_replace_all(".bed.gz","")

virtual_chip.dt <- opts$TFs %>% map(function(i) {
  fread(sprintf("%s/%s.bed.gz",io$virtual_chip.dir,i)) %>%
    setnames(c("chr","start","end","score")) %>%
    .[,idx:=sprintf("%s:%s-%s",chr,start,end)] %>%
    .[,c("chr","start","end"):=NULL] %>% .[,tf:=i] %>%
    return
}) %>% rbindlist

virtual_chip.mtx <- virtual_chip.dt %>% 
  dcast(idx~tf, value.var="score", fill=0) %>%
  matrix.please %>% Matrix::Matrix(.)

saveRDS(virtual_chip.mtx, file.path(io$outdir,"virtual_chip.mtx"))
virtual_chip.mtx <- readRDS(file.path(io$outdir,"virtual_chip_mtx.rds"))

###################################################################
## Update motifmatchr results using the virtual ChIP-seq library ##
###################################################################

opts$min_chip_score <- 0.1

# i <- "CLOCK"
for (i in colnames(virtual_chip.mtx)) {
  
  # peaks <- which(virtual_chip.mtx[,i]>0) %>% names
  peaks <- which(abs(virtual_chip.mtx[,i])>=opts$min_chip_score) %>% names
  # print(sprintf("%s: %d/%d",i,length(peaks),sum(assay(motifmatcher.se,"motifMatches")[,i])))
  
  assay(motifmatcher.se,"motifMatches")[,i][!rownames(assay(motifmatcher.se))%in%peaks] <- FALSE
  assay(motifmatcher.se,"motifScores")[,i][!rownames(assay(motifmatcher.se))%in%peaks] <- 0
  assay(motifmatcher.se,"motifCounts")[,i][!rownames(assay(motifmatcher.se))%in%peaks] <- 0
}

saveRDS(motifmatcher.se, file.path(io$outdir,"motifmatchr_virtual_chip.rds"))

##########
## Test ##
##########

# tf2peak_cor.dt[peak=="chr11:54220015-54220615"]
# peak_metadata.dt[idx=="chr8:90822247-90822847"]
# gene <- "TAL1"
# peak <- "chr4:100548071-100548671"
# 
# peak_metadata.dt[idx=="chr4:100548071-100548671"]
# 
# which(assay(motifmatcher.se[peak,])[1,])
# assay(tf2peak_cor.se[peak,gene])
# 
# to.plot <- data.table(
#   celltype = opts$celltypes,
#   # rna = logcounts(rna.sce.tf[gene])[1,],
#   rna = logcounts(rna.sce[str_to_title(gene)])[1,],
#   acc = assay(atac_pseudobulk_peakMatrix.se[peak,])[1,]
# )
# 
# 
# ggscatter(to.plot, x="rna", y="acc", fill="celltype", size=4, shape=21, 
#                 add="reg.line", add.params = list(color="black", fill="lightgray"), conf.int=TRUE) +
#   stat_cor(method = "pearson", label.x.npc = "middle", label.y.npc = "bottom") +
#   scale_fill_manual(values=opts$celltype.colors) +
#   labs(x=sprintf("%s expression",gene), y=sprintf("%s Peak accessibility",peak)) +
#   guides(fill=F) +
#   theme(
#     axis.text = element_text(size=rel(0.7))
#   )


##########
## TEST ##
##########

# tf.to.explore <- "GATA1"
# foo <- which( assay(tf2peak_cor.se[,tf.to.explore],"cor")[,1]!=0) %>% names
# bar <- which(assay(motifmatcher.se[,tf.to.explore],"motifMatches")[,1]==1) %>% names
# 
# peak.to.explore <- "chr1:170845679-170846279"
# foo <- which(assay(tf2peak_cor.se[peak.to.explore,],"cor")[1,]!=0) %>% names
# bar <- which(assay(motifmatcher.se[peak.to.explore],"motifMatches")[1,]==1) %>% names
# grep("GATA",bar)
# 
# assay(motifmatcher_new.se["chr2:153479346-153479946","GATA1"])
# assay(motifmatcher.se[predicted.binding.peaks,"GATA1"],"motifScores")[,1]
# assay(motifmatcher.se[predicted.binding.peaks,"GATA1"],"motifMatches")[,1]
# assay(motifmatcher.se[predicted.binding.peaks,"GATA1"],"motifCounts")[,1]
# 
# 
# peak.to.explore <- "chr3:87910471-87911071"
# assay(motifmatcher.se[peak.to.explore,"GATA1"],"motifScores")[,1]
# assay(motifmatcher.se[peak.to.explore,"GATA1"],"motifMatches")[,1]
# assay(motifmatcher.se[peak.to.explore,"GATA1"],"motifCounts")[,1]
# peak_metadata.dt[idx%in%peak.to.explore]