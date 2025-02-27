
#####################
## Define settings ##
#####################

# Load default settings
if (grepl("ricard",Sys.info()['nodename'])) {
  # source("/Users/ricard/scnmt_gastrulation/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/other/chip/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  # source("/homes/ricard/scnmt_gastrulation/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/other/chip/utils.R")
} else {
  stop("Computer not recognised")
}


# I/O

# Options

##############################
## Load genomic annotations ##
##############################

io$file <- "/Users/ricard/data/mm10_regulation/TF_ChIP/TAL1_blood_progenitors/GSM1692843_HB_Tal1.bed.gz"
granges <- io$file %>% 
  load_annotation_as_granges_list %>% .[[1]]

############################
## Load motif annotations ##
############################

io$pwm <- "/Users/ricard/data/gastrulation_multiome_10x/results/atac/archR/motif_seqlogo/PWMatrixList.rds"
pwms <- readRDS(io$pwm)

pwms2 <- load_motif_annotation(
  motifSet = "cisbp",
  species = "Mus musculus"
)[["motifs"]]

# Rename motifs
# names(pwms) <- names(pwms) %>% toupper %>% stringr::str_split(.,"_") %>% map_chr(1)
# names(pwms) <- gsub("TCFAP","TFAP",names(pwms))
# names(pwms) <- gsub("NKX2","NKX2-",names(pwms))
# names(pwms) <- gsub("NKX3","NKX3-",names(pwms))
# names(pwms) <- gsub("NKX6","NKX6-",names(pwms))

# Remove duplicated motifs
# pwms <- pwms[,!duplicated(names(pwms))]

###################
## Subset motifs ##
###################

pwms.filt <- pwms["TAL1"]

#################################
## Create motif overlap matrix ##
#################################

motifOverlap.se <- create_motif_overlap_matrix(
  pwms.filt, 
  granges, 
  genome = "BSgenome.Mmusculus.UCSC.mm10", 
  cutOff = 5e-05, 
  width = 7
)[["motifMatches"]]

rownames(motifOverlap.se) <- sprintf("%s:%s-%s",seqnames(granges), start(granges), end(granges))

mean(assay(motifOverlap.se)[,1])
names(which(assay(motifOverlap.se)[,1]==F))
names(which(assay(motifOverlap.se)[,1]==T))

assay(motifOverlap.se)[,1]["chr1:38767144-38767544"]

##########################################
## Plot fraction of peaks with TF motif ##
##########################################

#############
## Explore ##
#############

peak_metadata.dt <- fread(io$archR.peak.metadata) %>% .[,peak:=sprintf("%s:%s-%s",chr,start,end)]

peak_metadata.dt[peak=="chr18:47959154-47959554"]
tf2peak_cor.se <- readRDS(io$tf2peak_cor.se)[,"TAL1"]

assay(tf2peak_cor.se["chr18:47959154-47959554"],"cor")
pwms["TAL1"]$TAL1
source("/Users/ricard/gastrulation_multiome_10x/load_motifmatchR.R")
source("/Users/ricard/gastrulation_multiome_10x/load_motifmatchR_positions.R")
assay(motifmatcher.se["chr1:38766888-38767488",])[1,]["TAL1"]

library(GenomicRanges)
gr <- GRanges(
  seqnames = Rle(c("chr18"), c(1)),
  ranges = IRanges(47959154, 47959554),
  strand = Rle(strand(c("-")), c(1))
)

ov <- findOverlaps(query=gr, subject=motifmatcher_positions.se[["TAL1"]])
ov <- findOverlaps(query=gr, subject=motifmatcher_positions.se[["TAL1"]])
motifmatcher_positions.se[["TAL1"]][subjectHits(ov)]

########################################################
## Create BED file with the motif matches coordinates ##
########################################################

names(motifOverlap.se$motifPositions) %>% walk(function(i) {
  motif_matches.dt <- as.data.table(motifOverlap.se$motifPositions[[i]]) %>% 
    .[,c(1,2,3,5)] %>% 
    setnames(c("chr","start","end","strand")) %>%
    .[,anno:=sprintf("%s_%s",args$anno,i)] %>%
    .[,chr:=gsub("chr","",chr)] %>%
    .[,c("start","end"):=list(start-args$extend,end+args$extend)] %>%
    .[,id:=sprintf("%s:%s-%s",chr,start,end)] %>%
    .[,c("chr","start","end","strand","id","anno")] %>% 
    setorder(chr,start,end)
  
    fwrite(motif_matches.dt, sprintf("%s_%s.bed.gz",args$outprefix,i), sep="\t", quote=F, col.names = F)
})# %>% rbindlist %>% setorder(chr,start,end)


##########
## Save ##
##########

saveRDS(motifOverlap.se, paste0(args$outprefix,".rds"))



