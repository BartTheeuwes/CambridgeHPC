#####################
## Define settings ##
#####################

# load default setings
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
# io$metadata <- paste0(io$basedir,"/sample_metadata.txt.gz")
io$outdir <- paste0(io$basedir,"/results/atac/archR/sequence_conservation")

# Options

##########################
## Load peak annotation ##
##########################

# peakSet.gr <- readRDS(io$archR.peakSet.granges)
# io$archR.peak.metadata <- paste0(io$archR.directory,"/PeakCalls/all_peaks/peak_metadata.tsv.gz")
peakSet.dt <- fread(io$archR.peak.metadata) %>%
  .[,peak:=sprintf("%s_%s_%s",chr,start,end)] %>%
  .[,peakType:=factor(peakType,levels=c("Promoter","Intronic","Exonic","Distal"))]

# Filter
# peakSet.gr.filt <- peakSet.gr[peakSet.gr$score>opts$min.score,]

# Load peak stats
io$archR.peak.stats <- paste0(io$basedir,"/results/atac/archR/peak_calling/peak_stats_binarised.txt.gz")
peakStats.dt <- fread(io$archR.peak.stats)

################################
## Load sequence conservation ##
################################

io$phastcons <- "/Users/ricard/data/gastrulation_multiome_10x/results/atac/archR/sequence_conservation/phastcons_score_peaks.tsv.gz"
phastcons.dt <- fread(io$phastcons) %>% setnames("peak_name","peak") %>% 
  .[,peak:=sub(":","_",peak)] %>%
  
mean(phastcons.dt$peak%in%peakSet.dt$peak)

###############################
## Plot mean versus variance ##
###############################

# to.plot <- peakSet.dt[chr=="chr1",c("peak","peakType")] %>% 
#   merge(peakStats.dt,by="peak")
# 
# p <- ggscatter(to.plot, x="var_singlecell", y="var_pseudobulk", size=1) +
#   labs(x="Variance (single-cell)", y="Variance (pseudobulk)") +
#   theme(
#     axis.text = element_text(size=rel(0.7))
#   )
# 
# pdf(paste0(io$outdir,"/foo.pdf"), width = 6, height = 5)
# print(p)
# dev.off()
