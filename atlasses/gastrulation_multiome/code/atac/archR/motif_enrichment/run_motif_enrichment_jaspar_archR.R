# https://www.ArchRProject.com/bookdown/motif-and-feature-enrichment-with-archr.html

########################
## Load ArchR Project ##
########################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
} else {
  stop("Computer not recognised")
}

#####################
## Define settings ##
#####################

# I/O
# io$metadata <- paste0(io$basedir,"/processed/atac/archR/sample_metadata_after_archR.txt.gz")
io$outdir <- paste0(io$basedir,"/results/atac/archR/motif_enrichment")

# Options
opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E8.5_rep1",
  "E8.5_rep2"
)

#####################
## Update metadata ##
#####################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE] %>%
  .[sample%in%opts$samples]
stopifnot(sample_metadata$archR_cell %in% rownames(ArchRProject))

# subset celltypes with sufficient number of cells
opts$min.cells <- 100
sample_metadata <- sample_metadata %>%
  .[,N:=.N,by=c("celltype.predicted")] %>% .[N>opts$min.cells] %>% .[,N:=NULL]
table(sample_metadata$celltype.predicted)

##################
## Subset ArchR ##
##################

ArchRProject.filt <- ArchRProject[sample_metadata$archR_cell,]
table(getCellColData(ArchRProject.filt,"Sample")[[1]])
table(getCellColData(ArchRProject.filt,"celltype.predicted")[[1]])

##################
## Add Peak Set ##
##################

# Load Peak Set
# peaks.dt <- fread(io$archR.peakSet.bed) %>%
#   setnames(c("chr","start","end")) %>%
#   .[,id:=sprintf("%s:%s-%s",chr,start,end)]
# peaks.granges <- makeGRangesFromDataFrame(peaks.dt, keep.extra.columns = T)
# ArchRProject.filt <- addPeakSet(ArchRProject.filt, peaks.granges)

###################
## Sanity checks ##
###################

# Motif annotations require a PeakSet
getPeakSet(ArchRProject.filt) %>% head

###########################
## Add Motif annotations ##
###########################

opts$motif.pvalue.cutoff <- 5e-05  # default is 5e-05

# bar <- getPeakAnnotation(ArchRProject.filt, name="Motif_JASPAR")

# Load pre-computed motif annotation
# ArchRProject.filt@peakAnnotation <- readRDS(sprintf("%s/Annotations/JASPAR/peakAnnotation.rds",io$archR.directory))

# add motif set 
ArchRProject <- addMotifAnnotations(
  ArchRProject, 
  motifSet = "JASPAR",      
  collection = "CORE",  
  cutOff = opts$motif.pvalue.cutoff,   
  name = "Motif_JASPAR",
  force = TRUE
)
names(getPeakAnnotation(ArchRProject))

################################
## Calculate motif enrichment ##
################################

# The output of peakAnnoEnrichment() is a SummarizedExperiment object containing multiple assays 
# that store the results of enrichment testing with the hypergeometric test.

# upregulated TFs
opts$hypergeometric.cutoff <- "FDR <= 0.1 & Log2FC >= 0.5"

motifsUp <- peakAnnoEnrichment(ArchRProject,
  seMarker = markersPeaks,
  peakAnnotation = "Motif",
  cutOff = opts$hypergeometric.cutoff
)
names(assays(motifsUp))

# downregulated TFs
opts$hypergeometric.cutoff <- "FDR <= 0.1 & Log2FC <= -0.5"
motifsDown <- peakAnnoEnrichment(ArchRProject,
  seMarker = markersPeaks,
  peakAnnotation = "Motif",
  cutOff = opts$hypergeometric.cutoff
)
names(assays(motifsDown))

####################
## Prepare output ##
####################

dt.up <- assay(motifsUp) %>% as.data.frame %>% 
  as.data.table(keep.rownames = T) %>% setnames("rn","TF") %>%
  melt(id.vars=c("TF"), variable.name="celltype")

dt.down <- assay(motifsDown) %>% as.data.frame %>% 
  as.data.table(keep.rownames = T) %>% setnames("rn","TF") %>%
  melt(id.vars=c("TF"), variable.name="celltype")

##################
## Query output ##
##################

dt.up[celltype=="Erythroid"] %>% setorder(-value) %>% head(n=50) %>% View

##########
## Save ##
##########

outfile <- paste0(io$outdir,"/motif_enrichment_up.tsv.gz")
fwrite(dt.up, outfile, sep="\t")

outfile <- paste0(io$outdir,"/motif_enrichment_down.tsv.gz")
fwrite(dt.down, outfile, sep="\t")
