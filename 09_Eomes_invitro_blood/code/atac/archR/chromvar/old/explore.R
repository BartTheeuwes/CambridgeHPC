# https://www.archrproject.com/bookdown/chromvar-deviatons-enrichment-with-archr.html
library(GenomicRanges)

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
io$outdir <- paste0(io$basedir,"/results/atac/archR/TF_activity")

# Options
opts$samples <- c(
  # "E7.5_rep1",
  # "E7.5_rep2",
  # "E8.5_rep1",
  "E8.5_rep2"
)

#####################
## Update metadata ##
#####################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & pass_rnaQC==TRUE] %>%
  .[sample%in%opts$samples]
stopifnot(sample_metadata$archR_cell %in% rownames(ArchRProject))

##############################
## Load RNA expression data ##
##############################

# Load SingleCellExperiment
sce <- load_SingleCellExperiment(io$sce, cells = sample_metadata$cell, normalise = TRUE, remove_non_expressed_genes = TRUE)

# rename
all(colnames(sce) == sample_metadata$cell)
colnames(sce) <- sample_metadata$archR_cell

##################
## Subset ArchR ##
##################

ArchRProject.filt <- ArchRProject[sample_metadata$archR_cell,]
table(getCellColData(ArchRProject.filt,"Sample")[[1]])

#######################
## Fetch Peak Matrix ##
#######################

# Get Peak Matrix from ArchR
atac.peak.se <- getMatrixFromProject(ArchRProject.filt, useMatrix="PeakMatrix", binarize = T)
dim(atac.peak.se)

# Define peak names
peak_names <- rowRanges(atac.peak.se) %>% as.data.table %>% .[,id:=sprintf("%s:%s-%s",seqnames,start,end)] %>% .$id
rownames(atac.peak.se) <- peak_names

#############################
## Fetch motif2peak matrix ##
#############################

peakAnno <- getPeakAnnotation(ArchRProject.filt, name="Motif")
motif2peak.matches <- readRDS(peakAnno$Matches)

dim(motif2peak.matches)
motifs <- colnames(motif2peak.matches)

################################
## Load motif2gene annotation ##
################################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/atac/archR/load_motif_annotation.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/atac/archR/load_motif_annotation.R")
} else {
  stop("Computer not recognised")
}

motif2gene.dt <- motif2gene.dt %>%
  .[motif%in%motifs] %>%
  .[gene%in%rownames(sce)] %>%
  .[,N:=length(unique(motif)),by="gene"] %>% .[N==1] %>% .[,N:=NULL]

sce <- sce[motif2gene.dt$gene,]

# opts$min.motif.score <- 7

# motif2peak.positions <- readRDS(peakAnno$Positions)

# motif2peak.positions.dt <- motifs %>% map(function(i) {
#   motif2peak.positions[[i]] %>%
#     as.data.table() %>%
#     setnames("seqnames","chr") %>%
#     .[,motif:=i] %>%
#     return
# }) %>% rbindlist %>% 
#   # .[,idx:=sprintf("%s:%s-%s",chr,start,end)] %>%
#   .[,motif:=factor(motif,levels=motifs)] %>%
#   .[,chr:=factor(chr,levels=opts$chr)]
  
## 
  
motif2peak.matches.dt <- motifs %>% head %>% map(function(i) {
  print(i)
  tmp <- motif2peak.matches[,i]
  tmp <- tmp[which(assay(tmp)[,1])]
  dt <- rowRanges(tmp) %>% as.data.table() %>% 
    .[,c("seqnames","start","end","score")] %>%
    setnames(c("seqnames","score"),c("chr","peak_score")) %>%
    .[,motif:=i] %>%
    return
}) %>% rbindlist %>%
  # .[,idx:=sprintf("%s:%s-%s",chr,start,end)] %>%
  .[,motif:=factor(motif,levels=motifs)] %>%
  .[,chr:=factor(chr,levels=opts$chr)]

motif2peak.matches.dt[,.N,by="motif"]

rowRanges(matches)
colnames(matches)

i = motifs[1]
gene <- motif2gene.dt[motif==i,gene]

# for (i in )
tmp <- motif2peak.matches[,i][,1]
tmp <- tmp[which(assay(tmp)[,1])]
dt <- rowRanges(tmp) %>% as.data.table() %>% 
  .[,c("seqnames","start","end","score")] %>%
  setnames(c("seqnames","score"),c("chr","peak_score")) %>%
  .[,idx:=sprintf("%s:%s-%s",chr,start,end)]

peak.accessibility <- assay(atac.peak.se[dt$idx,]) %>% as.matrix
tf.expr <- logcounts(sce[gene,])[1,]
cells <- intersect(names(tf.expr),colnames(peak.accessibility))

cor(tf.expr,peak.accessibility)
# calculate matrix of correlation coefficients per TF and peak with the corresponding motif