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
io$outdir <- paste0(io$basedir,"/results/atac/archR/test")

# Options
opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2"
  # "E8.0_rep1",
  # "E8.0_rep2",
  # "E8.5_rep1",
  # "E8.5_rep2"
)

opts$motif_annotation <- "Motif_cisbp"

#####################
## Update metadata ##
#####################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & doublet_call==FALSE] %>%
  .[sample%in%opts$samples]

##################
## Subset ArchR ##
##################

ArchRProject.filt <- ArchRProject[sample_metadata$cell,]
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

peakAnno <- getPeakAnnotation(ArchRProject.filt, name=opts$motif_annotation)
motif2peak.matches <- readRDS(peakAnno$Matches)

# Rename motifs
colnames(motif2peak.matches) <- colnames(motif2peak.matches) %>% toupper %>% stringr::str_split(.,"_") %>% map_chr(1)

# Remove duplicated motifs
motif2peak.matches <- motif2peak.matches[,!duplicated(colnames(motif2peak.matches))]

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

motif2gene.dt <- motif2gene.dt %>% .[gene%in%colnames(motif2peak.matches)]

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
  
####################################################
## Load coexpression matrix between TFs and genes ##
####################################################

cor_rna.mtx <- readRDS("/Users/ricard/data/gastrulation_multiome_10x/results/rna/coexpression/correlation_matrix_tf2gene.rds")

# Restrict to marker genes
marker_genes.dt <- fread(io$rna.atlas.marker_genes)
cor_rna.mtx <- cor_rna.mtx[,colnames(cor_rna.mtx)%in%marker_genes.dt$gene]

##########
## Test ##
##########

i <- "ETV2"

motifs <- colnames(motif2peak.matches)
tmp <- motif2peak.matches[as.logical(assay(motif2peak.matches[,"ETV2"])==1),"ETV2"]

cor_dt <- data.table(
  gene = colnames(cor_rna.mtx),
  cor = cor_rna.mtx[i,]
) %>% .[,sign:=c("Down","Up")[as.numeric(cor>0)+1]]

dt <- rowData(tmp)[,c("score","distToGeneStart","distToTSS","nearestGene","peakType","idx")] %>% as.data.table %>%
  setnames("nearestGene","gene") %>%
  .[,id:=sprintf("%s:%s-%s",seqnames(rowRanges(tmp)), start(rowRanges(tmp)), end(rowRanges(tmp)))] %>%
  merge(cor_dt,by="gene")# %>%
  # merge(marker_genes.dt[,c("celltype","gene")],by="gene", allow.cartesian=T)
  # .[toupper(gene)%in%motifs]

cor_dt <- data.table(
  gene = colnames(cor_rna.mtx),
  cor = cor_rna.mtx[i,]
) %>% .[,sign:=c("Down","Up")[as.numeric(cor>0)+1]]

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

################
## Test ArchR ##
################

R.utils::sourceDirectory("/Users/ricard/git/ArchR/R/", verbose=T, modifiedOnly=FALSE)
# R.utils::sourceDirectory("/homes/ricard/git/ArchR/R/", verbose=T, modifiedOnly=FALSE)

# peaks.bed <- fread(io$archR.peakSet.bed) %>% .[,foo:=sprintf("%s_%s_%s",V1,V2,V3)]
# foo <- sprintf("%s_%s_%s",seqnames(ArchRProject.filt@peakSet),start(ArchRProject.filt@peakSet),end(ArchRProject.filt@peakSet))

ArchRProj = ArchRProject.filt
useMatrix = "GeneScoreMatrix"
useSeqnames = NULL
verbose = TRUE
binarize = TRUE
threads = getArchRThreads()
logFile = createLogFile("getMatrixFromProject")
