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
io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/trajectories/test")
io$trajectories.inputdir <- c(
  # "ectoderm" = paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/trajectories/ectoderm_trajectory"),
  # "endoderm" = paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/trajectories/endoderm_trajectory"),
  # "mesoderm" = paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/trajectories/mesoderm_trajectory")
  "blood" = paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/trajectories/blood_trajectory_knn50")
)

# Options
opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E8.0_rep1",
  "E8.0_rep2",
  "E8.5_rep1",
  "E8.5_rep2"
)

opts$motif_annotation <- "Motif_cisbp"

#####################
## Load metadata ##
#####################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & doublet_call==FALSE] %>%
  .[sample%in%opts$samples]

##################
## Subset ArchR ##
##################

ArchRProject.filt <- ArchRProject[sample_metadata$cell,]
table(getCellColData(ArchRProject.filt,"Sample")[[1]])

###############
## Load data ##
###############

# Load trajectory
io$pseudotime <- "/Users/ricard/data/gastrulation_multiome_10x/results/rna/trajectories/blood_trajectory/blood_trajectory.txt.gz"
trajectory.dt <- fread(io$pseudotime)[,c("cell","PC1")]
trajectory.mtx <- trajectory.dt %>% matrix.please

# Load data.tables with RNA + chromVAR estimates
# chromvar_rna_dt <- names(io$trajectories.inputdir) %>% 
#   map(function(x) fread(sprintf(sprintf("%s/chromvar_rna.txt.gz",io$trajectories.inputdir[[x]]))) %>%
#         .[,trajectory:=x]) %>%
#   rbindlist

# Load correlation results
cor_dt <- names(io$trajectories.inputdir) %>% 
  map(function(x) fread(sprintf(sprintf("%s/correlation_results.txt.gz",io$trajectories.inputdir[[x]]))) %>%
        .[,trajectory:=x]) %>%
  rbindlist

cor_dt[,sig:=padj_fdr<0.05 & abs(r)>=0.15]
  
#######################
## Fetch Peak Matrix ##
#######################

# Get Peak Matrix from ArchR
atac.peak.se <- getMatrixFromProject(ArchRProject.filt, useMatrix="PeakMatrix", binarize = T)
dim(atac.peak.se)

# Define peak names
peak_names <- rowRanges(atac.peak.se) %>% as.data.table %>% .[,id:=sprintf("%s:%s-%s",seqnames,start,end)] %>% .$id
rownames(atac.peak.se) <- peak_names

###################################################
## Load correlation matrix between TFs and genes ##
###################################################

cor_rna.mtx <- readRDS("/Users/ricard/data/gastrulation_multiome_10x/results/rna/coexpression/correlation_matrix_tf2gene.rds")

# Restrict to marker genes
marker_genes.dt <- fread(io$rna.atlas.marker_genes)
cor_rna.mtx <- cor_rna.mtx[,colnames(cor_rna.mtx)%in%marker_genes.dt$gene]


#############################
## Load motif2peak matrix ##
#############################

motif2peak.matches <- readRDS(getPeakAnnotation(ArchRProject.filt, name=opts$motif_annotation)$Matches)

# Rename motifs
colnames(motif2peak.matches) <- colnames(motif2peak.matches) %>% toupper %>% stringr::str_split(.,"_") %>% map_chr(1)

# Remove duplicated motifs
motif2peak.matches <- motif2peak.matches[,!duplicated(colnames(motif2peak.matches))]

##############################################################
## Find genes, TFs and peaks that correlate with pseudotime ##
##############################################################

TF.cor <- cor_dt[sig==T & abs(r)>0.30 & rna_sign=="Up" & chromvar_sign=="Up"]

cells <- intersect(trajectory.dt$cell, colnames(atac.peak.se))
peaks.to.keep <- names(which(rowSums(assay(atac.peak.se[,cells]))>100))

peak.cor <- cor(trajectory.mtx[cells,],t(as.matrix(assay(atac.peak.se[peaks.to.keep,cells]))))[1,]
# peak.cor <- names(which(peak.cor[1,]>0.25))

peak_cor.dt <- data.table(
  peak_id = names(peak.cor),
  atac_cor = peak.cor
)# %>% .[,sign:=c("Down","Up")[as.numeric(cor>0)+1]]

##########
## Test ##
##########

i <- "KLF1"

motifs <- colnames(motif2peak.matches)
tmp <- motif2peak.matches[as.logical(assay(motif2peak.matches[,i])==1),i]

rna_cor.dt <- data.table(
  gene = colnames(cor_rna.mtx),
  rna_cor = cor_rna.mtx[i,]
)# %>% .[,sign:=c("Down","Up")[as.numeric(cor>0)+1]]

dt <- rowData(tmp)[,c("score","distToGeneStart","distToTSS","nearestGene","peakType","idx")] %>% as.data.table %>%
  setnames("nearestGene","gene") %>%
  .[,peak_id:=sprintf("%s:%s-%s",seqnames(rowRanges(tmp)), start(rowRanges(tmp)), end(rowRanges(tmp)))] %>%
  merge(rna_cor.dt,by="gene") %>%
  merge(peak_cor.dt,by="peak_id")
# merge(marker_genes.dt[,c("celltype","gene")],by="gene", allow.cartesian=T)
# .[toupper(gene)%in%motifs]



################################
## Identify TF feedback loops ##
################################

dt <- rowData(tmp)[,c("score","distToGeneStart","distToTSS","nearestGene","peakType","idx")] %>% as.data.table %>%
  setnames("nearestGene","gene") %>%
  .[,peak_id:=sprintf("%s:%s-%s",seqnames(rowRanges(tmp)), start(rowRanges(tmp)), end(rowRanges(tmp)))]
