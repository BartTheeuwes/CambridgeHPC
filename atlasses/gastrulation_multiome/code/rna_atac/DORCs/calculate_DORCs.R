suppressMessages(library(furrr))
suppressMessages(library(argparse))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--samples',      type="character",    nargs="+",  help='Samples')
# p$add_argument('--celltypes',  type="character",    nargs="+",  help='Cell type')
p$add_argument('--gene_window',  type="integer",                  help='Gene window')
p$add_argument('--ncores',       type="integer",      default=1,  help='Number of cores')
p$add_argument('--denoised',     action="store_true",             help='Use denoised ATAC data?')
p$add_argument('--remove_ExE_celltypes', action="store_true",   help='Remove ExE cell types?')
p$add_argument('--test_mode',    action="store_true",             help='Test mode? subset number of cells')
p$add_argument('--outdir',       type="character",                help='Output file')
args <- p$parse_args(commandArgs(TRUE))

## START TEST
args$samples <- c("E7.5_rep1", "E7.5_rep2", "E8.0_rep1", "E8.0_rep2", "E8.5_rep1", "E8.5_rep2")
args$gene_window <- 1e4
args$remove_ExE_celltypes <- FALSE
args$ncores <- 1
args$test_mode <- TRUE
args$outdir <- "/Users/ricard/data/gastrulation_multiome_10x/results/rna_atac/DORCs"
## END TEST

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

# Define multiprocessing
options(future.globals.maxSize = 7 * 1024^3)
if (args$ncores>1) {
  plan(multiprocess, workers=args$ncores)
} else {
  plan(sequential)
}

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & pass_rnaQC==TRUE & doublet_call==FALSE] %>%
  .[sample%in%args$samples]
stopifnot(sample_metadata$cell %in% rownames(ArchRProject))

##################
## Subset ArchR ##
##################

ArchRProject.filt <- ArchRProject[sample_metadata$cell,]
table(ArchRProject.filt$sample)

##############################
## Load RNA expression data ##
##############################

cat("Loading SingleCellExperiment...\n")

# Load SingleCellExperiment
sce <- load_SingleCellExperiment(io$sce, cells = sample_metadata$cell, normalise = TRUE, remove_non_expressed_genes = TRUE)

# Filter genes
sce <- sce[!grepl("Rik|Gm",rownames(sce)),]

# Load gene metadata
gene.metadata <- fread(io$gene_metadata) %>%
  # .[,tss:=ifelse(strand=="+",start,end)] %>%
  .[symbol%in%rownames(sce)] %>%
  .[,c("chr","start","end","symbol","ens_id")] %>%
  setnames("symbol","gene")

#######################################
## Load chromatin accessibility data ##
#######################################

cat("Loading Peak Matrix...\n")

# Load peak matrix as a SummarizedExperiment object
# ArchRProject.filt <- addPeakMatrix(ArchRProject.filt, binarize = TRUE, force = TRUE)
atac.peaks.se <- getMatrixFromProject(ArchRProject.filt, binarize = TRUE, useMatrix = "PeakMatrix")[,sample_metadata$cell]

# Make sure that chr order is the same between @peakSet and the PeakMatrix
stopifnot(levels(seqnames(ArchRProject.filt@peakSet)) == levels(seqnames(rowRanges(atac.peaks.se))))

# Define peak names
tmp <- rowRanges(atac.peaks.se) %>% as.data.table %>% .[,peak:=sprintf("%s:%s-%s",seqnames,start,end)]
rownames(atac.peaks.se) <- tmp$peak

# Add peak metadata to RowRanges
peak_metadata_external <- fread(io$archR.peak.metadata) %>%
  .[,peak:=sprintf("%s:%s-%s",chr,start,end)] %>%
  .[,c("peak","idx","chr","start","end","width","strand","score","distToGeneStart","distToTSS","nearestGene","peakType","GC")] %>%
  GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T)
names(peak_metadata_external) <- peak_metadata_external$peak
peak_metadata_external <- peak_metadata_external[rownames(atac.peaks.se)]
rowRanges(atac.peaks.se)  <- peak_metadata_external


# saveRDS(atac.peaks.se, paste0(args$outdir,"/atac_SummarizedExperiment.rds"))

################################################
## Load denoised chromatin accessibility data ##
###############################################


##########################################
## Link peaks2genes by genomic distance ##
##########################################

cat("Overlapping genes with peaks...\n")

peak.metadata_tmp <- peak.metadata %>% copy %>%
  setkey(chr,start,end)

gene.metadata_tmp <- gene.metadata %>% 
  .[,c("start", "end") := list(start-args$gene_window, end+args$gene_window)] %>% 
  setkey(chr,start,end)

# Do the overlap
ov <- foverlaps(peak.metadata_tmp, gene.metadata_tmp, nomatch = 0) %>% 
  setnames(c("i.start","i.end"),c("peak.start","peak.end")) %>%
  .[,c("gene.start","gene.end") := list(start+args$gene_window, end-args$gene_window)] %>% 
  .[,c("start","end"):=NULL] %>%
  .[,c("start_dist","end_dist"):=list( gene.end-peak.start, gene.start-peak.end)] %>%
  .[,c("start_dist","end_dist"):=list( ifelse(end_dist<0 & start_dist>0,0,start_dist), ifelse(end_dist<0 & start_dist>0,0,end_dist) )] %>%
  .[,dist:=ifelse(abs(start_dist)<abs(end_dist),abs(start_dist),abs(end_dist))] %>% 
  .[,c("idx","strand","start_dist","end_dist","width"):=NULL]

################################################
## Calculate Spearman correlation coefficient ##
################################################

cat("Calculating correlation coefficients between gene-peak pairs...\n")

# Define genes
genes <- unique(ov$gene)
if (args$test_mode) { genes <- genes %>% head(n=10) }

r.dt <- genes %>% future_map(function(i) {
  
  foo <- as.matrix(logcounts(sce[i,])) %>% t %>% as.data.table(keep.rownames="cell") %>%
    melt(id.vars="cell", value.name="expr", variable.name="gene")
  
  peaks <- unique(ov[gene==i,peak])
  bar <- assay(atac.peaks.se[peaks,]) %>% as.matrix %>% t %>%
    as.data.table(keep.rownames="cell") %>%
    melt(id.vars="cell", value.name="accessibility", variable.name="peak")
  
  merge(foo,bar,by="cell") %>%
    .[,.(r=cor(accessibility, expr, method="spearman") %>% round(3)), by=c("gene","peak")] %>%
    return(.)
}) %>% rbindlist

#########################
## Estimate background ##
#########################

cat("Loading background peaks...\n")

# Compute background peaks
# ArchRProject.filt <- addBgdPeaks(ArchRProject.filt, method = "ArchR", force = TRUE)

# Load previously computed background peaks
background.peaks.se <- getBgdPeaks(ArchRProject.filt)  # file: metadata(getPeakSet(ArchRProject.filt))$bgdPeaks
stopifnot(seqlevels(rowRanges(background.peaks.se)) == seqlevels(rowRanges(atac.peaks.se)))

# Define peak names
rownames(background.peaks.se) <- rowRanges(background.peaks.se) %>% as.data.table %>% 
  setnames("seqnames","chr") %>% .[,peak:=sprintf("%s:%s-%s",chr,start,end)] %>% .$peak
stopifnot(rownames(background.peaks.se)%in%rownames(atac.peaks.se))

################################
## Generate null distribution ##
################################

cat("Calculating Null distribution using background peaks...\n")

# calculate the expected population mean and expected population standard deviation for the Spearman correlations. 

r_null.dt <- genes %>% future_map(function(i) {
  
  foo <- as.matrix(logcounts(sce[i,])) %>% t %>% as.data.table(keep.rownames="cell") %>%
    melt(id.vars="cell", value.name="expr", variable.name="gene")
  peaks <- unique(ov[gene==i,peak])
  peaks %>% map(function(j) {
    
    background.peaks <- rownames(background.peaks.se)[assay(background.peaks.se[j,])[1,]]
    bar <- assay(atac.peaks.se[background.peaks,]) %>% as.matrix %>% t %>%
      as.data.table(keep.rownames="cell") %>%
      melt(id.vars="cell", value.name="accessibility", variable.name="background_peak")
    
    merge(foo,bar,by="cell") %>%
      .[,.(r=cor(accessibility, expr, method="spearman")), by = c("gene","background_peak")] %>%
      .[,peak:=j] %>%
      return(.)
  }) %>% rbindlist }) %>% rbindlist %>% 
  .[,.(bgd_mean=round(mean(r),6), bgd_sd=round(sd(r),6)) ,by = c("gene","peak")] %>%
  merge(r.dt, by=c("gene","peak")) %>%
  .[,z:=round((r-bgd_mean)/bgd_sd,3)] %>%  # z_score = (obs-pop.mean)/pop.sd
  .[,p:=round(pnorm(abs(z),lower.tail=F),6)] %>%
  .[,padj_fdr:=round(p.adjust(p, method="fdr"),6)] %>%
  .[,sig:=padj_fdr<0.10]

# Save
outfile <- sprintf("%s/DORCs_correlation_estimates_%s_%s.tsv.gz",args$outdir, paste(args$samples,collapse="-"), args$gene_window)
cat(sprintf("Saving in %s...\n",outfile))
fwrite(r_null.dt, outfile, sep="\t", quote=F)

##################
## Define DORCs ##
##################

cat("Defining DORCs...\n")

# foo[,sum(sig,na.rm=T),by="gene"]
# To define DORCs (a set of nearby peaks per gene), we rank genes by the number of significantly associated peaks ( ± 50kb around TSSs, p < 0.05).
