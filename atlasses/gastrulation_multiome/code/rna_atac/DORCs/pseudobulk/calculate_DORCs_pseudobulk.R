suppressMessages(library(furrr))
suppressMessages(library(argparse))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--celltypes',      type="character",    nargs="+",  help='celltypes')
# p$add_argument('--celltypes',  type="character",    nargs="+",  help='Cell type')
p$add_argument('--gene_window',  type="integer",                  help='Gene window')
p$add_argument('--ncores',       type="integer",      default=1,  help='Number of cores')
p$add_argument('--remove_ExE_celltypes', action="store_true",   help='Remove ExE cell types?')
p$add_argument('--test_mode',    action="store_true",             help='Test mode? subset number of cells')
p$add_argument('--outdir',       type="character",                help='Output file')
args <- p$parse_args(commandArgs(TRUE))

## START TEST
args$celltypes = c("Epiblast", "Primitive_Streak")
args$gene_window <- 1e4
args$remove_ExE_celltypes <- FALSE
args$ncores <- 1
args$test_mode <- TRUE
args$outdir <- "/Users/ricard/data/gastrulation_multiome_10x/results/rna_atac/DORCs/pseudobulk"
## END TEST



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

############################
## Load RNA and ATAC data ##
############################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/rna_atac/rna_vs_chromvar/pseudobulk/load_rna_chromvar_pseudobulk.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/rna_atac/rna_vs_chromvar/pseudobulk/load_rna_chromvar_pseudobulk.R")
} else {
  stop("Computer not recognised")
}

##########################################
## Link peaks2genes by genomic distance ##
##########################################

cat("Overlapping genes with peaks...\n")

peak_metadata.tmp <- peak_metadata.dt %>% copy %>%
  setkey(chr,start,end)

gene_metadata.tmp <- gene_metadata.dt %>% 
  .[,c("start", "end") := list(start-args$gene_window, end+args$gene_window)] %>% 
  setkey(chr,start,end)

# Do the overlap
ov <- foverlaps(peak_metadata.tmp, gene_metadata.tmp, nomatch = 0) %>% 
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
fwrite(r_null.dt, outfile, sep="\t", quote=F)

