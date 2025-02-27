library(GenomicRanges)

#####################
## Define settings ##
#####################

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

# Options
opts$motif_annotation <- "CISBP"

# I/O
io$chip_dir.prefix <- "~/data/mm10_regulation/TF_ChIP"
io$outdir <- file.path(io$basedir,"results_new/rna_atac/virtual_chipseq/CISBP/individual_peaks"); dir.create(io$outdir, showWarnings = F)
io$virtual_chip.dir <- file.path(io$basedir,"results_new/rna_atac/virtual_chipseq/CISBP")
io$virtual_chip.mtx <- file.path(io$virtual_chip.dir,"virtual_chip.mtx")

###################################
## Load virtual ChIP-seq library ##
###################################

virtual_chip.mtx <- readRDS(io$virtual_chip.mtx)

# Fix peak names (temporary)
if (!all(grepl(":",rownames(virtual_chip.mtx)))) {
  rownames(virtual_chip.mtx) <- sub("-",":",rownames(virtual_chip.mtx))
}

######################
## Load motifmatchR ##
######################

io$motifmatcher.se <- sprintf("%s/Annotations/%s-Scores.rds",io$archR.directory,opts$motif_annotation)
motifmatcher.se <- readRDS(io$motifmatcher.se)

io$motifmatcher_positions.se <- sprintf("%s/Annotations/%s-Positions.rds",io$archR.directory,opts$motif_annotation)
motifmatcher_positions.se <- readRDS(io$motifmatcher_positions.se)

stopifnot(colnames(motifmatcher.se)==names(motifmatcher_positions.se))

################
## Subset TFs ##
################

source(here::here("atac/archR/load_motif_annotation.R"))

motifs <- intersect(colnames(motifmatcher.se),motif2gene.dt$motif)
genes <- intersect(colnames(virtual_chip.mtx),motif2gene.dt$gene)
motif2gene_filt.dt <- motif2gene.dt[motif%in%motifs & gene%in%genes]
motifs <- motif2gene_filt.dt$motif
genes <- motif2gene_filt.dt$gene

tmp <- genes; names(tmp) <- motifs

stopifnot(motif2gene_filt.dt$motif%in%colnames(motifmatcher.se))
stopifnot(motif2gene_filt.dt$motif%in%names(motifmatcher_positions.se))
stopifnot(motif2gene_filt.dt$gene%in%colnames(virtual_chip.mtx))

motifmatcher.se <- motifmatcher.se[,motifs]
motifmatcher_positions.se <- motifmatcher_positions.se[motifs]
colnames(motifmatcher.se) <- tmp[colnames(motifmatcher.se)]
names(motifmatcher_positions.se) <- tmp[names(motifmatcher_positions.se)]

virtual_chip.mtx <- virtual_chip.mtx[,genes]

########################
## Load peak metadata ##
########################

peak_metadata.dt <- fread(io$archR.peak.metadata) %>%
  .[,c("chr","start","end","score")] %>%
  .[,chr:=factor(chr)] %>%
  .[,idx:=sprintf("%s:%s-%s",chr,start,end)] %>%
  .[idx%in%rownames(virtual_chip.mtx)]

peak_metadata.gr <- makeGRangesFromDataFrame(peak_metadata.dt, keep.extra.columns = F)

###########################
## Plot individual peaks ##
###########################

# opts$min_motif_score <- 4

peaks.to.plot <- c("chr2:32154225-32154825", "chr9:61808496-61809096", "chr4:117630220-117630820", "chr17:46618359-46618959","chr5:23922696-23923296","chr2:28598060-28598660")
stopifnot(peaks.to.plot%in%rownames(virtual_chip.mtx))

# i <- "chr1:133329196-133329796"
for (i in peaks.to.plot) {
  
  # Identify motif locations + scores
  peak.gr <- rowRanges(motifmatcher.se[i,])
  TFs <- which(assay(motifmatcher.se[i,],"motifMatches")[1,]) %>% names
  tmp <- TFs %>% map(function(j) {
    hits <- findOverlaps(
      query = peak.gr,
      subject = motifmatcher_positions.se[[j]],
      ignore.strand=TRUE
    ) %>% subjectHits()
    if (length(hits)>0) {
      tmp <- motifmatcher_positions.se[[j]][hits]
      tmp$TF <- j
      return(tmp)
    } else {
      print(sprintf("No hits found for %s",j))
      return(NULL)
    }
  }) 
  motif_locations.dt <- tmp[!sapply(tmp,is.null)] %>% as(., "GRangesList") %>% unlist %>%
    as.data.table %>% setnames("score","motif_score") %>% setorder(start) 
  
  # Add in silico ChIP-seq scores
  stopifnot(motif_locations.dt$TF%in%colnames(virtual_chip.mtx))
  insilico_chip_score <- virtual_chip.mtx[i,motif_locations.dt$TF]
  
  # Prepare data.table for plotting
  to.plot <- motif_locations.dt %>% .[,insilico_chip_score:=insilico_chip_score]
  
  # Remove palindromic motif duplicates
  to.plot[,N:=.N,by=c("seqnames","start","end","TF")] 
  to.plot <- rbind(to.plot[N==1],to.plot[N==2 & strand=="+"]) %>% .[,N:=NULL]
  
  # Filter by mininum motif score
  # to.plot <- to.plot %>% .[score>=opts$min_motif_score]
  # to.plot <- motif_locations.dt %>% .[score<=opts$min_motif_score,score:=opts$min_motif_score]
  
  # Plot
  to.plot.text <- to.plot %>% .[insilico_chip_score>=0.25] %>% setorder(-insilico_chip_score) %>% head(n=6)
  
  p <- ggplot(to.plot, aes_string(x="start", y="motif_score")) +
    geom_jitter(aes(fill=insilico_chip_score, size=motif_score, shape=strand, alpha=insilico_chip_score), width = 10, height=0.15) +
    scale_x_continuous(limits=c(start(peak.gr),end(peak.gr))) +
    # scale_fill_gradientn(colours = rev(terrain.colors(10))) +
    # scale_fill_gradient2(low = "gray50", high = "gray10") +
    scale_fill_gradient2(low = "blue", mid="gray90", high = "red", limits=c(-1,1)) +
    scale_shape_manual(values=c(21,24)) +
    scale_size_continuous(range=c(0.5,3)) +
    scale_alpha_continuous(range=c(0.25,1)) +
    ggrepel::geom_text_repel(aes(label=TF), size=3, max.overlaps=Inf, data=to.plot.text) +
    # geom_segment(aes_string(xend="score"), size=0.5, yend=0) +
    labs(x="Genomic position (ATAC peak)", y="Motif score", title=i) +
    theme_classic() +
    guides(size="none", alpha="none") +
    theme(
      plot.title = element_text(hjust=0.5),
      legend.position = "top",
      axis.text = element_text(size=rel(0.8), color="black"),
      axis.title = element_text(size=rel(1.1), color="black")
    )
  
  pdf(file.path(io$outdir,sprintf("peak_%s_motif_location_virtual_chip.pdf",sub(":","-",i))), width=6, height=4)
  print(p)
  dev.off()
}
