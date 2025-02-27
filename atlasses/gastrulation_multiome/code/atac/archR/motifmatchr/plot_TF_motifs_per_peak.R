
#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

# Options
opts$motif_annotation <- "Motif_cisbp_lenient"


# I/O
io$outdir <- paste0(io$basedir,"/results/atac/archR/peaks/test")

######################
## Load motifmatchR ##
######################

source(here::here("load_motifmatchR.R"))
source(here::here("load_motifmatchR_positions.R"))

########################
## Load peak metadata ##
########################

peak_metadata.dt <- fread(io$archR.peak.metadata) %>% 
  .[,idx:=sprintf("%s:%s-%s",chr,start,end)]

opts$min.peak_score <- 50
peaks <- peak_metadata.dt[score>=opts$min.peak_score,idx]

motifmatcher.se <- motifmatcher.se[peaks,]

###################################
## Load virtual ChIP-seq library ##
###################################

# TFs <- assay(motifmatcher.se[i])[1,] %>% which %>% names
TFs.chip <- list.files(io$virtual_chip.dir, "*_filt.bed.gz") %>% str_replace_all("_filt.bed.gz","")

##########
## Plot ##
##########

number_motifs_per_peak <- rowSums(assay(motifmatcher.se,"motifMatches"))
opts$min_motif_score <- 3

peaks.to.plot <- number_motifs_per_peak[number_motifs_per_peak>=50 & number_motifs_per_peak<=100] %>% names %>% sample(size=25)

i <- "chrX:103822354-103822954"
for (i in peaks.to.plot) {
  
  peak.gr <- rowRanges(motifmatcher.se[i,])
  
  TFs <- which(assay(motifmatcher.se[i,],"motifMatches")[1,]) %>% names
  TFs <- TFs[TFs%in%TFs.chip]
  
  # Load virtual ChIP-seq
  virtual_chip.dt <- TFs %>% map(function(j) {
    file <- sprintf("%s/%s.bed.gz",io$virtual_chip.dir,j)
    if (file.exists(file)) {
      fread(file, select=c(1,2,3,8)) %>%
        setnames(c("chr","start","end","score")) %>%
        .[!is.na(score)] %>%
        .[,idx:=sprintf("%s:%s-%s",chr,start,end)] %>%
        .[idx==i] %>%
        .[,c("chr","start","end"):=NULL] %>%
        .[,tf:=j] %>%
        return
    }
  }) %>% rbindlist
  
  ## START TEST
  foo <- motifmatcher_positions.se[[j]][seqnames(motifmatcher_positions.se[[j]])==seqnames(peak.gr)]
  foo[(start(foo) >= start(peak.gr)) & (end(foo) <= end(peak.gr))]
  ## END TEST
  
  motif_locations.dt <- TFs %>% map(function(j) {
    hits <- findOverlaps(
      query = peak.gr,
      subject = motifmatcher_positions.se[[j]][seqnames(motifmatcher_positions.se[[j]])==seqnames(peak.gr)]
    ) %>% subjectHits()
    tmp <- motifmatcher_positions.se[[j]][hits]
    tmp$TF <- j
    return(tmp)
  }) %>% as(., "GRangesList") %>% unlist %>%
    as.data.table %>% setorder(start)
  
  # Plot
  to.plot <- motif_locations.dt %>% .[score<=opts$min_motif_score,score:=opts$min_motif_score]
  # to.plot.text <- to.plot %>% setorder(-score) %>% head(n=15)
  if (nrow(to.plot)>50) {
    to.plot.text <- to.plot[sample(.N,50)]
  } else {
    to.plot.text <- to.plot
  }
  
  
  p <- ggplot(to.plot, aes_string(x="start", y="score")) +
    geom_jitter(aes(fill=score, size=score, shape=strand), width = 0.05, height=0.05) +
    scale_x_continuous(limits=c(start(peak.gr),end(peak.gr))) +
    # scale_fill_gradientn(colours = terrain.colors(10)) +
    scale_fill_gradient2(low = "gray50", high = "gray10") +
    scale_shape_manual(values=c(21,24)) +
    scale_size_continuous(range=c(1,4)) +
    ggrepel::geom_text_repel(aes(label=TF), size=3, max.overlaps=Inf, data=to.plot.text) +
    # geom_segment(aes_string(xend="score"), size=0.5, yend=0) +
    labs(x="Genomic position (ATAC peak)", y="Motif score", title=i) +
    theme_classic() +
    guides(fill=F, size=F) +
    theme(
      plot.title = element_text(hjust=0.5),
      legend.position = "none",
      axis.text = element_text(size=rel(0.8), color="black"),
      axis.title = element_text(size=rel(1.1), color="black")
    )
  
  pdf(sprintf("%s/peak_%s_motifs_location.pdf",io$outdir,sub(":","-",i)), width=7, height=5)
  print(p)
  dev.off()
}



#############
## Explore ##
#############

# Plot legend
p <- ggplot(to.plot, aes_string(x="start", y="score")) +
  geom_point(aes(shape=strand), size=5, color="gray10") +
  # scale_shape_manual(values=c(1,2)) +
  # guides(fill=F, size=F) +
  theme_classic() +
  theme(
    legend.position = "top"
  )

pdf(sprintf("%s/legend.pdf",io$outdir), width=7, height=5)
print(p)
dev.off()