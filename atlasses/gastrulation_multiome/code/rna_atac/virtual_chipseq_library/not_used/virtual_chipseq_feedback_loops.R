
#####################
## Define settings ##
#####################

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

# I/O
io$virtual_chip.dir <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/virtual_chipseq_lenient")
io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/virtual_chipseq_lenient/feedback_loops")

# Options
# opts$TFs <- c("T", "ZIC2", "TAL1", "GATA1", "FOXA2", "GATA4")
opts$TFs <- list.files(io$virtual_chip.dir, ".bed.gz") %>% str_replace_all(".bed.gz","")
opts$max.distance <- 1e5

#############################
## Load peak2gene linkages ##
#############################

# peak2gene.dt <- fread(io$archR.peak2gene.all) %>%
peak2gene.dt <- fread(io$archR.peak2gene.nearest) %>%
  .[,gene:=toupper(gene)] %>%
  .[gene%in%opts$TFs & dist<opts$max.distance] %>%
  .[,peak:=sprintf("chr%s:%s-%s",chr,peak.start,peak.end)]

###################################
## Load virtual ChIP-seq library ##
###################################

virtual_chip.dt <- opts$TFs %>% map(function(i) {
  fread(sprintf("%s/%s.bed.gz",io$virtual_chip.dir,i)) %>%
    setnames(c("chr","start","end","correlation_score","max_accessibility_score","motif_score","motif_counts","score")) %>%
    .[,idx:=sprintf("%s:%s-%s",chr,start,end)] %>%
    .[,c("chr","start","end"):=NULL] %>%
    .[,tf:=i] %>%
    return
}) %>% rbindlist

#########################
## Find feedback loops ##
#########################

opts$min.score <- 0.25

virtual_chip_filt.dt <- virtual_chip.dt %>%
  .[score>=opts$min.score & idx%in%peak2gene.dt$peak] %>%
  .[,c("idx","tf","score")] %>%
  setnames("idx","peak") %>%
  merge(peak2gene.dt[,c("chr","peak.start","peak.end","peak","gene","dist")], by="peak") %>%
  .[gene==tf]

# Save
to.save <- virtual_chip_filt.dt[,c("chr","peak.start","peak.end","tf","score","dist")]
fwrite(to.save, file.path(io$outdir,"TF_predicted_feedback_loops.txt.gz"), sep="\t")

##########
## Plot ##
##########

to.plot <- virtual_chip_filt.dt[,c("tf","peak","score","dist")] %>%
  .[score>=0.5,score:=0.5] %>%
  .[,log_dist:=log10(dist+1)]

ggplot(to.plot, aes_string(x="log_dist", y="score")) +
  geom_jitter(aes(size=score), fill="gray70", shape=21, width = 0.25, height=0.01) +
  ggrepel::geom_text_repel(aes(label=tf), size=3, max.overlaps=Inf, data=to.plot) +
  # geom_segment(aes_string(xend="score"), size=0.5, yend=0) +
  labs(x="Genomic distance from the TSS (log10)", y="TF Binding score") +
  theme_classic() +
  guides(size="none") +
  theme(
    plot.title = element_text(hjust=0.5),
    legend.position = "none",
    axis.text = element_text(size=rel(0.8), color="black"),
    axis.title = element_text(size=rel(1.1), color="black")
  )

  

