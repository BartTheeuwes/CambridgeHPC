library(ggpubr)

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
} else {
  stop("Computer not recognised")
}

######################################
## Load RNA-based doublet inference ##
######################################

io$rna.doublet <- "/Users/ricard/data/gastrulation_multiome_10x/results/rna/doublets/sample_metadata_after_doublets.txt.gz"
rna_doublet.dt <- fread(io$rna.doublet) %>%
  .[,c("cell","sample","hybrid_call","hybrid_score")] %>%
  setnames(c("hybrid_call","hybrid_score"),c("RNA_doublet_call","RNA_doublet_score"))

######################################
## Load ATAC-based doublet inference ##
######################################

io$atac.doublet <- "/Users/ricard/data/gastrulation_multiome_10x/results/atac/archR/doublets/doublet_scores.txt.gz"
atac_doublet.dt <- fread(io$atac.doublet) %>%
  setnames(c("DoubletScore","DoubletEnrichment"),c("ATAC_doublet_score","ATAC_doublet_enrichment"))

###########
## merge ##
###########

doublet_dt <- merge(rna_doublet.dt, atac_doublet.dt, by=c("cell"))

#################
## Scatterplot ##
#################

ggscatter(doublet_dt, x="RNA_doublet_score","ATAC_doublet_enrichment", size=1,
          add="reg.line", add.params = list(color="blue", fill="lightgray"), conf.int=TRUE) +
  labs(x="RNA doublet score", y="ATAC doublet score") +
  stat_cor(method = "pearson") +
  facet_wrap(~sample) +
  theme(
    axis.text = element_text(size=rel(0.75))
  )


#######################
## Overlay with UMAP ##
#######################