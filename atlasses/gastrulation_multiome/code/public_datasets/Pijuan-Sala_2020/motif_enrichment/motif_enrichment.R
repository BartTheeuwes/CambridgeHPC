
########################
## Load ArchR project ##
########################

source("/Users/ricard/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/archR/load_archR_project.R")

#####################
## Define settings ##
#####################

io$outdir <- paste0(io$basedir,"/results/motif_enrichment/archR")

#########################
## Define marker peaks ##
#########################

# get marker features
markersPeaks <- getMarkerFeatures(
  ArchRProj = ArchRProject, 
  useMatrix = "PeakMatrix", 
  groupBy = "celltype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# io$markersPeaks <- "/Users/ricard/data/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/results/markers/markerPeaks.rds"
# saveRDS(markersPeaks,io$markersPeaks)
markersPeaks <- readRDS(io$markersPeaks)

# returns a SummarizedExperiment object
# markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
# names(markerList)
# head(markerList[["C1"]])

# returns a GRanges object
# markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
# dt <- names(markerList) %>% 
#   map(function(i) as.data.table(markerList[[i]]) %>% .[,celltype:=i]) %>% 
#   rbindlist %>%
#   setnames("seqnames","chr") %>%
#   .[,c("strand","width"):=NULL]
# head(dt)

###########################
## Add Motif annotations ##
###########################

# https://www.ArchRProject.com/bookdown/motif-and-feature-enrichment-with-archr.html

opts$motifSet <- "cisbp"     # [JASPAR2016, JASPAR2018, JASPAR2020, cisbp, encode, homer]
opts$collection <- "CORE"    # only for JASPAR motif sets.
opts$motif.pvalue.cutoff <- 5e-05  # default is 5e-05

# add motif set 
ArchRProject <- addMotifAnnotations(
  ArchRProject, 
  motifSet = opts$motifSet,      
  collection = opts$collection,  
  cutOff = opts$motif.pvalue.cutoff,   
  name = "Motif"
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

####################
## Sequence Logos ##
####################

library(ggseqlogo)

motifs <- getPeakAnnotation(ArchRProject)[["motifs"]]

for (i in names(motifs)) {
  m <- 0.25*exp(as.matrix(motifs[[i]]))
  # seqLogo::seqLogo(m)
  p <- ggseqlogo( m )
  pdf(sprintf("%s/seqlogo_%s.pdf",io$outdir,i), width=4, height=2)
  print(p)
  dev.off()
}

###############
## Line plot ##
###############

dt <- dt.up

for (i in unique(dt$celltype)) {
  to.plot <- dt[celltype==i] %>% 
    setorder(-value) %>%
    .[,TF:=stringr::str_split(TF,"_") %>% map_chr(1)] %>%
    .[,.(value=mean(value)), by=c("celltype","TF")] %>%
    head(n=25) %>%
    .[,TF:=factor(TF,levels=rev(TF))]
  
  p <- ggplot(to.plot, aes(x=TF, y=value)) +
    geom_point(size=2) +
    geom_segment(aes_string(xend="TF"), size=0.75, yend=0) +
    labs(x="q-value (hypergeometric test)") +
    coord_flip() +
    theme_bw()
  
  pdf(sprintf("%s/archr_motif_enrichment_lineplot_%s.pdf",io$outdir,i), width=4, height=10)
  print(p)
  dev.off()
}

##########
## Save ##
##########

outfile <- paste0(io$outdir,"/motif_enrichment_up.tsv.gz")
fwrite(dt.up, outfile, sep="\t")

outfile <- paste0(io$outdir,"/motif_enrichment_down.tsv.gz")
fwrite(dt.down, outfile, sep="\t")
