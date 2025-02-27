library(scran)

#####################
## Define settings ##
#####################

source("/Users/ricard/gastrulation_multiome_10x/settings.R")
io$outdir <- paste0(io$basedir,"/results/rna/doublets")

opts$max_doublet_score <- 5000

###############
## Load data ##
###############

# load sample metadata
sample_metadata <- fread(io$metadata) %>% .[pass_rnaQC==TRUE]

# load SingleCellExperiment
sce <- readRDS(io$sce)[,sample_metadata$cell]
dim(sce)

#############################
## Calculate doublet score ##
#############################

# scDblFinder::computeDoubletDensity
foo <- doubletCells(sce)

############################
## Update sample metadata ##
############################

sample_metadata <- fread(io$metadata) %>% 
    .[,doublet_score:=NULL] %>%
    merge(data.table(cell=colnames(sce), doublet_score=foo), by="cell",all.x=TRUE) %>%
    .[doublet_score>=opts$max_doublet_score,doublet_score:=opts$max_doublet_score]

head(sample_metadata)

# Save output
# fwrite(sample_metadata, io$metadata, sep="\t", na="NA", quote=F)

##########
## Plot ##
##########

to.plot <- sample_metadata %>% 
    melt(id.vars=c("cell"), measure.vars=c("doublet_score"))

p <- gghistogram(to.plot, x="value", bins=50) +
    geom_vline(xintercept=opts$max_doublet_score, linetype="dashed") +
    theme(
        axis.text =  element_text(size=rel(0.8)),
        legend.position = "right"
    )

# table(sample_metadata$doublet_score<opts$max_doublet_score)

pdf(sprintf("%s/doublet_score.pdf",io$outdir), width=9, height=5, useDingbats = F)
print(p)
dev.off()

