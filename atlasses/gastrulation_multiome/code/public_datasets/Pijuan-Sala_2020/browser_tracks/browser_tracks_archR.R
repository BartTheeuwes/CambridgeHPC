
#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/utils.R")
} else {
  stop("Computer not recognised")
}

# I/O
io$outdir <- paste0(io$basedir,"/results/BrowserTrack"); dir.create(io$outdir, showWarnings = F)

# Options
opts$celltypes = c(
  "Surface_ectoderm",
  "Notochord",
  "Gut",
  "Cardiomyocytes",
  "Mid_Hindbrain",
  "Endothelium",
  "Paraxial_mesoderm",
  "Spinal_cord",
  "Somitic_mesoderm",
  "Erythroid",
  "Neural_crest",
  "Mixed_mesoderm",
  "NMP",
  # "Forebrain",
  "ExE_endoderm",
  "Allantois",
  "Mesenchyme",
  "Pharyngeal_mesoderm"
)

########################
## Load ArchR project ##
########################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/load_archR_project.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/load_archR_project.R")
} else {
  stop("Computer not recognised")
}

##########################
## Update cell metadata ##
##########################

sample_metadata <- sample_metadata %>%
  .[celltype%in%opts$celltypes]

table(sample_metadata$celltype)

#########################
## Subset ArchR object ##
#########################

sample_metadata <- sample_metadata[cell%in%rownames(ArchRProject)]
ArchRProject.filt <- ArchRProject[sample_metadata$cell]

# add celltype to ArchR's CellColData
# tmp <- sample_metadata %>% 
#   .[cell%in%rownames(ArchRProject.filt)] %>% setkey(cell) %>% .[rownames(ArchRProject.filt)] %>%
#   as.data.frame() %>% tibble::column_to_rownames("cell")
# stopifnot(all(tmp$TSSEnrichment_atac == getCellColData(ArchRProject.filt, "TSSEnrichment")[[1]]))
# ArchRProject.filt <- addCellColData(
#   ArchRProject.filt,
#   data = tmp$celltype, 
#   name = "celltype",
#   cells = rownames(tmp),
#   force = TRUE
# )
stopifnot(sample_metadata$celltype == ArchRProject.filt$celltype)

#######################
## Load marker genes ##
#######################

marker_genes.dt <- fread(io$rna.atlas.marker_genes)

# GenomicRanges
genes.gr <- getGeneAnnotation(ArchRProject.filt)[["genes"]]
genes.gr <- genes.gr[genes.gr$symbol%in%unique(marker_genes.dt$gene)]

####################
## Browser tracks ##
####################

genes.to.plot <- genes.gr$symbol %>% head(n=5) %>% unname
# i <- "chr2:39483639-39484239"

# opts$extend.upstream <- 1e4
# opts$extend.downstream <- 1e4
# opts$tileSize <- 50

opts$tileSize <- 30

# Ugly hack
celltype.order = opts$celltypes
stopifnot(sort(celltype.order)==sort(unique(ArchRProject.filt$celltype)))
rename <- paste(1:length(celltype.order),celltype.order,sep="_")
names(rename) <- celltype.order
ArchRProject.filt$celltype2 <- stringr::str_replace_all(ArchRProject.filt$celltype,rename)
opts$celltype.colors2 <- opts$celltype.colors[names(opts$celltype.colors)%in%unique(ArchRProject.filt$celltype)]
names(opts$celltype.colors2) <- stringr::str_replace_all(names(opts$celltype.colors2),rename)

for (i in genes.to.plot) {
  
  to.plot.gr <- genes.gr[genes.gr$symbol==i]
  gene.length <- abs(end(to.plot.gr) - start(to.plot.gr))
  start(to.plot.gr) <- start(to.plot.gr) - round(gene.length/1.5)
  end(to.plot.gr) <- end(to.plot.gr) + round(gene.length/1.5)
  
  # Plot
  p <- plotBrowserTrack(
    ArchRProj = ArchRProject.filt, 
    region = to.plot.gr,
    geneSymbol = i,
    # useMatrix = "GeneScoreMatrix",
    groupBy = "celltype2", 
    tileSize = opts$tileSize,
    # upstream = opts$extend.upstream,
    # downstream = opts$extend.downstream,
    pal = opts$celltype.colors2,
    plotSummary = c("bulkTrack", "featureTrack", "geneTrack"),
    sizes = c(13, 1, 1),
  )
  
  pdf(sprintf("%s/%s_BrowserTrack.pdf",io$outdir,i), width = 9, height = 5)
  grid::grid.draw(p[[1]])
  dev.off()
}
