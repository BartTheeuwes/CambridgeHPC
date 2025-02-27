library(ggpubr)

matrix.please<-function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[[1]]
  m
}

##############
## Settings ##
##############

source("/Users/ricard/gastrulation_multiome_10x/settings.R")
io$atlas.marker_genes <- paste0(io$atlas.basedir,"/results/marker_genes/E8.5/marker_genes.txt.gz")
io$outdir <- paste0(io$basedir,"/results/rna/celltype_affinity"); dir.create(io$outdir, showWarnings = F)

###############
## Load data ##
###############

# Load metadata
sample_metadata <- sample_metadata %>%
  .[pass_rnaQC==T & stage=="E8.5"]
  # .[,c("id_rna","stage")]


# Load SingleCellExperiment
sce <- readRDS(io$sce)[,sample_metadata$cell]

# Load gene markers
marker_genes.dt <- fread(io$atlas.marker_genes) %>%
  .[gene%in%rownames(sce)]

################
## Parse data ##
################

# Subset genes
sce <- sce[unique(marker_genes.dt$gene),]
dim(sce)

# if (isTRUE(opts$scale)) {
#   stop("TO-DO")
# } else {
#   expr.matrix <- as.matrix(logcounts(sce))
# }


################################
## Computations per cell type ##
################################

# Calculate atlas-based cell type "score" for each predicted celltype in the query
# dt.celltype <- unique(sce$celltype.mapped) %>% map(function(i) {
#   logcounts(sce[,sce$celltype.mapped==i]) %>% 
#     rowMeans %>% as.data.table(keep.rownames = T) %>% 
#     setnames(c("ens_id","expr")) %>%
#     merge(marker_genes.dt[,c("celltype","ens_id")], by="ens_id", allow.cartesian=TRUE) %>%
#     .[,.(score=mean(expr)),by="celltype"] %>%
#     .[,c("celltype.mapped"):=i]
# }) %>% rbindlist 

###########################
## Computations per cell ##
###########################

# Calculate atlas-based cell type "score" for each cell in the query
dt.cell <- logcounts(sce) %>% as.data.table(keep.rownames = T) %>% 
  melt(id.vars="rn") %>%
  setnames(c("gene","cell","expr")) %>%
  merge(marker_genes.dt[,c("celltype","gene")], by="gene", allow.cartesian=TRUE) %>%
  .[,.(score=mean(expr)),by=c("cell","celltype")]

########################
## Plot per cell type ##
########################

colors <- opts$celltype.colors[names(opts$celltype.colors)%in%unique(marker_genes.dt$celltype)]

# names(colors) <- names(colors) %>% stringr::str_replace_all("_"," ")
# 
# # Plot number of marker genes per cell types
# for (i in unique(dt$celltype.mapped)) {
#   
#   to.plot <- dt.celltype[celltype.mapped==i] %>%
#     .[,celltype:=stringr::str_replace_all(celltype,"_"," ")] %>%
#     .[,celltype:=factor(celltype,levels=names(colors))]
#   
#   p <- ggbarplot(to.plot, x="celltype", y="score", fill="celltype") +
#     scale_fill_manual(values=colors) +
#     labs(x="", y="Cell type affinity") +
#     theme(
#       axis.text.y = element_text(size=rel(0.75)),
#       axis.text.x = element_text(colour="black",size=rel(0.8), angle=90, hjust=1, vjust=0.5),
#       axis.ticks.x = element_blank(),
#       legend.position = "none"
#   )
#   
#   pdf(sprintf("%s/%s_affinity.pdf",io$outdir,i), width = 9, height = 5)
#   print(p)
#   dev.off()
# }


###################
## Plot per cell ##
###################

# for (i in as.character(unique(dt.cell$cell))) {
for (i in head(as.character(unique(dt.cell$cell)),n=50)) {
  
  to.plot <- dt.cell[cell==i] %>%
    # .[,celltype:=stringr::str_replace_all(celltype,"_"," ")] %>%
    .[,celltype:=factor(celltype,levels=names(colors))]
  
  p <- ggbarplot(to.plot, x="celltype", y="score", fill="celltype") +
    scale_fill_manual(values=opts$celltype.colors) +
    labs(x="", y="Cell type affinity") +
    theme(
      axis.text.y = element_text(size=rel(0.75)),
      axis.text.x = element_text(colour="black",size=rel(0.8), angle=90, hjust=1, vjust=0.5),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    )
  
  pdf(sprintf("%s/%s.pdf",io$outdir,i), width = 9, height = 5)
  print(p)
  dev.off()
}


#################################################################
## Cell type prediction by correlation of cell type affinities ##
#################################################################

io$celltype_affinities <- "/Users/ricard/data/gastrulation10x/results/celltype_affinity/E8.0-E8.5/celltype_affinity_E8.5.txt.gz"
reference_celltype_affinities <- fread(io$celltype_affinities) %>%
  setnames(c("celltype.ref","score","celltype"))

# dt.cell.test <- dt.cell %>% copy %>% .[cell%in%head(unique(dt.cell$cell),n=1000)]# %>% setnames("celltype","celltype.ref")

celltype.predictions <- merge(reference_celltype_affinities, dt.cell, by="celltype",allow.cartesian=T) %>%
  setnames(c("score.x","score.y"),c("score_atlas","score_query")) %>%
  .[,.(cor=cor(score_query,score_atlas)),by=c("cell","celltype.ref")] %>%
  setorder(cell,-cor) %>% 
  .[,head(.SD,n=1), by="cell"] %>%
  setnames(c("cell","celltype.predicted","celltype.affinity"))

# to.plot <- unique(dt.cell$cell) %>% map(function(i) {
# to.plot <- head(unique(dt.cell$cell),n=10) %>% map(function(i) {
#   dt.cell[cell==i] %>%
#     setnames("celltype","celltype.ref") %>%
#     merge(reference_celltype_affinities,by="celltype.ref") %>%
#     setnames(c("score.x","score.y"),c("score_query","score_atlas")) %>%
#     .[,.(cor=cor(score_query,score_atlas)),by=c("cell","celltype")]
# }) %>% rbindlist

############################
## Update sample metadata ##
############################

sample_metadata <- fread(io$metadata) %>%
  merge(celltype.predictions, by="cell", all.x=TRUE)

fwrite(sample_metadata, io$metadata, sep="\t", na="NA", quote=F)

