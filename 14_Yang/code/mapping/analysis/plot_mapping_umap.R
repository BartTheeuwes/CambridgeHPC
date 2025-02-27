# NOTE THAT THIS IS CURRENTLY IMPLEMENTED ONLY FOR THE MNN MAPPING

suppressPackageStartupMessages(library(argparse))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--query_metadata',        type="character",                               help='Cell metadata (after mapping)')
p$add_argument('--samples',         type="character",       nargs="+",   help='Samples')
p$add_argument('--atlas_metadata',        type="character",                               help='Cell metadata (after mapping)')
p$add_argument('--outdir',          type="character",                               help='Output file')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$query_metadata <- "/Users/argelagr/data/gastrulation_histones/results/rna/mapping/sample_metadata_after_mapping.txt.gz"
# args$atlas_metadata <- "/Users/argelagr/data//gastrulation10x/sample_metadata.txt.gz"
# args$outdir <- "/Users/argelagr/data/gastrulation_histones/results/rna/mapping/pdf"
## END TEST ##

#####################
## Define settings ##
#####################

# load default setings
here::i_am("processing/1_create_seurat_rna.R")
source(here::here("settings.R"))

# Options

# Dot size
opts$size.mapped <- 0.18
opts$size.nomapped <- 0.1

# Transparency
opts$alpha.mapped <- 0.65
opts$alpha.nomapped <- 0.35

####################
## Load 10x atlas ##
####################

# Load atlas cell metadata
meta_atlas <- fread(args$atlas_metadata) 

# Extract precomputed dimensionality reduction coordinates
umap.dt <- meta_atlas %>%
  .[,c("cell","umapX","umapY","celltype", "celltype_extended")] %>%
  setnames(c("umapX","umapY"),c("V1","V2")) %>%
    .[, mapped:="scRNA-seq atlas"]

###################
## Load metadata ##
###################

sample_metadata <- fread(args$query_metadata) %>%
    .[pass_rnaQC==TRUE & sample%in%args$samples & !is.na(closest.cell_mnn)] %>%
    .[, c("closest.cell_mnn","celltype.mapped_mnn", "celltype_extended.mapped_mnn", "tdTom")] %>%
    merge(., umap.dt[,c("cell", "V1","V2")], by.x='closest.cell_mnn', by.y='cell') %>%
    setnames(c("closest.cell_mnn", "celltype.mapped_mnn","celltype_extended.mapped_mnn"),c("cell", "celltype","celltype_extended")) %>%
    .[, mapped:="Chimaera"] %>%
    .[, c('cell', 'V1', 'V2', 'celltype', 'celltype_extended', 'mapped', 'tdTom')]

to.plot = rbind(umap.dt, copy(sample_metadata)[,tdTom:=NULL])

##############################
## Define plotting function ##
##############################


plot.dimred <- function(plot_df, query.label, atlas.label = "Atlas") {
  
  # Define dot size  
  size.values <- c(opts$size.mapped, opts$size.nomapped)
  names(size.values) <- c(query.label, atlas.label)
  
  # Define dot alpha  
  alpha.values <- c(opts$alpha.mapped, opts$alpha.nomapped)
  names(alpha.values) <- c(query.label, atlas.label)
  
  # Define dot colours  
  colour.values <- c("red", "lightgrey")
  names(colour.values) <- c(query.label, atlas.label)
  
  # Plot
  ggplot(plot_df, aes(x=V1, y=V2)) +
    ggrastr::geom_point_rast(aes(size=mapped, alpha=mapped, colour=mapped)) +
    scale_size_manual(values = size.values) +
    scale_alpha_manual(values = alpha.values) +
    scale_colour_manual(values = colour.values) +
    # labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
    guides(colour = guide_legend(override.aes = list(size=6))) +
    theme_classic() +
    theme(
      legend.position = "top", 
      legend.title = element_blank(),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank()
    )
}


####################
## Plot all cells ##
####################

#plot mapping
to.plot = rbind(umap.dt, copy(sample_metadata)[,tdTom:=NULL])
p <- plot.dimred(to.plot, query.label = "Chimaera", atlas.label = "scRNA-seq atlas")

pdf(sprintf("%s/umap_mapped_allcells.pdf",args$outdir), width=8, height=6.5)
print(p)
dev.off()

# Plot celltypes
to.plot = rbind(copy(umap.dt)[,tdTom:=TRUE], copy(umap.dt)[,tdTom:=FALSE], sample_metadata) %>% # umap.dt[,tdTom:=TRUE],
    .[,tdTom:=ifelse(tdTom==TRUE, 'KO', 'WT')]

p =  ggplot(to.plot[mapped=='Chimaera'], aes(x=V1, y=V2)) +
    ggrastr::geom_point_rast(data = to.plot[mapped=='scRNA-seq atlas'], 
                             colour='lightgrey', size=opts$size.nomapped , alpha=opts$alpha.nomapped) +
    ggrastr::geom_point_rast(aes(colour=celltype), size=opts$size.mapped , alpha=opts$alpha.mapped) +
    scale_colour_manual(values = opts$celltype.colors) +
    facet_wrap(~tdTom) +
    # labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
    theme_classic() +
    theme(
      legend.position = "none", 
      text = element_text(size=15), 
      legend.title = element_blank(),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank()
    )

pdf(sprintf("%s/umap_mapped_celltype.pdf",args$outdir), width=8, height=5)
    print(p)
dev.off()

p =  ggplot(to.plot[mapped=='Chimaera'], aes(x=V1, y=V2)) +
    ggrastr::geom_point_rast(data = to.plot[mapped=='scRNA-seq atlas'], 
                             colour='lightgrey', size=opts$size.nomapped , alpha=opts$alpha.nomapped) +
    ggrastr::geom_point_rast(aes(colour=celltype_extended), size=opts$size.mapped , alpha=opts$alpha.mapped) +
    scale_colour_manual(values = opts$celltype_extended.colors) +
    facet_wrap(~tdTom) +
    # labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
    theme_classic() +
    theme(
      legend.position = "none", 
      text = element_text(size=15), 
      legend.title = element_blank(),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank()
    )

pdf(sprintf("%s/umap_mapped_celltype_extended.pdf",args$outdir), width=8, height=5)
    print(p)
dev.off()

###############################
## Plot one sample at a time ##
###############################
sample_metadata <- fread(args$query_metadata) %>%
    .[pass_rnaQC==TRUE & sample%in%args$samples & !is.na(closest.cell_mnn)] %>%
    .[, c("closest.cell_mnn","sample", 'celltype.mapped_mnn')] %>%
    setnames('celltype.mapped_mnn', 'celltype') %>% 
    merge(., umap.dt[,c("cell", "V1","V2")], by.x='closest.cell_mnn', by.y='cell') %>%
    setnames(c("closest.cell_mnn"),c("cell")) %>%
    .[, mapped:="Chimaera"] %>%
    .[, c('cell', 'V1', 'V2', 'mapped', 'sample', 'celltype')]

for (i in args$samples) {  
    to.plot = rbind(umap.dt[,c('cell', 'V1', 'V2', 'mapped')][,`:=`(sample='atlas', celltype=NA)], 
        copy(sample_metadata)[sample==i][,mapped:=sample])
    opts$size.mapped = 0.3

    p =  ggplot(to.plot[mapped==i], aes(x=V1, y=V2)) +
        ggrastr::geom_point_rast(data = to.plot[mapped=='scRNA-seq atlas'], 
                                 colour='lightgrey', size=opts$size.nomapped , alpha=opts$alpha.nomapped) +
        ggrastr::geom_point_rast(aes(colour=celltype), size=opts$size.mapped , alpha=opts$alpha.mapped) +
        scale_colour_manual(values = opts$celltype.colors) +
        theme_classic() +
        theme(
          legend.position = "none", 
          text = element_text(size=15), 
          legend.title = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank()
        )
    
  pdf(sprintf("%s/umap_mapped_%s.pdf",args$outdir,i), width=8, height=6.5)
      print(p)
  dev.off()
}
