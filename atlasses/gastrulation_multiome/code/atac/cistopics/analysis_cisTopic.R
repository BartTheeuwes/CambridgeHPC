suppressPackageStartupMessages(library(argparse))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--samples',         type="character",                nargs='+',     help='Samples')
p$add_argument('--nfeatures',       type="integer",    default=1000,               help='Number of features')
p$add_argument('--ntopics',         type="integer",      nargs="+",    help='Number of topics')
p$add_argument('--n_neighbors',     type="integer",    default=30,   nargs='+',     help='(UMAP) Number of neighbours')
p$add_argument('--min_dist',        type="double",     default=0.3,  nargs='+',     help='(UMAP) Minimum distance')
p$add_argument('--colour_by',       type="character",  default="celltype.mapped",  nargs='+',  help='Metadata columns to colour the UMAP by')
p$add_argument('--seed',            type="integer",    default=42,                  help='Random seed')
p$add_argument('--indir',           type="character",                               help='Input directory')
p$add_argument('--outdir',          type="character",                               help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

## START TEST
args$samples <- opts$samples
args$nfeatures <- 1e5
args$ntopics <- 50
args$n_neighbors <- 30
args$min_dist <- 0.3
args$seed <- 42
args$colour_by <- c("celltype.mapped","sample","stage","log_nFrags_atac")
args$indir <- paste0(io$basedir,"/results/atac/cistopic/all_cells")
args$outdir <- paste0(io$basedir,"/results/atac/cistopic/all_cells")
## END TEST


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

# I/O
io$pdfdir <- paste0(args$outdir,"/pdf")

# Options

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & doublet_call==FALSE & sample%in%args$samples] %>%
  .[,log_nFrags_atac:=log10(nFrags_atac)]

##########################
## Load cisTopic output ##
##########################

# outfile <- sprintf("%s/%s_lsi_features%d_dims%d_%sbatchcorrectionby%s.txt.gz",args$outdir, paste(args$samples,collapse="-"), args$nfeatures, args$ntopics, args$batch.method, paste(args$batch.variable,collapse="-"))
# E7.5_rep1-E7.5_rep2-E8.0_rep1-E8.0_rep2-E8.5_rep1-E8.5_rep2_cistopic_features100000_ntopics50.txt.gz
io$input.file <- sprintf("%s/%s_cistopic_features%d_ntopics%d.txt.gz",args$indir, paste(args$samples,collapse="-"), args$nfeatures, args$ntopics)
topics.mtx <- fread(io$input.file) %>% matrix.please

#########################
## Plot cisTopic stats ##
#########################

# Plot model selection
# pdf(sprintf("%s/%s_model_selection_features%d_ntopics%d.pdf",io$pdfdir, paste(args$samples,collapse="-"), args$nfeatures, ncol(topics.mtx)), width=7, height=3)
# ggline(lik.dt, x="ntopics", y="lik")
# dev.off()

# Plot correlation between factors
# pdf(sprintf("%s/%s_corrplot_features%d_ntopics%d.pdf",io$pdfdir, paste(args$samples,collapse="-"), args$nfeatures, ncol(topics.mtx)), width=7, height=3)
# corrplot::corrplot(cor(topics.mtx), tl.cex=0.2, tl.col="black")
# dev.off()


##########
## UMAP ##
##########

for (i in args$n_neighbors) {
  for (j in args$min_dist) {
    
    # Run UMAP
    set.seed(args$seed)
    umap.dt <- uwot::umap(topics.mtx, n_neighbors = i, min_dist = j, metric = "cosine") %>%
      as.data.table(keep.rownames=F) %>%
      setnames(c("umap1","umap2")) %>%
      .[,cell:=rownames(topics.mtx)]
    
    # Plotw
    to.plot <- umap.dt %>%
      merge(sample_metadata,by="cell")
    
    for (k in args$colour_by) {
      
      p <- ggplot(to.plot, aes_string(x="umap1", y="umap2", fill=k)) +
        ggrastr::geom_point_rast(size=1.5, shape=21, stroke=0.05) +
        theme_classic() +
        theme(
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()
        )
      
      if (k%in%c("celltype.mapped","celltype.predicted")) {
        p <- p + scale_fill_manual(values=opts$celltype.colors) +
          theme(
            legend.position="none",
            legend.title=element_blank()
          )
      }
      
      # Save UMAP plot
      outfile <- sprintf("%s/%s_umap_nfeatures%d_ntopics%d_neigh%d_dist%s_%s.pdf",io$pdfdir, paste(args$samples,collapse="-"), args$nfeatures, args$ntopics, i, j, k)
      pdf(outfile, width=7, height=5)
      print(p)
      dev.off()
    }
    
    # Save UMAP coordinates
    outfile <- sprintf("%s/%s_umap_nfeatures%d_ntopics%d_neigh%d_dist%s.txt.gz",args$outdir, paste(args$samples,collapse="-"), args$nfeatures, args$ntopics, i, j)
    fwrite(umap.dt, outfile)
  }
}



############################################
## Boxplots of topic values per cell type ##
############################################

# # Plot topic score per cell type
# topic.dt <- topics.mtx %>% as.data.table(keep.rownames = T) %>% 
#   setnames("rn","cell") %>%
#   melt(id.vars=c("cell","celltype.mapped"), variable.name="topic") %>%
#   setnames("celltype.mapped","celltype")

# for (i in unique(topic.dt$topic)) {
#   to.plot <- topic.dt[topic==i] %>% 
#     .[,celltype:=factor(celltype,levels=names(opts$celltype.colors))]
#   p <- ggboxplot(to.plot, x="celltype", y="value", fill="celltype", outlier.shape=NA) +
#     scale_fill_manual(values=opts$celltype.colors) +
#     geom_hline(yintercept=0, linetype="dashed") +
#     labs(x="", y=sprintf("%s z-score",i)) +
#     theme(
#       legend.position = "none",
#       axis.text.x = element_text(color="black", angle=40, hjust=1, size=rel(0.8)),
#       axis.text.y = element_text(color="black", size=rel(0.8))
#     )
#   png(sprintf("%s/%s_%s_%s_boxplot.png",args$outdir,paste(args$samples,collapse="-"),args$nfeatures,i), width = 800, height = 400)
#   print(p)
#   dev.off()
# }