suppressPackageStartupMessages(library(furrr))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(cisTopic))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--samples',         type="character",                nargs='+',     help='Samples')
p$add_argument('--nfeatures',       type="integer",    default=1000,               help='Number of features')
p$add_argument('--ntopics',      type="integer",      nargs="+",    help='Number of topics')
p$add_argument('--batch.variable',  type="character",                               help='Metadata column to apply batch correction on')
p$add_argument('--batch.method',    type="character",  default="MNN",               help='Batch correctin method ("Harmony" or "MNN")')
p$add_argument('--n_neighbors',     type="integer",    default=30,   nargs='+',     help='(UMAP) Number of neighbours')
p$add_argument('--min_dist',        type="double",     default=0.3,  nargs='+',     help='(UMAP) Minimum distance')
p$add_argument('--remove_ExE_celltypes', action="store_true",   help='Remove ExE cell types?')
p$add_argument('--colour_by',       type="character",  default="celltype.mapped",  nargs='+',  help='Metadata columns to colour the UMAP by')
p$add_argument('--ncores',       type="integer",      default=1,    help='Number of cores')
p$add_argument('--seed',            type="integer",    default=42,                  help='Random seed')
p$add_argument('--outdir',          type="character",                               help='Output directory')
p$add_argument('--test_mode',    action="store_true",               help='Test mode? subset number of cells')
args <- p$parse_args(commandArgs(TRUE))

## START TEST
# args$samples <- opts$samples[1:2]
# args$nfeatures <- 1e4
# args$ntopics <- c(10,20)
# # args$batch.variable <- "stage"
# # args$batch.method <- "Harmony"
# args$colour_by <- c("celltype.mapped","sample","log_nFrags_atac","doublet_call")
# args$ncores <- 1
# args$test_mode <- TRUE
# args$outdir <- paste0(io$basedir,"/results/atac/cistopic/test")

## END TEST

########################
## Load ArchR project ##
########################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
} else {
  stop("Computer not recognised")
}

#####################
## Define settings ##
#####################

# I/O
# io$metadata <- paste0(io$basedir, "/sample_metadata.txt.gz")
io$pdfdir <- sprintf("%s/pdf",args$outdir); dir.create(io$pdfdir,showWarnings = F)

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & doublet_call==FALSE & sample%in%args$samples] %>%
  .[,log_nFrags_atac:=log10(nFrags_atac)]

if (args$remove_ExE_celltypes) {
  sample_metadata <- sample_metadata %>% .[!celltype.predicted%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}


if (args$test_mode) {
  sample_metadata <- head(sample_metadata,n=100)
  args$ntopics <- c(10,20)
}

table(sample_metadata$sample)

###################
## Sanity checks ##
###################

stopifnot(args$colour_by %in% colnames(sample_metadata))
# stopifnot(unique(sample_metadata$celltype.mapped) %in% names(opts$celltype.colors))

if (length(args$batch.variable)>0) {
  stopifnot(args$batch.variable%in%colnames(sample_metadata))
  if (length(unique(sample_metadata[[args$batch.variable]]))==1) {
    message(sprintf("There is a single level for %s, no batch correction applied",args$batch.variable))
    args$batch.variable <- NULL
  } else {
    library(batchelor)
  }
}

##################
## Subset ArchR ##
##################

ArchRProject.filt <- ArchRProject[sample_metadata$cell]

# ArchRProject.filt@sampleColData <- ArchRProject.filt@sampleColData[args$samples,,drop=F]

# Update ArchR metadata
# stopifnot(sample_metadata$cell == rownames(getCellColData(ArchRProject.filt)))
# ArchRProject.filt <- addCellColData(ArchRProject.filt,
#   data = sample_metadata$stage, 
#   name = "stage",
#   cells = sample_metadata$cell,
#   force = TRUE
# )

######################
## Load peak matrix ##
######################

# Get Peak Matrix from ArchR
atac.peak.se <- getMatrixFromProject(ArchRProject.filt, useMatrix="PeakMatrix", binarize = TRUE)
dim(atac.peak.se)

# Define peak names
peak_names <- rowRanges(atac.peak.se) %>% as.data.table %>% .[,id:=sprintf("%s:%s-%s",seqnames,start,end)] %>% .$id
rownames(atac.peak.se) <- peak_names

#######################
## Feature selection ##
#######################

# Load peak variability estimates
peakStats.dt <- fread(io$archR.peak.stats)

# Define highly variable peaks
peaks <- peakStats.dt %>% 
  setorder(-var_pseudobulk) %>% 
  head(n=args$nfeatures) %>% 
  .[,peak:=stringr::str_replace(peak,"_",":")] %>%
  .[,peak:=stringr::str_replace(peak,"_","-")] %>%
  .$peak

# Subset
atac.peak.se.filt <- atac.peak.se[peaks,]
dim(atac.peak.se.filt)
# rm(atac.peak.se)

##############
## CisTopic ##
##############

# Create cisTopicObject
cisTopicObject <- createcisTopicObject(
  count.matrix = assay(atac.peak.se.filt,"PeakMatrix"), 
  project.name = sprintf('Multiome10x_%s',paste(args$samples,collapse="-")), 
  keepCountsMatrix = FALSE
)


# Run Latent Dirichlet Allocation
# cisTopicObject <- runCGSModels(cisTopicObject, topic, burnin = 125, iterations = 250, nCores = 2, seed = 123)

# Run Latent Dirichlet Allocation with WarpLDA
cisTopicObject <- runWarpLDAModels(
  object = cisTopicObject, 
  topic = args$ntopics, 
  iterations = 500, 
  nCores = args$ncores, 
  seed = args$seed,
  returnType = "allModels"
)

#####################
## Model selection ##
#####################

# Returns a cisTopic object (when the input is a cisTopic object) with the selected model 
# stored in object@selected.model, and the log likelihoods of the models in object@log.lik
cisTopicObject.best_model <- selectModel(
  cisTopicObject, 
  type = "maximum", 
  keepBinaryMatrix = F, 
  keepModels = F
)

# Create data.table 
lik <- cisTopicObject@models %>% map_dbl(function(x) x$log.likelihoods[2,ncol(x$log.likelihoods)])
lik.dt <- data.table(
  ntopics = names(lik),
  lik = lik
)

# Retrieve normalised topic-cell assignments
topics.mtx <- modelMatSelection(cisTopicObject.best_model, 'cell', "Z-score") %>% t 

######################
## Batch correction ##
######################

if (length(args$batch.variable)>0) {
  print(sprintf("Applying %s batch correction for variable: %s", args$batch.method, args$batch.variable))
  stop("Not implemented")
}

##########
## Save ##
##########

# outfile <- sprintf("%s/%s_lsi_features%d_dims%d_%sbatchcorrectionby%s.txt.gz",args$outdir, paste(args$samples,collapse="-"), args$nfeatures, args$ntopics, args$batch.method, paste(args$batch.variable,collapse="-"))
outfile <- sprintf("%s/%s_cistopic_features%d_ntopics%d.txt.gz",args$outdir, paste(args$samples,collapse="-"), args$nfeatures, args$ntopics)
fwrite(round(topics.mtx,3) %>% as.data.table(keep.rownames=T) %>% setnames("rn","cell"), outfile, sep="\t")

# save all models
# io$outfile.allmodels <- sprintf("%s/%s_%s_cistopic_warpLDA_allModels.rds",args$outdir,paste(args$samples,collapse="-"),args$nfeatures)
# saveRDS(cisTopicObject, io$outfile.allmodels)

# save best model
outfile <- sprintf("%s/%s_cistopic_bestModel_features%d_ntopics%d.rds",args$outdir,paste(args$samples,collapse="-"),args$nfeatures, args$ntopics)
saveRDS(cisTopicObject.best_model,outfile)



#########################
## Plot cisTopic stats ##
#########################


# Plot model selection
# pdf(sprintf("%s/%s_model_selection_features%d_ntopics%d.pdf",io$pdfdir, paste(args$samples,collapse="-"), args$nfeatures, ncol(topics.mtx)), width=7, height=3)
# ggline(lik.dt, x="ntopics", y="lik")
# dev.off()

# Plot correlation between topics
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

    # Plot
    to.plot <- umap.dt %>% merge(sample_metadata,by="cell")

    for (k in args$colour_by) {
      
      p <- ggplot(to.plot, aes_string(x="umap1", y="umap2", fill=k)) +
        geom_point(size=1.5, shape=21, stroke=0.05) +
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
      outfile <- sprintf("%s/%s_umap_nfeatures%d_ntopics%d_neigh%d_dist%s_%s.pdf",io$pdfdir, paste(args$samples,collapse="-"), args$nfeatures, ncol(topics.mtx), i, j, k)
      pdf(outfile, width=7, height=5)
      print(p)
      dev.off()
    }
    
    # Save UMAP coordinates
    outfile <- sprintf("%s/%s_umap_nfeatures%d_ntopics%d_neigh%d_dist%s.txt.gz",args$outdir, paste(args$samples,collapse="-"), args$nfeatures, ncol(topics.mtx), i, j)
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

