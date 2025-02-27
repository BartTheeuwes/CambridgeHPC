# TO DO:
# Include filtering for non-informative genes

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(argparse))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--sce',             type="character",                               help='SingleCellExperiment file')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--features',        type="integer",    default=1000,                help='Number of features')
p$add_argument('--npcs',            type="integer",    default=30,                  help='Number of PCs')
p$add_argument('--n_neighbors',     type="integer",    default=30,     help='(UMAP) Number of neighbours')
p$add_argument('--min_dist',        type="double",     default=0.3,     help='(UMAP) Minimum distance')
p$add_argument('--colour_by',       type="character",  default="celltype.mapped_mnn",  nargs='+',  help='Metadata columns to colour the UMAP by')
p$add_argument('--remove_ExE_cells', action="store_true",                                 help='Remove ExE cells?')
p$add_argument('--seed',            type="integer",    default=42,                  help='Random seed')
p$add_argument('--outdir',          type="character",                               help='Output file')
p$add_argument('--samples',         type="character",                nargs='+',     help='Samples')
p$add_argument('--vars_to_regress', type="character",                nargs='+',     help='Metadata columns to regress out')
p$add_argument('--batch_correction',type="character",                               help='Metadata column to apply batch correction on')
# p$add_argument('--test',      action = "store_true",                       help='Testing mode')


args <- p$parse_args(commandArgs(TRUE))
print(class(args$vars_to_regress))
print(args$vars_to_regress)

#####################
## Define settings ##
#####################
here::i_am("rna/processing/1_create_seurat_rna.R")
source(here::here("settings.R"))
source(here::here("utils.R"))

dir.create(args$outdir, recursive=TRUE, showWarnings = FALSE)

## START TEST ##
# args$samples <- opts$samples
# args$sce <- io$rna.sce
# args$metadata <- io$metadata
# # args$metadata <- paste0(io$basedir,"/results/rna/doublets/sample_metadata_after_doublets.txt.gz")
# args$features <- 1000
# args$npcs <- 30
# args$test <- FALSE
# args$colour_by <- c("celltype.mapped","sample")
# args$vars.to.regress <- c("nFeature_RNA","mitochondrial_percent_RNA")
# args$batch.correction <- c("sample")
# args$remove_ExE_cells <- FALSE
# args$outdir <- paste0(io$basedir,"/results/rna/dimensionality_reduction/test")
## END TEST ##

# if (isTRUE(args$test)) print("Test mode activated...")

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
   .[pass_rnaQC==TRUE & doublet_call==FALSE & sample%in%args$samples]


if (args$remove_ExE_cells) {
  print("Removing ExE cells...")
  sample_metadata <- sample_metadata %>%
    .[!celltype.mapped_mnn%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}


###################
## Sanity checks ##
###################

stopifnot(args$colour_by %in% colnames(sample_metadata))
# stopifnot(unique(sample_metadata$celltype.mapped) %in% names(opts$celltype.colors))

 if (args$batch_correction %in% c('day', 'sample', 'day_sample', 'replicate')) {
     source(here::here("rna/mapping/run/mnn/mapping_functions_extended.R"))
     library(batchelor)
   }

 if (length(args$vars_to_regress)>0) {
  stopifnot(args$vars_to_regress%in%colnames(sample_metadata))
 }


###############
## Load data ##
###############

# Load RNA expression data as SingleCellExperiment object
sce <- load_SingleCellExperiment(args$sce, cells=sample_metadata$cell, normalise = TRUE)

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

colData(sce)$day = factor(colData(sce)$day, levels=c(
                                        "D5",
                                        "D4.5", 
                                        "D4", 
                                        "D3.5", 
                                        "D3"
                                          )) # add factor levels to days

#######################
## Feature selection ##
#######################
## Find HVGs - detection on WT samples only
# Get gene metadata
gene_metadata <- fread(io$gene_metadata) %>% .[,c("chr","ens_id","symbol")] %>%
  .[symbol!="" & ens_id%in%rownames(sce)] %>%
  .[!duplicated(symbol)]

# Imprinted genes
imprint = gene_metadata[c(grep('maternally', gene_metadata$description),
                       grep('paternally', gene_metadata$description)), symbol]
# Other imprinted genes: 
#- Nnat (https://www.genecards.org/cgi-bin/carddisp.pl?gene=NNAT)
#- Grb10 (https://www.genecards.org/cgi-bin/carddisp.pl?gene=GRB10)

# decomp <- modelGeneVar(sce, block=colData(sce)$sample) # Only detect HVGs from WT samples
# decomp <- decomp[decomp$mean > 0.01,]
# hvgs <- decomp[order(decomp$FDR),] %>% 
#     as.data.table(., keep.rownames=T) %>% 
#     .[grep("^Rik|Rik$|^mt-|^Rps|^Rpl|^Gm",rn,invert=T)] %>% # filter out non-informative genes
#     .[grep("^Hbb|^Hba",rn,invert=T)] %>% # Removing Haem genes due to spilling over from erythroids to other cells 
#     .[!rn %in% c(imprint, 'Grb10', 'Nnat')] %>%  # remove imprinted genes
#     .[!rn %in% c("Xist", "Tsix")] %>%  # remove Xist & Tsix
#     .[!rn=="tomato-td"] %>% # remove tomato itself
#     .[!rn%in%gene_metadata[chr=="chrY",symbol]] %>%
#      head(n=args$features) %>% .$rn

# Filter some genes manually
genes_keep = rownames(sce)
genes_keep <- genes_keep[grep("^Rik|Rik$|^mt-|^Rps|^Rpl|^Gm",genes_keep,invert=T)] # filter out non-informative genes
genes_keep <- genes_keep[grep("^Hbb|^Hba",genes_keep,invert=T)] # test removing Haem genes 
genes_keep <- genes_keep[!genes_keep %in% c(imprint, 'Grb10', 'Nnat')] # remove imprinted genes
#genes_keep <- genes_keep[!genes_keep %in% c("Xist", "Tsix")] # remove Xist & Tsix
#genes_keep <- genes_keep[!genes_keep == "tomato-td"] # remove tomato itself
#genes_keep <- genes_keep[!genes_keep %in% gene_metadata[chr=="chrY",symbol]] # no genes on y-chr 

sce_wt <- sce[genes_keep, sample_metadata[genotype=='WT', cell]]

# Find variable genes using Seurat
# Also try with only using WT for HVG detection
hvgs = VariableFeatures(FindVariableFeatures(as.Seurat(sce_wt), nfeatures = args$features))
sce_filt <- sce[hvgs,]

############################
## Regress out covariates ##
############################

 if (length(args$vars_to_regress)>0) {
   print(sprintf("Regressing out variables: %s", paste(args$vars_to_regress,collapse=" ")))
   logcounts(sce_filt) <- RegressOutMatrix(
     mtx = logcounts(sce_filt), 
     covariates = colData(sce_filt)[,args$vars_to_regress,drop=F]
   )
 }

############################
## PCA + Batch correction ##
############################

# Not the most efficient way lol but it works
 if (args$batch_correction=='day') {
    print('batch correcting by day')
    pca <- multiBatchPCA(sce_filt, batch = colData(sce_filt)[[args$batch_correction]], d = args$npcs)
    pca.corrected <- reducedMNN(pca)$corrected
    colnames(pca.corrected) <- paste0("PC",1:ncol(pca.corrected))
    pca.corrected = pca.corrected[match(colnames(sce_filt), rownames(pca.corrected)),]
    reducedDim(sce_filt, "PCA") <- pca.corrected
 } else if (args$batch_correction=='sample') {
    print('batch correcting by sample')
    pca <- multiBatchPCA(sce_filt, batch = colData(sce_filt)[[args$batch_correction]], d = args$npcs)
    pca.corrected <- reducedMNN(pca)$corrected
    colnames(pca.corrected) <- paste0("PC",1:ncol(pca.corrected))
    pca.corrected = pca.corrected[match(colnames(sce_filt), rownames(pca.corrected)),]
    reducedDim(sce_filt, "PCA") <- pca.corrected
 } else if (args$batch_correction=='replicate') {
    print('batch correcting by replicate')
    pca <- multiBatchPCA(sce_filt, batch = colData(sce_filt)[[args$batch_correction]], d = args$npcs)
    pca.corrected <- reducedMNN(pca)$corrected
    colnames(pca.corrected) <- paste0("PC",1:ncol(pca.corrected))
    pca.corrected = pca.corrected[match(colnames(sce_filt), rownames(pca.corrected)),]
    reducedDim(sce_filt, "PCA") <- pca.corrected
 } else if (args$batch_correction=='day_sample') {
    print('batch correcting by day and sample')
    meta = colData(sce_filt)
    #get order: oldest to youngest; most cells to least cells
    order_df = meta[!duplicated(meta$sample), c("day", "sample")]
    order_df$ncells = sapply(order_df$sample, function(x) sum(meta$sample == x))
    order_df$day = factor(order_df$day, 
                            levels = rev(c(
                                        "D5",
                                        "D4.5", 
                                        "D4", 
                                        "D3.5", 
                                        "D3"
                                          )))
    order_df = order_df[order(order_df$day, order_df$ncells, decreasing = TRUE),]
    order_df$day = as.character(order_df$day)

    all_correct = doBatchCorrect(counts = logcounts(sce_filt), 
                                 timepoints = meta$day, 
                                 samples = meta$sample, 
                                 timepoint_order = order_df$day, 
                                 sample_order = order_df$sample, 
                                 npc = args$npcs)
                         
   colnames(all_correct) <- paste0("PC",1:ncol(all_correct))
   reducedDim(sce_filt, "PCA") <- all_correct
 } else {
    print('PCA without batch correction')
    sce_filt <- runPCA(sce_filt, ncomponents = args$npcs, ntop=args$features)  
 }

# Save PCA coordinates
 outfile <- sprintf("%s/pca_features%d_pcs%d.txt.gz",args$outdir, args$features, args$npcs)
 pca.dt <- reducedDim(sce_filt,"PCA") %>% round(3) %>% as.data.table(keep.rownames = T) %>% setnames("rn","cell")
 fwrite(pca.dt, outfile)

##########
## UMAP ##
##########

# Run
set.seed(args$seed)
sce_filt <- runUMAP(sce_filt, dimred="PCA", n_neighbors = args$n_neighbors, min_dist = args$min_dist)

# Fetch UMAP coordinates
umap.dt <- reducedDim(sce_filt,"UMAP") %>% as.data.table %>% 
  .[,cell:=colnames(sce_filt)] %>%
  setnames(c("UMAP1","UMAP2","cell"))

# Save UMAP coordinates
fwrite(umap.dt, sprintf("%s/umap_features%d_pcs%d_neigh%d_dist%s.txt.gz",args$outdir, args$features, args$npcs, args$n_neighbors, args$min_dist))

##########
## Plot ##
##########
to.plot <- reducedDim(sce_filt,"UMAP") %>% as.data.table %>% 
    .[,cell:=colnames(sce_filt)] %>%
    merge(sample_metadata, by="cell")
    
for (i in args$colour_by) {

  p <- ggplot(to.plot, aes_string(x="V1", y="V2", col=i)) +
    geom_point(size=0.2) +
    theme_classic() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )

  if (grepl("celltype.mapped",i)) {
    p <- p + scale_color_manual(values=opts$celltype.colors) +
      theme(
        legend.position="none",
        legend.title=element_blank()
      )
  }
  
  if (grepl("celltype_extended.mapped",i)) {
    p <- p + scale_color_manual(values=opts$celltype_extended.colors) +
      theme(
        legend.position="none",
        legend.title=element_blank()
      )
  }
    
  if (grepl("day",i)) {
    p <- p + scale_color_manual(values=opts$days.colors) +
      theme(
        legend.position="none",
        legend.title=element_blank()
      )
  }
    
  if (grepl("genotype",i)) {
    p <- p + scale_color_manual(values=opts$genotype.colors) +
      theme(
        legend.position="none",
        legend.title=element_blank()
      )
  }
    
  if (class(to.plot[[i]]) %in% c('integer', 'numeric')) {
    p <- p + scale_colour_gradientn(colors=rainbow(4)) +
      theme(
        legend.title=element_blank()
      )
  }      

  # Save UMAP plot
  outfile <- file.path(args$outdir,sprintf("umap_features%d_pcs%d_neigh%d_dist%s_%s.pdf",args$features, args$npcs, args$n_neighbors, args$min_dist, i))
  pdf(outfile, width=7, height=5)
  print(p)
  dev.off()
}

