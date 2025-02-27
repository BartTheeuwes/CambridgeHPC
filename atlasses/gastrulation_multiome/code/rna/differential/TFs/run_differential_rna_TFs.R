suppressMessages(library(SingleCellExperiment))
suppressMessages(library(scater))
suppressMessages(library(edgeR))
suppressMessages(library(argparse))

## Initialize argument parser ##
p <- ArgumentParser(description='')
p$add_argument('--groupA',    type="character",    help='group A')
p$add_argument('--groupB',    type="character",    help='group B')
# p$add_argument('--stages',    type="character",    default="all", nargs="+", help='Stages to use')
p$add_argument('--group_label',    type="character",    help='Group label')
p$add_argument('--test_mode', action="store_true", help='Test mode? subset number of cells')
p$add_argument('--outfile',   type="character",    help='Output file')
args <- p$parse_args(commandArgs(TRUE))

## START TEST
# args$groupA <- "Erythroid1"
# args$groupB <- "Neural_crest"
# args$group_label <- "celltype.mapped_mnn"
# args$test_mode <- TRUE
# args$stages <- c("E8.5")
# args$outfile <- c("/hps/nobackup2/research/stegle/users/ricard/gastrulation_multiome_10x/results/differential/foo.tsv.gz")
## END TEST

#####################
## Define settings ##
#####################

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))
source(here::here("rna/differential/utils.R"))

# I/O
io$TFs <- file.path(io$basedir,"processed_new/atac/archR/Annotations/CISBP_TFs.txt.gz")

# Define groups
opts$groups <- c(args$groupA,args$groupB)

# Define FDR threshold
opts$threshold_fdr <- 0.01

# Define minimum logFC for significance
opts$min.logFC <- 1

# For a given gene, the minimum fraction of cells that must express it in at least one group
# Note that this threshold has to be lower for TFs than for genes because they are less expressed
opts$min_detection_rate_per_group <- 0.25

# Sanity checks
# stopifnot(args$stages%in%opts$stages)
stopifnot(args$groupA%in%opts$celltypes)
stopifnot(args$groupB%in%opts$celltypes)

######################
## Load list of TFs ##
######################

TFs <- fread(io$TFs)[["gene"]]

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE]

stopifnot(args$group_label%in%colnames(sample_metadata))

sample_metadata <- sample_metadata %>%
  setnames(args$group_label,"group") %>%
  # .[,group:=eval(as.name(args$group_label))] %>%
  .[group%in%c(args$groupA,args$groupB)] %>%
  .[,group:=factor(group,levels=opts$groups)] %>% setorder(group) # Sort cells so that groupA comes before groupB

if (isTRUE(args$test_mode)) {
  print("Testing mode activated")
  sample_metadata <- sample_metadata %>% split(.,.$group) %>% map(~ head(.,n=250)) %>% rbindlist
}

table(sample_metadata$group)

#########################
## Load RNA expression ##
#########################

# Load SingleCellExperiment object
rna.sce <- load_SingleCellExperiment(
  file = io$rna.sce, 
  normalise = TRUE, 
  cells = sample_metadata$cell
)
rna.sce$group <- sample_metadata$group

# Load gene metadata
gene_metadata.dt <- fread(io$gene_metadata) %>%
  .[symbol%in%rownames(rna.sce)] %>%
  .[,c("symbol","ens_id")] %>%
  setnames("symbol","gene")

################
## Subset TFs ##
################

TFs <- intersect(TFs,toupper(rownames(rna.sce)))
rna_tf.sce <- rna.sce[stringr::str_to_title(TFs),]
rownames(rna_tf.sce) <- toupper(rownames(rna_tf.sce))

rm(rna.sce)

# Filter genes
rna_tf.sce <- rna_tf.sce[rownames(rna_tf.sce)%in%toupper(gene_metadata.dt$gene),]

gene_metadata_tf.dt <- gene_metadata.dt %>%
  .[,gene:=toupper(gene)] %>%
  .[gene%in%rownames(rna_tf.sce)]

################
## Parse data ##
################

# calculate detection rate per gene
cdr.dt <- data.table(
  rownames(rna_tf.sce),
  rowMeans(logcounts(rna_tf.sce[,rna_tf.sce$group==opts$groups[1]])>0) %>% round(2),
  rowMeans(logcounts(rna_tf.sce[,rna_tf.sce$group==opts$groups[2]])>0) %>% round(2)
) %>% setnames(c("gene",sprintf("detection_rate_%s",opts$groups[1]),sprintf("detection_rate_%s",opts$groups[2])))
# .[,cdr_diff:=abs(out[,(sprintf("detection_rate_%s",opts$groups[1])),with=F][[1]] - out[,(sprintf("detection_rate_%s",opts$groups[2])),with=F][[1]])] %>%

################################################
## Differential expression testing with edgeR ##
################################################

out <- doDiffExpr(rna_tf.sce, opts$groups, opts$min_detection_rate_per_group) %>%
  # Add sample statistics
  .[,c("groupA_N","groupB_N"):=list(table(sample_metadata$group)[1],table(sample_metadata$group)[2])]%>% 
  # setnames(c("groupA_N","groupB_N"),c(sprintf("N_%s",opts$groups[1]),sprintf("N_%s",opts$groups[2]))) %>%
  # Add gene statistics
  merge(cdr.dt, all.y=T, by="gene") %>%
  # merge(gene_metadata_tf.dt, all.y=T, by="gene") %>%
  # Calculate statistical significance
  # .[, sig := (padj_fdr<=opts$threshold_fdr & abs(logFC)>=opts$min.logFC)] %>%
  # .[is.na(sig),sig:=FALSE] %>%
  setorder(padj_fdr, na.last=T)

# Parse columns
out[,c("p.value","padj_fdr","logFC","log_padj_fdr"):=list(signif(p.value,digits=3), signif(padj_fdr,digits=3), round(logFC,3),round(log_padj_fdr,3))]

##################
## Save results ##
##################

fwrite(out, args$outfile, sep="\t", na="NA", quote=F)
