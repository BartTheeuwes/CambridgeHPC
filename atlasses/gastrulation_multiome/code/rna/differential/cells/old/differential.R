here::i_am("rna/differential/cells/differential.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

suppressMessages(library(edgeR))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',    type="character",    help='Cell metadata file')
p$add_argument('--sce',    type="character",    help='SingleCellExperiment file')
p$add_argument('--groupA',    type="character",    help='group A')
p$add_argument('--groupB',    type="character",    help='group B')
p$add_argument('--samples',    type="character",    default="all", nargs="+", help='Samples to use')
p$add_argument('--celltypes',    type="character",    default="all", nargs="+", help='Celltypes to use')
p$add_argument('--min_cells',       type="integer",       default=5,      help='Minimum number of cells per group')
p$add_argument('--group_variable',    type="character",    help='Group label')
p$add_argument('--test_mode', action="store_true", help='Test mode? subset number of cells')
p$add_argument('--outfile',   type="character",    help='Output file')
args <- p$parse_args(commandArgs(TRUE))

## START TEST
# args$metadata <- io$metadata
# args$sce <- io$rna.sce
# args$groupA <- "Nascent_mesoderm"
# args$groupB <- "Neural_crest"
# args$group_variable <- "celltype.mapped"
# args$celltypes <- "all"
# args$samples <- "all"
# args$min_cells <- 5
# args$test_mode <- FALSE
# args$outfile <- NULL
## END TEST

#####################
## Define settings ##
#####################

# Load utils
source(here::here("rna/differential/utils.R"))

# Sanity checks
# stopifnot(args$samples%in%opts$samples)
# stopifnot(args$groupA%in%opts$celltypes)
# stopifnot(args$groupB%in%opts$celltypes)

# Define cell types
if (args$celltypes[1]=="all") {
  args$celltypes <- opts$celltypes
} else {
  stopifnot(args$celltypes%in%c(opts$celltypes,"Erythroid","Blood_progenitors"))
}

# Define groups
opts$groups <- c(args$groupA,args$groupB)

# Define FDR threshold
opts$threshold_fdr <- 0.01

# Define minimum logFC for significance
opts$min.logFC <- 1

# For a given gene, the minimum fraction of cells that must express it in at least one group
opts$min_detection_rate_per_group <- 0.40

# Rename celltypes
# opts$rename_celltypes <- c(
#   "Erythroid3" = "Erythroid",
#   "Erythroid2" = "Erythroid",
#   "Erythroid1" = "Erythroid",
#   "Blood_progenitors_1" = "Blood_progenitors",
#   "Blood_progenitors_2" = "Blood_progenitors"
#   # "Intermediate_mesoderm" = "Mixed_mesoderm",
#   # "Paraxial_mesoderm" = "Mixed_mesoderm",
#   # "Nascent_mesoderm" = "Mixed_mesoderm",
#   # "Pharyngeal_mesoderm" = "Mixed_mesoderm"
#   # "Visceral_endoderm" = "ExE_endoderm"
# )

########################
## Load cell metadata ##
########################

sample_metadata <- fread(args$metadata) %>%
  # .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,opts$rename_celltypes)] %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE & celltype%in%args$celltypes]

stopifnot(args$group_variable%in%colnames(sample_metadata))

sample_metadata <- sample_metadata %>%
  setnames(args$group_variable,"group") %>%
  # .[,group:=eval(as.name(args$group_variable))] %>%
  .[group%in%opts$groups] %>%
  .[,group:=factor(group,levels=opts$groups)] %>% setorder(group) # Sort cells so that groupA comes before groupB

if (isTRUE(args$test_mode)) {
  print("Testing mode activated")
  sample_metadata <- sample_metadata %>% split(.,.$group) %>% map(~ head(.,n=250)) %>% rbindlist
}

table(sample_metadata$group)

# TO-DO: CHECK FOR MINIMUM NUMBER OF CELLS PER GROUP
# args$min_cells

#########################
## Load RNA expression ##
#########################

# Load SingleCellExperiment object
sce <- load_SingleCellExperiment(
  file = args$sce, 
  normalise = TRUE, 
  cells = sample_metadata$cell
)
sce$group <- sample_metadata$group

# Load gene metadata
# gene_metadata <- fread(io$gene_metadata) %>%
#   .[symbol%in%rownames(sce)] %>%
#   .[,c("symbol","ens_id")] %>%
#   setnames("symbol","gene")

################
## Parse data ##
################

# calculate detection rate per gene
cdr.dt <- data.table(
  rownames(sce),
  rowMeans(logcounts(sce[,sce$group==opts$groups[1]])>0) %>% round(2),
  rowMeans(logcounts(sce[,sce$group==opts$groups[2]])>0) %>% round(2)
) %>% setnames(c("gene",sprintf("detection_rate_%s",opts$groups[1]),sprintf("detection_rate_%s",opts$groups[2])))
# .[,cdr_diff:=abs(out[,(sprintf("detection_rate_%s",opts$groups[1])),with=F][[1]] - out[,(sprintf("detection_rate_%s",opts$groups[2])),with=F][[1]])] %>%

# Filter genes
# sce <- sce[rownames(sce)%in%gene_metadata$gene,]

################################################
## Differential expression testing with edgeR ##
################################################

out <- doDiffExpr(sce, opts$groups, opts$min_detection_rate_per_group) %>%
  .[,c("groupA_N","groupB_N"):=list(table(sample_metadata$group)[1],table(sample_metadata$group)[2])]%>% 
  merge(cdr.dt, all.y=T, by="gene") %>%
  setorder(padj_fdr, na.last=T)

# Parse columns
out[,c("padj_fdr","logFC"):=list(signif(padj_fdr,digits=3), round(logFC,3))]

##################
## Save results ##
##################

fwrite(out, args$outfile, sep="\t", na="NA", quote=F)
