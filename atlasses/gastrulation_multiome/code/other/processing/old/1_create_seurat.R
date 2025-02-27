suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(argparse))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--inputdir',        type="character",                    help='Input directory')
p$add_argument('--outputdir',       type="character",                    help='Output directory')
p$add_argument('--samples',         type="character",       nargs="+",   help='Samples')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
}

## START TEST ##
args <- list()
args$inputdir <- paste0(io$basedir,"/original")
args$outputdir <- paste0(io$basedir,"/processed/rna")
args$samples <- opts$samples
## END TEST ##

###############
## Load data ##
###############

stopifnot(args$samples%in%opts$samples)
files.10x <- sprintf("%s/%s/filtered_feature_bc_matrix",args$inputdir,args$samples)
names(files.10x) <- args$samples

counts <- Read10X(files.10x, strip.suffix = TRUE)
lapply(counts,dim)

###################
## Create Seurat ##
###################

seurat <- CreateSeuratObject(
  counts = counts["Gene Expression"][[1]],
  project = "Gastrulation Multiome 10x"
)

# Add ATAC modality
# seurat[["ATAC"]] <- CreateAssayObject(counts = counts["Peaks"][[1]])

seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^mt-")

seurat

#####################
## Create metadata ##
#####################

metadata <- seurat@meta.data %>% as.data.table %>% .[,orig.ident:=NULL] %>%
  .[,c("cell","barcode","sample","nFeature_RNA","nCount_RNA","percent.mt")]

head(metadata)

# metadata <- seurat@meta.data %>%
#   tibble::rownames_to_column("barcode") %>%
#   as.data.table %>%
#   .[,orig.ident:=NULL]
# fwrite(metadata, io$metadata, sep="\t", quote=F)


##########
## Save ##
##########

saveRDS(sratseurat, paste0(args$outputdir,"/seurat.rds"))
fwrite(metadata, paste0(args$outputdir,"/metadata.txt.gz"), quote=F, na="NA", sep="\t")