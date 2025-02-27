here::i_am("rna/processing/1_create_seurat_rna.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(Seurat))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--inputdir',        type="character",                    help='Input directory')
p$add_argument('--samples',         type="character",       nargs="+",   help='Samples')
p$add_argument('--min_reads',           type="integer",    default=1000,                  help='Minimum number of reads')
p$add_argument('--outdir',       type="character",                    help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

# These settings are important for consistency with ArchR, which provides little flexibility to edit cell names
opts$trim.barcode <- FALSE
opts$sample_cell_separator <- "#"

## START TEST ##
# args <- list()
# args$inputdir <- paste0(io$basedir,"/original")
# args$samples <- opts$samples[1:2]
# args$min_reads <- 1000
# args$outdir <- paste0(io$basedir,"/processed/rna")
## END TEST ##

##############################
## Load and merge data sets ##
##############################

stopifnot(args$samples%in%opts$samples)

count_mtx <- list()
cell.info <- list()

for (i in args$samples) {
  print(i)
    
  # Load cell metadata
  barcode.loc <- sprintf("%s/%s/outs/unfiltered_feature_bc_matrix/barcodes.tsv.gz",args$inputdir,i)
  cell.info[[i]] <- fread(barcode.loc, header=F) %>%
    setnames("barcode") %>%
    .[,barcode:=ifelse(rep(opts$trim.barcode,.N),gsub("-1","",barcode),barcode)] %>%
    .[,c("sample","cell"):=list(i,sprintf("%s%s%s",i,opts$sample_cell_separator,barcode))]
  dim(cell.info[[i]])
  
  # Load matrix  
  count_mtx[[i]] <- Read10X(sprintf("%s/%s/outs/unfiltered_feature_bc_matrix",args$inputdir,i))[["Gene Expression"]]

  # Basic filtering
  count_mtx[[i]] <- count_mtx[[i]][,colSums(count_mtx[[i]])>=args$min_reads]
  cell.info[[i]] <- cell.info[[i]][barcode%in%colnames(count_mtx[[i]])] %>% setkey(barcode) %>% .[colnames(count_mtx[[i]])]
  stopifnot(sum(is.na(cell.info[[i]]$cell))==0)
  stopifnot(sum(is.na(cell.info[[i]]$barcode))==0)

}

print(lapply(count_mtx,dim))

#######################
## Keep common genes ##
#######################

genes <- Reduce("intersect",lapply(count_mtx,rownames))
for (i in 1:length(count_mtx)) {
  count_mtx[[i]] <- count_mtx[[i]][genes,]
}

stopifnot(length(unique(lapply(count_mtx,nrow)))==1)
stopifnot(length(unique(lapply(count_mtx,rownames)))==1)

#################
## Concatenate ##
#################

# Concatenate cell metadata
cell.info <- rbindlist(cell.info)
rownames(cell.info) <- cell.info$cell

# Concatenate matrices
count_mtx <- do.call("cbind",count_mtx)
colnames(count_mtx) <- cell.info$cell

##################
## Filter genes ##
##################

# Remove duplicated genes
count_mtx <- count_mtx[!duplicated(rownames(count_mtx)),]

# Sanity checks
stopifnot(sum(duplicated(rownames(count_mtx)))==0)
stopifnot(sum(duplicated(colnames(count_mtx)))==0)

##########################
## Create Seurat object ##
##########################

cell.info.to.seurat <- cell.info[cell%in%colnames(count_mtx)] %>% setkey(cell) %>% .[colnames(count_mtx)] %>% as.data.frame
rownames(cell.info.to.seurat) <- cell.info.to.seurat$cell
stopifnot(rownames(cell.info.to.seurat)==colnames(count_mtx))
stopifnot(sum(is.na(rownames(cell.info.to.seurat$cell)))==0)

seurat <- CreateSeuratObject(count_mtx, meta.data = cell.info.to.seurat)

head(seurat@meta.data)

# Add mitochondrial percenatge
seurat[["mitochondrial_percent_RNA"]] <- PercentageFeatureSet(seurat, pattern = "^mt-") %>% round(2)

# Add ribosomal RNA content
ribo.genes <- grep(pattern = "^Rp[l|s]", x = rownames(seurat), value = TRUE)
seurat[["ribosomal_percent_RNA"]] <- PercentageFeatureSet(seurat, features = ribo.genes) %>% round(2)

#####################
## Create metadata ##
#####################

metadata <- seurat@meta.data %>% as.data.table %>% .[,orig.ident:=NULL] %>%
  .[,c("cell","barcode","sample","nFeature_RNA","nCount_RNA","mitochondrial_percent_RNA","ribosomal_percent_RNA")]

# Add manual metadata columns
stopifnot(unique(metadata$sample)%in%names(opts$sample2alias))
metadata[,alias:=stringr::str_replace_all(sample,opts$sample2alias)]

stopifnot(unique(metadata$sample)%in%names(opts$sample2day))
metadata[,day:=stringr::str_replace_all(sample,opts$sample2day)]

stopifnot(unique(metadata$sample)%in%names(opts$sample2genotype))
metadata[,genotype:=stringr::str_replace_all(sample,opts$sample2genotype)]

stopifnot(unique(metadata$sample)%in%names(opts$sample2rep))
metadata[,replicate:=stringr::str_replace_all(sample,opts$sample2rep)]

##########
## Save ##
##########

fwrite(metadata, file.path(args$outdir,"metadata.txt.gz"), quote=F, na="NA", sep="\t")
saveRDS(seurat, file.path(args$outdir,"seurat.rds"))

