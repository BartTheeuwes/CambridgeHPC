suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(argparse))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--sce',             type="character",                               help='SingleCellExperiment file')
p$add_argument('--TFs_file',             type="character",                               help='txt file with a list of TFs')
p$add_argument('--outdir',          type="character",                               help='Output file')

args <- p$parse_args(commandArgs(TRUE))


#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

## START TEST ##
# args$sce <- file.path(io$basedir,"results/rna/pseudobulk/SingleCellExperiment_pseudobulk_celltype.mapped_mnn.rds")
# args$TFs_file <- file.path(io$archR.directory,"Annotations/Motif_cisbp_TFs.txt.gz")
# args$outdir <- file.path(io$basedir,"results/rna/coexpression")
## END TEST ##

###################
## Load metadata ##
###################

# sample_metadata <- fread(io$metadata) %>%
#   .[pass_rnaQC==TRUE & doublet_call==FALSE] %>%
#   .[celltype.predicted%in%opts$celltypes] 

# if (opts$remove.ExE.celltypes) {
#   sample_metadata <- sample_metadata %>%
#     .[!celltype.mapped%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
#   opts$celltypes <- opts$celltypes[!opts$celltypes%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
# }

#########################
## Load pseudobulk RNA ##
#########################

# Load SingleCellExperiment
rna.sce <- readRDS(args$sce)#[,opts$celltypes]

##################################################
## Split RNA expression matrix into TF vs genes ##
##################################################

TFs <- fread(args$TFs_file)[["gene"]] %>% unique

print("The following TFs are not found in the SingleCellExperiment:")
print(TFs[!TFs%in%toupper(rownames(rna.sce))])

rna_tfs.sce <- rna.sce[toupper(rownames(rna.sce))%in%TFs,]
rna_genes.sce <- rna.sce
rownames(rna_tfs.sce) <- toupper(rownames(rna_tfs.sce))

print(sprintf("Number of TFs: %s",nrow(rna_tfs.sce)))
print(sprintf("Number of genes: %s",nrow(rna_genes.sce)))

##########################
## Correlation analysis ##
##########################

# tf2gene
tf2gene_cor.mtx <- cor(t(logcounts(rna_tfs.sce)),t(logcounts(rna_genes.sce))) %>% round(2)
saveRDS(tf2gene_cor.mtx, file.path(args$outdir,"correlation_matrix_tf2gene_pseudobulk.rds"))

# tf2tf
tf2tf_cor.mtx <- cor(t(logcounts(rna_tfs.sce)),t(logcounts(rna_tfs.sce))) %>% round(2)
saveRDS(tf2tf_cor.mtx, file.path(args$outdir,"correlation_matrix_tf2tf_pseudobulk.rds"))

##########
## Plot ##
##########

# tf2gene_cor.mtx <- readRDS(file.path(args$outdir,"correlation_matrix_tf2gene.rds"))

# i <- "FOXA2"
# j <- "Cab39l"

# to.plot <- data.table(
#   TF = logcounts(rna_tfs.sce[i,])[1,],
#   target_gene = logcounts(rna_genes.sce[j,])[1,],
#   celltype = colnames(rna_tfs.sce)
# )


# ggscatter(to.plot, x="TF", y="target_gene", fill="celltype", size=4, shape=21, 
#           add="reg.line", add.params = list(color="black", fill="lightgray"), conf.int=TRUE) +
#   stat_cor(method = "pearson") +
#   scale_fill_manual(values=opts$celltype.colors) +
#   labs(x=sprintf("%s expression",i), y=sprintf("%s expression",j)) +
#   guides(fill=F) +
#   theme(
#     plot.title = element_text(hjust = 0.5, size=rel(0.85)),
#     axis.text = element_text(size=rel(0.7))
#   )
