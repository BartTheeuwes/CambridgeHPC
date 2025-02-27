#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

# I/O
io$outdir <- file.path(io$basedir,"results/rna/coexpression")

# Options
opts$remove.ExE.celltypes <- FALSE

###################
## Load metadata ##
###################

sample_metadata <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE] %>%
  .[celltype.predicted%in%opts$celltypes] 

if (opts$remove.ExE.celltypes) {
  sample_metadata <- sample_metadata %>%
    .[!celltype.mapped%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
  opts$celltypes <- opts$celltypes[!opts$celltypes%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}

#########################
## Load pseudobulk RNA ##
#########################

# Load SingleCellExperiment
rna.sce <- readRDS(io$rna.pseudobulk.sce)[,opts$celltypes]

##################################################
## Split RNA expression matrix into TF vs genes ##
##################################################

TFs <- fread(io$TFs)[[1]]

rna.sce.tf <- rna.sce[toupper(rownames(rna.sce))%in%TFs,]
rna.sce.target <- rna.sce
rownames(rna.sce.tf) <- toupper(rownames(rna.sce.tf))

##########################
## Correlation analysis ##
##########################

# tf2gene
tf2gene_cor.mtx <- cor(t(logcounts(rna.sce.tf)),t(logcounts(rna.sce.target))) %>% round(2)
saveRDS(tf2gene_cor.mtx, file.path(io$outdir,"correlation_matrix_tf2gene.rds"))

# tf2tf
tf2tf_cor.mtx <- cor(t(logcounts(rna.sce.tf)),t(logcounts(rna.sce.tf))) %>% round(2)
saveRDS(tf2tf_cor.mtx, file.path(io$outdir,"correlation_matrix_tf2tf.rds"))

##########
## Plot ##
##########

tf2gene_cor.mtx <- readRDS(file.path(io$outdir,"correlation_matrix_tf2gene.rds"))

i <- "FOXA2"
j <- "Cab39l"

to.plot <- data.table(
  TF = logcounts(rna.sce.tf[i,])[1,],
  target_gene = logcounts(rna.sce.target[j,])[1,],
  celltype = colnames(rna.sce.tf)
)


ggscatter(to.plot, x="TF", y="target_gene", fill="celltype", size=4, shape=21, 
          add="reg.line", add.params = list(color="black", fill="lightgray"), conf.int=TRUE) +
  stat_cor(method = "pearson") +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(x=sprintf("%s expression",i), y=sprintf("%s expression",j)) +
  guides(fill=F) +
  theme(
    plot.title = element_text(hjust = 0.5, size=rel(0.85)),
    axis.text = element_text(size=rel(0.7))
  )
