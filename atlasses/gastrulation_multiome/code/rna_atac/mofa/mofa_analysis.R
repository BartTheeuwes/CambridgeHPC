
suppressMessages(library(MOFA2))
suppressPackageStartupMessages(library(argparse))

here::i_am("rna_atac/mofa/mofa_analysis.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--mofa_model',        type="character",                               help='MOFA model')
p$add_argument('--outdir',        type="character",                               help='Output directory')

args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

## START TEST ##
# args$metadata <- file.path(io$basedir,"results_new/rna_atac/mofa/all_cells/PeakMatrix/sample_metadata.txt.gz")
# args$mofa_model <- file.path(io$basedir,"results_new/rna_atac/mofa/all_cells/PeakMatrix/mofa.hdf5")
# args$outdir <- file.path(io$basedir, "results_new/rna_atac/mofa/all_cells/PeakMatrix/pdf")

args$metadata <- file.path(io$basedir,"results_new/rna_atac/mofa/all_cells/sample_metadata.txt.gz")
args$mofa_model <- file.path(io$basedir,"results_new/rna_atac/mofa/all_cells/mofa_allcells_pca_lsi.rds")
args$outdir <- file.path(io$basedir, "results_new/rna_atac/mofa/all_cells/pdf")
## END TEST ##

# I/O
dir.create(args$outdir, showWarnings=F)

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata)

###############
## Load MOFA ##
###############

# MOFAobject <- load_model(args$mofa_model, load_data = T)
MOFAobject <- readRDS(args$mofa_model)

# factors_names(MOFAobject) <- paste0("Factor ", 1:get_dimensions(MOFAobject)[["K"]])
#########################
## Add sample metadata ##
#########################

cells <- as.character(unname(unlist(MOFA2::samples_names(MOFAobject))))

sample_metadata_to_mofa <- copy(sample_metadata) %>%
  setnames("cell","sample") %>%
  .[sample%in%cells] %>% setkey(sample) %>% .[cells]
stopifnot(all(cells==sample_metadata_to_mofa$cell))

samples_metadata(MOFAobject) <- sample_metadata_to_mofa

#################################
## Correlation between factors ##
#################################

plot_factor_cor(MOFAobject)

####################
## Subset factors ##
####################

# r2 <- MOFAobject@cache$variance_explained$r2_per_factor
# factors <- sapply(r2, function(x) x[,"RNA"]>0.01)
# MOFAobject <- subset_factors(MOFAobject, which(apply(factors,1,sum)>=1))
# factors(MOFAobject) <- paste("Factor",1:get_dimensions(MOFAobject)[["K"]], sep=" ")

#############################
## Plot variance explained ##
#############################

p <- plot_variance_explained(MOFAobject, plot_total = T)[[2]]

pdf(sprintf("%s/mofa_var_explained_total.pdf",args$outdir), width=6, height=3)
print(p)
dev.off()

p <- plot_variance_explained(MOFAobject, factors = 1:25, x="view", y="factor", max_r2 = 5) +
  theme(legend.position = "top")

pdf(sprintf("%s/mofa_var_explained.pdf",args$outdir), width=6, height=8)
print(p)
dev.off()

##################
## Plot factors ##
##################

plot_factor(MOFAobject, factors = c(1), color_by = "celltype.predicted", group_by = "celltype.predicted", add_violin = T, add_boxplot = T, add_dots = F, dodge=T, legend=F) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )


plot_factor(MOFAobject, factors = 3, color_by = "celltype.predicted", group_by = "celltype.predicted", 
            add_boxplot = T, add_dots = F, add_violin = T, dodge=T, legend=F) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# plot_factors(MOFAobject, factors = c(1,3), color_by = "celltype.predicted")

p <- plot_factors(MOFAobject, factors = c(1,2), color_by = "celltype.predicted", dot_size = 1) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme(
    legend.position = "none"
  )

pdf(sprintf("%s/mofa_scatterplot_factor12.pdf",args$outdir), width=6, height=5)
print(p)
dev.off()

######################
## Sumarise factors ##
######################

levels_df <- MOFAobject@samples_metadata[,c("sample","celltype.predicted")] %>% setnames("celltype.predicted","level")
p <- summarise_factors(MOFAobject, levels_df, factors = 1:25, abs = F, return_data = F) +
  guides(x = guide_axis(angle = 90)) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(color="black", size=rel(0.75)),
    axis.text.y = element_text(color="black", size=rel(0.75)),
  )

pdf(sprintf("%s/summarise_factor_celltype.pdf",args$outdir), width=7, height=5)
print(p)
dev.off()

##################
## Plot weights ##
##################

plot_weights(MOFAobject, factor = 1, view="ATAC", nfeatures = 10, text_size = 3)
plot_weights(MOFAobject, factor = 1, view="RNA", nfeatures = 15, text_size = 4)

#######################################
## Correlate factors with covariates ##
#######################################

foo <- correlate_factors_with_covariates(MOFAobject, covariates = c("nFeature_RNA","nFrags_atac","ribosomal_percent_RNA","mitochondrial_percent_RNA"), return_data = T)

# pheatmap::pheatmap(foo)

######################
## Batch correction ##
######################

# library(batchelor)

# Select factors to use 
factors.to.use <- 1:get_dimensions(MOFAobject)[["K"]]
# factors.to.use <- factors.to.use[!factors.to.use%in%c("3")]
# factors.to.use <- c(1,2,3,6,7)

# Extract factors
Z <- get_factors(MOFAobject, factors=factors.to.use)[[1]]

# MNN correction
Z <- reducedMNN(Z, batch=MOFAobject@samples_metadata$stage)$corrected
colnames(Z) <- colnames(Z)

# UMAP
umap_embedding <- uwot::umap(Z, n_neighbors=25, min_dist=0.50, metric="cosine")
# umap_embedding <- uwot::umap(Z.corrected, n_neighbors=25, min_dist=0.40, metric="cosine")

# Plot
to.plot <- umap_embedding %>% as.data.table %>%
  .[,sample:=rownames(Z)] %>%
  merge(MOFAobject@samples_metadata[,c("sample","stage","celltype.predicted")] %>% as.data.table)

p <- ggplot(to.plot, aes(x=V1, y=V2, fill=celltype.predicted)) +
  geom_point(size=1, shape=21, stroke=0.1) +
  scale_fill_manual(values=opts$celltype.colors) +
  guides(fill = guide_legend(override.aes = list(size=2))) +
  theme_classic() +
  ggplot_theme_NoAxes() +
  theme(
    legend.position="none"
  )

pdf(sprintf("%s/mofa_umap_celltype.pdf",args$outdir), width=7, height=5)
print(p)
dev.off()

# to.plot %>% .[,stage:=factor(stage,levels=opts$stages)]
p <- ggplot(to.plot, aes(x=V1, y=V2, fill=stage)) +
  geom_point(size=1, shape=21, stroke=0.1) +
  scale_fill_manual(values=opts$stage.colors) +
  guides(fill = guide_legend(override.aes = list(size=2))) +
  theme_classic() +
  ggplot_theme_NoAxes() +
  theme(
    legend.position="none"
  )

pdf(sprintf("%s/mofa_umap_stage.pdf",args$outdir), width=7, height=5)
print(p)
dev.off()

# plot_dimred(MOFAobject, method="UMAP", color_by="stage")

#########################
## Contribution scores ##
#########################

factors.to.use <- "all"

r2.per.sample <- calculate_variance_explained_per_sample(MOFAobject)[[1]]
# foo <- rowSums(r2.per.sample)
# r2.per.sample["E8.0_rep1#GCGAAGTAGTACCGCA-1",]
# foo["E8.0_rep1#GCGAAGTAGTACCGCA-1"]
# sort(foo) %>% head
MOFAobject <- calculate_contribution_scores(MOFAobject, factors = factors.to.use, scale = TRUE)

cells <- intersect(names(which(r2.per.sample[,"RNA"]>=3)), names(which(r2.per.sample[,"ATAC"]>=3)))

to.plot <- MOFAobject@samples_metadata[MOFAobject@samples_metadata$sample%in%cells,c("sample","RNA_contribution","ATAC_contribution","celltype.predicted")] %>% as.data.table

order.celltypes <- to.plot[,median(RNA_contribution),by="celltype.predicted"] %>% setorder(-V1) %>% .$celltype.predicted
to.plot[,celltype.predicted:=factor(celltype.predicted, levels=order.celltypes)]


p <- ggplot(to.plot, aes(x=celltype.predicted, y=RNA_contribution)) +
  geom_boxplot(aes(fill = celltype.predicted), alpha=0.9, outlier.shape=NA, coef=1.5) +
  coord_flip(ylim = c(0.10,0.90)) +
  geom_hline(yintercept=0.5, linetype="dashed", size=0.5) +
  scale_fill_manual(values=opts$celltype.colors, drop=F) +
  theme_classic() +
  labs(y="RNA contribution score", x="") +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  )


pdf(sprintf("%s/contribution_scores.pdf",args$outdir), width=5, height=7)
print(p)
dev.off()


##################################################################
## Plot cumulative variance explained per view vs factor number ##
##################################################################

factors <- 1:25

r2.dt <- MOFAobject@cache$variance_explained$r2_per_factor[[1]][factors,] %>%
  as.data.table %>% .[,factor:=as.factor(factors)] %>%
  melt(id.vars="factor", variable.name="view", value.name = "r2") %>%
  .[,cum_r2:=cumsum(r2), by="view"]

# threshold.var <- 5
# max.factor <- max(which(apply(r2,1,sum) >= threshold.var))

p <- ggline(r2.dt, x="factor", y="cum_r2", color="view") +
  # scale_color_manual(values=opts$colors.views) +
  labs(x="Factor number", y="Cumulative variance explained (%)") +
  # geom_vline(xintercept = max.factor, linetype="dashed") +
  theme(
    legend.title = element_blank(), 
    legend.position = "top",
    axis.text = element_text(size=rel(0.8))
  )

pdf(paste0(args$outdir,"/r2_vs_factor.pdf"), width=8, height=5, useDingbats = F)
print(p)
dev.off()


