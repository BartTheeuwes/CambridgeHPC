suppressMessages(library(MOFA2))
suppressMessages(library(ggpubr))

#####################
## Define settings ##
#####################

source("/Users/ricard/gastrulation_multiome_10x/settings.R")

################
## Load mofa ##
################

# file <- paste0(io$basedir,"/mofa/hdf5/mofa.hdf5")
# mofa <- load_model(file)

io$mofa <- paste0(io$basedir, "/results/rna_atac/mofa/mofa_model_blood_rna.rds")
mofa <- readRDS(io$mofa)

#########################
## Add sample metadata ##
#########################

# cells <- as.character(unname(unlist(MOFA2::samples_names(mofa))))
# 
# sample_metadata.mofa <- copy(sample_metadata) %>%
#   setnames("cell","sample") %>%
#   setnames("batch","group") %>%
#   .[sample%in%cells] %>% setkey(sample) %>% .[cells]
# stopifnot(all(cells==sample_metadata.mofa$cell))
# 
# samples_metadata(mofa) <- sample_metadata.mofa

#################################
## Correlation between factors ##
#################################

plot_factor_cor(mofa)

####################
## Subset factors ##
####################

# r2 <- mofa@cache$variance_explained$r2_per_factor
# factors <- sapply(r2, function(x) x[,"RNA"]>0.01)
# mofa <- subset_factors(mofa, which(apply(factors,1,sum)>=1))
# factors(mofa) <- paste("Factor",1:get_dimensions(mofa)[["K"]], sep=" ")

#############################
## Plot variance explained ##
#############################

plot_variance_explained(mofa, x="view", y="factor", max_r2 = 10)

##################
## Plot factors ##
##################

plot_factor(mofa, factors = c(3), color_by = "celltype.predicted", group_by = "celltype.predicted", add_violin = T, add_boxplot = T, add_dots = F, dodge=T, legend=F) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )


plot_factor(mofa, factors = 7, color_by = "Gpc6", group_by = "celltype.predicted", dodge=T) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# plot_factors(mofa, factors = c(1,7), color_by = "Gpc6")

plot_factors(mofa, factors = c(3,8), color_by = "celltype.predicted") +
  scale_fill_manual(values=opts$celltype.colors)

######################
## Sumarise factors ##
######################

levels_df <- mofa@samples_metadata[,c("sample","celltype.predicted")] %>% setnames("celltype.predicted","level")
summarise_factors(mofa, levels_df, factors = "all", abs = F, return_data = F)

##################
## Plot weights ##
##################

plot_weights(mofa, factor = 8, view="ATAC_chromVAR", nfeatures = 10, text_size = 3)
plot_weights(mofa, factor = 5, view="RNA", nfeatures = 15, text_size = 4)

#######################################
## Correlate factors with covariates ##
#######################################

correlate_factors_with_covariates(mofa, covariates = c("nFeature_RNA","nFrags_atac","mtFraction_RNA"))

######################
## Batch correction ##
######################

library(batchelor)

# Select factors to use 
factors <- 1:get_dimensions(mofa)[["K"]]
factors.to.use <- factors[!factors%in%c("4")]
# factors.to.use <- c(1,2,3,6,7)

# Extract factors
Z <- get_factors(mofa, factors=factors.to.use)[[1]]

# MNN correction
mofa.corrected <- reducedMNN(Z, batch=mofa@samples_metadata$stage)$corrected
colnames(mofa.corrected) <- colnames(Z)

# mofa.corrected <- Z

# UMAP
umap_embedding <- uwot::umap(mofa.corrected, n_neighbors=15, min_dist=0.15, metric="cosine")

# Plot
to.plot <- umap_embedding %>% as.data.table %>%
  .[,sample:=rownames(Z)] %>%
  merge(mofa@samples_metadata[,c("sample","stage","celltype.predicted")] %>% as.data.table)

ggplot(to.plot, aes(x=V1, y=V2, fill=celltype.predicted)) +
  geom_point(size=1.75, shape=21, stroke=0.1) +
  scale_fill_manual(values=opts$celltype.colors) +
  guides(fill = guide_legend(override.aes = list(size=2))) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position="none",
    legend.title=element_blank()
  )

ggplot(to.plot, aes(x=V1, y=V2, fill=stage)) +
  geom_point(size=1.75, shape=21, stroke=0.1) +
  guides(fill = guide_legend(override.aes = list(size=2))) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position="right",
    legend.title=element_blank()
  )


##########
## UMAP ##
##########

# Select factors
factors <- 1:get_dimensions(mofa)[["K"]]
factors.to.use <- factors[!factors%in%c("4")]
# factors.to.use <- c(1,2,3,6,7)

# Run UMAP
mofa <- run_umap(mofa, factors=factors.to.use, min_dist = 0.3)

plot_dimred(mofa, method="UMAP", color_by = "stage")
plot_dimred(mofa, method="UMAP", color_by = "celltype.predicted", legend=F) +
  scale_fill_manual(values=opts$celltype.colors)
plot_dimred(mofa, method="UMAP", color_by = "Factor1")

#########################
## Contribution scores ##
#########################

mofa <- calculate_contribution_scores(mofa, factors = factors.to.use, scale = TRUE)
head(mofa@cache$contribution_scores)
head(mofa@samples_metadata)

to.plot <- mofa@samples_metadata[,c("sample","celltype.predicted","RNA_contribution")] %>% 
  as.data.table(keep.rownames = T)

# order.celltypes <- to.plot[,median(RNA_fraction),by="celltype"] %>% setorder(-V1) %>% .$celltype
# to.plot[,celltype:=factor(celltype, levels=order.celltypes)]

ggboxplot(to.plot, x="celltype.predicted", y="RNA_contribution", fill="celltype.predicted", outlier.shape=NA) +
  labs(x="", y="RNA contribution") +
  geom_hline(yintercept=0.5, linetype="dashed") +
  scale_fill_manual(values=opts$celltype.colors) +
  theme(
    axis.text.x = element_text(size=rel(0.7), color="black", angle=40, vjust=1, hjust=1),
    axis.ticks.x = element_blank(),
    legend.position = "none"
    # legend.position = "right", legend.text = element_text(size=rel(0.75))
  )

###############
## Plot data ##
###############

plot_factor(mofa, factors = 5, color_by = "LYZ")

plot_factor(mofa, factors = 5, color_by = "IGLL1", group_by="BioClassification", dodge=TRUE) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

plot_data_heatmap(mofa,
                  factor = 4,
                  features = 25,
                  view = "ATAC_chromVAR",
                  denoise = FALSE,
                  legend = TRUE,
                  # min.value = 0, max.value = 6,
                  cluster_rows = T, cluster_cols = F,
                  show_colnames = F, show_rownames = T,
                  scale="row"
                  # annotation_samples = "Category",  annotation_colors = list("Category"=opts$colors), annotation_legend = F
)

plot_data_scatter(mofa, factor=1, view=1, color_by = "lab", features = 8, dot_size = 2)


