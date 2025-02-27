suppressMessages(library(MOFA2))

#####################
## Define settings ##
#####################

source("/Users/ricard/gastrulation_multiome_10x/settings.R")

################
## Load mofa ##
################

# file <- paste0(io$basedir,"/mofa/hdf5/mofa.hdf5")
# mofa <- load_model(file)

object <- load_model("/tmp/mofa_20210205-144602.hdf5", remove_outliers = T)
io$mofa <- paste0(io$basedir, "/results/rna/mofa/mofa_model_blood_rna.rds")
mofa <- readRDS(io$mofa)

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

plot_factor(mofa, factors = 1, color_by = "celltype.predicted", group_by = "celltype.predicted", add_violin = T, add_boxplot = T, add_dots = F, dodge=T, legend=F) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )


plot_factor(mofa, factors = c(1,5,10,14), color_by = "batch", group_by = "batch", add_violin = T, add_boxplot = T, add_dots = F, dodge=T) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

plot_factors(mofa, factors = c(1,2), color_by = "batch")

plot_factors(mofa, factors = c(1,2), color_by = "celltype.predicted", legend=F) +
  scale_fill_manual(values=opts$celltype.colors)


######################
## Sumarise factors ##
######################

levels_df <- mofa@samples_metadata[,c("sample","celltype.predicted")] %>% setnames("celltype.predicted","level")
summarise_factors(mofa, levels_df, factors = "all", abs = F, return_data = F)

##################
## Plot weights ##
##################

plot_weights(mofa, factor = 2, view="RNA", nfeatures = 15, text_size = 4)

#######################################
## Correlate factors with covariates ##
#######################################

correlate_factors_with_covariates(mofa, covariates = c("nFeature_RNA","mtFraction_RNA"))

######################
## Batch correction ##
######################

library(batchelor)

# Select factors to use 
factors <- 1:get_dimensions(mofa)[["K"]]
factors.to.use <- factors
# factors.to.use <- factors[!factors%in%c("3")]

# Extract factors
Z <- get_factors(mofa, factors=factors.to.use)[[1]]

# Scale factors by variance explained
scaling.factor <- get_variance_explained(mofa)[["r2_per_factor"]][[1]] / sum(get_variance_explained(mofa)[["r2_per_factor"]][[1]])
Z <- sweep(Z, MARGIN = 2, STATS = scaling.factor, FUN = "*")
plot(apply(Z,2,var))

# MNN correction
mofa.corrected <- reducedMNN(Z, batch=mofa@samples_metadata$batch)$corrected
colnames(mofa.corrected) <- colnames(Z)

# mofa.corrected <- Z

# UMAP
umap_embedding <- uwot::umap(mofa.corrected, n_neighbors=25, min_dist=0.30, metric="cosine")

# Plot
to.plot <- umap_embedding %>% as.data.table %>%
  .[,sample:=rownames(Z)] %>%
  merge(mofa@samples_metadata[,c("sample","stage","celltype.predicted")] %>% as.data.table)

ggplot(to.plot, aes(x=V1, y=V2, fill=stage)) +
  geom_point(size=1.25, shape=21, stroke=0.1) +
  guides(fill = guide_legend(override.aes = list(size=2))) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position="right",
    legend.title=element_blank()
  )

ggplot(to.plot, aes(x=V1, y=V2, fill=celltype.predicted)) +
  geom_point(size=1.25, shape=21, stroke=0.1) +
  scale_fill_manual(values=opts$celltype.colors) +
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
# factors.to.use <- factors[!factors%in%c("3")]
factors.to.use <- factors

# Run UMAP
mofa <- run_umap(mofa, factors=factors.to.use, n_neighbors = 25, min_dist = 0.3)

plot_dimred(mofa, method="UMAP", color_by = "batch")
plot_dimred(mofa, method="UMAP", color_by = "celltype.predicted") +
  scale_fill_manual(values=opts$celltype.colors)

#########################
## Contribution scores ##
#########################

mofa <- calculate_contribution_scores(mofa, scale = TRUE)

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
                  factor = 12,
                  features = 25,
                  denoise = FALSE,
                  legend = TRUE,
                  cluster_rows = T, cluster_cols = F,
                  show_colnames = F, show_rownames = T,
                  scale="row"
)

plot_data_scatter(mofa, factor=1, view=1, color_by = "lab", features = 8, dot_size = 2)



####################
## Compare to PCs ##
####################

Z <- get_factors(mofa)[[1]]
W <- get_weights(mofa)[[1]]

apply(Z,2,mean)
sce_filt <- runPCA(sce_filt, ncomponents = 15, ntop=nrow(sce_filt))
cor(colMeans(assay(sce,"logcounts")), reducedDim(sce_filt, "PCA"))
r <- cor(Z, reducedDim(sce_filt, "PCA")) %>% abs
corrplot::corrplot(r)


plot(apply(reducedDim(sce_filt, "PCA"),2,var))
plot(apply(get_factors(mofa)[[1]],2,var))
plot(apply(get_weights(mofa)[[1]],2,var))
plot(apply(sweep(get_factors(mofa)[[1]], MARGIN = 2, STATS = get_variance_explained(mofa)[["r2_per_factor"]][[1]], FUN = "*"),2,var))

apply(reducedDim(sce_filt, "PCA"),2,max)

apply(W,2,max)
apply(W,2,var)

