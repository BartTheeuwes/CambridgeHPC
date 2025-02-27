suppressMessages(library(MOFA2))
suppressMessages(library(ggpubr))

#####################
## Define settings ##
#####################

source("/Users/ricard/gastrulation_multiome_10x/settings.R")

################
## Load MOFAobject ##
################

# file <- paste0(io$basedir,"/MOFAobject/hdf5/MOFAobject.hdf5")
# MOFAobject <- load_model(file)

io$MOFAobject <- paste0(io$basedir, "/results/rna_atac/mofa/blood/mofa_model.rds")
MOFAobject <- readRDS(io$MOFAobject)

#########################
## Add sample metadata ##
#########################

# cells <- as.character(unname(unlist(MOFA2::samples_names(MOFAobject))))
# 
# sample_metadata.MOFAobject <- copy(sample_metadata) %>%
#   setnames("cell","sample") %>%
#   setnames("batch","group") %>%
#   .[sample%in%cells] %>% setkey(sample) %>% .[cells]
# stopifnot(all(cells==sample_metadata.MOFAobject$cell))
# 
# samples_metadata(MOFAobject) <- sample_metadata.MOFAobject

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

plot_variance_explained(MOFAobject, x="view", y="factor", max_r2 = 10)

##################
## Plot factors ##
##################

plot_factor(MOFAobject, factors = 1:3, color_by = "celltype.mapped", group_by = "celltype.mapped", add_violin = T, add_boxplot = T, add_dots = F, dodge=T, legend=F) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )


plot_factor(MOFAobject, factors = 1, color_by = "Epha5", group_by = "celltype.mapped", dodge=T) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

plot_factors(MOFAobject, factors = c(1,3), color_by = "stage")

plot_factors(MOFAobject, factors = c(1,4), color_by = "celltype.mapped") +
  scale_fill_manual(values=opts$celltype.colors)

######################
## Sumarise factors ##
######################

levels_df <- MOFAobject@samples_metadata[,c("sample","celltype.mapped")] %>% setnames("celltype.mapped","level")
summarise_factors(MOFAobject, levels_df, factors = "all", abs = F, return_data = F)

##################
## Plot weights ##
##################

plot_weights(MOFAobject, factor = 3, view="ATAC_chromVAR", nfeatures = 10, text_size = 3)
plot_weights(MOFAobject, factor = 3, view="RNA", nfeatures = 15, text_size = 4)

#######################################
## Correlate factors with covariates ##
#######################################

correlate_factors_with_covariates(MOFAobject, covariates = c("nFeature_RNA","nFrags_atac","mtFraction_RNA"))

######################
## Batch correction ##
######################

library(batchelor)

# Select factors to use 
factors <- 1:get_dimensions(MOFAobject)[["K"]]
factors.to.use <- factors[!factors%in%c("3")]
# factors.to.use <- c(1,2,3,6,7)

# Extract factors
Z <- get_factors(MOFAobject, factors=factors.to.use)[[1]]

# MNN correction
MOFAobject.corrected <- reducedMNN(Z, batch=MOFAobject@samples_metadata$stage)$corrected
colnames(MOFAobject.corrected) <- colnames(Z)

# MOFAobject.corrected <- Z

# UMAP
umap_embedding <- uwot::umap(MOFAobject.corrected, n_neighbors=15, min_dist=0.15, metric="cosine")

# Plot
to.plot <- umap_embedding %>% as.data.table %>%
  .[,sample:=rownames(Z)] %>%
  merge(MOFAobject@samples_metadata[,c("sample","stage","celltype.mapped")] %>% as.data.table)

ggplot(to.plot, aes(x=V1, y=V2, fill=celltype.mapped)) +
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
factors <- 1:get_dimensions(MOFAobject)[["K"]]
factors.to.use <- factors[!factors%in%c("3")]
# factors.to.use <- c(1,2,3,6,7)

# Run UMAP
MOFAobject <- run_umap(MOFAobject, factors=factors.to.use, min_dist = 0.15)

plot_dimred(MOFAobject, method="UMAP", color_by = "stage")
plot_dimred(MOFAobject, method="UMAP", color_by = "celltype.mapped", legend=F) +
  scale_fill_manual(values=opts$celltype.colors)
plot_dimred(MOFAobject, method="UMAP", color_by = "Factor1")

#########################
## Contribution scores ##
#########################

factors.to.use <- "all" # c(1,2,3)
MOFAobject <- calculate_contribution_scores(MOFAobject, factors = factors.to.use, scale = TRUE)
head(MOFAobject@cache$contribution_scores)
head(MOFAobject@samples_metadata)

to.plot <- MOFAobject@samples_metadata[,c("sample","celltype.mapped","RNA_contribution")] %>% 
  as.data.table(keep.rownames = T)

# order.celltypes <- to.plot[,median(RNA_fraction),by="celltype"] %>% setorder(-V1) %>% .$celltype
# to.plot[,celltype:=factor(celltype, levels=order.celltypes)]

ggboxplot(to.plot, x="celltype.mapped", y="RNA_contribution", fill="celltype.mapped", outlier.shape=NA) +
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

plot_factor(MOFAobject, factors = 5, color_by = "LYZ")

plot_factor(MOFAobject, factors = 5, color_by = "IGLL1", group_by="BioClassification", dodge=TRUE) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

plot_data_heatmap(MOFAobject,
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

plot_data_scatter(MOFAobject, factor=1, view=1, color_by = "lab", features = 8, dot_size = 2)



###################
## MEFISTO plots ##
###################

plot_factors_vs_cov(MOFAobject, factors=1:5)

smoothness <- get_scales(MOFAobject)

plot_smoothness(MOFAobject)
