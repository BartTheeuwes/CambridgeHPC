suppressMessages(library(MOFA2))

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/utils.R")
} else {
  stop("Computer not recognised")
}

io$MOFAobject <- paste0(io$basedir, "/results/rna_atac/mofa/TFs/mofa_TFs.rds")
io$outdir <- paste0(io$basedir, "/results/rna_atac/mofa/TFs/pdf"); dir.create(io$outdir, showWarnings = F)

###############
## Load MOFA ##
###############

MOFAobject <- readRDS(io$MOFAobject)

MOFAobject@samples_metadata$celltype.predicted <- factor(MOFAobject@samples_metadata$celltype.predicted, levels=opts$celltypes)

factors_names(MOFAobject) <- paste0("Factor ", 1:get_dimensions(MOFAobject)[["K"]])
  
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

p <- plot_variance_explained(MOFAobject, plot_total = T)[[2]]

pdf(sprintf("%s/mofa_var_explained_total.pdf",io$outdir), width=6, height=3)
print(p)
dev.off()

p <- plot_variance_explained(MOFAobject, factors = 1:25, x="view", y="factor", max_r2 = 7) +
  theme(legend.position = "top")

pdf(sprintf("%s/mofa_var_explained.pdf",io$outdir), width=6, height=8)
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

plot_factors(MOFAobject, factors = c(1,2), color_by = "OTX2_TF_RNA")

p <- plot_factors(MOFAobject, factors = c(1,3), color_by = "celltype.predicted") +
  scale_fill_manual(values=opts$celltype.colors) +
  theme(
    legend.position = "none"
  )

pdf(sprintf("%s/mofa_scatterplot_factor12.pdf",io$outdir), width=6, height=5)
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

pdf(sprintf("%s/summarise_factor_celltype.pdf",io$outdir), width=7, height=5)
print(p)
dev.off()

##################
## Plot weights ##
##################

plot_weights(MOFAobject, factor = 8, view="ATAC_chromVAR", nfeatures = 10, text_size = 3)
plot_weights(MOFAobject, factor = 1, view="TF_RNA", nfeatures = 15, text_size = 4)

#######################################
## Correlate factors with covariates ##
#######################################

foo <- correlate_factors_with_covariates(MOFAobject, covariates = c("nFeature_RNA","nFrags_atac","ribosomal_percent_RNA","mitochondrial_percent_RNA"), return_data = T)

##########
## UMAP ##
##########

factors.to.use <- 1:get_dimensions(MOFAobject)[["K"]]

set.seed(42)
MOFAobject <- run_umap(MOFAobject, factors=factors.to.use, n_neighbors=25, min_dist=0.40, metric="cosine")

p <- plot_dimred(MOFAobject, method="UMAP", color_by="celltype.predicted", legend=F, stroke=0.1) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme(
    axis.line = element_blank()
  )

pdf(sprintf("%s/mofa_umap_celltype.pdf",io$outdir), width=7, height=5)
print(p)
dev.off()

p <- plot_dimred(MOFAobject, method="UMAP", color_by="stage", legend=F, stroke=0.1) +
  scale_fill_manual(values=opts$stage.colors) +
  theme(
    axis.line = element_blank()
  )

pdf(sprintf("%s/mofa_umap_stage.pdf",io$outdir), width=7, height=5)
print(p)
dev.off()

#########################
## Contribution scores ##
#########################

factors.to.use <- "all"
MOFAobject <- calculate_contribution_scores(MOFAobject, factors = factors.to.use, scale = TRUE)

to.plot <- MOFAobject@samples_metadata[,c("sample","celltype.predicted","TF_RNA_contribution")] %>% as.data.table %>%
  setnames("sample","cell")

order.celltypes <- to.plot[,median(TF_RNA_contribution),by="celltype.predicted"] %>% setorder(-V1) %>% .$celltype.predicted
to.plot[,celltype.predicted:=factor(celltype.predicted, levels=order.celltypes)]

p <- ggplot(to.plot, aes(x=celltype.predicted, y=TF_RNA_contribution)) +
  geom_boxplot(aes(fill = celltype.predicted), alpha=0.9, outlier.shape=NA, coef=1) +
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


pdf(sprintf("%s/contribution_scores.pdf",io$outdir), width=5, height=7)
print(p)
dev.off()

###############
## Plot data ##
###############

plot_dimred(MOFAobject, method="UMAP", color_by="OTX2_TF_RNA", legend=F, stroke=0.1)
plot_dimred(MOFAobject, method="UMAP", color_by="MEIS1_TF_RNA", legend=F, stroke=0.1)

plot_dimred(MOFAobject, method="UMAP", color_by="TAL1_TF_chromVAR", legend=F, stroke=0.1)

# plot_data_heatmap(MOFAobject,
#                   factor = 4,
#                   features = 25,
#                   view = "ATAC_chromVAR",
#                   denoise = FALSE,
#                   legend = TRUE,
#                   # min.value = 0, max.value = 6,
#                   cluster_rows = T, cluster_cols = F,
#                   show_colnames = F, show_rownames = T,
#                   scale="row"
#                   # annotation_samples = "Category",  annotation_colors = list("Category"=opts$colors), annotation_legend = F
# )

# plot_data_scatter(MOFAobject, factor=1, view=1, color_by = "lab", features = 8, dot_size = 2)



