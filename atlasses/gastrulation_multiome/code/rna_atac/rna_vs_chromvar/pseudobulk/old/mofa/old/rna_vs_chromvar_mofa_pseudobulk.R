library(MOFA2)

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
} else {
  stop("Computer not recognised")
}


# I/O
io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/mofa")
# io$mofa.output <- paste0(io$basedir, "/results/rna_atac/mofa/NMPs/mofa_model.rds")

# Options

##############################
## Load RNA expression data ##
##############################

# TO-DO: REPLACE BY MULTIOME DATA
expr.dt <- fread(io$rna.atlas.average_expression_per_celltype) %>% 
  .[,c("group","gene","mean_expr")] %>%
  setnames(c("celltype","gene","expr")) %>% 
  .[,celltype:=gsub(" ","_",celltype)] %>% .[,celltype:=gsub("/","_",celltype)]

#######################################
## Load chromatin accessibility data ##
#######################################

io$archR.chromvar.deviations.pseudobulk <- "/Users/ricard/data/gastrulation_multiome_10x/results/atac/archR/chromvar/cisbp/archr_motif_chromvar_scores_celltype.txt.gz"

acc.dt <- fread(io$archR.chromvar.deviations.pseudobulk) %>%
  setnames(c("celltype","gene","chromvar")) %>%
  .[chromvar<0,chromvar:=0]

###########
## Merge ##
###########

expr_acc.dt <- merge(acc.dt, expr.dt,by = c("celltype","gene"))

############
## Filter ##
############

to.plot <- expr_acc.dt[,.(mean_expr=mean(expr)),by=c("gene")]

expr_acc.dt <- expr_acc.dt[,mean_expr:=mean(expr),by=c("gene")] %>% .[mean_expr>0.001] %>% .[,mean_expr:=NULL]

data_to_mofa.dt <- expr_acc.dt %>% 
  melt(id.vars=c("celltype","gene")) %>%
  setnames(c("sample","feature","view","value")) %>%
  .[view=="chromvar",feature:=paste0(feature,"_chromvar")]

##############
## Run MOFA ##
##############

MOFAobject <- create_mofa_from_df(data_to_mofa.dt)

# Model options
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 10

# Prepare MOFA object
MOFAobject <- prepare_mofa(
  MOFAobject,
  # data_options = data_opts,
  model_options = model_opts
  # training_options = train_opts
)

# Train the model
MOFAobject <- run_mofa(MOFAobject)

# Save
saveRDS(MOFAobject, io$mofa.output)

##############
## Analysis ##
##############

MOFAobject <- readRDS(io$mofa)

plot_variance_explained(MOFAobject, x="view", y="factor")


plot_factor(MOFAobject, factors = 10, color_by = "sample", dodge=F, dot_size=7, legend = FALSE) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

plot_factor(MOFAobject, factors = 6, color_by = "Mecom_chromvar", dot_size=7)


plot_factors(MOFAobject, factors = c(2,4), dot_size=7, color_by = "sample", legend = FALSE) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

plot_factors(MOFAobject, factors = c(1,3), dot_size=7, color_by = "Mecom_chromvar")

##########
## UMAP ##
##########

MOFAobject <- run_umap(MOFAobject, factors = "all")
plot_dimred(MOFAobject, method="UMAP", color_by = "sample", stroke=0.1, dot_size=10, legend=F) +
  scale_fill_manual(values=opts$celltype.colors)


##########
## t-SNE ##
##########

MOFAobject <- run_tsne(MOFAobject, perplexity=5)
plot_dimred(MOFAobject, method="TSNE", color_by = "sample", stroke=0.1, dot_size=10, legend=F) +
  scale_fill_manual(values=opts$celltype.colors)


##################
## Plot weights ##
##################

plot_weights(MOFAobject, view="expr", factor = 1, nfeatures = 20, scale = F)
plot_top_weights(MOFAobject, view="chromvar", factor = 6, nfeatures = 20, scale = F)


###############
## Plot data ##
###############

plot_data_heatmap(
  MOFAobject,
  factor = 1,
  features = 25,
  view = "chromvar",
  denoise = FALSE,
  legend = TRUE,
  # min.value = 0, max.value = 6,
  cluster_rows = T, cluster_cols = F,
  show_colnames = F, show_rownames = T,
  scale="row"
)

plot_data_scatter(
  MOFAobject, 
  factor = 1, 
  view = 1, 
  color_by = "sample", 
  sign = "negative",
  features = 8, 
  dot_size = 3, 
  legend=FALSE
) + scale_fill_manual(values=opts$celltype.colors)


