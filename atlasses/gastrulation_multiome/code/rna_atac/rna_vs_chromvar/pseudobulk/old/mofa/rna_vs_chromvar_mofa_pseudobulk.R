library(MOFA2)
options(ggrepel.max.overlaps = Inf)

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
io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/pseudobulk/mofa")
# io$mofa.output <- paste0(io$basedir, "/results/rna_atac/mofa/NMPs/mofa_model.rds")

# Options
opts$celltypes <- c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  # "PGC",
  "Anterior_Primitive_Streak",
  "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  "Mixed_mesoderm",
  "Intermediate_mesoderm",
  "Caudal_Mesoderm",
  "Paraxial_mesoderm",
  "Somitic_mesoderm",
  "Pharyngeal_mesoderm",
  "Cardiomyocytes",
  "Allantois",
  "ExE_mesoderm",
  "Mesenchyme",
  "Haematoendothelial_progenitors",
  "Endothelium",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
  # "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm",
  "Visceral_endoderm",
  "ExE_endoderm",
  "ExE_ectoderm",
  "Parietal_endoderm"
)

##############################
## Load RNA expression data ##
##############################

sce.pseudobulk <- readRDS(io$rna.pseudobulk.sce)

rownames(sce.pseudobulk) <- toupper(rownames(sce.pseudobulk))

#######################################
## Load chromatin accessibility data ##
#######################################

opts$motif_annotation <- "Motif_cisbp" # Motif_JASPAR2020_human
chromvar.se.pseudobulk <- readRDS(sprintf("%s/pseudobulk_DeviationMatrix_%s_summarized_experiment.rds",io$archR.pseudobulk.chromvar,opts$motif_annotation))
chromvar.se.pseudobulk <- chromvar.se.pseudobulk[rowData(chromvar.se.pseudobulk)$seqnames=="z",]
rownames(chromvar.se.pseudobulk) <- rowData(chromvar.se.pseudobulk)$name %>% toupper %>% stringr::str_split(.,"_") %>% map_chr(1)

###########
## Merge ##
###########

samples <- opts$celltypes[opts$celltypes %in% intersect(colnames(chromvar.se.pseudobulk), colnames(sce.pseudobulk))]
chromvar.se.pseudobulk <- chromvar.se.pseudobulk[,samples]
sce.pseudobulk <- sce.pseudobulk[,samples]

################################
## Load motif2gene annotation ##
################################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/atac/archR/load_motif_annotation.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/atac/archR/load_motif_annotation.R")
} else {
  stop("Computer not recognised")
}

motif2gene.dt <- motif2gene.dt %>%
  .[gene%in%rownames(chromvar.se.pseudobulk) & gene%in%rownames(rna.sce.pseudobulk)]  %>%
  .[,N:=length(unique(motif)),by="gene"] %>% .[N==1] %>% .[,N:=NULL]

rna.sce.pseudobulk <- rna.sce.pseudobulk %>% .[motif2gene.dt$gene,]
chromvar.se.pseudobulk <- chromvar.se.pseudobulk %>% .[motif2gene.dt$gene,]

############
## Filter ##
############

##########################
## Tensor decomposition ##
##########################

##############
## Run MOFA ##
##############

# Prepare data
chromvar.mtx <- assay(chromvar.se.pseudobulk)
rna.mtx <- logcounts(sce.pseudobulk)
rownames(rna.mtx) <- tolower(rownames(rna.mtx))
rownames(rna.mtx) <- paste(toupper(substr(rownames(rna.mtx), 1, 1)), substr(rownames(rna.mtx), 2, nchar(rownames(rna.mtx))), sep="")

# Create MOFA object
MOFAobject <- create_mofa_from_matrix(list("RNA"=rna.mtx, "chromVAR"=chromvar.mtx))

# Model options
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 8

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

p <- plot_variance_explained(MOFAobject, x="view", y="factor", max_r2 = 30)

pdf(sprintf("%s/var_explained.pdf",io$outdir), width=7.5, height=6)
print(p)
dev.off()

p <- plot_factor(MOFAobject, factors = 1:6, color_by = "sample", dodge=F, dot_size=5, legend = FALSE) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

pdf(sprintf("%s/beeswarm_factors.pdf",io$outdir), width=8, height=5)
print(p)
dev.off()

plot_factor(MOFAobject, factors = 6, color_by = "Mecom_chromvar", dot_size=7)


p <- plot_factors(MOFAobject, factors = c(1,3), dot_size=5, color_by = "sample", legend = FALSE) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
pdf(sprintf("%s/factors_13.pdf",io$outdir), width=5.5, height=5)
print(p)
dev.off()

plot_factors(MOFAobject, factors = c(2,3), dot_size=7, color_by = "STAT3")

plot_factors(MOFAobject, factors = c(1,3), dot_size=5, color_by = "sample", legend = FALSE) +
  scale_fill_manual(values=opts$celltype.colors)

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

plot_weights(MOFAobject, view="RNA", factor = 1, nfeatures = 20, scale = F)
plot_top_weights(MOFAobject, view="RNA", factor = 1, nfeatures = 20, sign = "positive")
plot_top_weights(MOFAobject, view="chromVAR", factor = 1, nfeatures = 20, scale = F)


###############
## Plot data ##
###############

plot_data_heatmap(
  MOFAobject,
  factor = 1,
  features = 25,
  view = "RNA",
  denoise = FALSE,
  legend = TRUE,
  # min.value = 0, max.value = 6,
  cluster_rows = T, cluster_cols = F,
  show_colnames = F, show_rownames = T,
  scale="row"
)

for (factor in factors_names(MOFAobject)[1:5]) {
  for (view in views_names(MOFAobject)) {
    for (sign in c("positive","negative")) {
      p <- plot_data_scatter(
        MOFAobject, 
        factor = factor, 
        view = view, 
        color_by = "sample", 
        sign = sign,
        features = 9, 
        dot_size = 3, 
        legend = FALSE
      ) + scale_fill_manual(values=opts$celltype.colors)
      
      pdf(sprintf("%s/scatterplots/%s_%s_%s_scatterplot.pdf",io$outdir,factor,view,sign), width=7.5, height=6)
      print(p)
      dev.off()
    }
  }
}



#########################
## Contribution scores ##
#########################

MOFAobject <- calculate_contribution_scores(MOFAobject)

to.plot <- MOFAobject@samples_metadata[,c("sample","RNA_contribution")] %>% 
  as.data.table(keep.rownames = T) %>%
  setorder(-RNA_contribution)

ggbarplot(to.plot, x="sample", y="RNA_contribution", fill="sample", order=to.plot$sample) +
  labs(x="", y="RNA contribution") +
  geom_hline(yintercept=0.5, linetype="dashed") +
  scale_fill_manual(values=opts$celltype.colors) +
  theme(
    axis.text.x = element_text(size=rel(0.7), color="black", angle=40, vjust=1, hjust=1),
    axis.ticks.x = element_blank(),
    legend.position = "none"
    # legend.position = "right", legend.text = element_text(size=rel(0.75))
  )

#############################
## RNA vs chromVAR weights ##
#############################

w.rna <- get_weights(MOFAobject, views="RNA", factor="all", as.data.frame=T) %>% as.data.table %>%
  .[,feature:=gsub("_RNA","",feature) %>% toupper] %>%
  .[,value:=value/max(abs(value)),by=c("factor")]
w.acc <- get_weights(MOFAobject, views="chromVAR", factor="all", as.data.frame=T) %>% as.data.table %>%
  .[,feature:=gsub("_chromVAR","",feature)] %>%
  .[,value:=value/max(abs(value)),by=c("factor")]

# Merge loadings
w.dt <- merge(
  w.rna[,c("feature","factor","value")], 
  w.acc[,c("feature","factor","value")], 
  by = c("feature","factor")
)

# Scatterplots
for (i in unique(w.dt$factor)) {
  
  to.plot <- w.dt[factor==i]
  to.label <- w.dt %>% 
    .[factor==i] %>%
    .[,value:=abs(value.x)+abs(value.y)] %>% setorder(-value) %>% head(n=20)
  
  p <- ggpubr::ggscatter(to.plot, x="value.x", y="value.y", size=1.25, add="reg.line", add.params = list(color="blue", fill="lightgray"), conf.int=TRUE) +
    coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) +
    scale_x_continuous(breaks=c(-1,0,1)) +
    scale_y_continuous(breaks=c(-1,0,1)) +
    ggrepel::geom_text_repel(data=to.label, aes(x=value.x, y=value.y, label=feature), size=3,  max.overlaps=100) +
    geom_vline(xintercept=0, linetype="dashed") +
    geom_hline(yintercept=0, linetype="dashed") +
    stat_cor(method = "pearson") +
    labs(x=sprintf("RNA weights (%s)",i), y=sprintf("chromVAR weights (%s)",i))
  
  pdf(sprintf("%s/scatterplots_weights/%s_RNA_vs_chromVAR_weights.pdf",io$outdir,i), width=6.5, height=5)
  print(p)
  dev.off()
}

