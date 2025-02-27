source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# Options
opts$motif_annotation <- "CISBP"

# I/O
io$rna_sce_pseudobulk_file <- file.path(io$basedir,"results/rna/pseudobulk/celltype/SingleCellExperiment_pseudobulk.rds")
io$atac_chromvar_chip_pseudobulk_file <- file.path(io$basedir,sprintf("results/atac/archR/chromvar_chip/pseudobulk/chromVAR_chip_%s_archr.rds",opts$motif_annotation))
io$outdir <- file.path(io$basedir,"results/rna_atac/rna_vs_chromvar_chip/pseudobulk/per_gene/fig"); dir.create(io$outdir, showWarnings=F, recursive = T)

#######################################
## Load pseudobulk RNA and ATAC data ##
#######################################

# Load pseudobulk RNA expression
rna_pseudobulk.sce <- readRDS(io$rna_sce_pseudobulk_file)

# Load chromVAR-ChIP matrix
atac_pseudobulk_chromvar.se <- readRDS(io$atac_chromvar_chip_pseudobulk_file)

# Fetch TF RNA expression matrix
rna_pseudobulk_tf.se <- rna_pseudobulk.sce[str_to_title(rownames(atac_pseudobulk_chromvar.se)),]
rownames(rna_pseudobulk_tf.se) <- toupper(rownames(rna_pseudobulk_tf.se))

########################
## Prepare data table ##
########################

atac_chromvar_pseudobulk.dt <- assay(atac_pseudobulk_chromvar.se) %>% t %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","celltype") %>%
  melt(id.vars=c("celltype"), variable.name="gene", value.name="chromvar_zscore")

rna_tf_pseudobulk.dt <- logcounts(rna_pseudobulk_tf.se) %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","gene") %>%
  data.table::melt(id.vars="gene", variable.name="celltype", value.name="expr")

###########
## Merge ##
###########

rna_chromvar.dt <- merge(
  rna_tf_pseudobulk.dt,
  atac_chromvar_pseudobulk.dt,
  by = c("celltype","gene")
)

######################################
## Scatter plot of individual genes ##
######################################

genes.to.plot <- unique(rna_chromvar.dt$gene)

i <- "FOXA2"
for (i in genes.to.plot) {

  to.plot <- rna_chromvar.dt[gene==i]

  to.plot.text <- rbind(
    to.plot %>% setorder(-expr) %>% head(n=7),
    to.plot %>% setorder(-chromvar_zscore) %>% head(n=7)
  ) %>% unique
  
  p <- ggscatter(to.plot, x="expr", y="chromvar_zscore", fill="celltype", size=5, shape=21, 
                  add="reg.line", add.params = list(color="black", fill="lightgray"), conf.int=TRUE) +
    stat_cor(method = "pearson", label.x.npc = "middle", label.y.npc = "bottom") +
    ggrepel::geom_text_repel(data=to.plot.text, aes(label=celltype), size=3) +
    scale_fill_manual(values=opts$celltype.colors) +
    labs(x=sprintf("%s expression",i), y=sprintf("Accessibility of %s targets (z-score)",i)) +
    guides(fill="none") +
    theme(
      axis.text = element_text(size=rel(0.7))
    )
  
  pdf(file.path(io$outdir,sprintf("%s_%s_rna_vs_chromvar_chip_pseudobulk.pdf",i,opts$motif_annotation)), width = 8, height = 5)
  print(p)
  dev.off()
}

