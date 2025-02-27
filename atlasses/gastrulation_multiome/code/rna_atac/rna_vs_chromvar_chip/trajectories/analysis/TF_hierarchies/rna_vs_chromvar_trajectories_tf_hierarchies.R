
#####################
## Define settings ##
#####################

# Load default settings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
} else {
  stop("Computer not recognised")
}

# I/O
io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/trajectories/tf_hierarchies"); dir.create(io$outdir, showWarnings=F)
io$trajectories.inputdir <- c(
  # "ectoderm" = paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/trajectories/ectoderm_trajectory_knn50"),
  # "endoderm" = paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/trajectories/endoderm_trajectory_knn50"),
  # "mesoderm" = paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/trajectories/mesoderm_trajectory_knn50"),
  "blood" = paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/trajectories/blood_trajectory_knn50")
)

# Options

###############
## Load data ##
###############

# Load data.tables with RNA + chromVAR estimates
rna_chromvar.dt <- names(io$trajectories.inputdir) %>%
  map(function(x) fread(sprintf(sprintf("%s/chromvar_rna.txt.gz",io$trajectories.inputdir[[x]]))) %>%
        .[,trajectory:=x]) %>%
  rbindlist# %>%
  # .[chromvar_zscore<0,chromvar_zscore:=0]

# Load correlation results
cor.dt <- names(io$trajectories.inputdir) %>% 
  map(function(x) fread(sprintf(sprintf("%s/correlation_results.txt.gz",io$trajectories.inputdir[[x]]))) %>%
        .[,trajectory:=x]) %>%
  rbindlist %>% setorder(-trajectory,padj_fdr)

# Load trajectories
trajectory.dt <- names(io$trajectories.inputdir) %>% 
  map(function(x) fread(sprintf(sprintf("%s/%s_trajectory.txt.gz",io$trajectories.inputdir[[x]],x)), select=c(1,2)) %>%
        .[,trajectory:=x]) %>%
  rbindlist

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[cell%in%trajectory.dt$cell] %>%
  merge(trajectory.dt, by="cell") %>%
  setnames("celltype.predicted","celltype")

trajectory.dt <- trajectory.dt %>% 
  merge(sample_metadata[,c("cell","celltype")], by="cell") %>% 
  setorder(PC1)

#################
## Filter data ##
#################

# Filter by correlation
cor_filt.dt <- cor.dt[sig==T & abs(r)>0.30 & rna_sign=="Up"]

# cor_filt.dt <- cor.dt[sig==T & abs(r)>0.30]
rna_chromvar.dt <- rna_chromvar.dt[gene%in%cor_filt.dt$gene]

# Filter by coverage

cells <- intersect(sample_metadata$cell,rna_chromvar.dt$cell)
trajectory.dt <- trajectory.dt[cell%in%cells]
# sample_metadata <- sample_metadata[cell%in%cells]
rna_chromvar_filt.dt <- rna_chromvar.dt[cell%in%cells] %>% setkey(cell) %>% .[trajectory.dt$cell]

################
## Clustering ##
################

to.plot <- rna_chromvar_filt.dt %>% dcast(gene~cell,value.var="expr") %>% 
  tibble::column_to_rownames("gene") %>%
  .[,trajectory.dt$cell] %>%
  apply(1,minmax.normalisation) %>% t

kmeans.out <- kmeans(to.plot, centers=5)

clusters <- unique(kmeans.out$cluster) %>% as.character

for (i in clusters) {
  
  tfs.to.plot <- which(kmeans.out$cluster==i) %>% names
  
  to.plot <- rna_chromvar_filt.dt %>% 
    .[gene%in%tfs.to.plot] %>%
    merge(trajectory.dt,by="cell")# %>%
    # .[,expr:=expr/max(expr),by=c("gene")]
  
  p <- ggplot(to.plot, aes(x=PC1, y=expr)) +
    # ggrastr::geom_point_rast(aes(fill=celltype), size=1.25, shape=21, stroke=0.1) +
    geom_point(aes(fill=celltype), size=1, shape=21, stroke=0.1) +
    stat_smooth(method="loess", color="black", alpha=0.75, span=0.5) +
    geom_rug(aes(color=celltype), sides="b") +
    facet_wrap(~gene, scales="free_y") +
    scale_color_manual(values=opts$celltype.colors) +
    scale_fill_manual(values=opts$celltype.colors) +
    guides(fill=F, color=F) +
    # scale_fill_manual(values=opts$celltype.colors) +
    # scale_fill_brewer(palette="Dark2") +
    labs(x="Pseudotime", y="RNA expression") +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(size=rel(0.75), color="black"),
      axis.ticks.x = element_blank(),
      legend.title = element_blank(),
      legend.position="top"
    )
  
  pdf(sprintf("%s/rna_chromvar_vs_pseudotime_cluster%s.pdf",io$outdir,i))
  print(p)
  dev.off()
  
}
#############
## Heatmap ##
#############

annotation_col.df <- trajectory.dt[,c("cell","celltype")] %>% tibble::column_to_rownames("cell")
mycolors <- list(celltype = opts$celltype.colors[unique(sample_metadata$celltype)])

# to.plot <- rna_chromvar_filt.dt %>% dcast(gene~cell,value.var="chromvar_zscore") %>% 
to.plot <- rna_chromvar_filt.dt %>% dcast(gene~cell,value.var="expr") %>% 
  tibble::column_to_rownames("gene") %>%
  .[,trajectory.dt$cell]
stopifnot(all(rownames(annotation_col.df)==colnames(to.plot)))

pheatmap(
  mat = to.plot, 
  cluster_cols = F, cluster_rows = T,
  # color = colorRampPalette(c("white", "black"))(100),
  annotation_col = annotation_col.df,
  annotation_colors = mycolors,
  show_colnames = FALSE, show_rownames = TRUE,
  legend = F, annotation_legend = F,
  treeheight_row = 0,
  scale = "row",
  width = 10, height = 5,
  filename = sprintf("%s/foo.pdf",io$outdir)
)


#############
## Explore ##
#############



tfs.to.plot <- c("RUNX1", "GATA1", "IKZF1", "KLF1", "FOXO3", "JUN")

to.plot <- rna_chromvar_filt.dt %>% 
  .[gene%in%tfs.to.plot] %>%
  melt(id.vars=c("cell","gene"), measure.vars=c("chromvar_zscore","expr"), variable.name="modality") %>%
  merge(trajectory.dt,by="cell") %>%
  .[,value_scaled:=value/max(value), by=c("modality","gene")]

ggplot(to.plot, aes(x=PC1, y=value, group=gene)) +
  stat_smooth(aes(color=gene), method="loess", alpha=0.75, span=0.5, se=F) +
  scale_color_brewer(palette="Dark2") +
  facet_wrap(~modality, nrow=1, scales="free_y") +
  labs(x="Pseudotime", y="") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=rel(0.75), color="black"),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    legend.position="right"
  )

pdf(sprintf("%s/rna_chromvar_vs_pseudotime_cluster%s.pdf",io$outdir,i))
print(p)
dev.off()