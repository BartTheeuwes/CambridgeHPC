
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
io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/trajectories/rna_vs_chromvar_temporal_dynamics")
io$trajectories.inputdir <- c(
  "ectoderm" = paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/trajectories/ectoderm_trajectory_knn50"),
  "endoderm" = paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/trajectories/endoderm_trajectory_knn50"),
  "mesoderm" = paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/trajectories/mesoderm_trajectory_knn50"),
  "blood" = paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/trajectories/blood_trajectory_knn50")
)

# Options

###############
## Load data ##
###############

# Load data.tables with RNA + chromVAR estimates
chromvar_rna_dt <- names(io$trajectories.inputdir) %>%
  map(function(x) fread(sprintf(sprintf("%s/chromvar_rna.txt.gz",io$trajectories.inputdir[[x]]))) %>%
        .[,trajectory:=x]) %>%
  rbindlist# %>%
  # .[chromvar_zscore<0,chromvar_zscore:=0]

# Load correlation results
cor_dt <- names(io$trajectories.inputdir) %>% 
  map(function(x) fread(sprintf(sprintf("%s/correlation_results.txt.gz",io$trajectories.inputdir[[x]]))) %>%
        .[,trajectory:=x]) %>%
  rbindlist %>% setorder(-trajectory,padj_fdr)

unique(cor_dt$gene)

# cor_dt[sig==T & rna_sign=="Down" & chromvar_sign=="Down"] %>% View
# cor_dt[sig==T & rna_sign=="Up" & chromvar_sign=="Up"] %>% View

################
## Parse data ##
################

# Remove outliers
opts$cutoff <- 3 * 1.96
chromvar_rna_dt <- chromvar_rna_dt %>%
  .[,c("mean_expr","sd_expr"):=list(mean(expr,na.rm=T),sd(expr,na.rm=T)),by=c("trajectory","gene")] %>%
  .[,z_score:=abs(expr-mean_expr)/sd_expr, by=c("trajectory","gene","cell")] %>%
  .[z_score<=opts$cutoff] %>% .[,c("mean_expr"):=NULL]

#####################################################################################
## For each region of the pseudotime, plot the difference between RNA and chromvar ##
#####################################################################################

opts$breaks <- 8

chromvar_rna_dt.filt <- chromvar_rna_dt %>% copy %>%
  .[,expr:=minmax.normalisation(expr), by=c("trajectory","gene")] %>%
  .[,chromvar_zscore:=minmax.normalisation(chromvar_zscore), by=c("trajectory","gene")] %>%
  # .[,expr:=expr/max(expr), by=c("trajectory","gene")] %>%
  # .[,chromvar_zscore:=chromvar_zscore/max(chromvar_zscore), by=c("trajectory","gene")] %>%
  .[,diff:=expr-chromvar_zscore] %>%
  .[,pseudotime_group:=cut(PC1,breaks=opts$breaks), by=c("trajectory")] %>%
  .[,.(diff=mean(diff), diff_se=sd(diff)),by=c("trajectory","pseudotime_group","gene")]

for (i in unique(chromvar_rna_dt.filt$trajectory)) {
  
  to.plot <- chromvar_rna_dt.filt %>% 
    .[trajectory==i] %>% 
    .[gene%in%cor_dt[trajectory==i & sig==T & abs(r>0.15) & rna_sign=="Up",gene]]
    
  p <- ggplot(to.plot, aes(x=pseudotime_group, y=diff)) +
    geom_boxplot(outlier.shape=NA, coef=1) +
    # stat_smooth(method="loess", color="black", alpha=0.75, span=0.5) +
    # geom_rug(aes(color=celltype.predicted), sides="b") +
    geom_hline(yintercept=0, linetype="dashed") +
    # scale_color_manual(values=opts$celltype.colors) +
    # scale_fill_distiller(palette = "YlOrRd", direction=1) +
    # guides(fill=F, color=F) +
    # coord_cartesian(ylim=c(-0.45,0.45)) +
    labs(x="Pseudotime (discretised)", y="RNA - chromVAR difference") +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.title = element_blank(),
      legend.position="top"
    )
  pdf(sprintf("%s/diff_vs_pseudotime_boxplots_%s.pdf",io$outdir,i,i), width=7, height=4)
  print(p)
  dev.off()
}


#########################
## Plot individual TFs ##
#########################

# genes.to.plot <- chromvar_rna_dt[,list(var_expr=var(expr), var_chromvar=var(chromvar_zscore)),by="gene"] %>%
#   .[,var_combined:=var_expr*var_chromvar] %>% setorder(-var_combined) %>% .$gene %>% head(n=50)
chromvar_rna_dt.filt <- chromvar_rna_dt %>% copy %>%
  # .[gene%in%genes.to.plot] %>%
  .[,expr:=minmax.normalisation(expr), by=c("trajectory","gene")] %>%
  # .[,chromvar_zscore:=chromvar_zscore/max(chromvar_zscore), by=c("trajectory","gene")] %>%
  .[,chromvar_zscore:=minmax.normalisation(chromvar_zscore), by=c("trajectory","gene")] %>%
  .[,diff:=expr-chromvar_zscore]

# chromvar_rna_dt.filt[trajectory=="blood" & gene=="HOXB9"] %>% View

##########
## Plot ##
##########

for (i in unique(chromvar_rna_dt.filt$trajectory)) {
  outdir <- file.path(io$outdir,i); dir.create(outdir, showWarnings = F)
  genes.to.plot <- cor_dt[trajectory==i & sig==T & abs(r>0.25) & rna_sign=="Up",gene]
  for (j in genes.to.plot) {
    
    to.plot <- chromvar_rna_dt.filt %>%
      .[trajectory==i & gene==j] %>%
      melt(id.vars=c("cell","PC1","celltype.predicted"), measure.vars=c("chromvar_zscore","expr","diff"), variable.name="modality")# %>%
    
    p <- ggplot(to.plot, aes(x=PC1, y=value)) +
      geom_point(aes(fill=celltype.predicted), size=1.75, shape=21, stroke=0.1) +
      stat_smooth(method="loess", color="black", alpha=0.75, span=0.5) +
      geom_rug(aes(color=celltype.predicted), sides="b") +
      geom_hline(yintercept=0, linetype="dashed") +
      facet_wrap(~modality, nrow=3, scales="free_y") +
      scale_color_manual(values=opts$celltype.colors) +
      scale_fill_manual(values=opts$celltype.colors) +
      guides(fill=F, color=F) +
      labs(x="Pseudotime", y=j) +
      theme_classic() +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.position="top"
      )
   
    png(sprintf("%s/%s_rna_chromvar_diff_pseudotime.png",outdir,j), width = 600, height = 800)
    print(p)
    dev.off()
  }
}


