
#####################
## Define settings ##
#####################
# load 
# Load default settings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/utils.R")
} else if (grepl("Workstation",Sys.info()['nodename'])){
  source("/home/lijingyu/gastrulation/gastrulation_multiome_10x/settings.R")
  source("/home/lijingyu/gastrulation/gastrulation_multiome_10x/utils.R")
} else {
  stop("Computer not recognised")
}

# I/O
io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/trajectories/blood_trajectory")
io$trajectories.inputdir <- c(
  "blood" = paste0(io$basedir,"/results/rna_atac/rna_vs_acc/trajectories/blood_trajectory")
)

# Options

###############
## Load data ##
###############

# Load data.tables with RNA + chromVAR estimates
DORC_rna_dt <- names(io$trajectories.inputdir) %>%
  map(function(x) fread(sprintf(sprintf("%s/DORC_rna_raw.txt.gz",io$trajectories.inputdir[[x]]))) %>%
        .[,trajectory:=x]) %>%
  rbindlist# %>%
  # .[chromvar_zscore<0,chromvar_zscore:=0]

# Load correlation results
# cor_dt <- names(io$trajectories.inputdir) %>% 
#   map(function(x) fread(sprintf(sprintf("%s/correlation_results.txt.gz",io$trajectories.inputdir[[x]]))) %>%
#         .[,trajectory:=x]) %>%
#   rbindlist %>% setorder(-trajectory,padj_fdr)
# 
# unique(cor_dt$gene)

# cor_dt[sig==T & rna_sign=="Down" & chromvar_sign=="Down"] %>% View
# cor_dt[sig==T & rna_sign=="Up" & chromvar_sign=="Up"] %>% View

################
## Parse data ##
################

# Remove outliers
opts$cutoff <- 3 * 1.96
DORC_rna_dt <- DORC_rna_dt %>%
  .[,c("mean_expr","sd_expr"):=list(mean(expr,na.rm=T),sd(expr,na.rm=T)),by=c("trajectory","gene")] %>%
  .[,z_score:=abs(expr-mean_expr)/sd_expr, by=c("trajectory","gene","cell")] %>%
  .[z_score<=opts$cutoff] %>% .[,c("mean_expr"):=NULL]

#####################################################################################
## For each region of the pseudotime, plot the difference between RNA and chromvar ##
#####################################################################################

opts$breaks <- 8

DORC_rna_dt.filt <- DORC_rna_dt %>% copy %>%
  .[,expr:=minmax.normalisation(expr), by=c("trajectory","gene")] %>%
  .[,DORC_score:=minmax.normalisation(DORC_score), by=c("trajectory","gene")] %>%
  # .[,expr:=expr/max(expr), by=c("trajectory","gene")] %>%
  # .[,chromvar_zscore:=chromvar_zscore/max(chromvar_zscore), by=c("trajectory","gene")] %>%
  .[,diff:=DORC_score-expr] %>%
  .[,pseudotime_group:=cut(PC1,breaks=opts$breaks), by=c("trajectory")] 
  # .[,.(diff=mean(diff), diff_se=sd(diff)),by=c("trajectory","pseudotime_group","gene")]

genes.to.plot <- 'Mrap'
for (i in unique(DORC_rna_dt.filt$trajectory)) {
  
  to.plot <- DORC_rna_dt.filt %>% 
    .[trajectory==i] %>% 
    .[gene==genes.to.plot]
    
  p <- ggplot(to.plot, aes(x=pseudotime_group, y=diff)) +
    geom_boxplot(outlier.shape=NA, coef=1) +
    # stat_smooth(method="loess", color="black", alpha=0.75, span=0.5) +
    # geom_rug(aes(color=celltype.predicted), sides="b") +
    geom_hline(yintercept=0, linetype="dashed") +
    # scale_color_manual(values=opts$celltype.colors) +
    # scale_fill_distiller(palette = "YlOrRd", direction=1) +
    # guides(fill=F, color=F) +
    # coord_cartesian(ylim=c(-0.45,0.45)) +
    labs(x="Pseudotime (discretised)", y="DORC_score - RNA difference") +
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

# genes.to.plot <- DORC_rna_dt[,list(var_expr=var(expr), var_chromvar=var(chromvar_zscore)),by="gene"] %>%
#   .[,var_combined:=var_expr*var_chromvar] %>% setorder(-var_combined) %>% .$gene %>% head(n=50)
DORC_rna_dt.filt <- DORC_rna_dt %>% copy %>%
  # .[gene%in%genes.to.plot] %>%
  .[,expr:=minmax.normalisation(expr), by=c("trajectory","gene")] %>%
  # .[,chromvar_zscore:=chromvar_zscore/max(chromvar_zscore), by=c("trajectory","gene")] %>%
  .[,DORC_score:=minmax.normalisation(DORC_score), by=c("trajectory","gene")] %>%
  .[,diff:=DORC_score-expr]

# DORC_rna_dt.filt[trajectory=="blood" & gene=="HOXB9"] %>% View

##########
## Plot ##
##########
genes.to.plot <- c('Trim10')
for (i in unique(DORC_rna_dt.filt$trajectory)) {
  outdir <- file.path(io$outdir,i); dir.create(outdir, showWarnings = F)
  
  for (j in genes.to.plot) {
    
    to.plot <- DORC_rna_dt.filt %>%
      .[trajectory==i & gene==j] %>%
      melt(id.vars=c("cell","PC1","celltype" ), measure.vars=c("DORC_score","expr","diff"), variable.name="modality")# 
    to.plot$celltype <- factor(to.plot$celltype,levels =opts$celltypes )
    p1 <- ggplot(to.plot, aes(x=PC1, y=value)) +
      geom_point(aes(fill=celltype), size=1.75, shape=21, stroke=0.1) +
      stat_smooth(method="loess", color="black", alpha=0.75, span=0.5) +
      geom_rug(aes(color=celltype), sides="b") +
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
    my_comparisons <- list(c(opts$celltypes[1],opts$celltypes[2]),c(opts$celltypes[2],opts$celltypes[3]),
                       c(opts$celltypes[3],opts$celltypes[4]),c(opts$celltypes[4],opts$celltypes[5]),
                       c(opts$celltypes[5],opts$celltypes[6]))
    p2 <- ggboxplot(to.plot, x="celltype", y="value", fill="celltype", outlier.shape=NA) +
      scale_fill_manual(values=opts$celltype.colors) +
      labs(x="", y="") +
      # stat_compare_means(comparisons = list( c("Somitic_mesoderm", "NMP"), c("NMP", "Spinal_cord") ), label="p.signif", hide.ns=T) +
      facet_wrap(~modality, nrow=3, scales="free_y") +
      # coord_cartesian(ylim=c(0,5)) +
      theme_classic() +
      guides(x = guide_axis(angle = 90)) +
      theme(
        # axis.text.x = element_text(color="black", size=rel(0.75)),
        # legend.title = element_blank(),
        # legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
      stat_compare_means(comparisons = my_comparisons)
    p <- cowplot::plot_grid(plotlist=list(p1,p2), nrow = 1, rel_widths = c(1/2,1/2))
    
    png(sprintf("%s/%s_rna_DORC_diff_pseudotime.png",outdir,j), width = 800, height = 700)
    print(p)
    dev.off()
  }
}


