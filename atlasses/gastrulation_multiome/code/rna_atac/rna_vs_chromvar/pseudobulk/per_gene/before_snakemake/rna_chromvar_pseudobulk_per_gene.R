#####################
## Define settings ##
#####################

# load default setings
source(here::here("settings.R"))
source(here::here("utils.R"))

# Options
opts$celltypes = c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  "PGC",
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
  "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm",
  "Visceral_endoderm",
  "ExE_endoderm",
  "ExE_ectoderm",
  "Parietal_endoderm"
)

opts$motif_annotation <- "Motif_cisbp_lenient"

# I/O
# io$archR.pseudobulk.deviations.se <- sprintf("%s/results/atac/archR/chromvar/pseudobulk/chromRNA_deviations_%s_pseudobulk_archr.rds",io$basedir,opts$motif_annotation)
# io$archR.pseudobulk.deviations.se <- sprintf("%s/results/atac/archR/chromvar/pseudobulk/chromVAR_deviations_summarized_experiment_%s_pseudobulk_all_peaks_archr.rds",io$basedir,opts$motif_annotation)
io$outdir <- file.path(io$basedir,"results/rna_atac/rna_vs_chromvar/pseudobulk/per_gene")

#######################################
## Load pseudobulk RNA and ATAC data ##
#######################################

source(here::here("rna_atac/load_rna_atac_pseudobulk.R"))

###########
## Merge ##
###########

rna_chromvar.dt <- merge(
  rna_tf_pseudobulk.dt,
  atac_chromvar_pseudobulk.dt,
  by = c("celltype","gene")
)

length(unique(rna_chromvar.dt$gene))

##########################
## Correlation analysis ##
##########################

opts$threshold_fdr <- 0.10

cor_rna_vs_chromvar_per_gene.dt <- rna_chromvar.dt %>% copy %>%
  .[,c("chromvar_zscore","expr"):=list(chromvar_zscore + rnorm(n=.N,mean=0,sd=1e-5), expr + rnorm(n=.N,mean=0,sd=1e-5))] %>% # add some noise 
  .[, .(V1 = unlist(cor.test(chromvar_zscore, expr)[c("estimate", "p.value")])), by = c("gene")] %>%
  .[, para := rep(c("r","p"), .N/2)] %>% 
  data.table::dcast(gene ~ para, value.var = "V1") %>%
  .[,"padj_fdr" := list(p.adjust(p, method="fdr"))] %>%
  .[, sig := p<=opts$threshold_fdr] %>% 
  setorder(p, na.last = T)

# Save
fwrite(cor_rna_vs_chromvar_per_gene.dt, file.path(io$outdir,"cor_rna_vs_chromVAR+_pseudobulk.txt.gz"), sep="\t", quote=F)

#############################################
## Load pre-computed correlation estimates ##
#############################################

cor_rna_vs_chromvar_per_gene.dt <- fread(paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/pseudobulk/per_gene/cor_rna_vs_chromVAR+_pseudobulk.txt.gz")) %>%
  .[,cor_sign:=as.factor(c("Repressor","Activator")[(r>0)+1])]

##################
## Volcano plot ##
##################

to.plot <- cor_rna_vs_chromvar_per_gene.dt %>%
  .[padj_fdr<=1e-14,padj_fdr:=1e-14] %>%
  .[,log_pval:=-log10(padj_fdr+1e-100)] %>%
  .[,dot_size:=minmax.normalisation(log_pval)]

negative_hits <- to.plot[sig==TRUE & r<0,gene]
positive_hits <- to.plot[sig==TRUE & r>0,gene]
all <- nrow(to.plot)

xlim <- max(abs(to.plot$r), na.rm=T)
ylim <- max(-log10(to.plot$padj_fdr+1e-100), na.rm=T)

p <- ggplot(to.plot, aes(x=r, y=log_pval)) +
  geom_segment(aes(x=0, xend=0, y=0, yend=ylim-1), color="orange", size=0.5) +
  geom_jitter(aes(fill=sig, size=dot_size, alpha=dot_size), width=0.05, shape=21) + 
  # ggrepel::geom_text_repel(data=head(to.plot[sig==T & r<0],n=25), aes(x=r, y=log_pval, label=gene), size=3,  max.overlaps=100, segment.color = NA) +
  # ggrepel::geom_text_repel(data=head(to.plot[sig==T & r>0],n=25), aes(x=r, y=log_pval, label=gene), size=3,  max.overlaps=100, segment.color = NA) +
  ggrepel::geom_text_repel(data=to.plot[sig==T & r<(-0.50)][sample(.N,12)], aes(x=r, y=log_pval, label=gene), size=4,  max.overlaps=100, segment.color = NA) +
  ggrepel::geom_text_repel(data=to.plot[sig==T & r>0.75][sample(.N,12)], aes(x=r, y=log_pval, label=gene), size=4,  max.overlaps=100, segment.color = NA) +
  scale_fill_manual(values=c("black","red")) +
  # scale_size_manual(values=c(0.5,1)) +
  scale_size_continuous(range = c(0.2,2)) + 
  scale_alpha_continuous(range=c(0.25,1)) +
  scale_x_continuous(limits=c(-xlim-0.15,xlim+0.15)) +
  scale_y_continuous(limits=c(0,ylim+1)) +
  annotate("text", x=0, y=ylim+1, size=4, label=sprintf("(%d)", all)) +
  annotate("text", x=-xlim-0.15, y=ylim+1, size=4, label=sprintf("%d (-)",length(negative_hits))) +
  annotate("text", x=xlim+0.15, y=ylim+1, size=4, label=sprintf("%d (+)",length(positive_hits))) +
  labs(x="Pearson correlation (TF RNA expr vs TF chromRNA score)", y=expression(paste("-log"[10],"(p.value)"))) +
  theme_classic() +
  theme(
    axis.text = element_text(size=rel(0.75), color='black'),
    axis.title = element_text(size=rel(1.0), color='black'),
    legend.position="none"
  )

pdf(sprintf("%s/volcano_RNA_vs_chromVAR+.pdf",io$outdir), width = 7, height = 5)
print(p)
dev.off()


######################################
## Scatter plot of individual genes ##
######################################

genes.to.plot <- unique(rna_chromvar.dt$gene)
# genes.to.plot <- cor_rna_vs_chromvar_per_gene.dt[sig==T & abs(r)>0.25,gene]

for (i in genes.to.plot) {
  
  outfile <- sprintf("%s/individual_genes/%s_rna_vs_chromvar_correlated_peaks_pseudobulk.png",io$outdir,i)
  
  if (file.exists(outfile)) {
    print(sprintf("file or %s already exists...",i))
  } else {
    
    to.plot <- rna_chromvar.dt[gene==i]
    p1 <- ggscatter(to.plot, x="expr", y="chromvar_zscore", fill="celltype", size=4, shape=21, 
                    add="reg.line", add.params = list(color="black", fill="lightgray"), conf.int=TRUE) +
      stat_cor(method = "pearson", label.x.npc = "middle", label.y.npc = "bottom") +
      scale_fill_manual(values=opts$celltype.colors) +
      labs(x=sprintf("%s expression",i), y=sprintf("%s Motif accessibility (z-score)",i)) +
      guides(fill=F) +
      theme(
        axis.text = element_text(size=rel(0.7))
      )
    
    
    to.plot2 <- to.plot %>% melt(id.vars=c("celltype","gene"))
    p2 <- ggbarplot(to.plot2, x="celltype", y="value", fill="celltype") +
      facet_wrap(~variable, nrow=2, scales="free_y") +
      scale_fill_manual(values=opts$celltype.colors) +
      geom_hline(yintercept=0, linetype="dashed") +
      labs(x="", y="") +
      theme_classic() +
      guides(x = guide_axis(angle = 90)) +
      theme(
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()
      )
    
    
    p <- cowplot::plot_grid(plotlist=list(p1,p2), nrow = 1, rel_widths = c(1/2,1/2))
    
    # pdf(sprintf("%s/individual_genes/%s_rna_vs_chromvar.pdf",io$outdir,i), width=10, height=4)
    png(outfile, width = 1000, height = 500)
    print(p)
    dev.off()
  }
}

#######################################
## Stacked barplots per TF (all TFs) ##
#######################################

genes.to.plot <- unique(cor_rna_vs_chromvar_per_gene.dt$gene)

to.plot <- chromvar.dt %>% copy %>%
  .[gene%in%genes.to.plot] %>%
  .[chromvar_zscore<0,chromvar_zscore:=0] %>%
  .[,var:=var(chromvar_zscore),by="gene"] %>% .[var>0] %>% .[,var:=NULL] %>%
  .[,value:=minmax.normalisation(chromvar_zscore), by="gene"] %>%
  merge(cor_rna_vs_chromvar_per_gene.dt[,c("gene","cor_sign")])

# Sort by euclidean distance
# stopifnot(unique(to.plot$gene)%in%rownames(atac.chromvar.se))
# TF.order <- colnames(chromvar.mtx)[order(irlba::prcomp_irlba(t(chromvar.mtx), n=1)$x[,1])]
# TF.order <- TF.order[TF.order%in%to.plot$gene]
# to.plot[,gene:=factor(gene,levels=TF.order)]
                   
# to.plot <- to.plot[value>0.25] 

p <- ggplot(to.plot, aes(x=gene, y=value)) +
  geom_bar(aes(fill=celltype), stat="identity", color="black", position="fill") +
  # facet_wrap(~cor_sign, nrow=2, scales="free_x") +
  scale_fill_manual(values=opts$celltype.colors[names(opts$celltype.colors)%in%unique(to.plot$celltype)]) +
  theme_classic() +
  labs(x="", y="Motif accessibility (chromVAR+)") +
  # guides(x = guide_axis(angle = 90)) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    # axis.text.x = element_blank(),
    axis.text.x = element_text(color="black", size=rel(0.85)),
    # axis.text.y = element_text(color="black", size=rel(1.0)),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )

pdf(sprintf("%s/TF_chromVAR_celltype_stacked_barplots_neural_crest.pdf",io$outdir), width=6, height=6)
print(p)
dev.off()


###############################################
## Stacked barplots per TF (only marker TFs) ##
###############################################

marker_tfs.dt <- fread(io$rna.atlas.marker_TFs.up)

genes.to.plot <- c("SOX10", "FOXD3", "DLX2", "SOX9", "TFAP2A", "TFAP2B", "NR2F1", "PAX7","PKNOX1") # Neural crest
genes.to.plot <- c("SOX10", "FOXD3", "DLX2", "SOX9", "TFAP2A", "TFAP2B", "NR2F1", "PAX7","PKNOX1") # NMPs

to.plot <- chromvar.dt %>% copy %>%
  .[gene%in%genes.to.plot] %>%
  .[chromvar_zscore<0,chromvar_zscore:=0] %>%
  .[,var:=var(chromvar_zscore),by="gene"] %>% .[var>0] %>% .[,var:=NULL] %>%
  .[,value:=minmax.normalisation(chromvar_zscore), by="gene"] %>%
  merge(cor_rna_vs_chromvar_per_gene.dt[,c("gene","cor_sign")])

# Sort by euclidean distance
# stopifnot(unique(to.plot$gene)%in%rownames(atac.chromvar.se))
# TF.order <- colnames(chromvar.mtx)[order(irlba::prcomp_irlba(t(chromvar.mtx), n=1)$x[,1])]
# TF.order <- TF.order[TF.order%in%to.plot$gene]
# to.plot[,gene:=factor(gene,levels=TF.order)]

to.plot <- to.plot[value>0.25] 

p <- ggplot(to.plot, aes(x=gene, y=value)) +
  geom_bar(aes(fill=celltype), stat="identity", color="black", position="fill") +
  # facet_wrap(~cor_sign, nrow=2, scales="free_x") +
  scale_fill_manual(values=opts$celltype.colors[names(opts$celltype.colors)%in%unique(to.plot$celltype)]) +
  theme_classic() +
  labs(x="", y="Motif accessibility (chromVAR+)") +
  # guides(x = guide_axis(angle = 90)) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    # axis.text.x = element_blank(),
    axis.text.x = element_text(color="black", size=rel(0.85)),
    # axis.text.y = element_text(color="black", size=rel(1.0)),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )

pdf(sprintf("%s/TF_chromVAR_celltype_stacked_barplots_neural_crest.pdf",io$outdir), width=6, height=6)
print(p)
dev.off()

###################################
## Variability vs predictability ##
###################################

to.plot <- cor_rna_vs_chromvar_per_gene.dt %>% 
  merge(rna_chromvar.dt[,.(var_expr=var(expr)),by="gene"], by = "gene") %>% 
  .[,r2:=r**2] %>% setorder(-r2)

ggscatter(to.plot, x="r2", y="var_expr", size=1) + 
  # add="reg.line", add.params = list(color="black", fill="lightgray"), conf.int=TRUE) +
  # stat_cor(method = "pearson", label.x.npc = "middle", label.y.npc = "bottom") +
  ggrepel::geom_text_repel(data=head(to.plot[r2<0.15 & var_expr>2],n=10), aes(x=r2, y=var_expr, label=gene), size=3, color="darkred", max.overlaps=Inf) +
  ggrepel::geom_text_repel(data=head(to.plot[r2>0.60 & var_expr>4],n=15), aes(x=r2, y=var_expr, label=gene), size=3, color="darkgreen", max.overlaps=Inf) +
  labs(x="Predictability", y="Variability") +
  # scale_size_manual(values=c("TRUE"=2, "FALSE"=1)) + guides(size=F) +
  theme(
    axis.text = element_text(size=rel(0.7))
  )

###########################
## Analysis of residuals ##
###########################

lm_residuals.dt <- rna_chromvar.dt %>% copy %>%
  .[,c("chromvar_zscore","expr"):=list(chromvar_zscore + rnorm(n=.N,mean=0,sd=1e-5), expr + rnorm(n=.N,mean=0,sd=1e-5))] %>% # add some noise 
  .[,residual:=lm(formula=expr~chromvar_zscore)[["residuals"]], by=c("gene")]

# Extract high positive residuals (higher expression than predicted by chromatin accessibility)
high_lm_residuals.dt <- lm_residuals.dt[residual>2]
# rna_dt <- rna_dt %>% merge(foo, by="id_rna") %>%
#   .[,covariate:=NULL]

genes.highVar.lowPred <- to.plot[r2<0.25 & var_expr>2.5,gene]
foo <- lm_residuals.dt#[gene%in%genes.highVar.lowPred] # %>% .[,abs_residual:=abs(residual)]

order.celltypes <- foo[,median(residual),by="celltype"] %>% setorder(-V1) %>% .$celltype
foo[,celltype:=factor(celltype, levels=order.celltypes)]

ggboxplot(foo, x="celltype", y="residual", fill="celltype", coef=0.1, outlier.shape = NA) +
  scale_fill_manual(values=opts$celltype.colors) +
  geom_hline(yintercept=0, linetype="dashed") +
  guides(x = guide_axis(angle = 90)) +
  coord_cartesian(ylim=c(-5,5)) +
  labs(x="", y="Residuals") +
  theme(
    legend.position = "none",
    # axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text = element_text(size=rel(0.7))
  )


#############
## Explore ##
#############

# genes.to.plot <- c("RARA")
# 
# for (i in genes.to.plot) {
#   
#   to.plot <- rna_chromvar.dt[gene==i] %>%
#     # .[,dot_size:=minmax.normalisation(abs(expr)*minmax.normalisation(chromvar_zscore))] %>% 
#     .[,dot_size:=minmax.normalisation(abs(expr))] %>% 
#     setorder(-dot_size)
#   
#   # p <- ggscatter(to.plot, x="chromvar_zscore", y="expr", fill="celltype", size="dot_size", shape=21,
#   p <- ggscatter(to.plot, x="chromvar_zscore", y="expr", fill="celltype", size=4, shape=21,
#                   add="reg.line", add.params = list(color="black", fill="lightgray"), conf.int=TRUE) +
#     stat_cor(method = "pearson", label.x.npc = "middle", label.y.npc = "bottom") +
#     # ggrepel::geom_text_repel(data=head(to.plot,n=5), aes(x=chromvar_zscore, y=expr, label=celltype), size=3) +
#     scale_fill_manual(values=opts$celltype.colors) +
#     # scale_size_continuous(range = c(1,7)) +
#     labs(y="RNA expression", x=sprintf("Motif accessibility (chromRNA z-score)")) +
#     guides(fill=F, size=F) +
#     theme(
#       axis.text = element_text(size=rel(0.7)),
#       axis.title = element_text(size=rel(0.8))
#     )
#   
#   # pdf(sprintf("%s/test/%s_rna_vs_chromvar_correlated_peaks_pseudobulk.pdf",io$outdir,i), width = 6.5, height = 4)
#   pdf(sprintf("%s/test/%s_scatterplot_rna_vs_chromRNA_pseudobulk.pdf",io$outdir,i), width = 5, height = 5)
#   print(p)
#   dev.off()
# }
