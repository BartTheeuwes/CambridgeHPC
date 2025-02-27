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


#####################
## Define settings ##
#####################

# Options
opts$stage <- c(
  # "E7.5"
  # "E8.0"
  "E8.5"
)

opts$max.cells <- 9000

# opts$motif_annotation <- "Motif_JASPAR2020_human"
opts$motif_annotation <- "Motif_cisbp"

# Denoising options
opts$denoise <- TRUE
opts$knn <- 25

# I/O
io$outdir <- sprintf("%s/results/rna_atac/rna_vs_chromvar/all_cells_per_stage/%s",io$basedir,opts$stage); dir.create(io$outdir, showWarnings = F)

###################
## Load metadata ##
###################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & pass_rnaQC==TRUE & doublet_call==FALSE] %>%
  .[stage%in%opts$stage] %>%
  .[,celltype.predicted:=factor(celltype.predicted,levels=opts$celltypes)] 

# Downsample number of cells
if (nrow(sample_metadata)>opts$max.cells) {
  set.seed(42)
  sample_metadata <- sample_metadata[sample(1:nrow(sample_metadata),opts$max.cells)]
}

table(sample_metadata$sample)

################################
## Load RNA and chromVAR data ##
################################

# I/O
io$pca.rna <- sprintf("%s/results/rna/dimensionality_reduction/%s/%s_pca_features2500_pcs30.txt.gz",io$basedir,opts$stage,paste(sort(unique(sample_metadata$sample)), collapse="-"))
io$pca.atac <- sprintf("%s/results/atac/archR/dimensionality_reduction/PeakMatrix/%s/%s_lsi_features50000_ndims30.txt.gz",io$basedir,opts$stage,paste(sort(unique(sample_metadata$sample)), collapse="-"))

args <- opts
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/rna_atac/rna_vs_chromvar/load_rna_chromvar_single_cells.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/rna_atac/rna_vs_chromvar/load_rna_chromvar_single_cells.R")
} else {
  stop("Computer not recognised")
}

# Merge
chromvar_rna_dt <- merge(
  rna_dt,
  chromvar_dt %>% 
    merge(sample_metadata[,c("cell","celltype.predicted")], by="cell"),
  by = c("cell","gene")
)

length(unique(chromvar_rna_dt$gene))
length(unique(chromvar_rna_dt$cell))

# fwrite(chromvar_rna_dt, sprintf("%s/chromvar_rna.txt.gz",args$outdir), quote=F, sep="\t", na="NA")
# chromvar_rna_dt <- fread(sprintf("%s/chromvar_rna_blood.txt.gz",args$outdir))

##########################
## Correlation analysis ##
##########################

opts$threshold_fdr <- 0.10

cor.dt <- chromvar_rna_dt %>% copy %>%
  .[,c("chromvar_zscore","expr"):=list(chromvar_zscore + rnorm(n=.N,mean=0,sd=1e-5), expr + rnorm(n=.N,mean=0,sd=1e-5))] %>% # add some noise 
  .[, .(V1 = unlist(cor.test(chromvar_zscore, expr)[c("estimate", "p.value")])), by = c("gene")] %>%
  .[, para := rep(c("r","p"), .N/2)] %>% 
  data.table::dcast(gene ~ para, value.var = "V1") %>%
  .[,"padj_fdr" := list(p.adjust(p, method="fdr"))] %>%
  # .[,"log_padj_fdr" := list(-log10(padj_fdr))] %>%
  .[, sig := padj_fdr <= opts$threshold_fdr] %>% 
  setorder(padj_fdr, na.last = T)

head(cor.dt,n=3)

# Save
# fwrite(cor.dt, sprintf("%s/cor_rna_vs_chromvar_%s_%s.txt.gz",io$outdir,opts$stage,opts$knn), sep="\t", quote=F)
fwrite(cor.dt, sprintf("%s/cor_rna_vs_chromvar.txt.gz",io$outdir), sep="\t", quote=F)

##################
## Volcano plot ##
##################

to.plot <- cor.dt

negative_hits <- to.plot[sig==TRUE & r<0,gene]
positive_hits <- to.plot[sig==TRUE & r>0,gene]
all <- nrow(to.plot)

xlim <- max(abs(to.plot$r), na.rm=T)
ylim <- max(-log10(to.plot$padj_fdr+1e-100), na.rm=T)

p <- ggplot(to.plot, aes(x=r, y=-log10(padj_fdr+1e-100))) +
  geom_segment(aes(x=0, xend=0, y=0, yend=ylim-1), color="orange", size=0.5) +
  geom_point(aes(color=sig, size=sig)) +
  # ggrastr::geom_point_rast(aes(color=sig, size=sig)) +
  ggrepel::geom_text_repel(data=head(to.plot[sig==T],n=50), aes(x=r, y=-log10(padj_fdr+1e-100), label=gene), size=3,  max.overlaps=100) +
  scale_color_manual(values=c("black","red")) +
  scale_size_manual(values=c(0.75,1.25)) +
  scale_x_continuous(limits=c(-xlim-0.5,xlim+0.5)) +
  scale_y_continuous(limits=c(0,ylim+6)) +
  annotate("text", x=0, y=ylim+6, size=4, label=sprintf("(%d)", all)) +
  annotate("text", x=-xlim-0.5, y=ylim+6, size=4, label=sprintf("%d (-)",length(negative_hits))) +
  annotate("text", x=xlim+0.5, y=ylim+6, size=4, label=sprintf("%d (+)",length(positive_hits))) +
  labs(x="Pearson correlation", y=expression(paste("-log"[10],"(p.value)"))) +
  theme_classic() +
  theme(
    axis.text = element_text(size=rel(0.75), color='black'),
    axis.title = element_text(size=rel(1.0), color='black'),
    legend.position="none"
  )

pdf(sprintf("%s/volcano_pearson_correlation.pdf",io$outdir), width = 9, height = 6)
# png(sprintf("%s/volcano_plots/volcano_pearson_correlation.png",io$outdir), width = 800, height = 500)
print(p)
dev.off()


#####################################
## Scatter plot of individual hits ##
#####################################

# # genes.to.plot <- unique(chromvar_rna_dt$gene)
# # genes.to.plot <- cor.dt[sig==T & abs(r)>0.25,gene]
# 
# for (i in genes.to.plot) {
#   
#   outfile <- sprintf("%s/individual_genes/%s_rna_vs_chromvar_knn%s_allcells.png",io$outdir,i,opts$knn)
#   
#   if (file.exists(outfile)) {
#     print(sprintf("file or %s already exists...",i))
#   } else {
#     
#     to.plot <- chromvar_rna_dt[gene==i]
#     
#     p1 <- ggscatter(to.plot, x="chromvar_zscore", y="expr", color="celltype.predicted", size=1, 
#                     add="reg.line", add.params = list(color="black", fill="lightgray"), conf.int=TRUE) +
#       stat_cor(method = "pearson", label.x.npc = "middle", label.y.npc = "bottom") +
#       scale_colour_manual(values=opts$celltype.colors) +
#       labs(y=sprintf("%s expression",i), x=sprintf("%s Motif accessibility (z-score)",i)) +
#       guides(color=F) +
#       theme(
#         axis.text = element_text(size=rel(0.7))
#       )
#     
#     
#     to.plot2 <- to.plot %>% melt(id.vars=c("cell","celltype.predicted","gene"))
#     p2 <- ggboxplot(to.plot2, x="celltype.predicted", y="value", fill="celltype.predicted", outlier.shape=NA) +
#       facet_wrap(~variable, nrow=2, scales="free_y") +
#       scale_fill_manual(values=opts$celltype.colors) +
#       geom_hline(yintercept=0, linetype="dashed") +
#       labs(x="", y="") +
#       theme_classic() +
#       guides(x = guide_axis(angle = 90)) +
#       theme(
#         legend.position = "none",
#         axis.text.x = element_blank(),
#         axis.title.x = element_blank(),
#         axis.ticks.x = element_blank()
#       )
#     
#     
#     p <- cowplot::plot_grid(plotlist=list(p1,p2), nrow = 1, rel_widths = c(1/2,1/2))
#     
#     # pdf(sprintf("%s/individual_genes/%s_rna_vs_chromvar.pdf",io$outdir,i), width=10, height=4)
#     png(outfile, width = 1000, height = 500)
#     print(p)
#     dev.off()
#   }
# }
