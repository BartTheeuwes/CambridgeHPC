#####################
## Define settings ##
#####################

# load default setings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/utils.R")
} else {
  stop("Computer not recognised")
}

# I/O
io$rna.metacells.sce <- paste0(io$basedir,"/results/rna/metacells/SingleCellExperiment_1000metacells.rds")
io$outdir <- paste0(io$basedir,"/results/rna/metacells/pdf"); dir.create(io$outdir, showWarnings = F)

# Options
# opts$number_metacells <- 1000

#######################
## Load marker genes ##
#######################

marker_genes.dt <- fread(io$rna.atlas.marker_genes)

#################################################
## Load RNA expression at the pseudobulk level ##
#################################################

rna_pseudobulk.sce <- readRDS(io$rna.pseudobulk.sce)[,]
rna_pseudobulk.sce <- rna_pseudobulk.sce[rownames(rna_pseudobulk.sce) %in% unique(marker_genes.dt$gene)]

# Create data.table
rna_pseudobulk.dt <- logcounts(rna_pseudobulk.sce) %>% as.data.table(keep.rownames = T) %>%
  melt(id.vars="rn", variable.name="celltype", value.name="expr") %>% setnames("rn","gene") #%>%
  # .[,dataset:=as.factor("pseudobulk")]


###############################################
## Load RNA expression at the metacell level ##
###############################################

rna_metacells.sce <- readRDS(io$rna.metacells.sce)
rna_metacells.sce <- rna_metacells.sce[rownames(rna_metacells.sce) %in% unique(marker_genes.dt$gene)]


rna_metacell.dt <- logcounts(rna_metacells.sce) %>% as.data.table(keep.rownames = T) %>%
  melt(id.vars="rn", variable.name="metacell", value.name="expr") %>% setnames("rn","gene") %>%
  merge(colData(rna_metacells.sce)[,c("celltype.mapped"),drop=F] %>% as.data.table(keep.rownames = T) %>% setnames(c("metacell","celltype")), by="metacell")# %>%
  # .[,dataset:=as.factor("pseudobulk")] 

###########
## Merge ##
###########

# expr.dt <- merge(
#   rna_metacell.dt[,c("celltype","gene","expr")],
#   rna_pseudobulk.dt[,c("celltype","gene","expr")] %>% merge(),
#   by = c("celltype","gene"),
#   suffixes=c(".atlas",".query")
# )
# 
# length(unique(expr.dt$gene))


##########
## Plot ##
##########

# i <- "Hoxa10"
genes.to.plot <- rownames(rna_metacells.sce) %>% head(n=100)

for (i in genes.to.plot)  {
  
  to.plot <- merge(
    rna_pseudobulk.dt[gene==i],
    rna_metacell.dt[gene==i],
    by=c("gene","celltype"),
    suffixes=c(".pseudobulk",".metacell")
    # allow.cartesian=TRUE
  )
  
  
  p <- ggplot(to.plot, aes(x=expr.pseudobulk, y=expr.metacell)) +
    geom_point(aes(fill=celltype), color="black", shape=21, size=3, stroke=0.5) +
    stat_cor(method = "pearson") +
    # geom_abline(slope=1, intercept=0, linetype="dashed", color="black") +
    stat_smooth(method="lm", color="black", alpha=0.1, size=0.2) +
    scale_fill_manual(values=opts$celltype.colors) +
    labs(x="RNA expression (pseudobulk)", y="RNA expression (metacell)", title=i) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(color="black"),
      legend.position = "none"
    )
  
  pdf(sprintf("%s/individual_genes/%s_scatterplot_pseudobulk_vs_metacell.pdf",io$outdir,i), width = 7.5, height = 5)
  print(p)
  dev.off()
}

#####################################
## Calculate correlations per gene ##
#####################################

cor.dt <- expr.dt %>% 
  .[, cor.test(x=expr.atlas, y=expr.query)[c("estimate","p.value")], by = c("gene")] %>%
  # .[, para := rep(c("r","p"), .N/2)] %>% data.table::dcast(gene ~ para, value.var = "V1") %>%
  .[, padj_fdr := p.adjust(p.value, method="fdr")] %>%
  .[, sig := padj_fdr <= 0.1] %>%
  setorder(padj_fdr)

to.plot <- cor.dt

negative_hits <- to.plot[sig==TRUE & estimate<0,gene]
positive_hits <- to.plot[sig==TRUE & estimate>0,gene]
all <- nrow(to.plot)

xlim <- max(abs(to.plot$estimate), na.rm=T)
ylim <- max(-log10(to.plot$p.value), na.rm=T)

p <- ggplot(to.plot, aes(x=estimate, y=-log10(padj_fdr))) +
  labs(title="", x="Pearson correlation", y=expression(paste("-log"[10],"(q.value)"))) +
  # geom_hline(yintercept = -log10(opts$threshold_fdr), color="blue") +
  geom_segment(aes(x=0, xend=0, y=0, yend=ylim-1), color="orange") +
  ggrastr::geom_point_rast(aes(color=sig), size=1) +
  scale_color_manual(values=c("black","red")) +
  scale_x_continuous(limits=c(-xlim-0.05,xlim+0.05)) +
  scale_y_continuous(limits=c(0,ylim+1)) +
  annotate("text", x=0, y=ylim+1, size=5, label=sprintf("(%d)", all)) +
  annotate("text", x=-0.75, y=ylim+1, size=5, label=sprintf("%d (-)",length(negative_hits))) +
  annotate("text", x=0.75, y=ylim+1, size=5, label=sprintf("%d (+)",length(positive_hits))) +
  # ggrepel::geom_text_repel(data=head(to.plot[sig==T],n=top_genes), aes(x=estimate, y=-log10(p.value), label=symbol), size=5) +
  theme_classic() +
  theme(
    # plot.title=element_text(size=28, face='bold', margin=margin(0,0,10,0), hjust=0.5),
    # axis.text=element_text(size=rel(1.75), color='black'),
    # axis.title=element_text(size=rel(1.95), color='black'),
    legend.position="none"
  )

pdf(sprintf("%s/PijuanSala2019_comparison_volcano_plot.pdf",io$outdir), width = 6, height = 6)
print(p)
dev.off()

