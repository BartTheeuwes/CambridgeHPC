# here::i_am("rna/differential/wt_vs_ko/analysis/barplots.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

# Load utils
# source(here::here("differential/analysis/utils.R"))

##############
## Settings ##
##############

# I/O
io$basedir <- file.path(io$basedir,"test")
io$indir <- file.path(io$basedir,"results/rna/differential/old/wt_vs_ko")
io$outdir <- file.path(io$basedir,"results/rna/differential/old/wt_vs_ko/pdf"); dir.create(io$outdir, showWarnings = F)

# Options
opts$min.cells <- 50

###############
## Load data ##
###############

source(here::here("rna/differential/wt_vs_ko/analysis/load_data.R"))

####################
## Filter results ##
####################

# Filter
diff.dt <- diff.dt[groupA_N>=opts$min.cells & groupB_N>=opts$min.cells]

# Remove sex-specific differences
# dt.filt <- dt[gene!="Xist"]

# Filter out genes manually
# dt.filt <- dt.filt[!grep("[^Rik|^Gm|^Rpl]",gene)]
diff.dt <- diff.dt[!grep("^Hb",gene)]

# Subset to lineage markers
marker_genes.dt <- fread(io$rna.atlas.marker_genes.up)
diff_markers.dt <- diff.dt[gene%in%unique(marker_genes.dt$gene)]

################
## Polar plot ##
################

to.plot <- diff_markers.dt %>% copy %>%
  .[,.(N=sum(sig,na.rm=T)) ,by=c("celltype","sign")]

p <- ggplot(to.plot, aes(x=celltype, y=N)) +
  geom_bar(aes(fill = celltype), color="black", stat = 'identity') + 
  facet_wrap(~sign, nrow=1) +
  scale_fill_manual(values=opts$celltype.colors, drop=F) +
  coord_polar() +
  # guides(colour = guide_legend(override.aes = list(size=2), ncol=1)) +
  theme_bw() +
  theme(
    legend.position = "none",
    legend.text = element_text(size=rel(0.75)),
    legend.title = element_blank(),
    axis.title=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks=element_blank(),
    axis.line=element_blank(),
    axis.text.x = element_blank()
    # axis.text.x = element_text(angle= -76 - 360 / length(unique(to.plot$celltype)) * seq_along(to.plot$celltype))
  )

pdf(file.path(io$outdir,"DE_polar_plots_marker_genes.pdf"), width=11, height=8)
print(p)
dev.off()

##############
## Bar plot ##
##############

to.plot <- diff.dt[sig==T] %>% .[,.N, by=c("celltype","sign")]
to.plot <- diff_markers.dt[sig==T] %>% .[,.N, by=c("celltype","sign")]

p <- ggplot(to.plot, aes(x=celltype, y=N, fill=sign)) +
  geom_bar(color="black", stat = 'identity') + 
  # scale_fill_manual(values=opts$celltype.colors, drop=F) +
  labs(x="", y="Number of DE genes") +
  guides(x = guide_axis(angle = 90)) +
  theme_classic() +
  theme(
    legend.position = "none",
    # axis.line = element_blank(),
    axis.text.y = element_text(color="black", size=rel(1)),
    axis.text.x = element_text(color="black", size=rel(0.9))
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank()
  )


pdf(file.path(io$outdir,"DE_barplots_marker_genes.pdf"), width=8, height=4)
print(p)
dev.off()

##########################################
## Barplot of fraction of genes up/down ##
##########################################

to.plot <- diff_markers.dt[sig==T] %>%
  .[,direction:=c("Downregulated in KO","Upregulated in KO")[as.numeric(logFC>0)+1]] %>%
  .[,.N, by=c("celltype","direction")]

p <- ggplot(to.plot, aes(x=factor(celltype), y=N)) +
  geom_bar(aes(fill = direction), color="black", stat="identity") + 
  # facet_wrap(~class, scales="free_y") +
  labs(x="", y="Number of DE genes") +
  guides(x = guide_axis(angle = 90)) +
  theme_classic() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    # axis.line = element_blank(),
    axis.text.x = element_text(color="black", size=rel(0.75))
  )

pdf(file.path(io$outdir,"DE_barplots_marker_genes_direction.pdf"), width=8, height=4.5)
print(p)
dev.off()

#############################
## Cell fate bias barplots ##
#############################

to.plot <- diff_markers.dt %>%
  merge(
    marker_genes.dt[,c("gene","celltype")] %>% setnames("celltype","celltype_marker"), by = "gene", allow.cartesian=TRUE
  ) %>% .[,sum(sig), by=c("celltype","celltype_marker","sign")]

p <- ggplot(to.plot, aes(x=celltype, y=V1)) +
  geom_bar(aes(fill = celltype_marker), color="black", stat="identity") + 
  facet_wrap(~sign, scales="fixed", nrow=2) +
  # facet_wrap(~class, scales="fixed") +
  scale_fill_manual(values=opts$celltype.colors[names(opts$celltype.colors)%in%unique(to.plot$celltype_marker)], drop=F) +
  labs(x="", y="Number of DE genes") +
  guides(x = guide_axis(angle = 90)) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.line = element_blank(),
    axis.text.x = element_text(color="black", size=rel(0.75))
  )

pdf(sprintf("%s/DE_barplots_marker_genes_fate_bias.pdf",io$outdir), width=9, height=5)
print(p)
dev.off()

################################
## Cell fate bias polar plots ##
################################

to.plot <- diff_markers.dt %>%
  merge(
    marker_genes.dt[,c("gene","celltype")] %>% setnames("celltype","celltype_marker"), by = "gene", allow.cartesian=TRUE
  ) %>% .[,sum(sig), by=c("celltype","celltype_marker","sign")]

to.plot[V1>=7,V1:=7]

for (i in unique(to.plot$celltype)) {
  
  p <- ggplot(to.plot[celltype==i], aes(x=celltype_marker, y=V1)) +
    geom_bar(aes(fill = celltype_marker), color="black", stat = 'identity') + 
    facet_wrap(~sign, nrow=1) +
    scale_fill_manual(values=opts$celltype.colors, drop=F) +
    coord_polar()
    theme_bw() +
    theme(
      legend.position = "none",
      legend.text = element_text(size=rel(0.75)),
      legend.title = element_blank(),
      axis.title=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      axis.line=element_blank(),
      axis.text.x = element_blank()
      # axis.text.x = element_text(angle= -76 - 360 / length(unique(to.plot$celltype)) * seq_along(to.plot$celltype))
    )
  
  pdf(file.path(io$outdir,sprintf("%s_DE_polar_plot_fate_bias.pdf",i)), width=6, height=3)
  print(p)
  dev.off()
}

##########
## TEST ##
##########

diff_markers.dt[celltype=="NMP" & sig==TRUE & sign=="Downregulated in WT"] %>% View
