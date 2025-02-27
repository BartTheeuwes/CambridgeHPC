#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))


# I/O
io$housekeeping.genes <-"/Users/ricard/data/genesets/manual_genesets/housekeeping/housekeeping.tsv"
io$rna.pseudobulk.sce <- paste0(io$basedir,"/results/rna/pseudobulk/SingleCellExperiment_velocyto.rds")
io$outdir <- paste0(io$basedir,"/results_new/rna/marker_genes/nascent"); dir.create(io$outdir, showWarnings = F)


opts$celltypes <- c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  # "PGC",
  # "Anterior_Primitive_Streak",
  "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  # "Mixed_mesoderm",
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
  "Surface_ectoderm"
  # "Visceral_endoderm"
  # "ExE_endoderm",
  # "ExE_ectoderm"
  # "Parietal_endoderm"
)

#######################
## Load marker genes ##
#######################

# io$marker_genes_rna <- io$rna.atlas.marker_genes
io$marker_genes_rna <- file.path(io$basedir,"results_new/rna/differential/marker_genes/marker_genes_up.txt.gz")
marker_genes.dt <- fread(io$rna.atlas.marker_genes) %>%
  .[celltype%in%opts$celltypes]

# marker_genes.dt <- marker_genes.dt[score>=0.90]
# marker_genes.dt <- marker_genes.dt %>% .[grep("Rik",gene,invert = T)]

housekeeping.dt <- fread(io$housekeeping.genes, header=F) %>%
  setnames(c("ens_id","gene"))

pluripotency_genes <- c("Dppa2","Dppa3","Dppa4","Morc1","Tex19.1","Zfp981","Dppa5a","Zfp42")

#########################
## Load pseudobulk RNA ##
#########################

# Load SingleCellExperiment
rna_velocyto_pseudobulk.sce <- readRDS(io$rna.pseudobulk.sce)[,opts$celltypes]
assays(rna_velocyto_pseudobulk.sce)

# Subset features
# rna_velocyto_pseudobulk.sce <- rna_velocyto_pseudobulk.sce[rownames(rna_velocyto_pseudobulk.sce) %in% unique(marker_genes.dt$gene),]
# rna_velocyto_pseudobulk.sce <- rna_velocyto_pseudobulk.sce[!rownames(rna_velocyto_pseudobulk.sce) %in% unique(marker_genes.dt$gene),]
# rna_velocyto_pseudobulk.sce <- rna_velocyto_pseudobulk.sce[grep("Olfr",rownames(rna_velocyto_pseudobulk.sce)),]

# Prepare data.table
rna_unspliced.dt <- assay(rna_velocyto_pseudobulk.sce,"unspliced_log") %>%
  as.data.table(keep.rownames="gene") %>% 
  melt(id.vars="gene", variable.name="celltype", value.name="unspliced_expr")

rna_spliced.dt <- assay(rna_velocyto_pseudobulk.sce,"spliced_log") %>%
  as.data.table(keep.rownames="gene") %>% 
  melt(id.vars="gene", variable.name="celltype", value.name="spliced_expr")

###########
## Merge ##
###########

# Merge
rna.dt <- merge(rna_unspliced.dt, rna_spliced.dt, by = c("gene","celltype"))
length(unique(rna.dt$gene))

####################
## Quantification ##
####################

rna_olfactory.dt <- rna.dt[grep("Olfr",gene)] %>% 
  .[,class:="Negative control (olfactory receptors)"]

rna_housekeeping.dt <- rna.dt[gene%in%housekeeping.dt$gene] %>%
  .[,class:="Positive control (housekeeping genes)"]

rna_pluripotency.dt <- rna.dt[gene%in%pluripotency_genes] %>%
  .[,class:="Pluripotency genes"]

# rna_markers.dt <- rna.dt %>%
#   .[gene%in%unique(marker_genes.dt$gene)] %>%
#   .[,celltype:=factor(celltype,levels=opts$celltypes)] %>%
#   .[,class:=sprintf("%s markers",celltype)]

rna_markers.dt <- rna.dt %>%
  merge(marker_genes.dt[,c("celltype","gene")] %>% .[,class:=sprintf("%s markers",celltype)] %>% .[,celltype:=NULL], by="gene")
  
  
###############################################################
## Scatterplot of unspliced RNA expr versus spliced RNA expr ##
###############################################################
  
celltypes.to.plot <- c("Gut","Neural_crest","Endothelium","Cardiomyocytes","Notochord","Allantois")

to.plot <- rbindlist(list(rna_markers.dt, rna_pluripotency.dt, rna_olfactory.dt, rna_housekeeping.dt)) %>%
  # .[celltype%in%celltypes.to.plot] %>%
  .[,.(unspliced_expr=mean(unspliced_expr), spliced_expr=mean(spliced_expr)), by=c("class","celltype")] %>%
  .[,c("unspliced_expr","spliced_expr"):=list(minmax.normalisation(unspliced_expr),minmax.normalisation(spliced_expr))]
  
p1 <- ggscatter(to.plot[class%in%c("Positive control (housekeeping genes)","Negative control (olfactory receptors)","Pluripotency genes")], 
                x="unspliced_expr", y="spliced_expr", fill="celltype", shape=21) +
  geom_hline(yintercept=0.5, linetype="dashed") +
  geom_vline(xintercept=0.5, linetype="dashed") +
  facet_wrap(~class, nrow=1, scales="fixed") +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(x="Unspliced RNA expression (scaled)", y="Spliced RNA expression (scaled)") +
  coord_cartesian(ylim=c(0,1), xlim=c(0,1)) +
  scale_x_continuous(breaks=c(0,0.5,1)) + scale_y_continuous(breaks=c(0,0.5,1)) +
  theme(
    strip.text = element_text(size=rel(0.65)),
    axis.text = element_text(color="black", size=rel(0.50)),
    axis.title = element_text(color="black", size=rel(0.80)),
    legend.position = "none"
  )

p2 <- ggscatter(to.plot[!class%in%c("Positive control (housekeeping genes)","Negative control (olfactory receptors)","Pluripotency genes")], 
                x="unspliced_expr", y="spliced_expr", fill="celltype", shape=21) +
  geom_hline(yintercept=0.5, linetype="dashed") +
  geom_vline(xintercept=0.5, linetype="dashed") +
  facet_wrap(~class, nrow=3, scales="fixed") +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(x="Unspliced RNA expression (scaled)", y="Spliced RNA expression (scaled)") +
  coord_cartesian(ylim=c(0,1), xlim=c(0,1)) +
  scale_x_continuous(breaks=c(0,0.5,1)) + scale_y_continuous(breaks=c(0,0.5,1)) +
  theme(
    strip.text = element_text(size=rel(0.85)),
    axis.text = element_text(color="black", size=rel(0.50)),
    axis.title = element_text(color="black", size=rel(0.80)),
    legend.position = "none"
  )

p <- cowplot::plot_grid(plotlist=list(p1,p2), nrow=2, rel_heights=c(1/3,2/3))
pdf(sprintf("%s/marker_genes_rna_vs_accessibility.pdf",io$outdir), width = 7, height = 8)
print(p)
dev.off()



#############
## Explore ##
#############

# celltypes.to.plot <- c("Gut","Neural_crest","Endothelium","Cardiomyocytes","Notochord","Allantois")

to.plot <- rna_markers.dt[class=="Epiblast markers"] %>%
  .[,ratio:=unspliced_expr/spliced_expr] %>%
  .[,.(ratio=mean(ratio)), by=c("class","celltype")]

p <- ggbarplot(to.plot, x="celltype", y="ratio", fill="celltype") +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(x="", y="Ratio unspliced/spliced RNA expression") +
  theme(
    axis.text = element_text(color="black", size=rel(0.50)),
    axis.title = element_text(color="black", size=rel(0.80)),
    legend.position = "none"
  )

p <- ggscatter(to.plot, x="unspliced_expr", y="spliced_expr", fill="celltype", shape=21) +
  geom_hline(yintercept=0.5, linetype="dashed") +
  geom_vline(xintercept=0.5, linetype="dashed") +
  facet_wrap(~class, nrow=1, scales="fixed") +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(x="Unspliced RNA expression (scaled)", y="Spliced RNA expression (scaled)") +
  coord_cartesian(ylim=c(0,1), xlim=c(0,1)) +
  scale_x_continuous(breaks=c(0,0.5,1)) + scale_y_continuous(breaks=c(0,0.5,1)) +
  theme(
    strip.text = element_text(size=rel(0.65)),
    axis.text = element_text(color="black", size=rel(0.50)),
    axis.title = element_text(color="black", size=rel(0.80)),
    legend.position = "none"
  )



#############
## Explore ##
#############

io$outdir <- "/Users/ricard/data/gastrulation_multiome_10x/results/rna/velocity/pseudobulk"

genes.to.plot <- unique(marker_genes.dt$gene)

for (i in genes.to.plot) {
  
  to.plot <- rna.dt[gene==i]
  
  p <- ggscatter(to.plot, x="spliced_expr", y="unspliced_expr", fill="celltype", shape=21, size=5) +
    # geom_abline(slope=1, intercept=0) +
    scale_fill_manual(values=opts$celltype.colors) +
    labs(x="Spliced RNA expression (scaled)", y="Unspliced RNA expression (scaled)") +
    theme(
      axis.text = element_text(color="black", size=rel(0.75)),
      axis.title = element_text(color="black", size=rel(0.80)),
      legend.position = "none"
    )
  
  pdf(sprintf("%s/%s_spliced_vs_unspliced_pseudobulk.pdf",io$outdir,i), width = 6, height = 5)
  print(p)
  dev.off()
  
}


#############
## Explore ##
#############

to.plot <- rna.dt %>%
  .[gene%in%unique(marker_genes.dt$gene)] %>%
  .[,.(unspliced_expr=mean(unspliced_expr), spliced_expr=mean(spliced_expr)),by="celltype"] %>%
  merge(sample_metadata[,.N,by=c("celltype.mapped")] %>% setnames("celltype.mapped","celltype"), by="celltype") %>%
  .[,log_N:=log(N)]
  # melt(id.vars=c("celltype","N"))

p1 <- ggscatter(to.plot, x="spliced_expr", y="unspliced_expr", fill="celltype", shape=21, size=5) +
                # add="reg.line", add.params = list(color="black", fill="lightgray"), conf.int=TRUE) +
  scale_fill_manual(values=opts$celltype.colors) +

  labs(x="Spliced RNA expression (scaled)", y="Unspliced RNA expression (scaled)") +
  theme(
    axis.text = element_text(color="black", size=rel(0.75)),
    axis.title = element_text(color="black", size=rel(0.80)),
    legend.position = "none"
  )


p2 <- ggscatter(to.plot, x="spliced_expr", y="unspliced_expr", fill="log_N", shape=21, size=5) +
  # geom_abline(slope=1, intercept=0) +
  scale_fill_gradientn(colours = terrain.colors(10)) +
  labs(x="Spliced RNA expression (scaled)", y="Unspliced RNA expression (scaled)") +
  theme(
    axis.text = element_text(color="black", size=rel(0.75)),
    axis.title = element_text(color="black", size=rel(0.80)),
    legend.position = "right"
  )



p <- cowplot::plot_grid(plotlist=list(p1,p2), nrow=1, rel_widths=c(1/2,1/2))
pdf(sprintf("%s/rna_content_per_celltype_pseudobulk.pdf",io$outdir), width = 10, height = 4)
print(p)
dev.off()

