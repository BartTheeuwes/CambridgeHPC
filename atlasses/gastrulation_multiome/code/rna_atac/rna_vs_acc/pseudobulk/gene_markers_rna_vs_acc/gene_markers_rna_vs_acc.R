#####################
## Define settings ##
#####################

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))


# I/O
io$basedir <- file.path(io$basedir,"test")
io$housekeeping.genes <-"/Users/argelagr/data/genesets/manual_genesets/housekeeping/housekeeping.tsv"
io$marker_genes_rna <- file.path(io$basedir,"results/rna/differential/pseudobulk/celltype/parsed/marker_genes_filtered.txt.gz")
io$rna.pseudobulk.sce <- file.path(io$basedir,"results/rna/pseudobulk/celltype/SingleCellExperiment_pseudobulk.rds")
io$archR.pseudobulk.GeneMatrix.se <- file.path(io$basedir,"results/atac/archR/pseudobulk/celltype/GeneScoreMatrix_TSS/pseudobulk_GeneScoreMatrix_TSS_summarized_experiment.rds")
io$outdir <- file.path(io$basedir,"results/rna_atac/rna_vs_acc/pseudobulk/gene_markers_rna_vs_acc"); dir.create(io$outdir, showWarnings = F)

# Options
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

####################
## Load gene sets ##
####################

marker_genes.dt <- fread(io$marker_genes_rna) %>% .[celltype%in%opts$celltypes]

houseekeping.dt <- fread(io$housekeeping.genes, header=F) %>% setnames(c("ens_id","gene"))

pluripotency_genes <- c("Dppa2","Dppa3","Dppa4","Morc1","Tex19.1","Zfp981","Dppa5a","Zfp42")

##########################
## Load pseudobulk ATAC ##
##########################

# Load SummarizedExperiment
atac_GeneScores.se <- readRDS(io$archR.pseudobulk.GeneMatrix.se)[,opts$celltypes]

# Normalise ATAC data
assay(atac_GeneScores.se,"logcounts") <- log(1e6*(sweep(assay(atac_GeneScores.se),2,colSums(assay(atac_GeneScores.se),na.rm=T),"/"))+1)

# Prepare data.table
atac.dt <- assay(atac_GeneScores.se,"logcounts") %>% as.data.table(keep.rownames = "gene") %>%
  melt(id.vars="gene", variable.name="celltype", value.name="acc")

#########################
## Load pseudobulk RNA ##
#########################

# Load SingleCellExperiment
rna_pseudobulk.sce <- readRDS(io$rna.pseudobulk.sce)[,opts$celltypes]

# Prepare data.table
rna.dt <- logcounts(rna_pseudobulk.sce) %>%
  as.data.table(keep.rownames="gene") %>% 
  melt(id.vars="gene", variable.name="celltype", value.name="expr")

###########
## Merge ##
###########

# Sanity checks
foo <- unique(atac.dt$celltype)
bar <- unique(rna.dt$celltype)
length(intersect(foo,bar))


foo <- unique(atac.dt$gene)
bar <- unique(rna.dt$gene)
length(intersect(foo,bar))


# Merge
rna_atac.dt <- merge(rna.dt, atac.dt, by = c("gene","celltype"))
length(unique(rna_atac.dt$gene))

# split gene sets
rna_atac_olfactory.dt <- rna_atac.dt[grep("Olfr",gene)] %>% 
  .[,class:="Negative control (olfactory receptors)"]

rna_atac_housekeeping.dt <- rna_atac.dt[gene%in%houseekeping.dt$gene] %>%
  .[,class:="Positive control (housekeeping genes)"]

rna_atac_pluripotency.dt <- rna_atac.dt[gene%in%pluripotency_genes] %>%
  .[,class:="Pluripotency genes"]

rna_atac_markers.dt <- rna_atac.dt %>%
  .[gene%in%marker_genes.dt[,gene]] %>%
  .[,celltype:=factor(celltype,levels=opts$celltypes)] %>%
  .[,class:="Marker genes"]

############################################################################
## Boxplot of accessibility levels per cell type (using all marker genes) ##
############################################################################

# to.plot <- rbindlist(list(rna_atac_markers.dt,rna_atac_olfactory.dt,rna_atac_housekeeping.dt))# %>% .[acc>=opts$max.atac,acc:=opts$max.atac]
# 
# p <- ggboxplot(to.plot, x="celltype", y="acc", fill="celltype", outlier.shape=NA) +
#   facet_wrap(~class, nrow=2, scales="free_y") +
#   scale_fill_manual(values=opts$celltype.colors) +
#   # coord_cartesian(ylim=c(0,opts$max.atac)) +
#   labs(x="", y="Chromatin accessibility levels") +
#   theme(
#     axis.text.y = element_text(color="black", size=rel(0.75)),
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     legend.position = "none"
#   )
# 
# pdf(file.path(io$outdir,"marker_genes_accessibility_boxplots.pdf"), width = 10, height = 6)
# print(p)
# dev.off()

##################################################################################
## Scatterplot of average RNA expression versus average chromatin accessibility ##
##################################################################################
  
celltypes.to.plot <- c("Gut","Endothelium","Cardiomyocytes","Neural_crest")

rna_atac_markers.dt2 <- rna_atac.dt %>%
  merge(marker_genes.dt[celltype%in%celltypes.to.plot,c("celltype","gene")] %>% setnames("celltype","class"), by="gene") %>%
  .[,class:=sprintf("%s markers",class)]

# to.plot <- rbindlist(list(rna_atac_markers.dt2,rna_atac_pluripotency.dt,rna_atac_olfactory.dt,rna_atac_housekeeping.dt)) %>%
to.plot <- rbindlist(list(rna_atac_markers.dt2,rna_atac_olfactory.dt,rna_atac_housekeeping.dt)) %>%
  # .[celltype%in%celltypes.to.plot] %>%
  .[acc>=opts$max.atac,acc:=opts$max.atac] %>%
  # .[,.(expr=median(expr), acc=median(acc)), by=c("class","celltype")] %>%
  .[,.(expr=mean(expr), acc=mean(acc)), by=c("class","celltype")] %>%
  .[,c("expr","acc"):=list(minmax.normalisation(expr),minmax.normalisation(acc))]
  
p1 <- ggplot(to.plot[class%in%c("Positive control (housekeeping genes)","Negative control (olfactory receptors)")], 
             aes_string(x="acc", y="expr", fill="celltype", shape="class")) +
  geom_jitter(size=4, width=0.05, height=0.05) +
# p1 <- ggscatter(to.plot[class%in%c("Positive control (housekeeping genes)","Negative control (olfactory receptors)","Pluripotency genes")], "acc", y="expr", color="celltype", shape="class") +
  geom_hline(yintercept=0.5, linetype="dashed") +
  geom_vline(xintercept=0.5, linetype="dashed") +
  # facet_wrap(~class, nrow=1, scales="fixed") +
  scale_shape_manual(values=c(22,24)) +
  scale_fill_manual(values=opts$celltype.colors) + guides(fill="none") +
  labs(x="Chromatin accessibility (scaled)", y="RNA expression (scaled)") +
  coord_cartesian(ylim=c(0,1), xlim=c(0,1)) +
  scale_x_continuous(breaks=c(0,0.5,1)) + scale_y_continuous(breaks=c(0,0.5,1)) +
  guides(shape = guide_legend(override.aes = list(fill="gray70"))) +
  theme_classic() +
  theme(
    axis.text = element_text(color="black", size=rel(0.75)),
    axis.title = element_text(color="black", size=rel(1)),
    legend.position = "none",
    legend.title = element_blank()
  )

# p2 <- ggscatter(to.plot[!class%in%c("Positive control (housekeeping genes)","Negative control (olfactory receptors)","Pluripotency genes")], 
#                 x="acc", y="expr", fill="celltype", shape=21) +
p2 <- ggplot(to.plot[!class%in%c("Positive control (housekeeping genes)","Negative control (olfactory receptors)")], 
             aes_string(x="acc", y="expr", fill="celltype")) +
  geom_jitter(shape=21, size=3.5, width=0.025, height=0.025) +
  geom_hline(yintercept=0.5, linetype="dashed") +
  geom_vline(xintercept=0.5, linetype="dashed") +
  scale_fill_manual(values=opts$celltype.colors) +
  facet_wrap(~class, nrow=2, scales="fixed") +
  labs(x="Chromatin accessibility (scaled)", y="RNA expression (scaled)") +
  coord_cartesian(ylim=c(0,1), xlim=c(0,1)) +
  scale_x_continuous(breaks=c(0,0.5,1)) + scale_y_continuous(breaks=c(0,0.5,1)) +
  theme_classic() +
  theme(
    strip.text = element_text(size=rel(1)),
    axis.text = element_text(color="black", size=rel(0.75)),
    axis.title = element_text(color="black", size=rel(1)),
    legend.position = "none"
    # legend.direction = "none"
  )

p <- cowplot::plot_grid(plotlist=list(p1,p2), nrow=1, rel_widths=c(2/5,3/5))
pdf(file.path(io$outdir,"marker_genes_rna_vs_accessibility.pdf"), width = 7, height = 4)
print(p)
dev.off()

################
## PAGA graph ##
################

# source(here::here("load_paga_graph.R"))
# 
# # celltypes.to.plot <- c("Gut","Forebrain_Midbrain_Hindbrain","Pharyngeal_mesoderm","Cardiomyocytes","NMP","Spinal_cord")
# 
# to.plot <- rbindlist(list(rna_atac_markers.dt,rna_atac_pluripotency.dt,rna_atac_olfactory.dt,rna_atac_housekeeping.dt)) %>%
#   .[acc>=opts$max.atac,acc:=opts$max.atac] %>%
#   .[,.(expr=mean(expr), acc=mean(acc)), by=c("class","celltype")] %>%
#   .[,c("expr","acc"):=list(minmax.normalisation(expr),minmax.normalisation(acc))]
# 
# # Define colors
# rna.col.seq <- round(seq(0,1,0.1), 2)
# atac.col.seq <- round(seq(0,1,0.1), 2)
# rna.colors <- colorRampPalette(c("gray92", "darkgreen"))(length(rna.col.seq))
# atac.colors <- colorRampPalette(c("gray92", "darkblue"))(length(atac.col.seq)) 
# 
# 
# p.net <- ggnet2(
#   net = net.paga,
#   mode = c("x", "y"),
#   node.size = 0,
#   edge.size = 0.15,
#   edge.color = "grey",
#   label = FALSE,
#   label.size = 2.3
# )
# 
# classes.to.plot <- unique(to.plot$class)#[[1]]
# 
# for (i in classes.to.plot) {
#   
#   expr.values <- to.plot[class==i,c("celltype","expr")] %>% matrix.please %>% .[opts$celltypes,]
#   expr.colors <- round(expr.values,1) %>% map(~ rna.colors[which(rna.col.seq == .)]) %>% unlist
#   
#   # rna
#   p.rna <- p.net + geom_text(label = "\u25D0", aes(x=x, y=y), color=expr.colors, size=20, family = "Arial Unicode MS",
#                       data = p.net$data[,c("x","y")] %>% dplyr::mutate(expr=expr.colors)) +
#     scale_colour_manual(values=expr.colors) + 
#     labs(title=sprintf("RNA expression\n(%s)",i)) +
#     theme(
#       plot.title = element_text(hjust = 0.5)
#     )
#   
#   # accessibility
#   acc.values <- to.plot[class==i,c("celltype","acc")] %>% matrix.please %>% .[opts$celltypes,]
#   acc.colors <- round(acc.values,1) %>% map(~ atac.colors[which(atac.col.seq == .)]) %>% unlist
#   
#   p.acc <- p.net + geom_text(label = "\u25D1", aes(x=x, y=y), color=acc.colors, size=20, family = "Arial Unicode MS",
#                       data = p.net$data[,c("x","y")] %>% dplyr::mutate(acc=acc.colors)) +
#     scale_fill_manual(values=acc.colors) + 
#     labs(title=sprintf("Chromatin accessibility\n(%s)",i)) +
#     theme(
#       plot.title = element_text(hjust = 0.5)
#     )
# 
#   p <- cowplot::plot_grid(plotlist=list(p.rna,p.acc), nrow=1, scale = 0.9) + 
#     theme(
#       plot.title = element_text(hjust = 0.5, size=rel(1.25))
#     )
#   
#   outfile <- sprintf("%s/paga_%s_marker_genes_rna_vs_accessibility.png",io$outdir,str_replace_all(i," ","-"))  
#   png(outfile, width = 550, height = 400)
#   # pdf(outfile, width = 6, height = 4)
#   print(p)
#   dev.off()
# }


###########################
## Plot individual genes ##
###########################

opts$max.expr <- 12
opts$max.atac <- 5.5

tmp <- rbindlist(list(rna_atac_markers.dt,rna_atac_pluripotency.dt,rna_atac_olfactory.dt,rna_atac_housekeeping.dt)) %>%
  .[acc>=opts$max.atac,acc:=opts$max.atac] %>%
  .[expr>=opts$max.expr,expr:=opts$max.expr] %>%
  .[,c("expr_scaled","acc_scaled"):=list(minmax.normalisation(expr),minmax.normalisation(acc))]

genes.to.plot <- c("Actb","Hsp90ab1","Tex19.1","Dppa3","Olfr1416","Olfr1415","Ryr2","Pax3","Tbx20","Epha5","Gypa","Smim1")

to.plot <- tmp %>% .[gene%in%genes.to.plot] %>% .[,gene:=factor(gene,levels=genes.to.plot)]

p <- ggscatter(to.plot, x="acc_scaled", y="expr_scaled", fill="celltype", shape=21, size=3) +
  geom_hline(yintercept=0.5, linetype="dashed") +
  geom_vline(xintercept=0.5, linetype="dashed") +
  facet_wrap(~gene, nrow=2, scales="fixed") +
  coord_cartesian(ylim=c(0,1), xlim=c(0,1)) +
  scale_fill_manual(values=opts$celltype.colors) +
  scale_x_continuous(breaks=c(0,0.5,1)) + scale_y_continuous(breaks=c(0,0.5,1)) +
  labs(x="Chromatin accessibility (scaled)", y="RNA expression (scaled)") +
  theme(
    # plot.title = element_text(hjust = 0.75),
    strip.text = element_text(size=rel(1)),
    strip.background = element_rect(fill = NA, color = "black"),
    axis.text = element_text(color="black", size=rel(0.65)),
    axis.title = element_text(color="black", size=rel(0.9)),
    legend.position = "none",
    panel.spacing = unit(1, "lines")
  )
      
# pdf(file.path(io$outdir,sprintf("%s_rna_vs_accessibility.pdf",i)), width = 5, height = 4)
pdf(file.path(io$outdir,"rna_vs_accessibility_selected_genes.pdf"), width = 10, height = 6)
print(p)
dev.off()
# }


#######################
## Celltype distance ##
#######################

celltypes.to.plot <- c("Epiblast")

rna_atac_subset_markers.dt <- rna_atac.dt %>%
  merge(marker_genes.dt[celltype%in%celltypes.to.plot,c("celltype","gene")] %>% setnames("celltype","class"), by="gene") %>%
  .[,class:=sprintf("%s markers",class)]

to.plot <- rbindlist(list(rna_atac_subset_markers.dt,rna_atac_pluripotency.dt,rna_atac_olfactory.dt,rna_atac_housekeeping.dt)) %>%
  # .[class=="Erythroid3 markers"] %>%
  # .[acc>=opts$max.atac,acc:=opts$max.atac] %>%
  # .[,.(expr=median(expr), acc=median(acc)), by=c("class","celltype")] %>%
  .[,.(expr=mean(expr), acc=mean(acc)), by=c("class","celltype")] %>%
  .[,c("expr","acc"):=list(minmax.normalisation(expr),minmax.normalisation(acc))]

# ggscatter(to.plot[class=="Epiblast markers"], x="acc", y="expr", fill="celltype", shape=21) +
#   geom_hline(yintercept=0.5, linetype="dashed") +
#   geom_vline(xintercept=0.5, linetype="dashed") +
#   facet_wrap(~class, nrow=3, scales="fixed") +
#   scale_fill_manual(values=opts$celltype.colors) +
#   labs(x="Chromatin accessibility (scaled)", y="RNA expression (scaled)") +
#   # coord_cartesian(ylim=c(0,1), xlim=c(0,1)) +
#   # scale_x_continuous(breaks=c(0,0.5,1)) + scale_y_continuous(breaks=c(0,0.5,1)) +
#   theme(
#     strip.text = element_text(size=rel(0.85)),
#     axis.text = element_text(color="black", size=rel(0.50)),
#     axis.title = element_text(color="black", size=rel(0.80)),
#     legend.position = "none"
#   )

# Plot distance in PCA space versus expression
reducedDim(rna_pseudobulk.sce, "PCA") <- irlba::prcomp_irlba(t(logcounts(rna_pseudobulk.sce)), n=10)$x
# saveRDS(reducedDim(rna_pseudobulk.sce, "PCA"), "/Users/ricard/data/gastrulation_multiome_10x/results/rna/pseudobulk/rna_pca_pseudobulk.rds")

celltype_distance.dt <- dist(reducedDim(rna_pseudobulk.sce, "PCA")) %>% as.matrix %>% as.data.table(keep.rownames = T) %>%
  setnames("rn","celltypeA") %>%
  melt(id.vars="celltypeA", variable.name="celltypeB", value.name="distance")

to.plot2 <- to.plot[class=="Epiblast markers"] %>% 
  merge(
    celltype_distance.dt[celltypeA=="Epiblast"] %>% .[,celltypeA:=NULL] %>% setnames("celltypeB","celltype"), by="celltype"
  ) %>% melt(id.vars=c("celltype","class","distance"), variable.name="modality")

facet.labels <- c(expr = "Epiblast markers (RNA expr)", acc = "Epiblast markers (Gene acc.)")

p <- ggscatter(to.plot2, x="distance", y="value", fill="celltype", size=4, shape=21,
          add="reg.line", add.params = list(color="gray50", fill="lightgray"), conf.int=TRUE) +
  stat_cor(method = "pearson") +
  facet_wrap(~modality, nrow=1, scales="free",  labeller = as_labeller(facet.labels)) +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(x="Distance between celltypes", y="Value") +
  theme_classic() +
  theme(
    strip.text = element_text(size=rel(0.85)),
    axis.text = element_text(color="black", size=rel(0.60)),
    axis.title = element_text(color="black", size=rel(0.80)),
    legend.position = "none"
  )


pdf(sprintf("%s/epiblast_markers_rna_acc_vs_celltype_distance.pdf",io$outdir), width = 5, height = 3.5)
print(p)
dev.off()
