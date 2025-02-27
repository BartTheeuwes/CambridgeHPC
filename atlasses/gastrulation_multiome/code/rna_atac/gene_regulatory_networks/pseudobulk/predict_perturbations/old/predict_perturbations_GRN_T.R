
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
# io$archR.pseudobulk.deviations.se <- sprintf("%s/pseudobulk/pseudobulk_DeviationMatrix_%s_summarized_experiment.rds",io$archR.directory,opts$motif_annotation)
# io$archR.pseudobulk.deviations.se <- sprintf("%s/results/atac/archR/chromvar/pseudobulk/chromVAR_deviations_summarized_experiment_%s_pseudobulk_correlated_peaks.rds",io$basedir,opts$motif_annotation)
# io$pca.rna <- paste0(io$basedir,"/results/rna/dimensionality_reduction/all_cells/E7.5_rep1-E7.5_rep2-E8.0_rep1-E8.0_rep2-E8.5_rep1-E8.5_rep2_pca_features2500_pcs30_batchcorrectionbysample.txt.gz")

# Options
opts$celltypes = c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  # "PGC",
  # "Anterior_Primitive_Streak",
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
  # "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm"
  # "Visceral_endoderm",
  # "ExE_endoderm",
  # "ExE_ectoderm",
  # "Parietal_endoderm"
)

opts$aggregated.celltypes <- c(
  "Erythroid1" = "Erythroid",
  "Erythroid2" = "Erythroid",
  "Erythroid3" = "Erythroid",
  "Blood_progenitors_1" = "Blood_progenitors",
  "Blood_progenitors_2" = "Blood_progenitors",
  "Rostral_neurectoderm" = "Neurectoderm",
  "Caudal_neurectoderm" = "Neurectoderm",
  "Anterior_Primitive_Streak" = "Primitive_Streak"
)

##################################
## Load gene regulatory network ##
##################################

marker_genes.dt <- fread(io$rna.atlas.marker_genes) %>%
  .[celltype%in%opts$celltypes] %>%
  .[,celltype:=stringr::str_replace_all(celltype,opts$aggregated.celltypes)] %>%
  .[,c("celltype","gene")] %>% unique

tf2gene.mtx <- readRDS("/Users/ricard/data/gastrulation_multiome_10x/results/rna_atac/gene_regulatory_networks/pseudobulk/per_TF/T/corr_matrix_stringent.rds")
tf2gene.mtx <- tf2gene.mtx[1,abs(tf2gene.mtx)>0.5]

################################################
## Load pseudobulk RNA and chromVAR estimates ##
################################################

# opts$motif_annotation <- "Motif_cisbp"
# io$archR.pseudobulk.deviations.se <- sprintf("%s/pseudobulk/pseudobulk_DeviationMatrix_%s_summarized_experiment.rds",io$archR.directory,opts$motif_annotation)
# 
# if (grepl("ricard",Sys.info()['nodename'])) {
#   source("/Users/ricard/gastrulation_multiome_10x/rna_atac/load_rna_atac_pseudobulk.R")
# } else if (grepl("ebi",Sys.info()['nodename'])) {
#   source("/homes/ricard/gastrulation_multiome_10x/rna_atac/load_rna_atac_pseudobulk.R")
# } else {
#   stop("Computer not recognised")
# }

###############################
## Load motifmatcher results ##
###############################

# if (grepl("ricard",Sys.info()['nodename'])) {
#   source("/Users/ricard/gastrulation_multiome_10x/load_motifmatchR.R")
# } else if (grepl("ebi",Sys.info()['nodename'])) {
#   source("/homes/ricard/gastrulation_multiome_10x/load_motifmatchR.R")
# } else {
#   stop("Computer not recognised")
# }

#############
## Explore ##
#############

# Calculate average RNA expression across all TF target genes for each celltype
predicted_downregulated.genes <- names(tf2gene.mtx)[tf2gene.mtx>0]
predicted_upregulated.genes <- names(tf2gene.mtx)[tf2gene.mtx<0]
# predictions <- logcounts(rna.sce[predicted_upregulated.genes,]) %>% apply(2, minmax.normalisation) %>% colMeans


## Load DE genes from T ko
diff_rna.dt <- fread("/Users/ricard/data/10x_gastrulation_T_Chimera/results/differential/pseudobulk/DE_pseudobulk.txt.gz")
diff_rna_filt.dt <- diff_rna.dt %>% .[abs(logFC)>=0.5 & padj_fdr<=0.05]
# tmp <- diff_rna.dt[celltype=="Haematoendothelial_progenitors",gene]

# diff_rna.dt$gene %in% predicted_upregulated.genes
# diff_rna.dt$gene %in% predicted_downregulated.genes

############################################
## Plot number of well-predicted DE genes ##
############################################

foo <- data.table(
  class = "T_network",
  sign = c("downregulated","upregulated"),
  value = c(mean(predicted_downregulated.genes%in%diff_rna_filt.dt$gene), mean(predicted_upregulated.genes%in%diff_rna_filt.dt$gene))
)

ggbarplot(foo, x="sign", y="value", fill="sign", stat="identity") +
  labs(x="", y="Fraction of correctly predicted DE genes") +
  coord_cartesian(ylim=c(0,1)) +
  theme(
    axis.text.y = element_text(size=rel(0.75)),
    legend.position = "none"
  )


##########################################
## Plot network with the logFC per gene ##
##########################################

tf2gene_filt.mtx <- tf2gene.mtx[tf2gene.mtx>=0]

# Prepare data
node_list.dt <- data.table(node_id=1, node_name="T")
target_list.dt <- data.table(target_id=1:length(tf2gene_filt.mtx), target_name=names(tf2gene_filt.mtx))

edge_list.dt <- data.table(from="T", to=names(tf2gene_filt.mtx), weight=tf2gene_filt.mtx)

node_list_metadata.dt <- data.table(
  label = c(node_list.dt$node_name, target_list.dt$target_name),
  class = c(rep("TF",nrow(node_list.dt)), rep("gene",nrow(target_list.dt)))
)

net <- graph_from_data_frame(d = edge_list.dt, vertices = node_list_metadata.dt)

# Define groups
V(net)$group <- factor(V(net)$class, levels=c("TF","gene"))

# Define colors
pal <- grDevices::colorRamp(c("gray80", "purple"))( (1:100)/100 )
tmp <- diff_rna.dt[gene%in%c("T",names(tf2gene_filt.mtx))] %>% setkey(gene) %>% .[str_to_title(names(V(net)))]
tmp[,abs_logFC:=abs(logFC)] %>% .[abs_logFC>=8,abs_logFC:=8]
V(net)$color <- colourvalues::colour_values(tmp$abs_logFC, palette = pal)

# Create network layout
set.seed(42)
layout <- layout.sphere(net)

# Plot
ggnet2(
  mode = layout, 
  net = net,
  color = V(net)$color,
  node.size = c(8,20)[factor(V(net)$class)],
  label = TRUE,
  label.color = "black",
  label.size = c(3,5)[factor(V(net)$class)],
  arrow.size = 0.1,
  legend.position = "none"
) 


###################
## Cell fate bias ##
###################

to.plot <- data.table(
  TF = "T",
  gene = names(tf2gene.mtx),
  cor = tf2gene.mtx,
  sign = c("-","+")[as.numeric(tf2gene.mtx>0)+1]
) %>% merge(
  marker_genes.dt[,c("gene","celltype")] %>% setnames("celltype","celltype_marker"),
  by = c("gene"), allow.cartesian=TRUE
) %>%
  .[,.N, by=c("celltype_marker","sign","TF")] %>% 
  .[,value:=minmax.normalisation(N),by=c("TF","sign")]

ggplot(to.plot[sign=="+"], aes(x=factor(TF), y=N)) +
  geom_bar(aes(fill = celltype_marker), color="black", stat="identity") + 
  # facet_wrap(~sign, scales="free_y") +
  scale_fill_manual(values=opts$celltype.colors, drop=F) +
  labs(x="", y="Number of DE genes (scaled)") +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
  )

ggplot(to.plot[sign=="+"], aes(x=factor(celltype_marker), y=N)) +
  geom_bar(aes(fill = celltype_marker), color="black", stat="identity") + 
  facet_wrap(~sign, scales="free_y",nrow=2) +
  scale_fill_manual(values=opts$celltype.colors, drop=F) +
  labs(x="", y="Number of predicted DE genes") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right"
  )
