library(igraph)
library(GGally)
library(network)
library(sna)

#####################
## Define settings ##
#####################

# Load default settings
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
io$outdir <- paste0(io$basedir,"/results/atac/archR/chromvar/pseudobulk/network/per_celltype")


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
  # "ExE_ectoderm",
  "Parietal_endoderm"
)

opts$motif_annotation <- "Motif_cisbp"

###################################
## Load pseudobulk TF activities ##
###################################

io$archR.pseudobulk.deviations.se <- sprintf("%s/results/atac/archR/chromvar/pseudobulk/chromVAR_deviations_summarized_experiment_%s_pseudobulk_correlated_peaks.rds",io$basedir,opts$motif_annotation)
source("/Users/ricard/gastrulation_multiome_10x/rna_atac/rna_vs_chromvar/pseudobulk/load_rna_chromvar_pseudobulk.R")

chromvar.dt <- chromvar.dt[celltype%in%opts$celltypes]

#################################
## Calculate similarity matrix ##
#################################

chromvar.mtx <- chromvar.dt %>% 
  dcast(celltype~gene,value.var="chromvar_zscore") %>%
  matrix.please

chromvar.mtx[chromvar.mtx>100] <- 100

chromvar.cor <- cor(t(chromvar.mtx))

# Parse correlation results
diag(chromvar.cor) <- 0
chromvar.cor[chromvar.cor<0.30] <- 0

# Filter correlation results
# chromvar.cor <- chromvar.cor[rowSums(chromvar.cor>0)>3,rowSums(chromvar.cor>0)>3]
# chromvar.mtx <- chromvar.mtx[,rownames(chromvar.cor)]

####################
## Create network ##
####################

network <- network(chromvar.cor, mode="undirected")
# network_unweighted <- graph_from_adjacency_matrix(chromvar.cor, mode="undirected", weighted = TRUE, diag = TRUE)

# celltype.max <- rownames(chromvar.mtx)[apply(chromvar.mtx,2,which.max)]
# col <- opts$celltype.colors[celltype.max]

##############################################
## Plot network, colour by cell type labels ##
##############################################

set.seed(42)
p <- ggnet2(
  net = network,
  color = opts$celltype.colors[rownames(chromvar.cor)],
  mode = "fruchtermanreingold",
  node.size = 7,
  edge.size = 0.01,
  edge.color = "grey",
  label = TRUE,
  label.size = 3.1
)

pdf(paste0(io$outdir,"/TF_network_coloured_by_celltype_labelled.pdf"), width=10, height=8)
print(p)
dev.off()

#########################
## Community detection ##
#########################

foo <- cluster_louvain(network_unweighted)
foo$membership
# foo <- cluster_walktrap(network_unweighted)

set.seed(42)
p <- ggnet2(
  net = network,
  color = foo$membership,
  mode = "fruchtermanreingold",
  node.size = 3,
  edge.size = 0.01,
  edge.color = "grey",
  label = TRUE,
  label.size = 2.1
)

pdf(paste0(io$outdir,"/graph/TF_network_coloured_by_celltype.pdf"), width=10, height=8)
print(p)
dev.off()

#######################################
## Plot number of connections per TF ##
#######################################

to.plot <- data.table(
  gene = rownames(chromvar.cor),
  nedges = rowSums(chromvar.cor>0)
) %>% setorder(-nedges) %>% head(n=80) %>%
  .[,gene:=factor(gene,levels=rev(gene))]

p <- ggplot(to.plot, aes_string(x="gene", y="nedges"), fill="gray70") +
  geom_point(size=2) +
  geom_segment(aes_string(xend="gene"), size=0.5, yend=0) +
  coord_flip() +
  labs(x="", "Number of edges") +
  theme_classic() +
  theme(
    axis.ticks.y = element_blank(),
    axis.text = element_text(size=rel(0.75), color="black")
  )

pdf(paste0(io$outdir,"/nedges_per_TF.pdf"), width = 5, height = 11)
print(p)
dev.off()

#################################################
## For each TF, plot number of correlating TFs ##
#################################################

genes.to.plot <- rownames(chromvar.cor)

for (i in genes.to.plot) {
  
  to.plot <- data.table(
    gene = rownames(chromvar.cor),
    cor = chromvar.cor[i,]
  ) %>% setorder(-cor) %>% head(n=50) %>%
    .[,gene:=factor(gene,levels=rev(gene))]
  
  p <- ggplot(to.plot, aes_string(x="gene", y="cor"), fill="gray70") +
    geom_point(size=2) +
    geom_segment(aes_string(xend="gene"), size=0.5, yend=0) +
    coord_flip(ylim=c(0,1)) +
    labs(x="",y="Correlation coefficient") +
    theme_classic() +
    theme(
      axis.ticks.y = element_blank(),
      axis.text = element_text(size=rel(0.75), color="black")
    )
  
  pdf(sprintf("%s/correlations/%s_top_correlating_genes.pdf",io$outdir,i), width = 5, height = 9)
  print(p)
  dev.off()
}


##########
## Test ##
##########

to.plot <- data.table(
  chromvar.mtx["Erythroid1",],
  chromvar.mtx["Parietal_endoderm",],
  chromvar.mtx["ExE_ectoderm",],
  chromvar.mtx["Mixed_mesoderm",],
  gene = colnames(chromvar.mtx)
)

ggscatter(to.plot, x="V1", y="V3", size=1.5, add="reg.line", add.params = list(color="black", fill="lightgray"), conf.int=TRUE) +
  stat_cor(method = "pearson", label.x.npc = "middle", label.y.npc = "bottom") +
  ggrepel::geom_text_repel(data=to.plot[V1>0.40 & V3>0.40], aes(x=V1, y=V2, label=gene), size=3,  max.overlaps=100) +
  labs(x="TF activity in erythroid", y="TF activity in ExE ectoderm") +
  theme(
    axis.text = element_text(size=rel(0.7))
  )


