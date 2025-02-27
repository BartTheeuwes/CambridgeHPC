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
io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/pseudobulk/TF_activities/per_celltype")


# Options
opts$celltypes = c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  # "PGC",
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

###################################
## Load pseudobulk TF activities ##
###################################

tf_activities.dt <- fread(io$tf.activities.pseudobulk) %>%
  .[celltype%in%opts$celltypes]

opts$celltypes <-opts$celltypes[opts$celltypes %in% unique(tf_activities.dt$celltype)]

########################################################
## Correlate TF activities per cell type (across TFs) ##
########################################################

tf_activities.mtx <- tf_activities.dt %>% 
  dcast(celltype~gene,value.var="activity") %>%
  matrix.please %>% t

tf_activities.cor <- cor(tf_activities.mtx)

# Parse correlation results
diag(tf_activities.cor) <- 0
tf_activities.cor[tf_activities.cor<0.25] <- 0

# Filter correlation results
# tf_activities.cor <- tf_activities.cor[rowSums(tf_activities.cor>0)>3,rowSums(tf_activities.cor>0)>3]
# tf_activities.mtx <- tf_activities.mtx[,rownames(tf_activities.cor)]

####################
## Create network ##
####################

network <- network(tf_activities.cor, mode="undirected")
network_unweighted <- graph_from_adjacency_matrix(tf_activities.cor, mode="undirected", weighted = TRUE, diag = TRUE)

##############################################
## Plot network, colour by cell type labels ##
##############################################

set.seed(42)

p <- ggnet2(
  net = network,
  color = opts$celltype.colors[rownames(tf_activities.cor)],
  mode = "fruchtermanreingold",
  node.size = 6,
  edge.size = 0.15,
  edge.color = "grey",
  label = TRUE,
  label.size = 2.3
)

pdf(paste0(io$outdir,"/TF_activitiy_network_coloured_by_celltype.pdf"), width=10, height=8)
print(p)
dev.off()

###################################
## Scatterplots of TF activities ##
###################################

dir.create(paste0(io$outdir,"/scatterplots"), showWarnings = F)

for (i in opts$celltypes) {
  for (j in opts$celltypes) {
    
    to.plot <- data.table(
      V1 = tf_activities.mtx[,i],
      V2 = tf_activities.mtx[,j],
      gene = rownames(tf_activities.mtx)
    )
    
    p <- ggscatter(to.plot, x="V1", y="V2", size=1.5, add="reg.line", add.params = list(color="black", fill="lightgray", alpha=0.5), conf.int=TRUE) +
      stat_cor(method = "pearson", label.x.npc = "middle", label.y.npc = "bottom") +
      ggrepel::geom_text_repel(data=to.plot[V1>0.40 & V2>0.40], aes(x=V1, y=V2, label=gene), size=3,  max.overlaps=100) +
      labs(x=sprintf("TF activity in %s",i), y=sprintf("TF activity in %s",j)) +
      theme(
        axis.text = element_text(size=rel(0.7))
      )
    
    pdf(sprintf("%s/scatterplots/%s_vs_%s_tf_activities.pdf",io$outdir,i,j), width = 8, height = 8)
    print(p)
    dev.off()
  }
}

##########
## Test ##
##########

tmp <- data.table(
  blood = rowMeans(tf_activities.mtx[,c("Erythroid1","Erythroid2","Erythroid3")]),
  exe = rowMeans(tf_activities.mtx[,c("Parietal_endoderm","ExE_ectoderm")]),
  blood_exe = rowMeans(tf_activities.mtx[,c("Erythroid1","Erythroid2","Erythroid3","Parietal_endoderm","ExE_ectoderm")]),
  embryonic = rowMeans(tf_activities.mtx[,!colnames(tf_activities.mtx)%in%c("Erythroid1","Erythroid2","Erythroid3","Parietal_endoderm","ExE_ectoderm")]),
  gene = rownames(tf_activities.mtx)
) 

foo <- c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
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
  "NMP",
  "Rostral_neurectoderm",
  "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm"
)

bar <- c(
  "Erythroid1",
  "Erythroid2",
  "Erythroid3",
  "ExE_ectoderm",
  "Parietal_endoderm"
)

tmp <- tf_activities.dt %>%
  .[celltype%in%c(foo,bar)] %>%
  .[,blood_exe:=celltype%in%bar]

to.plot <- tmp[,.(
  diff = mean(.SD[blood_exe==TRUE,activity]) - mean(.SD[blood_exe==FALSE,activity]),
  p.value = t.test(x=.SD[blood_exe==TRUE,activity], y=.SD[blood_exe==FALSE,activity])[["p.value"]]
), by=c("gene")] %>% .[,sig:=p.value<0.001] %>% setorder(p.value)

negative_hits <- to.plot[sig==TRUE & diff<0,gene]
positive_hits <- to.plot[sig==TRUE & diff>0,gene]
all <- nrow(to.plot)

xlim <- max(abs(to.plot$diff), na.rm=T)
ylim <- max(-log10(to.plot$p.value+1e-100), na.rm=T)

p <- ggplot(to.plot, aes(x=diff, y=-log10(p.value+1e-100))) +
  geom_segment(aes(x=0, xend=0, y=0, yend=ylim-1), color="orange", size=0.5) +
  geom_point(aes(color=sig, size=sig)) +
  # ggrastr::geom_point_rast(aes(color=sig, size=sig)) +
  ggrepel::geom_text_repel(data=head(to.plot[sig==T & diff>0],n=10), aes(x=diff, y=-log10(p.value+1e-100), label=gene), size=3,  max.overlaps=100) +
  ggrepel::geom_text_repel(data=head(to.plot[sig==T & diff<0],n=10), aes(x=diff, y=-log10(p.value+1e-100), label=gene), size=3,  max.overlaps=100) +
  scale_color_manual(values=c("black","red")) +
  scale_size_manual(values=c(0.75,1.25)) +
  scale_x_continuous(limits=c(-xlim-0.1,xlim+0.1)) +
  scale_y_continuous(limits=c(0,ylim+3)) +
  annotate("text", x=0, y=ylim+3, size=4, label=sprintf("(%d)", all)) +
  annotate("text", x=-xlim-0.1, y=ylim+3, size=4, label=sprintf("%d (-)",length(negative_hits))) +
  annotate("text", x=xlim+0.1, y=ylim+3, size=4, label=sprintf("%d (+)",length(positive_hits))) +
  labs(x="Difference in activity between community 1 and 2", y=expression(paste("-log"[10],"(p.value)"))) +
  theme_classic() +
  theme(
    axis.text = element_text(size=rel(0.75), color='black'),
    axis.title = element_text(size=rel(1.0), color='black'),
    legend.position="none"
  )

# pdf(sprintf("%s/volcano_plots/volcano_pearson_correlation.pdf",io$outdir), width = 9, height = 6)
# png(sprintf("%s/volcano_plots/volcano_pearson_correlation.png",io$outdir), width = 800, height = 500)
print(p)
# dev.off()
