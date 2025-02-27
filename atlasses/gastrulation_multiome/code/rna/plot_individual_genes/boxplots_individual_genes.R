#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

# I/O
io$outdir <- paste0(io$basedir,"/results/rna/individual_genes")

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
	"Parietal_endoderm",
	"ExE_ectoderm"
)


############################
## Update sample metadata ##
############################

sample_metadata <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE] %>%
  .[celltype.predicted%in%opts$celltypes]
table(sample_metadata$celltype.predicted)

###############
## Load data ##
###############

# Load SingleCellExperiment object
sce <- load_SingleCellExperiment(io$rna.sce, normalise = TRUE, cells = sample_metadata$cell)

# Add sample metadata to the colData of the SingleCellExperiment
colData(sce) <- sample_metadata %>% as.data.frame %>% tibble::column_to_rownames("cell") %>%
  .[colnames(sce),] %>% DataFrame()

# Load gene metadata
gene_metadata <- fread(io$gene_metadata) %>%
  .[symbol%in%rownames(sce)]

################
## Parse data ##
################

##########
## Plot ##
##########

celltype.order <- opts$celltypes

# genes.to.plot <- rownames(sce)[grep("Gata",rownames(sce))]
genes.to.plot <- c("Creg1")


for (i in 1:length(genes.to.plot)) {
  gene <- genes.to.plot[i]
  print(sprintf("%s/%s: %s",i,length(genes.to.plot),gene))
  
  # Create data.table to plot
  to.plot <- data.table(
    cell = colnames(sce),
    expr = logcounts(sce[gene,])[1,]
  ) %>% merge(sample_metadata, by="cell") %>%
    .[,celltype.predicted:=factor(celltype.predicted,levels=celltype.order)]

  # Plot
  p <- ggplot(to.plot, aes(x=celltype.predicted, y=expr, fill=celltype.predicted)) +
    geom_violin(scale = "width", alpha=0.8) +
    geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.8) +
    # geom_jitter(size=2, shape=21, stroke=0.2, alpha=0.5) +
    scale_fill_manual(values=opts$celltype.colors) +
    theme_classic() +
    # labs(title=gene, x="",y=sprintf("%s expression",gene)) +
    guides(x = guide_axis(angle = 90)) +
    labs(title=gene, x="",y="RNA expression") +
    theme(
      strip.text = element_text(size=rel(1.0)),
      plot.title = element_text(hjust = 0.5, size=rel(1.1), color="black"),
      # plot.title = element_blank(),
      axis.text.x = element_text(colour="black",size=rel(0.8)),
      # axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(colour="black",size=rel(1.0)),
      axis.title.y = element_text(colour="black",size=rel(1.2)),
      legend.position="none",
      legend.title = element_blank(),
      legend.text = element_text(size=rel(1.1))
    )
    
  # pdf(sprintf("%s/%s.pdf",io$outdir,i), width=5, height=3.5)
  # jpeg(sprintf("%s/%s.jpeg",io$outdir,gene), width = 1400, height = 700)
  print(p)
  # dev.off()
}



##########
## Test ##
##########


gene <- "Foxa2"

to.plot <- data.table(
  cell = colnames(sce),
  expr = logcounts(sce[gene,])[1,]
) %>% merge(sample_metadata[,c("cell","celltype.predicted")], by="cell") %>%
  .[,celltype.predicted:=factor(celltype.predicted,levels=rev(celltype.order))]

p <- ggplot(to.plot, aes(x=celltype.predicted, y=expr, fill=celltype.predicted)) +
  # geom_violin(width=0.4, alpha=0.8) +
  geom_boxplot(width=0.4, outlier.shape=NA, alpha=0.8, coef=0.5) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme_classic() +
  coord_flip() +
  labs(x="",y="RNA expression") +
  theme(
    axis.line = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(colour="black",size=rel(1.0)),
    axis.title.x = element_text(colour="black",size=rel(1.2)),
    legend.position="none"
  )

pdf(sprintf("%s/test.pdf",io$outdir), width=2.5, height=8)
print(p)
dev.off()