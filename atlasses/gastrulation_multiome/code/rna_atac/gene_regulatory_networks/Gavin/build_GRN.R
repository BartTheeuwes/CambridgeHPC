
# Load default settings
library(tidyr)
library(ggplot2)
library(BuenColors)
library(ggrepel)
library(anndata)
library(plyr)
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/utils.R")
  source("/Users/ricard/gastrulation/gastrulation_multiome_10x/rna_atac/gene_regulatory_networks/GRN/functions/DORC_functions.R")
  source("/Users/ricard/gastrulation/gastrulation_multiome_10x/rna_atac/gene_regulatory_networks/GRN/functions/GRN_functions.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/utils.R")
  source("/homes/ricard/gastrulation/gastrulation_multiome_10x/rna_atac/gene_regulatory_networks/GRN/functions/DORC_functions.R")
  source("/homes/ricard/gastrulation/gastrulation_multiome_10x/rna_atac/gene_regulatory_networks/GRN/functions/GRN_functions.R")
} else if (grepl("Workstation",Sys.info()['nodename'])){
  source("/home/lijingyu/gastrulation/gastrulation_multiome_10x/settings.R")
  source("/home/lijingyu/gastrulation/gastrulation_multiome_10x/utils.R")
  source("/home/lijingyu/gastrulation/gastrulation_multiome_10x/rna_atac/gene_regulatory_networks/GRN/functions/DORC_functions.R")
  source("/home/lijingyu/gastrulation/gastrulation_multiome_10x/rna_atac/gene_regulatory_networks/GRN/functions/GRN_functions.R")
  
} else {
  stop("Computer not recognised")
}

################
## Define I/O ##
################

io$atac.peaks.bulk.se <- paste0(io$basedir,"/processed/atac/archR/pseudobulk/pseudobulk_PeakMatrix_summarized_experiment.rds")
io$rna.bulk.se <- paste0(io$basedir,"/processed/rna/pseudobulk/SingleCellExperiment.rds")
io$outdir <- paste0(io$basedir,"/results/rna_atac/GRN")
io$vritual_chip.mtx <- paste0(io$basedir,'/results/rna_atac/virtual_chipseq/virtual_chip_mtx.rds')
io$FigR_GRN <- paste0(io$outdir,"/FigR_GRN.csv")
io$DORC <- paste0(io$outdir,'/DORC.rds')
io$chip_GRN.mtx <- paste0(io$outdir,'/chip_GRN_mtx.csv')
io$chip_GRN.df <- paste0(io$outdir,'/chip_GRN_df.csv')


##############
## Load data##
##############

ATAC.se <- readRDS(io$atac.peaks.bulk.se)
assayNames(ATAC.se) <- 'counts'
peaks <- ATAC.se@elementMetadata
peaks <- paste0(peaks$seqnames,':',peaks$start,'-',peaks$end)
rownames(ATAC.se) <- 1:nrow(ATAC.se)

# Load paired RNA
RNA.se <- readRDS(io$rna.bulk.se)
marker_genes.dt <- fread(io$rna.atlas.marker_genes) 
RNA.se <- RNA.se[intersect(rownames(RNA.se),unique(marker_genes.dt$gene)),]
rnaMat <- assays(RNA.se)$counts

##############
##Run DORC ###
##############
# Use 8 cores here
DORC <- runGenePeakcorr(ATAC.se = ATAC.se,
                        RNAmat = rnaMat,
                        genome = "mm10",
                        windowPadSize = 50000,
                        nCores = 8,
                        n_bg = 100,
                        p.cut = NULL)

DORC$Peak <-  DORC$Peak %>% as.character(.) %>% mapvalues(., from=1:nrow(ATAC.se), to=peaks)

saveRDS(DORC,io$DORC)
DORC <- readRDS(io$DORC)

# Filter associations using correlation p-value 
DORC.filt <- DORC %>% filter(pvalZ <= 0.2)


##############
##Build GRN ##
##############
# There are two approaches, FigR or virtual chip

### Using virtual chip
virtual_chip.mtx <- readRDS(io$vritual_chip.mtx)
TF_target.df <- build_GRN_df(virtual_chip.mtx,DORC.filt)
TF_target.mtx <- build_GRN_matrix(virtual_chip.mtx,DORC.filt)
write.csv(TF_target_matrix,file=chip_GRN.mtx)
write.csv(TF_target_df,file=io$chip_GRN.df,row.names = FALSE)

### Using FigR
# Make J plot to rank genes by # significant peaks
dorcGenes <- dorcJplot(dorcTab = DORC.filt,
                       cutoff = 1,
                       returnGeneList = TRUE)

dorcGenes <- dorcGenes[dorcGenes%in%unique(marker_genes.dt$gene)]


DORCs.granges <-DORC.filt[DORC.filt$Gene %in% dorcGenes,]

# Compute gene score
rownames(rnaMat) <- gsub(x=rownames(rnaMat),pattern = "-",replacement = "",fixed = TRUE)
genes <- unique(DORCs.granges$Gene) %>% as.character
DORCs.mtx <- matrix(NA, nrow=length(genes), ncol=length(opts$celltypes))
rownames(DORCs.mtx) <- genes; colnames(DORCs.mtx) <- opts$celltypes

for (i in genes) {
  peaks <- DORCs.granges[DORCs.granges$Gene==i,]$Peak %>% as.character 
  DORCs.mtx[i,] <- colMeans(assay(ATAC.se[peaks,]))
}

## runFigR could compute 
#(1) correlation between TF expression and DORC accessibility
#(2) TF's motif enrichment in the DORC peaks

GRN_FigR <- runFigR(ATAC.se = ATAC.se,
                    dorcK = 30,
                    dorcTab = DORCs.granges,
                    genome = 'mm10', 
                    dorcMat = DORCs.mtx,
                    rnaMat = rnaMat,
                    n_bg = 50,
                    nCores = 8)
# transform and save
df <- GRN_FigR[GRN_FigR$Enrichment.log10P>1&GRN_FigR$Corr.log10P>1,c(1,2)]
df$value <- 1
TF_DORC_matrix <- spread(df, Motif, value)
TF_DORC_matrix[is.na(TF_DORC_matrix)] = 0
write.csv(TF_DORC_matrix,file=io$FigR_GRN)


## Visualisation

gAll <- ggplot(GRN_FigR,aes(Corr.log10P,Enrichment.log10P,color=Score)) + 
  geom_point(size=0.01,shape=16) + 
  theme_classic() + 
  labs(x="-sign*Corr.log10P",y="-sign(Enrichment.Z)*Enrichment.log10P")+
  scale_color_gradientn(colours = jdb_palette("solar_extra"),limits=c(-4,4),oob = scales::squish)

gAll
# Mean plot

rankDrivers(figR.d = GRN_FigR)
plotDrivers(figR.d = GRN_FigR,marker = "Foxa2",)

# Heatmap 
figR_heat <- plotfigRHeatmap(figR.d = GRN_FigR,
                             score.cut = 1.5,
                             DORCs = DORC_toplot,
                             column_names_gp=gpar(fontsize=5),
                             show_row_dend = FALSE,
                             row_names_side = "left")

# D3Network
DORC_toplot <- c(as.character(unique(GRN_FigR$DORC)[1:10]))
plotfigRNetwork(GRN_FigR,score.cut = 1.5,DORCs = DORC_toplot)


# Using ggnet2

DORC_toplot <- GRN_FigR[GRN_FigR$Motif=='Tal1'&abs(GRN_FigR$Score)>=1.5,]$DORC
lildat <- GRN_FigR %>% filter(DORC %in% DORC_toplot) %>% filter(abs(Score) >= 1.5)
lildat$Motif <- paste0(lildat$Motif, ".")
lildat$DORC <- paste0(lildat$DORC)
lildat$weight <- scales::rescale(abs(lildat$Score),to=c(0.5,3))

net <- network(lildat[,c("Motif","DORC")],
               directed = TRUE,
               matrix.type="edgelist")

set.edge.value(x = net,attrname = "weight",value = lildat$weight)


# Access node name, and assign as either DORC/Motif
net %v% "class" = ifelse(sapply(net$val,"[[",2) %in% lildat$DORC, "DORC", "Motif")

set.vertex.attribute(net,"vertex.names",gsub(pattern = ".",replacement = "",x = get.vertex.attribute(net,"vertex.names"),fixed = TRUE))

# Assign edge color based on positive or negative correlation
set.edge.attribute(net, "color", ifelse(lildat$Corr > 0, "firebrick", "steelblue"))

# set colors for each mode
col = c("Motif" = "gray60", "DORC" = "sandybrown")

ggnet2(net = net,color.legend = "C",
       size.legend = "Degree",
       label = TRUE,
       size="class",
       #size="outdegree",
       edge.size = "weight",
       edge.color = "color",
       #node.label=NA,
       edge.alpha = 0.75,
       alpha = 0.7,
       color="class",
       size.palette = c("DORC" = 3, "Motif" = 1),
       label.size = 3,
       #palette = "Set2",
       palette=col,
       #arrow.size = 4,
       #arrow.gap = 0.025,
       legend.position = "none") + 
  coord_equal()








