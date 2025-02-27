library(pheatmap)

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
io$chromvar_markers <- paste0(io$basedir,"/results/atac/archR/chromvar/differential/markers/marker_TFs.txt.gz")
io$outdir <- paste0(io$basedir,"/results/atac/archR/chromvar/differential/markers/pdf"); dir.create(io$outdir, showWarnings = F)

# Options
opts$motif_annotation <- "Motif_cisbp"
opts$rename.genes <- c(
  "TCFAP" = "TFAP",
  "NKX2" = "NKX2-",
  "NKX3" = "NKX3-",
  "NKX6" = "NKX6-"
)
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

###############
## Load data ##
###############

chromvar_markers.dt <- fread(io$chromvar_markers) %>%
  .[,celltype%in%opts$celltypes] %>%
  .[,gene:=idx %>% toupper %>% stringr::str_split(.,"_") %>% map_chr(1)] %>%
  .[,gene:=stringr::str_replace_all(gene,opts$rename.genes)] %>%
  .[,N:=.N,by=c("gene","celltype")] %>% .[N==1] %>% .[,N:=NULL]

chromvar_markers.dt <- chromvar_markers.dt %>%
  .[,celltype:=stringr::str_replace_all(celltype,opts$aggregated.celltypes)] %>%
  .[,.(score=mean(score)),by=c("celltype","gene")]

# Select the top N markers per cell type
chromvar_markers.dt <- chromvar_markers.dt %>% setorder(-score) %>% .[,.SD[1:3],by=c("celltype")]

###############################
## Load pseudobulk estimates ##
###############################

io$archR.pseudobulk.deviations.se <- sprintf("%s/results/atac/archR/chromvar/pseudobulk/chromVAR_deviations_summarized_experiment_%s_pseudobulk_correlated_peaks.rds",io$basedir,opts$motif_annotation)
source("/Users/ricard/gastrulation_multiome_10x/rna_atac/rna_vs_chromvar/pseudobulk/load_rna_chromvar_pseudobulk.R")

chromvar.dt <- chromvar.dt %>%
  .[gene%in%chromvar_markers.dt$gene] %>% 
  .[celltype%in%opts$celltypes] %>% 
  droplevels

chromvar.dt <- chromvar.dt %>%
  .[,celltype:=stringr::str_replace_all(celltype,opts$aggregated.celltypes)] %>%
  .[,.(chromvar_zscore=mean(chromvar_zscore)),by=c("celltype","gene")]

unique(chromvar.dt$celltype)

#############
## Heatmap ##
#############

chromvar.mtx <- chromvar.dt %>% 
  dcast(celltype~gene,value.var="chromvar_zscore") %>%
  matrix.please

chromvar.mtx[chromvar.mtx>=150] <- 150
# chromvar.mtx[chromvar.mtx>=10] <- 10
chromvar.mtx[chromvar.mtx<0] <- 0

# chromvar.mtx <- minmax.normalisation(chromvar.mtx)
chromvar.mtx <- log(chromvar.mtx+0.01)

pheatmap(
  mat = chromvar.mtx, 
  cluster_cols = T, cluster_rows = T, 
  color = colorRampPalette(c("#F2F2F2", "#CD0000"))(100),
  fontsize_row = 6,
  fontsize_col = 5,
  legend = FALSE,
  # angle_col = 0,
  # scale = "column",
  filename = sprintf("%s/heatmap_marker_TFs_chromVAR.pdf",io$outdir),
  width = 8, height = 4
)


##########
## Plot ##
##########

# Plot number of marker peaks per cell types

to.plot <- chromvar_markers.dt %>% .[,.N,by=c("celltype")]

p <- ggbarplot(to.plot, x="celltype", y="N", fill="celltype") +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(x="", y="Number of marker TFs") +
  theme(
    axis.text.y = element_text(size=rel(0.75)),
    axis.text.x = element_text(colour="black",size=rel(0.8), angle=90, hjust=1, vjust=0.5),
    axis.ticks.x = element_blank(),
    legend.position = "none"
)

pdf(sprintf("%s/barplot_number_marker_TFs.pdf",io$outdir), width = 9, height = 5)
print(p)
dev.off()

# Plot exclusivity of cell types
to.plot <- chromvar_markers.dt %>% .[,N:=.N,by="idx"]

p <- ggboxplot(to.plot, x="celltype", y="N", fill="celltype", color="black", outlier.shape=NA) +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(x="", y="Exclusivity of TF markers\n(the smaller the more exclusive)") +
  theme(
    axis.text.y = element_text(size=rel(0.75)),
    axis.title.y = element_text(size=rel(0.85)),
    axis.text.x = element_text(colour="black",size=rel(0.7), angle=90, hjust=1, vjust=0.5),
    legend.position = "none"
  )

pdf(sprintf("%s/boxplot_exclusivity.pdf",io$outdir), width = 9, height = 5)
print(p)
dev.off()

# Plot exclusivity of markers
to.plot <- chromvar_markers.dt %>%
  .[,.(Nx=.N),by="idx"] %>%
  .[,Nx:=factor(Nx)] %>%
  .[,.(Ny=.N),by="Nx"]

p <- ggbarplot(to.plot, x="Nx", y="Ny", fill="gray70") +
  labs(x="Number of different cell types per marker TF", y="") +
  theme(
    axis.text = element_text(size=rel(0.75)),
  )
pdf(sprintf("%s/boxplot_exclusivity2.pdf",io$outdir), width = 7, height = 5)
print(p)
dev.off()
