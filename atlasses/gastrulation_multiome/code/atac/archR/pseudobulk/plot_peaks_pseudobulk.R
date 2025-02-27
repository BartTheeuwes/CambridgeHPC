
########################
## Load ArchR project ##
########################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
} else {
  stop("Computer not recognised")
}

#####################
## Define settings ##
#####################

# I/O
io$outdir <- paste0(io$basedir,"/results/atac/archR/peaks")

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
  "ExE_ectoderm",
  "Parietal_endoderm"
)

#####################
## Load metadata ##
#####################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & doublet_call==FALSE] %>%
  .[sample%in%opts$samples & celltype.predicted%in%opts$celltypes] %>%
  .[,celltype.predicted:=factor(celltype.predicted,levels=opts$celltypes)]

##################
## Subset ArchR ##
##################

ArchRProject.filt <- ArchRProject[sample_metadata$cell,]
table(getCellColData(ArchRProject.filt,"Sample")[[1]])

#####################
## Load PeakMatrix ##
#####################

peakMatrix.se <- getMatrixFromProject(ArchRProject.filt, useMatrix = "PeakMatrix", binarize = TRUE)

# Load peak metadata
peak_metadata.dt <- fread(io$archR.peak.metadata) %>% 
  .[,idx:=sprintf("%s_%s-%s",chr,start,end)]

# Define peak names
peak_names <- rowRanges(peakMatrix.se) %>% as.data.table %>% .[,idx:=sprintf("%s_%s-%s",seqnames,start,end)] %>% .$id
rownames(peakMatrix.se) <- peak_names

##############################
## Load denoised PeakMatrix ##
##############################

################################
## Load pseudobulk PeakMatrix ##
################################

peakMatrix.pseudobulk.se <- readRDS(io$archR.pseudobulk.peakMatrix.se)[,opts$celltypes]

# Define peak names
peak_names <- rowData(peakMatrix.pseudobulk.se) %>% as.data.table %>% .[,idx:=sprintf("%s_%s-%s",seqnames,start,end)] %>% .$id
rownames(peakMatrix.pseudobulk.se) <- peak_names

##########
## Plot ##
##########

stopifnot(rownames(peakMatrix.se) == rownames(peakMatrix.pseudobulk.se))

peaks.to.plot <- rownames(peakMatrix.se)# %>% head(n=15)

for (i in 1:length(peaks.to.plot)) {
  peak <- peaks.to.plot[i]
  print(sprintf("%s/%s: %s",i,length(peaks.to.plot),peak))
  
  # Create data.table with pseudobulk estimates
  to.plot.pseudobulk <- data.table(
    celltype.predicted = colnames(peakMatrix.pseudobulk.se),
    value = assay(peakMatrix.pseudobulk.se[peak,])[1,]
  ) %>% .[,celltype.predicted:=factor(celltype.predicted,levels=opts$celltypes)] %>%
    .[,type:="pseudobulk"]
  
  # Create data.table with single-cells
  to.plot.singlecell <- data.table(
    cell = colnames(peakMatrix.se),
    value = assay(peakMatrix.se[peak,])[1,]
  ) %>% merge(sample_metadata[,c("cell","celltype.predicted")], by="cell") %>%
    .[,.(value=mean(value)),by="celltype.predicted"] %>%
    .[,type:="single_cell"]
  
  to.plot <- rbind(to.plot.pseudobulk, to.plot.singlecell)
  
  # Plot
  p <- ggbarplot(to.plot, x="celltype.predicted", y="value", fill="celltype.predicted") +
    facet_wrap(~type) +
    scale_fill_manual(values=opts$celltype.colors) +
    # theme_classic() +
    # labs(title=peak, x="",y=sprintf("%s expression",peak)) +
    guides(x = guide_axis(angle = 90)) +
    labs(title=peak, x="",y="Chromatin accessibility") +
    theme(
      # strip.text = element_text(size=rel(1.0)),
      plot.title = element_text(hjust = 0.5, size=rel(0.9), color="black"),
      # # plot.title = element_blank(),
      axis.text.y = element_text(colour="black",size=rel(0.8)),
      axis.text.x = element_text(colour="black",size=rel(0.4)),
      # axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      # axis.text.y = element_text(colour="black",size=rel(1.0)),
      # axis.title.y = element_text(colour="black",size=rel(1.2)),
      legend.position="none"
      # legend.title = element_blank(),
      # legend.text = element_text(size=rel(1.1))
    )
  
  pdf(sprintf("%s/%s.pdf",io$outdir,peak), width=7.5, height=6)
  # jpeg(sprintf("%s/%s.jpeg",io$outdir,peak), width = 1400, height = 700)
  print(p)
  dev.off()
}
