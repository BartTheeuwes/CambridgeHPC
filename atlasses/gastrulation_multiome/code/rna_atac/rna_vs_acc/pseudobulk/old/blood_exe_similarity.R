
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
io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/exploration"); dir.create(io$outdir, showWarnings = F)

# Options
opts$celltypes = c(
  # "Epiblast",
  # "Primitive_Streak",
  # "Caudal_epiblast",
  # "PGC",
  # "Anterior_Primitive_Streak",
  # "Notochord",
  # "Def._endoderm",
  # "Gut",
  # "Nascent_mesoderm",
  # "Mixed_mesoderm",
  # "Intermediate_mesoderm",
  # "Caudal_Mesoderm",
  # "Paraxial_mesoderm",
  # "Somitic_mesoderm",
  # "Pharyngeal_mesoderm",
  # "Cardiomyocytes",
  # "Allantois",
  # "ExE_mesoderm",
  # "Mesenchyme",
  # "Haematoendothelial_progenitors",
  # "Endothelium",
  # "Blood_progenitors_1",
  # "Blood_progenitors_2",
  # "Erythroid1",
  # "Erythroid2",
  "Erythroid3",
  # "NMP",
  # "Rostral_neurectoderm",
  # "Caudal_neurectoderm",
  # "Neural_crest",
  # "Forebrain_Midbrain_Hindbrain",
  # "Spinal_cord",
  # "Surface_ectoderm",
  # "Visceral_endoderm",
  # "ExE_endoderm",
  # "ExE_ectoderm",
  "Parietal_endoderm"
)

opts$motif_annotation <- "Motif_cisbp"

##############################################
## Load pseudobulk RNA expression estimates ##
##############################################

# Load SingleCellExperiment
rna.sce <- readRDS(io$rna.pseudobulk.sce)[,opts$celltypes]

# Rename genes
rownames(rna.sce) <- toupper(rownames(rna.sce))

################################
## Load TF activity estimates ##
################################

tf_activities <- fread("/Users/ricard/data/gastrulation_multiome_10x/results/rna_atac/rna_vs_chromvar/pseudobulk/TF_activities/tf_activities.txt.gz") %>%
  .[celltype%in%opts$celltypes]

tf_activities_filt <- tf_activities %>% copy %>%
  # .[,celltype2:="blood"] %>% .[celltype%in%c("Parietal_endoderm"),celltype2:="ExE"] %>%
  # .[,.(activity=mean(activity)),by=c("gene","celltype2")] %>%
  dcast(gene~celltype) %>%
  .[Parietal_endoderm>0.5 & Erythroid3>0.5]
  
tf_activities_filt <- tf_activities %>% copy %>%
  dcast(gene~celltype) %>%
  .[,diff:=Parietal_endoderm-Erythroid3]

##################
## Venn Diagram ##
##################

tmp <- tf_activities %>% .[activity>0.5]


p <- venn.diagram(
  x = split(tmp$gene,tmp$celltype),
  filename=NULL
  # fill=brewer.pal(n=6,"Dark2")
)

# pdf(file=sprintf("%s/venn_atacQC_vs_rnaQC.pdf",io$outdir))
grid.draw(p)
# dev.off()

#############################
## Load chromVAR estimates ##
#############################

# Find TFs that have high activity in Erythroid and Parietal endoderm

################################
## Load pseudobulk PeakMatrix ##
################################

peakMatrix.se <- readRDS(io$archR.pseudobulk.peakMatrix.se)[,opts$celltypes]

# Load peak metadata
peak_metadata.dt <- fread(io$archR.peak.metadata) %>% 
  .[,idx:=sprintf("%s:%s-%s",chr,start,end)]

# Define peak names
peak_names <- rowData(peakMatrix.se) %>% as.data.table %>% .[,idx:=sprintf("%s:%s-%s",seqnames,start,end)] %>% .$id
rownames(peakMatrix.se) <- peak_names

# Subset peaks
print(mean(rownames(peakMatrix.se) %in% peak_metadata.dt$idx))
peakMatrix.se <- peakMatrix.se[rownames(peakMatrix.se) %in% peak_metadata.dt$idx]

###############################
## Load motifmatcher results ##
###############################

motifmatcher.se <- readRDS(sprintf("%s/Annotations/%s-Matches-In-Peaks.rds",io$archR.directory,opts$motif_annotation))

colnames(motifmatcher.se) <- colnames(motifmatcher.se) %>% toupper %>% stringr::str_split(.,"_") %>% map_chr(1)

motifmatcher.se <- motifmatcher.se[,!duplicated(colnames(motifmatcher.se))]

tmp <- rowRanges(motifmatcher.se)
rownames(motifmatcher.se) <- sprintf("%s:%s-%s",seqnames(tmp), start(tmp), end(tmp))

motifmatcher.se <- motifmatcher.se[rownames(peakMatrix.se),]

##############################
## Load coexpression matrix ##
##############################

cor_rna.mtx <- readRDS("/Users/ricard/data/gastrulation_multiome_10x/results/rna/coexpression/correlation_matrix_tf2gene.rds")

#########################################################
## ##
#########################################################

# TFs <- intersect(colnames(motifmatcher.se),rownames(rna.sce))
TFs <- tf_activities_filt$gene %>% head(n=3)

cor.dt <- TFs %>% map(function(i) {
  print(i)
  active_peaks_i <- rownames(motifmatcher.se)[which(assay(motifmatcher.se[,i])==1)]
  # TF_peak.cor <- cor(assay(rna.sce[i,])[1,], t(assay(peakMatrix.se[active_peaks_i,])))[1,]
  corr_output <- psych::corr.test(t(assay(rna.sce[i,])), t(assay(peakMatrix.se[active_peaks_i,])))
  sig_peaks <- which(corr_output$p[1,]<0.10 & abs(corr_output$r[1,])>0.25)
  data.table(
    TF = i, 
    peak = names(sig_peaks), 
    cor = round(corr_output$r[1,][sig_peaks],2),  
    p = round(corr_output$p[1,][sig_peaks],5)
  ) %>% return
}) %>% rbindlist

###############
TFs <- fread("/Users/ricard/data/gastrulation_multiome_10x/results/atac/archR/chromvar/differential/Motif_cisbp_Parietal_endoderm_vs_Erythroid3.txt.gz") %>%
  .[MeanDiff>8]

TFs <- c("SOX9","SOX7")

for (i in TFs) {
  all_peaks_i <- rownames(motifmatcher.se)[which(assay(motifmatcher.se[,i])==1)]
  peaks_open_blood <- names(which(assay(peakMatrix.se[all_peaks_i,"Erythroid3"])[,1]>0.3))
  peaks_open_exe <- names(which(assay(peakMatrix.se[all_peaks_i,"Parietal_endoderm"])[,1]>0.3))
  
  print(length(peaks_open_blood))
  print(length(peaks_open_exe))
  print(length(intersect(peaks_open_blood,peaks_open_exe)))
  print("\n")
}
