
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
# io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc")

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
  "Haematoendothelial_progenitors",
  "Endothelium",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3"
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
  # "Parietal_endoderm"
)


################################################
## Load pseudobulk RNA and chromVAR estimates ##
################################################

opts$motif_annotation <- "Motif_cisbp"
io$archR.pseudobulk.deviations.se <- sprintf("%s/pseudobulk/pseudobulk_DeviationMatrix_%s_summarized_experiment.rds",io$archR.directory,opts$motif_annotation)

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/rna_atac/load_rna_atac_pseudobulk.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/rna_atac/load_rna_atac_pseudobulk.R")
} else {
  stop("Computer not recognised")
}

#############################
## Load peak2gene linkages ##
#############################

peak2gene.dt <- fread(io$archR.peak2gene.all) %>% 
  # .[gene%in%opts$genes] %>%
  .[,peak:=sprintf("chr%s:%s-%s",chr,peak.start,peak.end)]

###################################################
## Load TF2peak correlation results (pseudobulk) ##
###################################################

tf2peak_cor.dt <- fread(io$tf2peak_cor.dt) %>%
  # .[!is.na(cor) & peak%in%unique(peak2gene.dt$peak) & TF%in%opts$TFs] %>%
  .[,cor_sign:=c("-","+")[(cor>0)+1]]

#######################################################
## Load RNA vs chromVAR correlation results per gene ##
#######################################################

io$file <- paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/pseudobulk/per_gene/cor_rna_vs_chromvar_correlated_peaks_pseudobulk.txt.gz")
cor_rna_vs_chromvar_per_gene.dt <- fread(io$file) %>%
  .[,cor_sign:=as.factor(c("Repressor","Activator")[(r>0)+1])]

###############################
## Load motifmatcher results ##
###############################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/load_motifmatchR.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/load_motifmatchR.R")
} else {
  stop("Computer not recognised")
}


#######################################
## Load pseudobulk RNA and ATAC data ##
#######################################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/rna_atac/load_rna_atac_pseudobulk.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/rna_atac/load_rna_atac_pseudobulk.R")
} else {
  stop("Computer not recognised")
}

# - Identify epigenetic priming for one of the two lineages
# - Explore why so many repressors in the erythroid lineage
# - Distinguish between upregulated and downregulated peaksº
# repression of endothelium fate in the erythroid lineage?

####################################
## Load differnetial ATAC results ##
####################################

# atac.diff <- fread(paste0(io$archR.peak.differential.dir,"/PeakMatrix_Endothelium_vs_Erythroid3.txt.gz")) %>%
# atac.diff <- fread(paste0(io$archR.peak.differential.dir,"/PeakMatrix_Blood_progenitors_2_vs_Endothelium.txt.gz")) %>%
atac_END_vs_ERI.diff <-  fread(paste0(io$archR.peak.differential.dir,"/PeakMatrix_Endothelium_vs_Erythroid3.txt.gz")) %>%
  .[,sig:=FDR<0.01 & abs(MeanDiff)>=0.20] %>% 
  .[,sign:="Up in Endothelium"] %>% 
  .[MeanDiff<0,sign:=c("Up in Erythroid")] %>%
  .[sig==T]

atac_BP_vs_END.diff <-  fread(paste0(io$archR.peak.differential.dir,"/PeakMatrix_Blood_progenitors_2_vs_Endothelium.txt.gz")) %>%
  .[,sig:=FDR<0.01 & abs(MeanDiff)>=0.20] %>% 
  .[,sign:="Up in Blood_progenitors"] %>% 
  .[MeanDiff<0,sign:=c("Up in Endothelium")] %>%
  .[sig==T]

atac_HEP_vs_ERI.diff <-  fread(paste0(io$archR.peak.differential.dir,"/PeakMatrix_Haematoendothelial_progenitors_vs_Erythroid3.txt.gz")) %>%
  .[,sig:=FDR<0.01 & abs(MeanDiff)>=0.20] %>% 
  .[,sign:="Up in Haematoendothelial_progenitors"] %>% 
  .[MeanDiff<0,sign:=c("Up in Erythroid")] %>%
  .[sig==T]
  
atac_HEP_vs_END.diff <-  fread(paste0(io$archR.peak.differential.dir,"/PeakMatrix_Haematoendothelial_progenitors_vs_Endothelium.txt.gz")) %>%
  .[,sig:=FDR<0.01 & abs(MeanDiff)>=0.20] %>% 
  .[,sign:="Up in Haematoendothelial_progenitors"] %>% 
  .[MeanDiff<0,sign:=c("Up in Endothelium")] %>%
  .[sig==T]

atac_HEP_vs_BP.diff <- fread(paste0(io$archR.peak.differential.dir,"/PeakMatrix_Haematoendothelial_progenitors_vs_Blood_progenitors_2.txt.gz")) %>%
  .[,sig:=FDR<0.01 & abs(MeanDiff)>=0.20] %>% 
  .[,sign:="Up in Haematoendothelial_progenitors"] %>% 
  .[MeanDiff<0,sign:=c("Up in Blood_progenitors")] %>%
  .[sig==T]

#############
## Explore ##
#############

endothelium_up_from_hep_peaks <- atac_HEP_vs_END.diff[sign=="Up in Endothelium",idx]
hep_up_vs_endothelium_peaks <- atac_HEP_vs_END.diff[sign=="Up in Haematoendothelial_progenitors",idx]
hep_down_vs_endothelium_peaks <- atac_HEP_vs_END.diff[sign=="Up in Endothelium",idx]
blood_up_from_hep_peaks <- atac_HEP_vs_BP.diff[sign=="Up in Blood_progenitors",idx]
blood_up_vs_endothelium_peaks <- atac_BP_vs_END.diff[sign=="Up in Blood_progenitors",idx]
endothelium_up_vs_blood_peaks <- atac_BP_vs_END.diff[sign=="Up in Endothelium",idx]


to.plot <- atac.dt %>% 
  # .[peak%in%endothelium_up_vs_blood_peaks] 
  .[peak%in%blood_up_vs_endothelium_peaks] 

p <- ggboxplot(to.plot, x="celltype", y="acc", outlier.shape=NA) +
  # coord_cartesian(ylim=c(0,1)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(color="black", angle=30, hjust=1)
  )

# pdf(sprintf("%s/archr_umap_celltype_%s_LSIiter%s_nfeatures%s_neighb%s_mindist%s.pdf",io$outdir,opts$matrix,opts$lsi.iterations, opts$lsi.varFeatures, opts$umap.neighbours,opts$umap.minDist))
print(p)
# dev.off()


# ARE HEP derived from the endothelium?
######################
## Motif enrichment ##
######################

features <- hep_down_vs_endothelium_peaks# atac.diff[sig==T & sign=="Up in Haematoendothelial_progenitors",idx]
background <- peak_names# hep_up_vs_endothelium_peaks# atac.diff[sig==T & sign=="Up in Endothelium",idx]
motifmatcher_filt1.se <- motifmatcher.se[features,]
motifmatcher_filt2.se <- motifmatcher.se[background,]
# features <- Matrix::colSums(assay(motifmatcher_filt1.se))
# background <- Matrix::colSums(assay(motifmatcher_filt2.se))


query.motifs <- assay(motifmatcher.se[features, ])
background.motifs <- assay(motifmatcher.se[background, ])
query.counts <- colSums(query.motifs)
background.counts <- colSums(background.motifs)
percent.observed <- query.counts / length(features) * 100
percent.background <- background.counts / length(background) * 100
fold.enrichment <- percent.observed / percent.background
p.list <- vector(mode = "numeric")
for (i in seq_along(query.counts)) {
  p.list[[i]] <- phyper(
    q = query.counts[[i]] - 1,
    m = background.counts[[i]],
    n = nrow(x = background.motifs) - background.counts[[i]],
    k = length(x = features),
    lower.tail = FALSE
  )
}

fold.enrichment[fold.enrichment>=3] %>% names 
