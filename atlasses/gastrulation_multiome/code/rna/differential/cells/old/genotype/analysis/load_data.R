
##############
## Settings ##
##############

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
  "Blood_progenitors",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
  "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm"
  # "Visceral_endoderm",
  # "ExE_endoderm"
  # "ExE_ectoderm"
  # "Parietal_endoderm"
)

# Options
opts$wt.class <- "WT"
opts$ko.class <- "T_KO"

# Define statistical significance
opts$threshold_fdr <- 0.01
opts$min.logFC <- 1

#############################################
## Load results from differential analysis ##
#############################################

# i <- opts$groupA[1]; j <- opts$celltypes[1]
diff.dt <- opts$celltypes %>% map(function(j) {
  file <- sprintf("%s/%s_%s_vs_%s.txt.gz", io$indir,j,opts$wt.class,opts$ko.class)
  if (file.exists(file)) {
    fread(file, select=c(1,2,4,6,7)) %>% .[,c("celltype"):=list(j)]
  }
}) %>% rbindlist %>%
  .[,celltype:=factor(celltype,levels=opts$celltypes)] %>%
  .[, sign := ifelse(logFC>0,"Upregulated in KO","Downregulated in WT")] %>%
  .[, sig := (padj_fdr<=opts$threshold_fdr & abs(logFC)>=opts$min.logFC)]


# Print stats
print(sprintf("Number of celltypes: %s",length(unique(diff.dt$celltype))))
print(sprintf("Number of genes: %s",length(unique(diff.dt$gene))))
