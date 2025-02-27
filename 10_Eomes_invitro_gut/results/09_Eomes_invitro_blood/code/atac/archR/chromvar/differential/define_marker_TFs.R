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
io$archR.diff.dir <- paste0(io$basedir,"/results/atac/archR/chromvar/differential")
io$outdir <- paste0(io$basedir,"/results/atac/archR/chromvar/differential/markers"); dir.create(io$outdir, showWarnings = F)

# Options
# opts$groups <- strsplit(list.files(io$diff.dir, pattern="*.gz"),"_vs_") %>% map(~ .[[1]]) %>% unlist %>% unique
opts$celltypes <- c(
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

opts$min.MeanDiff <- 1.5
opts$fdr <- 0.01

# Minimum fraction of significant differential pairwise comparisons
opts$score <- 0.75

opts$motif_annotation <- "Motif_cisbp"


##################
## Load results ##
##################

dt <- opts$celltypes %>% map(function(i) { opts$celltypes %>% map(function(j) {
  file <- sprintf("%s/%s_%s_vs_%s.txt.gz", io$archR.diff.dir,opts$motif_annotation,i,j)
  if (file.exists(file)) {
    fread(file) %>% .[,c("celltypeA","celltypeB"):=list(as.factor(i),as.factor(j))] %>%
      return
  }
}) %>% rbindlist }) %>% rbindlist %>% 
  .[,sig:=FALSE] %>% .[abs(MeanDiff)>opts$min.MeanDiff & FDR<opts$fdr,sig:=TRUE] %>%
  .[,direction:=c("up","down")[as.numeric(MeanDiff<0)+1]]  # up = higher accessibility in celltype A

dt %>% setnames("name","idx")

ncelltypes <- length(intersect(unique(dt$celltypeA),unique(dt$celltypeB)))

####################
## Define markers ##
####################

foo <- dt[,.(score=sum(sig==T & direction=="up")), by=c("celltypeA","idx")] %>% setnames("celltypeA","celltype")
bar <- dt[,.(score=sum(sig==T & direction=="down")), by=c("celltypeB","idx")] %>% setnames("celltypeB","celltype")
  
dt.filt <- merge(foo,bar,by=c("celltype","idx"), all=TRUE) %>% .[,score:=score.x+score.y] %>%
  .[,c("score.x","score.y"):=NULL] %>%
  .[,score:=round(score/(ncelltypes+1),2)] %>%
  .[score>=opts$score] %>%
  setorder(celltype,-score)
rm(foo,bar)

# dt.filt[celltype=="Erythroid1"]
# dt[celltypeA=="Erythroid2" | celltypeB=="Erythroid2"]
# foo[celltypeA=="Mesoderm"] %>% View

##########
## Save ##
##########

fwrite(dt.filt, paste0(io$outdir,"/marker_TFs.txt.gz"))
