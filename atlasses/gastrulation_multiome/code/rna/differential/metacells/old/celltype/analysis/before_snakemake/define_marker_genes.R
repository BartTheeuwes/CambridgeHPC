#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

io$differential_results <- file.path(io$basedir,"results/rna/differential/metacells/celltype")
io$outdir <- file.path(io$basedir,"results/rna/differential/metacells/celltype/marker_genes"); dir.create(io$outdir, showWarnings = F)

#############
## Options ##
#############

opts$celltypes <- c(
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
)# %>% head(n=4)

# Minimum fraction of significant differential pairwise comparisons
opts$score <- 0.75

# Significance thresholds
opts$min.log2FC <- 1
opts$fdr <- 0.01

##################
## Load results ##
##################

# i <- "Gut"; j <- "NMP"
dt <- opts$celltypes %>% map(function(i) { opts$celltypes %>% map(function(j) {
  file <- file.path(io$differential_results,sprintf("%s_vs_%s.txt.gz",i,j))
  if (file.exists(file)) {
    fread(file, select = c(1,2,4)) %>% 
      .[abs(logFC)>=opts$min.log2FC & padj_fdr<=opts$fdr] %>%
      .[,c("celltypeA","celltypeB"):=list(i,j)] %>%
      return
  } }) %>% rbindlist }) %>% rbindlist %>%
  .[,celltypeA:=factor(celltypeA,levels=opts$celltypes)] %>%
  .[,celltypeB:=factor(celltypeB,levels=opts$celltypes)] %>%
  .[,direction:=as.factor(c("down","up"))[as.numeric(logFC<0)+1]]

ncelltypes <- unique(c(as.character(unique(dt$celltypeA)),as.character(unique(dt$celltypeB)))) %>% length

#########################
## Define marker genes ##
#########################

foo <- dt[,.(score=sum(direction=="up")), by=c("celltypeA","gene")] %>% setnames("celltypeA","celltype")
bar <- dt[,.(score=sum(direction=="down")), by=c("celltypeB","gene")] %>% setnames("celltypeB","celltype")

markers_genes.dt <- merge(foo,bar,by=c("celltype","gene"), all=TRUE) %>% 
  .[is.na(score.x),score.x:=0] %>% .[is.na(score.y),score.y:=0] %>%
  .[,score:=score.x+score.y] %>%
  .[,c("score.x","score.y"):=NULL] %>%
  .[,score:=round(score/(ncelltypes-1),2)] %>%
  # .[score>=opts$score] %>%
  setorder(celltype,-score)
rm(foo,bar)

##########
## Save ##
##########

# Save marker score for all combination of genes and cell types
length(unique(markers_genes.dt$gene))
length(unique(markers_genes.dt$celltype))
fwrite(markers_genes.dt, file.path(io$outdir,"marker_genes_upregulated_all.txt.gz"), sep="\t")

# Save marker score for strong markers
# markers_genes_filt.dt <- markers_genes.dt %>% .[score>=opts$score & logFC>=2]
markers_genes_filt.dt <- markers_genes.dt %>% .[score>=opts$score]
length(unique(markers_genes_filt.dt$gene))
length(unique(markers_genes_filt.dt$celltype))
fwrite(markers_genes_filt.dt, file.path(io$outdir,"marker_genes_upregulated_filtered.txt.gz"), sep="\t")

