########################
## Load ArchR Project ##
########################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/utils.R")
} else {
  stop("Computer not recognised")
}

#####################
## Define settings ##
#####################

# I/O
io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/pseudobulk/per_gene")

opts$motif_annotation <- "Motif_cisbp" # Motif_JASPAR2020_human

###################################
## Load precomputed correlations ##
###################################

io$cor_rna_chromvar.pseudobulk <- "/Users/ricard/data/gastrulation_multiome_10x/results/rna_atac/rna_vs_chromvar/pseudobulk/per_gene/cor_rna_vs_chromvar_pseudobulk.txt.gz"

cor_rna_chromvar.dt <- fread(io$cor_rna_chromvar.pseudobulk) 

# Divide into activators and repressors
cor_rna_chromvar.dt %>%
  .[,class:="activator"] %>%
  .[r<0,class:="repressor"]
table(cor_rna_chromvar.dt$class)

######################
## Load motifmatchr ##
######################

motifmatcher.se <- readRDS(sprintf("%s/Annotations/%s-Matches-In-Peaks.rds",io$archR.directory,opts$motif_annotation))
colnames(motifmatcher.se) <- colnames(motifmatcher.se) %>% toupper %>% stringr::str_split(.,"_") %>% map_chr(1)
motifmatcher.se <- motifmatcher.se[,!duplicated(colnames(motifmatcher.se))]
# tmp <- rowRanges(motifmatcher.se)
# rownames(motifmatcher.se) <- sprintf("%s:%s-%s",seqnames(tmp), start(tmp), end(tmp))

####################################
## Compare expression of regulons ##
####################################

foo <- fread("/Users/ricard/data/gastrulation_multiome_10x/results/rna/pyscenic/auc_gastrulation_all_add_cor_new.csv.gz", header=T) %>%
 .[,V1:=NULL] %>% .[,TF:=toupper(TF)]

i <- "ZEB1"
i <- "ZEB2"
i <- "ETV1"

bar <- foo[TF==i]



asd <- readRDS("/Users/ricard/data/gastrulation_multiome_10x/results/rna/coexpression/correlation_matrix_tf2gene.rds")
asdd <- data.table(
  target = colnames(asd),
  cor = asd[i,]
) %>% merge(bar,by="target")


foo[importance>50,mean(regulation),by="TF"] %>% View


motifmatcher.se.filt <- motifmatcher.se[assay(motifmatcher.se[,i])[,1]==1,]
as.data.table(rowRanges(motifmatcher.se.filt)) %>% .[nearestGene%in%asdd[importance>15,target]] %>% View
