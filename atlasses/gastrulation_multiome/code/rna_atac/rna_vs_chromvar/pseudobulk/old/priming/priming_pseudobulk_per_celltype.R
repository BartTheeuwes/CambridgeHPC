library(network)
library(sna)

#####################
## Define settings ##
#####################

# Load default settings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
} else {
  stop("Computer not recognised")
}

# I/O
io$paga.connectivity <- paste0(io$atlas.basedir,"/results/PAGA/PAGA_connectivity.csv")
io$paga.coordinates <- paste0(io$atlas.basedir,"/results/PAGA/PAGA_coordinates.csv")
io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/pseudobulk/priming")

# Options
opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E8.0_rep1",
  "E8.0_rep2",
  "E8.5_rep1",
  "E8.5_rep2"
)

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

opts$motif_annotation <- "Motif_cisbp" # Motif_JASPAR2020_human

###################
## Load metadata ##
###################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & pass_rnaQC==TRUE & doublet_call==FALSE] %>%
  .[sample%in%opts$samples & celltype.mapped%in%opts$celltypes] %>%
  .[,celltype.mapped:=factor(celltype.mapped,levels=opts$celltypes)] 

# subset celltypes with sufficient number of cells
opts$min.cells <- 25
sample_metadata <- sample_metadata %>%
  .[,N:=.N,by=c("celltype.mapped")] %>% .[N>opts$min.cells] %>% .[,N:=NULL] %>% droplevels
table(sample_metadata$celltype.mapped)

opts$celltypes <- unique(sample_metadata$celltype.mapped) %>% as.character

##############################################
## Load pseudobulk RNA expression estimates ##
##############################################

rna.sce.pseudobulk <- readRDS(io$rna.pseudobulk.sce)[,opts$celltypes]
rownames(rna.sce.pseudobulk) <- toupper(rownames(rna.sce.pseudobulk))

#####################################
## Load pseudobulk chromVAR scores ##
#####################################

# Precomputed
io$archR.pseudobulk.deviations.se <- sprintf("%s/pseudobulk/pseudobulk_DeviationMatrix_%s_summarized_experiment.rds",io$archR.directory,opts$motif_annotation)
chromvar.se.pseudobulk <- readRDS(io$archR.pseudobulk.deviations.se)

chromvar.se.pseudobulk <- chromvar.se.pseudobulk[rowData(chromvar.se.pseudobulk)$seqnames=="z",]
rownames(chromvar.se.pseudobulk) <- rowData(chromvar.se.pseudobulk)$name %>% toupper %>% stringr::str_split(.,"_") %>% map_chr(1)

# Recompute
# chromvar.se.pseudobulk <- getGroupSE(ArchRProject.filt, groupBy = "celltype.mapped", useMatrix = "DeviationMatrix")
# chromvar.se.pseudobulk <- chromvar.se.pseudobulk[rowData(chromvar.se.pseudobulk)$seqnames=="z",]
# rownames(chromvar.se.pseudobulk) <- rowData(chromvar.se.pseudobulk)$name
# chromvar.se.pseudobulk <- chromvar.se.pseudobulk[,opts$celltypes]
# atac.deviation.mtx.pseudobulk <- assay(chromvar.se.pseudobulk)

################################
## Load motif2gene annotation ##
################################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/atac/archR/load_motif_annotation.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/atac/archR/load_motif_annotation.R")
} else {
  stop("Computer not recognised")
}

motif2gene.dt <- motif2gene.dt %>%
  .[gene%in%rownames(chromvar.se.pseudobulk) & gene%in%rownames(rna.sce.pseudobulk)]  %>%
  .[,N:=length(unique(motif)),by="gene"] %>% .[N==1] %>% .[,N:=NULL]

rna.sce.pseudobulk <- rna.sce.pseudobulk %>% .[motif2gene.dt$gene,]
chromvar.se.pseudobulk <- chromvar.se.pseudobulk %>% .[motif2gene.dt$gene,]

################
## Parse data ##
################

chromvar_dt <- assay(chromvar.se.pseudobulk) %>%
  as.matrix %>% t %>% as.data.table(keep.rownames = T) %>%
  setnames("rn","celltype") %>% 
  melt(id.vars=c("celltype"), variable.name="gene", value.name="chromvar_zscore")

rna_dt <- logcounts(rna.sce.pseudobulk) %>% as.matrix %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","gene") %>%
  melt(id.vars="gene", variable.name="celltype", value.name="expr")

###########
## Merge ##
###########

rna_chromvar.dt <- merge(
  rna_dt,
  chromvar_dt,
  by = c("celltype","gene")
)

# chromvar_rna_dt[chromvar_zscore<0,chromvar_zscore:=0]

###########
## Blood ##
###########

celltypes.subset = c(
  "Haematoendothelial_progenitors",
  # "Endothelium",
  # "Blood_progenitors_1",
  # "Blood_progenitors_2",
  # "Erythroid1",
  # "Erythroid2",
  "Erythroid3"
)

foo <- 
bar <- rna_chromvar.dt[celltype==""]


# Find TFs which are expressed and display signatures of active chromatin in the most differentiated lineages
TFs <- rna_chromvar.dt[celltype%in%"Erythroid3"] %>% .[expr>6 & chromvar_zscore>2.5,gene]
rna_chromvar.dt[celltype%in%"Blood_progenitors_1"] %>% .[gene%in%TFs & expr<4 & chromvar_zscore>2.5,gene]

to.plot <- rna_chromvar.dt %>%
  .[celltype%in%celltypes.subset] %>% 
  dcast(gene~celltype, value.var=c("expr","chromvar_zscore")) %>%
  .[,chromvar_diff:=chromvar_zscore_Erythroid3-chromvar_zscore_Haematoendothelial_progenitors] %>%
  .[,expr_diff:=expr_Erythroid3-expr_Haematoendothelial_progenitors]

ggscatter(to.plot, x="chromvar_diff", y="expr_diff", size=1) +
  labs(y="RNA expression difference", x="Motif accessibility difference") +
  ggrepel::geom_text_repel(data=head(to.plot[expr_diff>0 & chromvar_zscore_Haematoendothelial_progenitors>2.5],n=15), aes(x=chromvar_diff, y=expr_diff, label=gene), size=3,  max.overlaps=100) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_vline(xintercept=0, linetype="dashed") +
  theme(
    axis.text = element_text(size=rel(0.7))
  )
  