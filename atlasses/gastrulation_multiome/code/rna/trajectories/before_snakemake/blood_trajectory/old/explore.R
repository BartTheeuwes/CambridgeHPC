library(ggpubr)

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
} else {
  stop("Computer not recognised")
}

# I/O
io$metadata <- paste0(io$basedir, "/sample_metadata.txt.gz")
io$outdir <- paste0(io$basedir, "/results/rna_atac/test")
io$mofa.output <- paste0(io$basedir, "/results/rna_atac/mofa/mofa_model.rds")

# Options
opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E8.5_rep1",
  "E8.5_rep2"
)

opts$celltypes <- c(
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
  # "Endothelium",
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

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & pass_rnaQC==TRUE] %>%
  .[sample%in%opts$samples & celltype.predicted%in%opts$celltypes]
table(sample_metadata$celltype.predicted)

###############################
## Load pseudobulk ATAC peaks ##
###############################

io$atac.pseudobulk.peaks <- paste0(io$archR.directory,"/pseudobulk/pseudobulk_PeakMatrix_summarized_experiment.rds")

pseudobulk.atac.peak.se <- readRDS(io$atac.pseudobulk.peaks)[,opts$celltypes]
row.ranges.dt <- rowData(pseudobulk.atac.peak.se) %>% as.data.table %>% 
  setnames("seqnames","chr") %>%
  .[,idx:=sprintf("%s_%s_%s",chr,start,end)]
rownames(pseudobulk.atac.peak.se) <- row.ranges.dt$idx

###################################
## Load pseudobulk chromVAR data ##
###################################

io$atac.pseudobulk.chromvar <- paste0(io$archR.directory,"/pseudobulk/pseudobulk_DeviationMatrix_JASPAR_summarized_experiment.rds")

pseudobulk.atac.chromvar.se <- readRDS(io$atac.pseudobulk.chromvar)[,opts$celltypes]

# Filter
opts$metric <- "z" # "deviations" or "z"
pseudobulk.atac.chromvar.se <- pseudobulk.atac.chromvar.se[rowData(pseudobulk.atac.chromvar.se)$seqnames == opts$metric,]
pseudobulk.atac.chromvar.se <- pseudobulk.atac.chromvar.se[!grepl("var",rowData(pseudobulk.atac.chromvar.se)$name),]

# Rename 
row.ranges.dt <- rowData(pseudobulk.atac.chromvar.se) %>% as.data.table %>% 
  .[,TF:=strsplit(toupper(name),"_") %>% map_chr(c(1))]
rownames(pseudobulk.atac.chromvar.se) <- row.ranges.dt$TF

chromvar.mtx <- assay(pseudobulk.atac.chromvar.se, "DeviationMatrix_JASPAR")
dim(chromvar.mtx)

# Create data.table
pseudobulk.chromvar.dt <- chromvar.mtx %>% as.data.table(keep.rownames = T) %>%
  setnames("rn","gene") %>% melt(id.vars="gene", variable.name="celltype", value.name="chromvar_score")

#########################################
## Load pseudobulk RNA expression data ##
#########################################

pseudobulk.rna.dt <- fread(io$rna.atlas.average_expression_per_celltype) %>%
  .[,gene:=toupper(gene)] %>%
  .[gene%in%rownames(chromvar.mtx)] %>%
  setnames("group","celltype") %>%
  .[,celltype:=gsub(" ","_",celltype)] %>%
  .[celltype%in%opts$celltypes]

rna.mtx <- pseudobulk.rna.dt %>% 
  dcast(gene~celltype,value.var="mean_expr") %>% matrix.please
dim(rna.mtx)

###########
## Merge ##
###########

to.plot <- merge(
  pseudobulk.rna.dt[,c("gene","celltype","mean_expr")],
  pseudobulk.chromvar.dt,
  by = c("celltype","gene")
) %>% 
  melt(id.vars=c("celltype","gene"), variable.name="modality", value.name="value") %>%
  .[,celltype:=factor(celltype,levels=opts$celltypes)]

##########
## Plot ##
##########

genes.to.plot <- c("GATA1","TCF7","SOX18")

# to.plot2 <- to.plot[gene%in%genes.to.plot & modality=="chromvar_score"]

genes.to.plot <- to.plot[modality=="mean_expr",mean(value),by="gene"] %>% .[V1>0.5] %>% .$gene

for (i in genes.to.plot) {
  p <- ggline(to.plot[gene==i], x="celltype", y="value") +
    facet_wrap(~modality, scales = "free_y") +
    geom_hline(yintercept=0, linetype="dashed") +
    scale_color_brewer(palette = "Dark2") +
    labs(x="", y="") +
    guides(x = guide_axis(angle = 90)) +
    theme(
      axis.text = element_text(size=rel(0.6)),
      legend.title = element_blank()
    )
  
  pdf(sprintf("%s/%s.pdf",io$outdir,i), width=6, height=4)
  print(p)
  dev.off()
}
