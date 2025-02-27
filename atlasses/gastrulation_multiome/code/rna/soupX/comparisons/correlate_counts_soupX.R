library(ggpointdensity)

#####################
## Define settings ##
#####################

# Load default settings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
}

io$rna.sce.soupX <- paste0(io$basedir,"/processed/rna_soupX/SingleCellExperiment.rds")
io$outdir <- paste0(io$basedir,"/results/rna/soupX/correlations")


##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE]

io$metadata <- paste0(io$basedir,"/results/rna_soupX/qc/sample_metadata_after_qc.txt.gz")
sample_metadata.soupX <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE]

cells <- intersect(sample_metadata$cell,sample_metadata.soupX$cell)

###############################
## Load SingleCellExperiment ##
###############################

sce <- load_SingleCellExperiment(
  file = io$rna.sce, 
  normalise = TRUE, 
  cells = cells,
  remove_non_expressed_genes = TRUE
)
dim(sce)

sce.soupX <- load_SingleCellExperiment(
  file = io$rna.sce.soupX, 
  normalise = TRUE, 
  cells = cells, 
  remove_non_expressed_genes = TRUE
)
dim(sce.soupX)

############################################
## Load soupX feature abundance estimates ##
############################################

genes.to.correlate <- intersect(rownames(sce),rownames(sce.soupX))

soup.dt <- opts$samples %>% map(function(i) {
  file <- sprintf("%s/original/%s/filtered_feature_bc_matrix/soupX/soup.tsv.gz",io$basedir,i)
  fread(file) %>% .[,sample:=i]
}) %>% rbindlist

###############################
## Correlate counts per gene ##
###############################

r <- genes.to.correlate %>% head %>% map(function(i) {
  cor(as.numeric(logcounts(sce[i,])),as.numeric(logcounts(sce.soupX[i,])))
}) %>% unlist

######################################
## Plot genes that do not correlate ##
######################################

##################################################
## Plot genes that are not enriched in the soup ##
##################################################

##############################################
## Plot genes that are enriched in the soup ##
##############################################

genes.to.plot <- soup.dt %>%
  .[,sum(est),by="gene"] %>% setorder(-V1) %>% .$gene %>% head(n=15)

for (i in genes.to.plot) {
  
  # to.plot <- rbind(
  #   data.table(cell=colnames(sce), expr=as.numeric(logcounts(sce[i,])), class="standard"),
  #   data.table(cell=colnames(sce.soupX), expr=as.numeric(logcounts(sce.soupX[i,])), class="soupX")
  # ) %>% dcast(cell~class, value.var="expr")
  to.plot <- data.table(
      cell = colnames(sce), 
      standard_expr = as.numeric(counts(sce[i,])),
      soupX_expr = as.numeric(counts(sce.soupX[i,]))
    ) %>% .[standard_expr>0 | soupX_expr>0]
  
  p <- ggplot(to.plot, aes(x = standard_expr, y = soupX_expr)) +
    # geom_pointdensity() + viridis::scale_color_viridis() +
    geom_point(size=1.5) +
    stat_smooth(method="lm") +
    geom_abline(slope=1,intercept=0, color="gray") +
    labs(x=sprintf("%s expr (standard)",i), y=sprintf("%s expr (soupX)",i)) +
    theme_bw() 
  
  # pdf(sprintf("%s/%s.pdf",io$outdir,i), width=7, height=7)
  png(sprintf("%s/%s.png",io$outdir,i), width=600, height=600)
  print(p)
  dev.off()
}

