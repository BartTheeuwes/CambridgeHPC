
library(SingleCellExperiment)
library(scran)
library(scater)
library(TFBSTools)
library(glmnet)
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
# io$outdir <- paste0(io$basedir, "/results/atac/signac/dimensionality_reduction")
io$mofa.output <- paste0(io$basedir, "/results/rna_atac/mofa/mofa_model.rds")

# Options
opts$remove.small.celltypes <- TRUE
opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E8.5_rep1",
  "E8.5_rep2"
)

opts$remove.small.celltypes <- TRUE

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & pass_rnaQC==TRUE] %>%
  .[sample%in%opts$samples]

# subset celltypes with sufficient number of cells
if (opts$remove.small.celltypes) {
  opts$min.cells <- 100
  sample_metadata <- sample_metadata %>%
    .[,N:=.N,by=c("celltype.predicted")] %>% .[N>opts$min.cells] %>% .[,N:=NULL]
}
opts$celltypes <- unique(sample_metadata$celltype.predicted)

###############################
## Load single-cell ATAC data ##
###############################

# if (grepl("ricard",Sys.info()['nodename'])) {
#   source("/Users/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
# } else if (grepl("ebi",Sys.info()['nodename'])) {
#   source("/homes/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
# } else {
#   stop("Computer not recognised")
# }
# stopifnot(sample_metadata$archR_cell %in% rownames(ArchRProject))
# ArchRProject.filt <- ArchRProject[sample_metadata$archR_cell]

# Get Peak Matrix from ArchR
# atac.peak.se <- getMatrixFromProject(ArchRProject.filt, useMatrix="PeakMatrix")
# dim(atac.peak.se)

# Rename cells
# colnames(atac.peak.se) <- sample_metadata %>% 
#   .[archR_cell%in%colnames(atac.peak.se)] %>%
#   setkey(archR_cell) %>% .[colnames(atac.peak.se)] %>% .$cell

# Define peak names
# peak_names <- rowRanges(atac.peak.se) %>% as.data.table %>% .[,id:=sprintf("%s_%s_%s",seqnames,start,end)] %>% .$id
# rownames(atac.peak.se) <- peak_names

###############################
## Load pseudobulk ATAC data ##
###############################

io$atac.pseudobulk.peaks <- paste0(io$archR.directory,"/pseudobulk/pseudobulk_PeakMatrix_summarized_experiment.rds")

pseudobulk.atac.peak.se <- readRDS(io$atac.pseudobulk.peaks)[,opts$celltypes]
row.ranges.dt <- rowData(pseudobulk.atac.peak.se) %>% as.data.table %>% 
  setnames("seqnames","chr") %>%
  .[,idx:=sprintf("%s_%s_%s",chr,start,end)]
rownames(pseudobulk.atac.peak.se) <- row.ranges.dt$idx

##################################
## Select highly variable peaks ##
##################################

opts$npeaks <- 1000

# calculate variability per peak
# note: use pseudobulk data, not single-cell data
peak.var <- data.table(
  id = rownames(pseudobulk.atac.peak.se),
  # var = assay(atac.peak.se,"PeakMatrix") %>% sparseMatrixStats::rowVars(.)
  var = assay(pseudobulk.atac.peak.se,"PeakMatrix") %>% sparseMatrixStats::rowVars(.)
) %>% setorder(-var)

# Extract highly variable peaks (use pseudobulk data)
peak.mtx <- pseudobulk.atac.peak.se %>% .[head(peak.var$id,opts$npeaks),] %>% assay(.,"PeakMatrix")
colnames(peak.mtx) <- colnames(pseudobulk.atac.peak.se)

# Binarise
# peak.mtx[peak.mtx>1] <- 1

#########################
## Load TF information ##
#########################

# Load
io$atac.peak.annotation <- sprintf("%s/Annotations/JASPAR/peakAnnotation.rds",io$archR.directory)
tf.dt <- readRDS(io$atac.peak.annotation)[["Motif_JASPAR"]][["motifSummary"]] %>%
  as.data.table

# Filter
tf.dt <- tf.dt %>% .[grep("::",name, invert = T)] %>% .[grep("(var.+)",name, invert = T)]
# unique(tf.dt$name)

# Rename
tf.dt[name=="TBXT",name:="T"]
tf.dt[name=="RAX2",name:="ZFP58"]
tf.dt[name=="VENTX",name:="OBOX5"]

# tf.dt$name[grep("FOXF",tf.dt$name)]

##############################
## Load RNA expression data ##
##############################

# Load SingleCellExperiment object
# sce <- load_SingleCellExperiment(io$sce, normalise = TRUE, cells = sample_metadata$cell)

# Add sample metadata to the colData of the SingleCellExperiment
# colData(sce) <- sample_metadata %>% as.data.frame %>% tibble::column_to_rownames("cell") %>%
#   .[colnames(sce),] %>% DataFrame()

# Filter
# rownames(sce) <- toupper(rownames(sce))
# sce_filt <- sce[rownames(sce)%in%tf.dt$name,]

# Create matrix
# rna.mtx <- assay(sce_filt,"logcounts")
# dim(rna.mtx)


#########################################
## Load pseudobulk RNA expression data ##
#########################################

pseudobulk.rna.dt <- fread(io$rna.atlas.average_expression_per_celltype) %>%
  .[,gene:=toupper(gene)] %>%
  setnames("group","celltype") %>%
  .[,celltype:=gsub(" ","_",celltype)] %>%
  .[gene%in%tf.dt$name & celltype%in%opts$celltypes]

rna.mtx <- pseudobulk.rna.dt %>% 
  dcast(gene~celltype,value.var="mean_expr") %>% matrix.please
dim(rna.mtx)

#########################
## Subset common cells ##
#########################

samples <- intersect(colnames(rna.mtx),colnames(peak.mtx))
rna.mtx <- rna.mtx[,samples]
peak.mtx <- peak.mtx[,samples]

#############################
## Summarise per cell type ##
#############################

# rna.mtx2 <- rna.mtx
# colnames(rna.mtx2) <- sample_metadata[cell%in%colnames(rna.mtx)] %>% setkey(cell) %>% .[colnames(rna.mtx)] %>% .$celltype.predicted
# rna.mtx2 <- rna.mtx2 %>% as.matrix %>% t %>% 
#   as.data.table(keep.rownames = TRUE) %>% .[, lapply(.SD, mean), by = "rn"]

# peak.mtx2 <- peak.mtx
# colnames(peak.mtx2) <- sample_metadata[cell%in%colnames(peak.mtx)] %>% setkey(cell) %>% .[colnames(peak.mtx)] %>% .$celltype.predicted
# peak.mtx2 <- peak.mtx2 %>% as.matrix %>% t %>% 
#   as.data.table(keep.rownames = TRUE) %>% .[, lapply(.SD, mean), by = "rn"]

####################################
## Simple correlation coefficient ##
####################################

cor.dt <- data.table(
  TF = rownames(rna.mtx),
  cor(peak.mtx[i,],t(rna.mtx))[1,]
)

plot(x=rna.mtx["GATA1",], y=peak.mtx[i,])

############################################
## Fit penalised linear regression model per peak ##
############################################

stopifnot(colnames(rna.mtx)==colnames(peak.mtx))

for (i in rownames(peak.mtx)) {
  
  # dt <- merge(peak.mtx[,c(i,"rn"),with=F], rna.mtx, by="celltype") %>% .[,celltype:=NULL]
  # y <- dt[,i,with=F][[1]]
  # x <- dt[,-i,with=F] %>% as.matrix
  
  x <- rna.mtx %>% t
  y <- peak.mtx[i,]
  
  # Fit standard linear model
  # df <- data.frame(cbind(y,x))
  # lm <- lm(y~., data=df)
  
  # Fit penalised linear regression model
  x.test <- x[,c("GATA1","CLOCK","CREB1")]
  
  glmnet.fit <- glmnet(x.test, y, family="gaussian", alpha=0, lower.limits=0, standardize = F)
  plot(glmnet.fit)
  
  cvfit <- cv.glmnet(x, y, alpha=0)
  plot(cvfit)
  
  # r2
  # r2 <- glmnet.fit$dev.ratio[which(cvfit$glmnet.fit$lambda == cvfit$lambda.1se)] 
  r2 <- glmnet.fit$dev.ratio[which(cvfit$glmnet.fit$lambda == cvfit$lambda.min)]
  
  # get coefficients
  weights.dt <- data.table(
    weight = coef(glmnet.fit, s=cvfit$lambda.min) %>% as.numeric,
    TF = rownames(coef(glmnet.fit))
  )
  
  # plot predicted vs true
  to.plot <- data.table(
    predicted_lambda_min = predict(glmnet.fit, newx = x.test, s = cvfit$lambda.min, type = "response")[,1],
    predicted_lambda_1se = predict(glmnet.fit, newx = x.test, s = cvfit$lambda.1se, type = "response")[,1],
    true = y,
    celltype = names(y)
  ) 
  
  
  ggplot(to.plot, aes(x=predicted_lambda_min,y=true)) +
    geom_point(aes(color=celltype), size=2.5) +
    scale_color_manual(values=opts$celltype.colors) +
    stat_cor(method = "pearson") +
    geom_abline(intercept=0, slope=1) +
    # stat_smooth(method="lm", color="black", alpha=0.25, size=0.2) +
    ggrepel::geom_text_repel(aes(label=celltype), size=3, data=to.plot[true>0.5 ]) +
    labs(x="Predicted accessibiltiy", y="True accessibility") +
    theme_classic() + 
    theme(
      axis.text = element_text(color="black"),
      legend.position = "none"
    )
}


##########
## Test ##
##########
# 
# sort(assay(pseudobulk.atac.peak.se)["chr11_97187762_97188362",])
# # pseudobulk.peak.var <- rowVars(assay(pseudobulk.atac.peak.se))
# 
# 
# foo <- data.table(
#   id = rownames(pseudobulk.atac.peak.se),
#   var = assay(pseudobulk.atac.peak.se,"PeakMatrix") %>% sparseMatrixStats::rowVars(.)
# ) 
# 
# bar <- data.table(
#   id = peak_names,
#   var = assay(atac.peak.se,"PeakMatrix") %>% sparseMatrixStats::rowVars(.)
# ) 
# 
# foobar <- merge(foo,bar,by="id")
# foobar %>% setnames(c("peak","variance_pseudobulk","variance_singlecell"))
# foobar[var.x>0.01 | var.y>0.01]
# 
# foobar[,c("rank_pseudobulk","rank_singlecell"):=list(rank(-variance_pseudobulk), rank(-variance_singlecell))]
# # plot(foobar$var.x,foobar$var.y)
# 
# cor(foobar$variance_pseudobulk,foobar$variance_singlecell)
# foobar[peak=="chr5_119483908_119484508"]
# 
# ggplot(foobar, mapping = aes(x = variance_pseudobulk, y = variance_singlecell)) +
#   ggpointdensity::geom_pointdensity() +
#   viridis::scale_color_viridis() +
#   theme_classic()
# 
# 
