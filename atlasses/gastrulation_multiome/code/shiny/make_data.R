library(HDF5Array)

################
## Define I/O ##
################

source(here::here("settings.R"))
source(here::here("utils.R"))

io$outdir <- paste0(io$basedir,"/shiny")
opts$motif_annotation <- "Motif_cisbp_lenient"

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE] %>%
  .[,celltype.mapped:=factor(celltype.mapped,levels=opts$celltypes)]

##############################
## Dimensionality reduction ##
##############################

io$umap.rna <- paste0(io$basedir,"/results/rna/dimensionality_reduction/all_cells/E7.5_rep1-E7.5_rep2-E8.0_rep1-E8.0_rep2-E8.5_rep1-E8.5_rep2_umap_features2500_pcs30_neigh25_dist0.3.txt.gz")
umap_rna.dt <- fread(io$umap.rna) %>% setnames(c("UMAP1","UMAP2"),c("V1","V2"))

io$umap.atac <- paste0(io$basedir,"/results/atac/archR/dimensionality_reduction/PeakMatrix/all_cells/E7.5_rep1-E7.5_rep2-E8.0_rep1-E8.0_rep2-E8.5_rep1-E8.5_rep2_umap_nfeatures50000_ndims50_neigh45_dist0.45.txt.gz")
umap_atac.dt <- fread(io$umap.atac) %>% setnames(c("cell","V1","V2"))

# io$umap.mofa <- paste0(io$basedir,"/results/atac/archR/dimensionality_reduction/PeakMatrix/all_cells/E7.5_rep1-E7.5_rep2-E8.0_rep1-E8.0_rep2-E8.5_rep1-E8.5_rep2_umap_nfeatures50000_ndims50_neigh45_dist0.45.txt.gz")
# umap_mofa.dt <- fread(io$umap.mofa) %>% setnames(c("cell","V1","V2"))

# fwrite(umap_rna.dt, paste0(io$outdir,"/umap_rna.txt.gz"))
# fwrite(umap_atac.dt, paste0(io$outdir,"/umap_atac.txt.gz"))

########################
## Load ArchR Project ##
########################

source(here::here("atac/archR/load_archR_project.R"))

##############################
## single-cell RNA and ATAC ##
##############################

opts$rna.cells <- sample_metadata %>% .[pass_rnaQC==TRUE & doublet_call==FALSE & !is.na(celltype.predicted),cell]
opts$atac.cells <- sample_metadata %>% .[pass_atacQC==TRUE & doublet_call==FALSE & !is.na(celltype.predicted),cell]

source(here::here("rna_atac/load_rna_atac_single_cells.R"))
rm(list = c("atac.peakMatrix.se"))

atac.GeneScoreMatrix.se <- readRDS(io$archR.GeneScoreMatrix.se)

rna.sce <- load_SingleCellExperiment(
  file = io$rna.sce, 
  cells = sample_metadata$cell, 
  normalise = TRUE, 
  remove_non_expressed_genes = TRUE
)

genes.to.use <- intersect(rownames(rna.sce), rownames(atac.GeneScoreMatrix.se))
genes.to.use <- genes.to.use[!grepl("*Rik|^Gm|^mt-|^Rps|^Rpl|^Olf",genes.to.use)]
atac.GeneScoreMatrix.se <- atac.GeneScoreMatrix.se[genes.to.use,]
rna.sce <- rna.sce[genes.to.use,]
# write.table(genes.to.use, paste0(io$outdir,"/genes.txt"), row.names = F, col.names = F, quote=F)

## RNA ##

rna_cells <- colnames(rna.sce)
# write.table(rna_cells, paste0(io$outdir,"/cells_rna.txt"), row.names = F, col.names = F, quote=F)
# io$hdf5.outfile = paste0(io$outdir,"/rna_expr.hdf5")
# if(file.exists(io$hdf5.outfile)) { file.remove(io$hdf5.outfile) }
# writeHDF5Array(x = DelayedArray(round(logcounts(rna.sce),2)), file = io$hdf5.outfile, name = "rna_expr_logcounts", verbose = TRUE)

## ATAC ##
atac.GeneScoreMatrix.se <- atac.GeneScoreMatrix.se[,colnames(atac.GeneScoreMatrix.se)%in%opts$atac.cells]
atac_cells <- colnames(atac.GeneScoreMatrix.se)
# write.table(atac_cells, paste0(io$outdir,"/cells_atac.txt"), row.names = F, col.names = F, quote=F)

# Denoise ATAC
# pca.atac <- fread(io$pca.atac) %>% matrix.please %>% .[atac_cells,]
# assay(atac.GeneScoreMatrix.se,paste0(assayNames(atac.GeneScoreMatrix.se)[1],"_denoised")) <- smoother_aggregate_nearest_nb(mat=as.matrix(assay(atac.GeneScoreMatrix.se)), D=pdist(pca.atac), k=25)

# Save ATAC
# io$hdf5.outfile = paste0(io$outdir,"/atac_GeneScoreMatrix.hdf5")
# if(file.exists(io$hdf5.outfile)) { file.remove(io$hdf5.outfile) }
# writeHDF5Array(x = DelayedArray(round(assay(atac.GeneScoreMatrix.se,"denoised"),2)), file = io$hdf5.outfile, name = "atac_GeneScoreMatrix", verbose = TRUE)

#############################
## Pseudobulk RNA and ATAC ##
#############################

io$archR.pseudobulk.deviations.se <- sprintf("%s/results/atac/archR/chromvar/pseudobulk/chromVAR_deviations_summarized_experiment_%s_pseudobulk_correlated_peaks.rds",io$basedir,opts$motif_annotation)

source(here::here("rna_atac/load_rna_atac_pseudobulk.R"))

rna_pseudobulk.dt <- rna_pseudobulk.dt[gene%in%genes.to.use]
atac_gene_scores_pseudobulk.dt <- atac_gene_scores_pseudobulk.dt[gene%in%genes.to.use]
rna_pseudobulk.sce <- rna_pseudobulk.sce[genes.to.use,]
atac_pseudobulk_GeneScoreMatrix.se <- atac_pseudobulk_GeneScoreMatrix.se[genes.to.use,]

rna_chromvar_pseudobulk.dt <- merge(
  rna_tf_pseudobulk.dt,
  atac_chromvar_pseudobulk.dt,
  by = c("celltype","gene")
) %>% .[,c("expr","chromvar_zscore"):=list(round(expr,1), round(chromvar_zscore,1))]
# fwrite(rna_chromvar.dt, paste0(io$outdir,"/rna_vs_chromvar_pseudobulk.txt.gz"))

TFs <- unique(rna_chromvar_pseudobulk.dt$gene)
# write.table(TFs, paste0(io$outdir,"/TFs.txt"), row.names = F, col.names = F, quote=F)

atac_gene_scores_pseudobulk.dt <- atac_gene_scores_pseudobulk.dt[gene%in%genes.to.use]

gene_rna_atac_pseudobulk.dt <- merge(
  rna_pseudobulk.dt,
  atac_gene_scores_pseudobulk.dt,
  by = c("celltype","gene")
) %>% .[,c("expr","acc"):=list(round(expr,1), round(acc,1))]
fwrite(gene_rna_atac_pseudobulk.dt, paste0(io$outdir,"/gene_rna_vs_acc_pseudobulk.txt.gz"))

##########
## PAGA ##
##########

source(here::here("load_paga_graph.R"))

saveRDS(net.paga, paste0(io$outdir,"/paga_network.rds"))

########################
## In silico ChIP-seq ##
########################

TFs <- list.files(io$virtual_chip.dir, ".bed.gz") %>% str_replace_all(".bed.gz","")

# TFs <- c("T", "ZIC2", "TAL1", "GATA1", "FOXA2", "GATA4", "CDX2")

virtual_chip.dt <- TFs %>% map(function(i) {
  fread(sprintf("%s/%s.txt.gz",io$virtual_chip.dir,i)) %>%
    setnames(c("chr","start","end","score","correlation_score","max_accessibility_score","motif_score","motif_counts")) %>%
    .[,idx:=sprintf("%s:%s-%s",chr,start,end)] %>%
    .[,c("chr","start","end"):=NULL] %>%
    .[,tf:=i] %>%
    return
}) %>% rbindlist

seq.ranges <- seq(0,0.75,by=0.01)

to.plot <- seq.ranges %>% map(function(i) {
  virtual_chip.dt[score>=i,log2(.N),by="tf"] %>% .[,min_score:=i] %>% return
}) %>% rbindlist %>% setnames("V1","log2_N") %>% .[,log2_N:=round(log2_N,2)]
fwrite(to.plot, paste0(io$outdir,"/insilico_chip_stats.txt.gz"))


#################################
## Differential RNA expression ##
#################################

# single-cell and pseudobulk

##########################################
## Differential chromatin accessibility ##
##########################################

# single-cell and pseudobulk

################################
## Differential TF activities ##
################################

# single-cell and pseudobulk


################################################
## RNA vs chromVAR per cell type (pseudobulk) ##
################################################

rna_chromvar.dt <- merge(
  rna_tf.dt,
  chromvar.dt,
  by = c("celltype","gene")
)

length(unique(rna_chromvar.dt$gene))
# fwrite(rna_chromvar.dt, paste0(io$outdir,"/rna_vs_chromvar_pseudobulk.txt.gz"))

############################################
## Gene regulatory networks per cell type ##
############################################








###########
## IGNORE BELOW ##
###########

save(tsnes,
     umaps,
     meta,
     genes,
     markers_celltype,
     file = "/nfs/research1/marioni/jonny/embryos/scripts/shiny/data.RData")


saveRDS(genes, file = "/nfs/research1/marioni/jonny/embryos/scripts/shiny/genes.rds")
saveRDS(meta, file = "/nfs/research1/marioni/jonny/embryos/scripts/shiny/meta.rds")

