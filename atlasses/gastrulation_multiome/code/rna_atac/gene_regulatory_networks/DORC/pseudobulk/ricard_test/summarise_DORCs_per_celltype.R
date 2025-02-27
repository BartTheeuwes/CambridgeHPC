
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

# Options
opts$motif_annotation <- "Motif_cisbp"

# I/O
io$outdir <- paste0(io$basedir,"/results/rna_atac/DORCs/pseudobulk")
io$archR.pseudobulk.deviations.se <- sprintf("%s/pseudobulk/pseudobulk_DeviationMatrix_%s_summarized_experiment.rds",io$archR.directory,opts$motif_annotation)

#######################################
## Load pseudobulk RNA and ATAC data ##
#######################################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/rna_atac/rna_vs_acc/pseudobulk/load_rna_atac_pseudobulk.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/rna_atac/rna_vs_acc/pseudobulk/load_rna_atac_pseudobulk.R")
} else {
  stop("Computer not recognised")
}

rna.sce <- rna.sce[,opts$celltypes]
rna.sce.tf <- rna.sce.tf[,opts$celltypes]
atac.chromvar.se <- atac.chromvar.se[,opts$celltypes]
atac.peakMatrix.se <- atac.peakMatrix.se[,opts$celltypes]

####################################
## Load pre-computed DORC results ##
####################################

DORCs.granges <- readRDS("/Users/ricard/data/gastrulation_multiome_10x/results/rna_atac/DORCs/pseudobulk/DORCs_samples_distance10000.rds")

##########
## Test ##
##########

genes <- unique(DORCs.granges$gene) %>% as.character
DORCs.mtx <- matrix(NA, nrow=length(genes), ncol=length(opts$celltypes))
rownames(DORCs.mtx) <- genes; colnames(DORCs.mtx) <- opts$celltypes

for (i in genes) {
  peaks <- DORCs.granges[DORCs.granges$gene==i]$peak %>% as.character
  DORCs.mtx[i,] <- colMeans(assay(atac.peakMatrix.se[peaks,]))
}

# save
saveRDS(DORCs.mtx, paste0(io$outdir,"/DORCs_score_matrix.rds"))

###################################################
## Plot RNA expression versus DORC accessibility ##
###################################################

genes.to.plot <- DORCs.granges %>% as.data.table() %>% .[,.N,by=c("gene")] %>% .[N>=5,gene] %>% as.character
genes.to.plot <- genes# %>% head(n=5)
# genes.to.plot <- grep("Hox",genes,value=T)

for (i in genes.to.plot) {
  outfile <- sprintf("%s/per_gene/%s_rna_vs_DORC_acc_pseudobulk.pdf",io$outdir,i)
  if (!file.exists(outfile)) {
    
    to.plot <- data.table(
      acc = DORCs.mtx[i,] %>% minmax.normalisation(),
      expr = logcounts(rna.sce[i,])[1,] %>% minmax.normalisation(),
      celltype = opts$celltypes
    )
    
    p <- ggplot(to.plot, aes(x=acc, y=expr)) +
      ggrastr::geom_point_rast(aes(fill=celltype), color="black", shape=21, size=3, stroke=0.5) +
      # geom_smooth(method="lm", color="black", fill="black", alpha=0.25) + stat_cor(method = "pearson") +
      scale_fill_manual(values=opts$celltype.colors) +
      # coord_cartesian(xlim=c(opts$min.expr,opts$max.expr+0.01), ylim=c(opts$min.acc,opts$max.acc+0.01)) +
      # ggrepel::geom_text_repel(aes(label=celltype), size=3, data=to.plot[expr>2 | atac>2]) +
      labs(x="DORC accessibility", y="RNA expression", title=i) +
      theme_classic() + 
      theme(
        plot.title = element_text(hjust=0.5),
        axis.text = element_text(color="black"),
        legend.position = "none"
      )
    
    pdf(outfile, width = 6, height = 5)
    # png(sprintf("%s/per_gene/%s_acc_vs_rna_pseudobulk.png",io$outdir,i), width = 500, height = 400)
    print(p)
    dev.off()
  }
}


DORCs.granges[DORCs.granges$gene=="Morn5"]
