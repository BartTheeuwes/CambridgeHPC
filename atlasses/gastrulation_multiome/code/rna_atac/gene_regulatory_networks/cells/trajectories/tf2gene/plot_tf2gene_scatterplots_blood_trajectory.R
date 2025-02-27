###############################
## Plot TF-gene scatterplots ##
###############################

tf.to.plot <- "GATA1"
gene.to.plot <- "Tmem9"

# tail(sort(tmp[tf.to.plot,][!is.na(tmp[tf.to.plot,])]), n=10)
# head(sort(tmp[tf.to.plot,][!is.na(tmp[tf.to.plot,])]), n=10)

# target_peaks_i <- names(which(virtual_chip.mtx[,tf.to.plot]>=opts$min.chip.threshold))
# peak2gene.dt[gene==gene.to.plot & peak%in%target_peaks_i]
corr.mtx.filt$r[tf.to.plot,gene.to.plot]

to.plot <- data.table(
  cell = colnames(rna_tf.sce),
  tf = logcounts(rna_tf.sce[tf.to.plot,])[1,],
  gene = logcounts(rna_genes.sce[gene.to.plot,])[1,]
) %>% merge(sample_metadata[,c("cell","celltype.predicted")], by="cell")


ggscatter(to.plot, x="tf", y="gene", color="celltype.predicted", size=1,
          add="reg.line", add.params = list(color="gray60"), conf.int=FALSE) + 
  stat_cor(method = "pearson") +
  scale_color_manual(values=opts$celltype.colors) + 
  theme(
    legend.position="none"
  )
