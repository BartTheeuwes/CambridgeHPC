#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

io$diff_expr_tf <- file.path(io$basedir,"results/rna/differential/celltype.mapped/TFs/diff_expr_TFs.txt.gz")
# io$differential_results <- file.path(io$basedir,"results/rna/differential/celltype.mapped/TFs")
io$outdir <- file.path(io$basedir,"results/rna/differential/celltype.mapped/TFs/marker_TFs"); dir.create(io$outdir, showWarnings = F)

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

##########################################
## Load differential expression results ##
##########################################

# # i <- "Gut"; j <- "NMP"
# diff.dt <- opts$celltypes %>% map(function(i) { opts$celltypes %>% map(function(j) {
#   file <- file.path(io$differential_results,sprintf("%s_vs_%s.txt.gz",i,j))
#   if (file.exists(file)) {
#     fread(file, select = c(1,2,4)) %>% 
#       # .[abs(logFC)>=opts$min.log2FC & padj_fdr<=opts$fdr] %>%
#       .[,sig:=abs(logFC)>=opts$min.log2FC & padj_fdr<=opts$fdr] %>%
#       .[,c("celltypeA","celltypeB"):=list(i,j)] %>%
#       return
#   } }) %>% rbindlist }) %>% rbindlist %>%
#   .[,celltypeA:=factor(celltypeA,levels=opts$celltypes)] %>%
#   .[,celltypeB:=factor(celltypeB,levels=opts$celltypes)] %>%
#   .[,direction:=as.factor(c("down","up"))[as.numeric(logFC<0)+1]]
# 
# ncelltypes <- unique(c(as.character(unique(diff.dt$celltypeA)),as.character(unique(diff.dt$celltypeB)))) %>% length
# 
# diff.dt[,.N,by="gene"]

##########################################
## Load differential expression results ##
##########################################

# i <- "Gut"; j <- "NMP"
diff.dt <- fread(io$diff_expr_tf) %>%
  .[,sig:=abs(logFC)>=opts$min.log2FC & padj_fdr<=opts$fdr] %>%
  .[is.na(sig),sig:=FALSE] %>%
  .[,direction:=as.factor(c("down","up"))[as.numeric(logFC<0)+1]]

ncelltypes <- unique(c(as.character(unique(diff.dt$celltypeA)),as.character(unique(diff.dt$celltypeB)))) %>% length

# diff.dt[,.N,by="gene"] %>% View

#########################
## Define marker genes ##
#########################

# This only takes into account upregulated markers
foo <- diff.dt[,.(score=sum(direction=="up" & sig==TRUE,na.rm=T)), by=c("celltypeA","gene")] %>% setnames("celltypeA","celltype")
bar <- diff.dt[,.(score=sum(direction=="down" & sig==TRUE,na.rm=T)), by=c("celltypeB","gene")] %>% setnames("celltypeB","celltype")

markers_genes.dt <- merge(foo,bar,by=c("celltype","gene"), all=TRUE) %>% 
  .[is.na(score.x),score.x:=0] %>% .[is.na(score.y),score.y:=0] %>%
  .[,score:=score.x+score.y] %>%
  .[,c("score.x","score.y"):=NULL] %>%
  .[,score:=round(score/(ncelltypes-1),2)] %>%
  setorder(celltype,-score)
rm(foo,bar)

# diff.dt[gene=="MESP2"] %>% View
# foo[gene=="MESP2"]
# bar[gene=="MESP2"]

###################################################
## Add average logFC values from pseudobulk data ##
###################################################

# io$diff_expr_tf_pseudobulk <- file.path(io$basedir,"results/rna/differential/celltype.mapped/pseudobulk/TFs/diff_expr_TFs_pseudobuk.txt.gz")
# 
# diff_expr_pseudobulk.dt <- fread(io$diff_expr_tf_pseudobulk)
# 
# foo <- diff_expr_pseudobulk.dt[,.(diff=mean(diff), N=.N),by=c("gene","celltypeA")] %>% setnames("celltypeA","celltype")
# bar <- diff_expr_pseudobulk.dt[,.(diff=mean(diff), N=.N),by=c("gene","celltypeB")] %>% setnames("celltypeB","celltype")
# foobar <- merge(foo,bar,by=c("gene","celltype"))
# 
# markers_genes.dt <- markers_genes.dt %>% merge(diff_expr_pseudobulk.dt, by=c("celltype","gene"))

##########
## Save ##
##########
  
# Save marker score for all combination of TFs and cell types
length(unique(markers_genes.dt$gene))
length(unique(markers_genes.dt$celltype))
markers_genes.dt[,.N,by="celltype"]
markers_genes.dt[,.N,by="gene"]
fwrite(markers_genes.dt, file.path(io$outdir,"marker_TFs_upregulated_all.txt.gz"), sep="\t")

# Save marker score for strong markers
markers_genes_filt.dt <- markers_genes.dt %>% .[score>=opts$score]
length(unique(markers_genes_filt.dt$gene))
length(unique(markers_genes_filt.dt$celltype))
markers_genes_filt.dt[,.N,by="celltype"]
markers_genes_filt.dt[,.N,by="gene"]
fwrite(markers_genes_filt.dt, file.path(io$outdir,"marker_TFs_upregulated_filtered.txt.gz"), sep="\t")

########################################################
## Plot score vs number of marker genes per cell type ##
#########################3##############################

seq.ranges <- seq(0.10,1,by=0.10)
names(seq.ranges) <- as.character(1:length(seq.ranges))

to.plot <- seq.ranges %>% map(function(i) {
  markers_genes.dt %>% .[score>=i] %>% .[,.N,by=c("celltype")] %>% .[,min_score:=i] %>% return
}) %>% rbindlist

p <- ggline(to.plot, x="min_score", y="N", color="celltype", plot_type = "l") +
  scale_color_manual(values=opts$celltype.colors) +
  labs(x="Minimum marker score", y="Number of marker TFs") +
  theme(
    axis.text = element_text(size=rel(0.65)),
    legend.position = "none"
  )

pdf(file.path(io$outdir,"lineplot_number_marker_TFs.pdf"), width = 7, height = 4)
print(p)
dev.off()

################################################
## Plot number of marker genes per cell types ##
################################################

to.plot <- markers_genes_filt.dt %>% .[,.N,by=c("celltype")]

p <- ggbarplot(to.plot, x="celltype", y="N", fill="celltype") +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(x="", y="Number of marker TFs") +
  theme(
    axis.text.y = element_text(size=rel(0.65)),
    axis.text.x = element_text(colour="black",size=rel(0.7), angle=90, hjust=1, vjust=0.5),
    axis.title = element_text(colour="black",size=rel(0.75)),
    axis.ticks.x = element_blank(),
    legend.position = "none"
)

pdf(file.path(io$outdir,"barplot_number_marker_TFs.pdf"), width = 8, height = 4)
print(p)
dev.off()

##################################
## Plot gene marker exclusivity ##
##################################

to.plot <- markers_genes_filt.dt %>%
  .[,.(Nx=.N),by="gene"] %>%
  .[,Nx:=factor(Nx)] %>%
  .[,.(Ny=.N),by="Nx"]

p <- ggbarplot(to.plot, x="Nx", y="Ny", fill="gray70") +
  labs(x="Number of different cell types per marker TF", y="") +
  theme(
    axis.text = element_text(size=rel(0.75)),
  )
pdf(file.path(io$outdir,"boxplot_exclusivity_per_TF.pdf"), width = 7, height = 5)
print(p)
dev.off()

################################################
## Plot gene marker exclusivity per cell type ##
################################################

to.plot <- markers_genes_filt.dt %>% .[,N:=.N,by="gene"]

p <- ggboxplot(to.plot, x="celltype", y="N", fill="celltype", color="black") +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(x="", y="Exclusivity of TF markers\n(the smaller the more exclusive)") +
  theme(
    axis.text.y = element_text(size=rel(0.75)),
    axis.title.y = element_text(size=rel(0.85)),
    axis.text.x = element_text(colour="black",size=rel(0.7), angle=90, hjust=1, vjust=0.5),
    legend.position = "none"
  )

pdf(file.path(io$outdir,"boxplot_exclusivity_per_celltype.pdf"), width = 9, height = 5)
print(p)
dev.off()


