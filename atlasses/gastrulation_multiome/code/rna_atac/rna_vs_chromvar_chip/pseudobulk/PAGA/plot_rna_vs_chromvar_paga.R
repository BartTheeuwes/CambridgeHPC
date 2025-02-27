
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

# I/O
io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/pseudobulk/per_gene/paga")

# Options
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

opts$motif_annotation <- "Motif_cisbp"

################################################
## Load pseudobulk RNA and chromVAR estimates ##
################################################

io$archR.pseudobulk.deviations.se <- sprintf("%s/results/atac/archR/chromvar/pseudobulk/chromVAR_deviations_summarized_experiment_%s_pseudobulk_correlated_peaks_archr.rds",io$basedir,opts$motif_annotation)

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/rna_atac/load_rna_atac_pseudobulk.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/rna_atac/load_rna_atac_pseudobulk.R")
} else {
  stop("Computer not recognised")
}

################
## Parse data ##
################

# Filter genes with low variability
genes.to.keep.rna <- rna_tf.dt[,.(var(expr)),by="gene"] %>% .[V1>0.1,gene] %>% as.character
genes.to.keep.chromvar <- chromvar.dt[,.(var(chromvar_zscore)),by="gene"] %>% .[V1>0.1,gene] %>% as.character
genes.to.keep <- intersect(genes.to.keep.rna,genes.to.keep.chromvar)
rna_tf.dt <- rna_tf.dt[gene%in%genes.to.keep]
chromvar.dt <- chromvar.dt[gene%in%genes.to.keep]

###########
## Merge ##
###########

rna_chromvar.dt <- merge(
  rna_tf.dt,
  chromvar.dt,
  by = c("celltype","gene")
)

#####################
## Load PAGA graph ##
#####################

opts$celltypes <- unique(rna_chromvar.dt$celltype)

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/load_paga_graph.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/load_paga_graph.R")
} else {
  stop("Computer not recognised")
}


##############################################
## Plot network, colour by cell type labels ##
##############################################

# p <- ggnet2(
#   net = net, 
#   mode = c("x", "y"), 
#   color = opts$celltype.colors[celltypes],
#   node.size = 6, 
#   edge.size = 0.15, 
#   edge.color = "grey", 
#   label = TRUE, 
#   label.size = 2.3
# )
# 
# pdf(paste0(io$outdir,"/paga_coloured_by_celltype.pdf"), width=6, height=4)
# print(p)
# dev.off()

#####################################################################
## Plot network, colour by gene expression and motif accessibility ##
#####################################################################

opts$scale <- TRUE

if (opts$scale) {
  rna.col.seq <- chromvar.col.seq <- round(seq(0,1,0.1), 2)
  rna_chromvar.dt[,expr:=minmax.normalisation(expr), by="gene"]
  rna_chromvar.dt[,chromvar_zscore:=minmax.normalisation(chromvar_zscore), by="gene"]
} else {
  opts$max.expr <- 9; opts$min.expr <- 2
  rna.col.seq <- round(seq(opts$min.expr,opts$max.expr,0.1), 2)
  rna_chromvar.dt[expr>=opts$max.expr,expr:=opts$max.expr]
  rna_chromvar.dt[expr<=opts$min.expr,expr:=opts$min.expr]
  
  opts$max.chromvar <- 6; opts$min.chromvar <- 0
  chromvar.col.seq <- round(seq(opts$min.chromvar,opts$max.chromvar,0.1), 2)
  rna_chromvar.dt[chromvar_zscore>=opts$max.chromvar,chromvar_zscore:=opts$max.chromvar]
  rna_chromvar.dt[chromvar_zscore<=opts$min.chromvar,chromvar_zscore:=opts$min.chromvar]
  # activity.col.seq <- round(seq(0,max(rna_chromvar.dt$activity,na.rm=T),0.05), 2)
}

# Define colors
rna.colors <- colorRampPalette(c("gray92", "darkgreen"))(length(rna.col.seq))
chromvar.colors <- colorRampPalette(c("gray92", "purple"))(length(chromvar.col.seq)) 
# activity.colors <- colorRampPalette(c("gray92", "#EE9A00"))(length(activity.col.seq)) 

# Define genes to plot
genes.to.plot <- unique(rna_chromvar.dt$gene)# %>% head(n=10)
# genes.to.plot <- c("TAL1","T","GLIS2","STAT3")

# genes.to.plot <- "Prdm1"

p <- ggnet2(
  net = net.paga,
  mode = c("x", "y"),
  node.size = 0,
  edge.size = 0.15,
  edge.color = "grey",
  label = FALSE,
  label.size = 2.3
)

for (i in genes.to.plot) {
  
  if (opts$scale) {
    outfile <- sprintf("%s/individual_genes/%s_paga_rna_vs_chromvar_correlated_peaks_scaled.png",io$outdir,i)
  } else {
    outfile <- sprintf("%s/individual_genes/%s_paga_rna_vs_chromvar_correlated_peaks.png",io$outdir,i)
  }
  
  if (file.exists(outfile)) {
    print(sprintf("file or %s already exists...",i))
  } else {
    
    expr.values <- rna_chromvar.dt[gene==i,c("celltype","expr")] %>% matrix.please %>% .[opts$celltypes,]
    expr.colors <- round(expr.values,1) %>% map(~ rna.colors[which(rna.col.seq == .)]) %>% unlist
  
    p1 <- p + geom_text(label = "\u25D0", aes(x=x, y=y), color=expr.colors, size=20, family = "Arial Unicode MS",
                  data = p$data[,c("x","y")] %>% dplyr::mutate(expr=expr.colors)) +
      scale_colour_manual(values=expr.colors) + 
      labs(title="RNA expression") +
      theme(
        plot.title = element_text(hjust = 0.5)
      )
  
    # motif accessibility
    acc.values <- rna_chromvar.dt[gene==i,c("celltype","chromvar_zscore")] %>% matrix.please %>% .[,1]# %>% .[opts$celltypes,]
    acc.colors <- round(acc.values,1) %>% map(~ chromvar.colors[which(chromvar.col.seq == .)]) %>% unlist
    
    p2 <- p + geom_text(label = "\u25D1", aes(x=x, y=y), color=acc.colors, size=20, family = "Arial Unicode MS",
                  data = p$data[,c("x","y")] %>% dplyr::mutate(acc=acc.colors)) +
      scale_fill_manual(values=acc.colors) + 
      labs(title="Motif accessibility") +
      theme(
        plot.title = element_text(hjust = 0.5)
      )
    
    # TF activity
    # activity.values <- rna_chromvar.dt[gene==i,c("celltype","activity")] %>% matrix.please %>% .[,1]# %>% .[opts$celltypes,]
    # activity_colors.i <- round(activity.values,1) %>% map(~ activity.colors[which(activity.col.seq == .)]) %>% unlist
    # 
    # p3 <- p + geom_text(label = "\u25D1", aes(x=x, y=y), color=activity_colors.i, size=20, family = "Arial Unicode MS",
    #                     data = p$data[,c("x","y")] %>% dplyr::mutate(activity=activity_colors.i)) +
    #   scale_fill_manual(values=activity_colors.i) + 
    #   labs(title="TF activity") +
    #   theme(
    #     plot.title = element_text(hjust = 0.5)
    #   )
    
    p.all <- cowplot::plot_grid(plotlist=list(p1,p2), nrow=1, scale = 0.95) + labs(title=i) + 
      theme(
        plot.title = element_text(hjust = 0.5, size=rel(1.25))
      )
    
    png(outfile, width = 650, height = 400)
    print(p.all)
    dev.off()
  }
}
  