library(GGally)
library(network)
library(sna)

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
} else if (grepl("Workstation",Sys.info()['nodename'])){
  source("/home/lijingyu/gastrulation/gastrulation_multiome_10x/settings.R")
  source("/home/lijingyu/gastrulation/gastrulation_multiome_10x/utils.R")
} else {
  stop("Computer not recognised")
}
# I/O
io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/pseudobulk/paga")
io$DORCs.mtx <- paste0(io$basedir,"/results/rna_atac/DORCs/pseudobulk/DORCs_score_matrix_pseudobulk.rds")
# Options


opts$motif_annotation <- "Motif_cisbp"


#################################################
## Load pseudobulk RNA and chromVAR estimates ##
################################################

# io$archR.pseudobulk.deviations.se <- sprintf("%s/pseudobulk/pseudobulk_DeviationMatrix_%s_summarized_experiment.rds",io$archR.directory,opts$motif_annotation)
"/home/lijingyu/gastrulation/data/gastrulation_multiome_10x/results/atac/archR/chromvar/pseudobulk/chromVAR_deviations_summarized_experiment_Motif_cisbp_pseudobulk_correlated_peaks.rds"
io$archR.pseudobulk.deviations.se <- sprintf("%s/results/atac/archR/chromvar/pseudobulk/chromVAR_deviations_summarized_experiment_%s_pseudobulk_correlated_peaks.rds",io$basedir,opts$motif_annotation)

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/rna_atac/rna_vs_chromvar/pseudobulk/load_rna_chromvar_pseudobulk.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/rna_atac/rna_vs_chromvar/pseudobulk/load_rna_chromvar_pseudobulk.R")
} else if (grepl("Workstation",Sys.info()['nodename'])){
  source("/home/lijingyu/gastrulation/gastrulation_multiome_10x/rna_atac/rna_vs_chromvar/pseudobulk/load_rna_chromvar_pseudobulk.R")
} else {
  stop("Computer not recognised")
}

################
## Parse data ##
################

# Filter genes with low variability
genes.to.keep.rna <- rna.dt[,.(var(expr)),by="gene"] %>% .[V1>0.05,gene] %>% as.character
genes.to.keep.chromvar <- chromvar.dt[,.(var(chromvar_zscore)),by="gene"] %>% .[V1>0.05,gene] %>% as.character
genes.to.keep <- intersect(genes.to.keep.rna,genes.to.keep.chromvar)
rna.dt <- rna.dt[gene%in%genes.to.keep]
chromvar.dt <- chromvar.dt[gene%in%genes.to.keep]

#####################
## Load DORC score ##
#####################

DORCs.mtx <- readRDS(io$DORCs.mtx)
DORC_score.dt <- DORCs.mtx %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","gene") %>%
  melt(id.vars="gene", variable.name="celltype", value.name="DORC_score")
DORC_score.dt$gene<- toupper(DORC_score.dt$gene)

###########
## Merge ##
###########

rna_chromvar_DORC_score.dt <- merge(
  rna.dt,
  chromvar.dt,
  by = c("celltype","gene")
) %>% merge(.,DORC_score.dt,by = c("celltype","gene"))

#############################
## Calculate TF activities ##
#############################

rna_chromvar_DORC_score.dt %>%
  .[,activity:=minmax.normalisation(expr)*minmax.normalisation(chromvar_zscore),by="gene"] %>%
  .[,activity:=minmax.normalisation(activity),by="gene"]


#####################
## Load PAGA graph ##
#####################

celltypes <- unique(rna_chromvar_DORC_score.dt$celltype)
source("/home/lijingyu/gastrulation/gastrulation_multiome_10x/rna_atac/rna_vs_chromvar/pseudobulk/PAGA/load_paga_graph.R")

##############################################
## Plot network, colour by cell type labels ##
##############################################

p <- ggnet2(
  net = net.paga,
  mode = c("x", "y"),
  color = opts$celltype.colors[celltypes],
  node.size = 6,
  edge.size = 0.15,
  edge.color = "grey",
  label = TRUE,
  label.size = 2.3
)

pdf(paste0(io$outdir,"/paga_coloured_by_celltype.pdf"), width=6, height=4)
print(p)
dev.off()

#####################################################################
## Plot network, colour by gene expression and motif accessibility ##
#####################################################################

opts$scale <- TRUE

if (opts$scale) {
  DORC_score.col.seq  <- activity.col.seq <-  rna.col.seq <- chromvar.col.seq <- round(seq(0,1,0.1), 2)
  rna_chromvar_DORC_score.dt[,expr:=minmax.normalisation(expr),by="gene"]
  rna_chromvar_DORC_score.dt[,chromvar_zscore:=minmax.normalisation(chromvar_zscore),by="gene"]
  rna_chromvar_DORC_score.dt[,DORC_score:=minmax.normalisation(DORC_score),by="gene"]
  
} else {
  opts$max.expr <- 9
  opts$min.expr <- 2
  rna.col.seq <- round(seq(opts$min.expr,opts$max.expr,0.1), 2)
  rna_chromvar_DORC_score.dt[expr>=opts$max.expr,expr:=opts$max.expr]
  rna_chromvar_DORC_score.dt[expr<=opts$min.expr,expr:=opts$min.expr]
  
  opts$max.chromvar <- 6
  opts$min.chromvar <- 0
  chromvar.col.seq <- round(seq(opts$min.chromvar,opts$max.chromvar,0.1), 2)
  rna_chromvar_DORC_score.dt[chromvar_zscore>=opts$max.chromvar,chromvar_zscore:=opts$max.chromvar]
  rna_chromvar_DORC_score.dt[chromvar_zscore<=opts$min.chromvar,chromvar_zscore:=opts$min.chromvar]
  
  opts$max.DORC_score <- 1
  opts$min.DORC_score <- 0
  chromvar.col.seq <- round(seq(opts$min.DORC_score,opts$max.DORC_score,0.1), 2)
  rna_chromvar_DORC_score.dt[DORC_score>=opts$max.DORC_score,DORC_score:=opts$max.DORC_score]
  rna_chromvar_DORC_score.dt[DORC_score<=opts$min.DORC_score,DORC_score:=opts$min.DORC_score]
  
  
  activity.col.seq <- round(seq(0,max(rna_chromvar_DORC_score.dt$activity,na.rm=T),0.05), 2)
}

# Define colors
rna.colors <- colorRampPalette(c("gray92", "darkgreen"))(length(rna.col.seq))
chromvar.colors <- colorRampPalette(c("gray92", "purple"))(length(chromvar.col.seq)) 
activity.colors <- colorRampPalette(c("gray92", "#EE9A00"))(length(activity.col.seq)) 
DORC_score.colors <- colorRampPalette(c("gray92", "red"))(length(DORC_score.col.seq)) 

# Define genes to plot
genes.to.plot <- unique(rna_chromvar_DORC_score.dt$gene)# %>% head(n=10)
# genes.to.plot <- c("EOMES","T")

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
  
  outfile <- sprintf("%s/individual_genes/%s_paga_rna_vs_chromvar.png",io$outdir,i)
  if (file.exists(outfile)) {
    print(sprintf("file or %s already exists...",i))
  } else {
    
    expr.values <- rna_chromvar_DORC_score.dt[gene==i,c("celltype","expr")] %>% matrix.please %>% .[celltypes,]
    expr.colors <- round(expr.values,1) %>% map(~ rna.colors[which(rna.col.seq == .)]) %>% unlist
    
    p1 <- p + geom_text(label = "\u25D0", aes(x=x, y=y), color=expr.colors, size=20, family = "Arial Unicode MS",
                        data = p$data[,c("x","y")] %>% dplyr::mutate(expr=expr.colors)) +
      scale_colour_manual(values=expr.colors) + 
      labs(title="RNA expression") +
      theme(
        plot.title = element_text(hjust = 0.5)
      )
    
    # motif accessibility
    acc.values <- rna_chromvar_DORC_score.dt[gene==i,c("celltype","chromvar_zscore")] %>% matrix.please %>% .[,1]# %>% .[celltypes,]
    acc.colors <- round(acc.values,1) %>% map(~ chromvar.colors[which(chromvar.col.seq == .)]) %>% unlist
    
    p2 <- p + geom_text(label = "\u25D1", aes(x=x, y=y), color=acc.colors, size=20, family = "Arial Unicode MS",
                        data = p$data[,c("x","y")] %>% dplyr::mutate(acc=acc.colors)) +
      scale_fill_manual(values=acc.colors) + 
      labs(title="Motif accessibility") +
      theme(
        plot.title = element_text(hjust = 0.5)
      )
    
    # TF activity
    activity.values <- rna_chromvar_DORC_score.dt[gene==i,c("celltype","activity")] %>% matrix.please %>% .[,1]# %>% .[celltypes,]
    activity_colors.i <- round(activity.values,1) %>% map(~ activity.colors[which(activity.col.seq == .)]) %>% unlist
    
    p3 <- p + geom_text(label = "\u25D1", aes(x=x, y=y), color=activity_colors.i, size=20, family = "Arial Unicode MS",
                        data = p$data[,c("x","y")] %>% dplyr::mutate(activity=activity_colors.i)) +
      scale_fill_manual(values=activity_colors.i) + 
      labs(title="TF activity") +
      theme(
        plot.title = element_text(hjust = 0.5)
      )
    # DORC score
    DORC_score.values <- rna_chromvar_DORC_score.dt[gene==i,c("celltype","DORC_score")] %>% matrix.please %>% .[,1]# %>% .[celltypes,]
    DORC_score_colors.i <- round(DORC_score.values,1) %>% map(~ DORC_score.colors[which(DORC_score.col.seq == .)]) %>% unlist
    
    p4 <- p + geom_text(label = "\u25D1", aes(x=x, y=y), color=DORC_score_colors.i, size=20, family = "Arial Unicode MS",
                        data = p$data[,c("x","y")] %>% dplyr::mutate(DORC_score=DORC_score_colors.i)) +
      scale_fill_manual(values=DORC_score_colors.i) + 
      labs(title="DORC score") +
      theme(
        plot.title = element_text(hjust = 0.5)
      )
    
    
    p1234 <- cowplot::plot_grid(plotlist=list(p1,p2,p3,p4), nrow=1, scale = 0.9) + labs(title=i) + 
      theme(
        plot.title = element_text(hjust = 0.5, size=rel(1.25))
      )
    
    png(outfile, width = 1800, height = 600)
    print(p1234)
    dev.off()
  }
}
