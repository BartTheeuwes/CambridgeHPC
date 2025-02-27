#####################
## Define settings ##
#####################
library(anndata)
library(GGally)
library(network)
library(sna)
library(SCENIC)
library(glmnet)
library(dplyr)
library(igraph)
library(intergraph)
library(visNetwork)
source(here::here("settings.R"))
source(here::here("utils.R"))
# Options

opts$celltypes <- setdiff(opts$celltypes, c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm"))
io$sce <- paste0(io$basedir, '/processed/rna/pseudobulk/SingleCellExperiment.rds')
io$outdir <-  paste0(io$basedir,"/results/rna_atac/GRN/pseudobulk")
io$chip_GRN.df <- paste0(io$outdir,'/chip_GRN_df.csv')
io$global_chip_GRN_coef.df <- paste0(io$outdir,'/global_chip_GRN_coef.csv')
# io$FigR_GRN <- paste0(io$basedir,"/results/rna_atac/GRN/FigR_GRN.csv")
dir.create(io$outdir)
setwd(io$outdir)


##############
## Load GRN ##
##############

GRN <- read.csv(io$chip_GRN.df)

##########################
## Load pseudobulk data ##
##########################

# Filter celltypes
sce <- readRDS(io$sce)
sce<- sce[,opts$celltypes]

marker_genes.dt <- fread(io$rna.atlas.marker_genes) 
sce <- sce[rownames(sce)%in%unique(marker_genes.dt$gene),]
exprMat <- counts(sce)
####################
## run regression ##
####################

allTFs <- unfactor(unique(GRN$TF))

# Select TFs
inputTFs <- allTFs[allTFs %in% rownames(exprMat)] 
GRN <- GRN[GRN$TF%in%inputTFs,]
target <- unfactor(unique(GRN$target))
percMatched <- length(inputTFs)/length(allTFs)
if(percMatched < .40) warning("Only ", length(inputTFs) ," (", round(percMatched*100), "%) of the ", length(allTFs)," TFs were found in the dataset")

# Selecting the features for regression

GRN_coef.df <- data.frame()
for (i in target){
  print(i)
  TF <- unfactor(GRN[GRN$target==i,]$TF)
  x <- t(exprMat[TF,])
  y <- exprMat[i,]
  if (length(TF)>1){
    cvfit <- cv.glmnet(x, y, alpha=0)
    opt_ridge <- glmnet(x,y , alpha = 0, lambda  = cvfit$lambda.min)
    df <- data.frame(target=i,TF=rownames(opt_ridge$beta),beta=opt_ridge$beta@x)}
  else{df <- data.frame(target=i,TF=TF,beta=mean(y)/mean(x))}
  GRN_coef.df <- rbind(df,GRN_coef.df)
}
GRN_coef.df <- as.data.table(GRN_coef.df)
save(GRN_coef.df,file=io$global_chip_GRN_coef.df)

# remove negative links
GRN_coef.df <-  GRN_coef.df[beta>0]

#  build coef matrix
genes <- rownames(exprMat)
GRN_coef.mtx <- matrix(data=0,nrow = length(genes),ncol = length(genes) )
rownames(GRN_coef.mtx) <- genes
colnames(GRN_coef.mtx) <- genes

## for TFs
TF <- unfactor(unique(GRN_coef.df$TF))
for (i in TF){
  targets <- unfactor(GRN_coef.df[TF==i,]$target)
  GRN_coef.mtx[i,targets] <- GRN_coef.df[TF==i]$beta
  GRN_coef.mtx[i,i] <- 1
}

## for targets
for (i in setdiff(genes,TF)){
  GRN_coef.mtx[i,i] <- 1
}

# simulate
gene_KO <- 'Tal1'
delta_X <- KO_simulate(gene_KO = gene_KO,GRN_coef.mtx = GRN_coef.mtx,exprMat = exprMat,n_propagation = 2)

# plot
source(here::here("load_paga_graph.R"))


opts$min.delta_X <- -9
delta_X.col.seq <- round(seq(opts$min.delta_X,0,1), 0)

colors <- colorRampPalette(c("red", "grey"))(length(delta_X.col.seq))

gene2plot <- 'Klf1'
delta_X.values <- delta_X[,gene2plot]
delta_X.values[delta_X.values<opts$min.delta_X] <- opts$min.delta_X
delta_X.colors <- round(delta_X.values,0) %>% map(~ colors[which(delta_X.col.seq == .)]) %>% unlist

p <- ggnet2(
  net = net.paga,
  mode = c("x", "y"),
  node.size = 0,
  edge.size = 0.15,
  edge.color = "grey",
  label = FALSE,
  label.size = 2.3
)
p1 <- p + geom_text(label = "\u25D0", aes(x=x, y=y), color=delta_X.colors[p$data$label], size=5) +
  scale_colour_manual(values=delta_X.colors) + 
  labs(title=paste(gene2plot,'change after',gene_KO,'KO') )+
  theme(
    plot.title = element_text(hjust = 0.5)
  )


alphas <- rep(0.5,length(opts$celltypes)); names(alphas) <- opts$celltypes
sizes <- rep(6,length(opts$celltypes)); names(sizes) <- opts$celltypes



p2 <- ggnet2(
  net = net.paga,
  mode = c("x", "y"),
  color = opts$celltype.colors[opts$celltypes],
  node.alpha = alphas,
  node.size = sizes,
  edge.size = 0.15,
  edge.color = "grey",
  label = TRUE,
  label.size = 3.5,
  legend.position = "none"
)

cowplot::plot_grid(plotlist=list(p1,p2), nrow=1)
# Load the n graph.
ig <- graph.data.frame(d = GRN_coef.df[1:100,], directed = FALSE)
