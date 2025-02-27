here::i_am("rna_atac/rna_vs_acc/cells/gene_expr_vs_acc/plot_gene_expr_vs_acc_cells.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--sce',  type="character",              help='RNA SingleCellExperiment (single cells)') 
p$add_argument('--gene_score_matrix',  type="character",              help='ATAC Gene score matrix (single cells)') 
p$add_argument('--outdir',          type="character",                help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
args <- list()
args$sce <- io$rna.sce # file.path(io$basedir,"results_new/rna/pseudobulk/SingleCellExperiment_pseudobulk_celltype.mapped_mnn.rds") # io$rna.pseudobulk.sce
args$gene_score_matrix <- file.path(io$basedir,"results_new/atac/archR/gene_scores/GeneScoreMatrix_tss.rds")
args$outdir <- file.path(io$basedir,"results_new/rna_atac/rna_vs_acc/per_gene")
## END TEST ##

# Options

# I/O
dir.create(args$outdir, showWarnings=F)
dir.create(file.path(args$outdir,"individual_genes"), showWarnings=F)

###################
## Load metadata ##
###################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & pass_rnaQC==TRUE & doublet_call==FALSE] %>%
  .[!is.na(celltype.mapped_mnn)]


############################
## Load RNA and ATAC data ##
############################

# Load SingleCellExperiment
rna.sce <- load_SingleCellExperiment(
  file = io$rna.sce, 
  cells = sample_metadata$cell, 
  normalise = TRUE, 
  remove_non_expressed_genes = TRUE
)

# Load ATAC SummarizedExperiment
atac_GeneScoreMatrix.se <- readRDS(args$gene_score_matrix)[,sample_metadata$cell]

# Cap accessibility values
assay(atac_GeneScoreMatrix.se)[assay(atac_GeneScoreMatrix.se)>=15] <- 15

##################
## Filter genes ##
##################

genes <- intersect(rownames(rna.sce), rownames(atac_GeneScoreMatrix.se))

genes <- genes[grep("Rik",genes,invert = T)]
genes <- genes[grep("mt-",genes,invert = T)]
genes <- genes[grep("Rps|Rpl",genes,invert = T)]
genes <- genes[grep("Olfr",genes,invert = T)]

genes <- fread(io$rna.atlas.marker_genes.up)[,gene] %>% unique

##########################
## Correlation analysis ##
##########################

cor.dt <- genes[1:25] %>% map(function(i) {
  print(sprintf("%s (%d/%d)",i,match(i,genes),length(genes)))
  tmp <- cor.test(logcounts(rna.sce[i,])[1,], assay(atac_GeneScoreMatrix.se[i,])[1,])
  data.table(
    gene = i,
    r = tmp[["estimate"]],
    p = tmp[["p.value"]]
  )
}) %>% rbindlist

# Save
fwrite(cor.dt, file.path(args$outdir,"cor_rna_vs_acc_per_gene_cells.txt.gz"), sep="\t", quote=F)

#############################################
## Load pre-computed correlation estimates ##
#############################################

# cor.dt <- fread(paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/per_gene/cor_rna_vs_acc_per_gene_pseudobulk.txt.gz")) %>%
#   .[,cor_sign:=as.factor(c("Repressor","Activator")[(r>0)+1])]

opts$threshold_fdr <- 0.10
cor.dt %>% .[, sig:=p<=opts$threshold_fdr]

# Plot

to.plot <- cor.dt %>%
  .[p<=1e-14,p:=1e-14] %>%
  .[,log_pval:=-log10(p+1e-100)] %>%
  .[,dot_size:=minmax.normalisation(log_pval)]

negative_hits <- to.plot[sig==TRUE & r<0,gene]
positive_hits <- to.plot[sig==TRUE & r>0,gene]
all <- nrow(to.plot)

xlim <- max(abs(to.plot$r), na.rm=T)
ylim <- max(-log10(to.plot$p+1e-100), na.rm=T)

p <- ggplot(to.plot, aes(x=r, y=log_pval)) +
  geom_segment(aes(x=0, xend=0, y=0, yend=ylim-1), color="orange", size=0.5) +
  geom_jitter(aes(fill=sig, size=dot_size, alpha=dot_size), width=0.05, shape=21) + 
  ggrepel::geom_text_repel(data=head(to.plot[sig==T & r<(-0.30)],n=25), aes(x=r, y=log_pval, label=gene), size=3,  max.overlaps=100, segment.color = NA) +
  ggrepel::geom_text_repel(data=head(to.plot[sig==T & r>0.30],n=25), aes(x=r, y=log_pval, label=gene), size=3,  max.overlaps=100, segment.color = NA) +
  # ggrepel::geom_text_repel(data=to.plot[sig==T & r<(-0.50)][sample(.N,12)], aes(x=r, y=log_pval, label=gene), size=4,  max.overlaps=100, segment.color = NA) +
  # ggrepel::geom_text_repel(data=to.plot[sig==T & r>0.75][sample(.N,12)], aes(x=r, y=log_pval, label=gene), size=4,  max.overlaps=100, segment.color = NA) +
  scale_fill_manual(values=c("black","red")) +
  # scale_size_manual(values=c(0.5,1)) +
  scale_size_continuous(range = c(0.2,2)) + 
  scale_alpha_continuous(range=c(0.25,1)) +
  scale_x_continuous(limits=c(-xlim-0.15,xlim+0.15)) +
  scale_y_continuous(limits=c(0,ylim+1)) +
  annotate("text", x=0, y=ylim+1, size=4, label=sprintf("(%d)", all)) +
  annotate("text", x=-xlim-0.15, y=ylim+1, size=4, label=sprintf("%d (-)",length(negative_hits))) +
  annotate("text", x=xlim+0.15, y=ylim+1, size=4, label=sprintf("%d (+)",length(positive_hits))) +
  labs(x="Pearson correlation (RNA expression vs gene accessibility)", y=expression(paste("-log"[10],"(p.value)"))) +
  theme_classic() +
  theme(
    axis.text = element_text(size=rel(0.75), color='black'),
    axis.title = element_text(size=rel(1.0), color='black'),
    legend.position="none"
  )

pdf(file.path(args$outdir,"volcano_cor_rna_vs_acc_per_gene_pseudobulk.pdf"), width = 7, height = 5)
print(p)
dev.off()

##########################
## Scatterplot per gene ##
##########################

facet.labels <- c(expr = "RNA expression", acc = "Gene accessibility")

# genes.to.plot <- c("Ubb", "Top2a", "Actg1", "Gapdh")
genes.to.plot <- cor.dt$gene

# i <- "Ubb"
for (i in genes.to.plot) {
  outfile <- file.path(args$outdir,sprintf("individual_genes/%s_rna_vs_acc_cells.png",i))
  if (!file.exists(outfile)) {
    
    to.plot <- data.table(
        cell = sample_metadata$cell,
        gene = i,
        celltype = sample_metadata$celltype.mapped_mnn,
        expr = logcounts(rna.sce[i,])[1,],
        acc = assay(atac_GeneScoreMatrix.se[i,])[1,]
      )
      
    p <- ggscatter(to.plot, x="acc", y="expr", color="celltype", size=1, 
                    add="reg.line", add.params = list(color="black", fill="lightgray"), conf.int=TRUE) +
      stat_cor(method = "pearson", label.x.npc = "middle", label.y.npc = "bottom") +
      scale_color_manual(values=opts$celltype.colors) +
      labs(y=sprintf("%s expression",i), x=sprintf("%s Gene accessibility",i)) +
      guides(color="none") +
      theme(
        axis.text = element_text(size=rel(0.8))
      )

    # p <- cowplot::plot_grid(plotlist=list(p1,p2), nrow = 1, rel_widths = c(1/2,1/2))
    png(outfile, width = 500, height = 400)
    # pdf(outfile, height=8, width=11)
    print(p)
    dev.off()

  }
}
