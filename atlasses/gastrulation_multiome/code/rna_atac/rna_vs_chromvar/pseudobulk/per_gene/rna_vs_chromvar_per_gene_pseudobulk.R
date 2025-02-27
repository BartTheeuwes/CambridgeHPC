suppressPackageStartupMessages(library(argparse))

here::i_am("rna_atac/rna_vs_chromvar/pseudobulk/per_gene/rna_vs_chromvar_per_gene_pseudobulk.R")

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--motif_annotation',  type="character",              help='Motif annotation') 
p$add_argument('--outdir',          type="character",                help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

## START TEST ##
# args <- list()
# args$motif_annotation <- "CISBP"
# args$outdir <- file.path(io$basedir,"results_new2/rna_atac/rna_vs_chromvar/pseudobulk/per_gene")
## END TEST ##

# I/O
dir.create(args$outdir, showWarnings=F)

# Options
opts$motif_annotation <- args$motif_annotation

#######################################
## Load pseudobulk RNA and ATAC data ##
#######################################

source(here::here("rna_atac/load_rna_atac_pseudobulk.R"))

###########
## Merge ##
###########

rna_chromvar.dt <- merge(
  rna_tf_pseudobulk.dt,
  atac_chromvar_pseudobulk.dt,
  by = c("celltype","gene")
)

length(unique(rna_chromvar.dt$gene))
length(unique(rna_chromvar.dt$motif))

##########################
## Correlation analysis ##
##########################

# opts$threshold_fdr <- 0.10

cor.dt <- rna_chromvar.dt %>% copy %>%
  .[,c("chromvar_zscore","expr"):=list(chromvar_zscore + rnorm(n=.N,mean=0,sd=1e-5), expr + rnorm(n=.N,mean=0,sd=1e-5))] %>% # add some noise 
  .[, .(V1 = unlist(cor.test(chromvar_zscore, expr)[c("estimate", "p.value")])), by = c("gene","motif")] %>%
  .[, para := rep(c("r","p"), .N/2)] %>% 
  data.table::dcast(gene+motif ~ para, value.var = "V1") %>%
  .[,"padj_fdr" := list(p.adjust(p, method="fdr"))] %>%
  # .[, sig := p<=opts$threshold_fdr] %>% 
  setorder(p, na.last = T)

# Save
fwrite(cor.dt, file.path(args$outdir,sprintf("cor_rna_vs_chromvar_%s_per_gene_pseudobulk.txt.gz",args$motif_annotation)), sep="\t", quote=F)

#############################################
## Load pre-computed correlation estimates ##
#############################################

# cor.dt <- fread(paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/pseudobulk/per_gene/cor_rna_vs_chromVAR+_pseudobulk.txt.gz")) %>%
#   .[,cor_sign:=as.factor(c("Repressor","Activator")[(r>0)+1])]

opts$threshold_fdr <- 0.10
cor.dt %>% .[, sig:=p<=opts$threshold_fdr]

##################
## Volcano plot ##
##################

to.plot <- cor.dt %>%
  .[padj_fdr<=1e-14,padj_fdr:=1e-14] %>%
  .[,log_pval:=-log10(padj_fdr+1e-100)] %>%
  .[,dot_size:=minmax.normalisation(log_pval)]

negative_hits <- to.plot[sig==TRUE & r<0,gene]
positive_hits <- to.plot[sig==TRUE & r>0,gene]
all <- nrow(to.plot)

xlim <- max(abs(to.plot$r), na.rm=T)
ylim <- max(-log10(to.plot$padj_fdr+1e-100), na.rm=T)

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
  labs(x="Pearson correlation (TF RNA expr vs TF chromVAR score)", y=expression(paste("-log"[10],"(p.value)"))) +
  theme_classic() +
  theme(
    axis.text = element_text(size=rel(0.75), color='black'),
    axis.title = element_text(size=rel(1.0), color='black'),
    legend.position="none"
  )

pdf(file.path(args$outdir,sprintf("volcano_cor_rna_vs_chromvar_%s_per_gene_pseudobulk.pdf",args$motif_annotation)), width = 7, height = 5)
print(p)
dev.off()


###############
## Load PAGA ##
###############



######################################
## Scatter plot of individual genes ##
######################################

facet.labels <- c(expr = "RNA expression", chromvar_zscore = "Motif accessibility (z-score)")

for (i in unique(rna_chromvar.dt$gene)) {
  for (j in unique(rna_chromvar.dt[gene==i,motif])) {

    outfile <- file.path(args$outdir,sprintf("%s_%s_%s_rna_vs_chromvar_pseudobulk.png",i,j,args$motif_annotation))

    if (file.exists(outfile)) {
      print(sprintf("%s already exists...",outfile))
    } else {
      to.plot <- rna_chromvar.dt[gene==i & motif==j]
      to.plot2 <- to.plot %>% melt(id.vars=c("celltype","gene","motif"))

      p1 <- ggscatter(to.plot, x="expr", y="chromvar_zscore", fill="celltype", size=4, shape=21, 
                      add="reg.line", add.params = list(color="black", fill="lightgray"), conf.int=TRUE) +
        stat_cor(method = "pearson", label.x.npc = "middle", label.y.npc = "bottom") +
        scale_fill_manual(values=opts$celltype.colors) +
        labs(x=sprintf("%s expression",i), y=sprintf("%s Motif accessibility (z-score)",i)) +
        guides(fill="none") +
        theme(
          axis.text = element_text(size=rel(0.7))
        )
      
      p2 <- ggbarplot(to.plot2, x="celltype", y="value", fill="celltype") +
        facet_wrap(~variable, nrow=2, scales="free_y",  labeller = as_labeller(facet.labels)) +
        scale_fill_manual(values=opts$celltype.colors) +
        geom_hline(yintercept=0, linetype="dashed") +
        labs(x="", y="") +
        theme_classic() +
        guides(x = guide_axis(angle = 90)) +
        theme(
          legend.position = "none",
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank()
        )
      
      p <- cowplot::plot_grid(plotlist=list(p1,p2), nrow = 1, rel_widths = c(1/2,1/2))
      
      # pdf(sprintf("%s/individual_genes/%s_rna_vs_chromvar.pdf",args$outdir,i), width=10, height=4)
      png(outfile, width = 1000, height = 500)
      print(p)
      dev.off()
    }
  }
}

