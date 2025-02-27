# here::i_am("")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$outdir <- file.path(io$basedir,"results/rna_atac/wt_vs_ko/diff_expr_vs_motif_enrichment"); dir.create(io$outdir, showWarnings=F, recursive=T)

# Options
opts$motif_annotation <- "CISBP"

##############################
## Load TF motif enrichment ##
##############################

motif2gene.dt <- fread(sprintf("%s/Annotations/%s_TFs.txt.gz",io$archR.directory,opts$motif_annotation))

motif_enrichment.dt <- fread(file.path(io$basedir,"results/atac/archR/differential/wt_vs_ko/PeakMatrix/pdf/motif_enrichment_DA_peaks_wt_vs_ko.txt.gz")) %>% 
  merge(motif2gene.dt[,c("motif","gene")], by="motif") %>%
  .[,log_pval:=-log10(pval)] %>%
  setnames("log_pval","motif_enrichment_log_pval")

######################################
## Load differential RNA expression ##
######################################

diff_rna_tf.dt <- fread(file.path(io$basedir,"results/rna/differential/wt_vs_ko/pseudobulk/WT_vs_T_KO_DE_pseudobulk.txt.gz")) %>%
  .[,gene:=toupper(gene)] %>% .[gene%in%unique(motif_enrichment.dt$gene)] %>%
  setnames("diff","diff_gene_expr")

###########
## Merge ##
###########

to.plot <- merge(
  motif_enrichment.dt[,c("gene","celltype","motif_enrichment_log_pval")],
  diff_rna_tf.dt,
  by = c("celltype","gene")
)

##########
## Plot ##
##########

celltypes.to.plot <- c("Cardiomyocytes","Neural_crest","Gut","Endothelium","NMP","Spinal_cord")

to.plot2 <- to.plot %>%
  .[celltype%in%celltypes.to.plot] %>%
  .[,dot_size:=minmax.normalisation(abs(diff_gene_expr)*motif_enrichment_log_pval)]

# to.plot.text <- to.plot2[motif_enrichment_log_pval>=2 & abs(diff_gene_expr)>=0.5]
to.plot.text <- rbind(
  to.plot2[celltype=="NMP" & motif_enrichment_log_pval>=5],
  to.plot2[celltype=="NMP" & abs(diff_gene_expr)>1.5],
  to.plot2[celltype=="NMP" & motif_enrichment_log_pval>=2 & abs(diff_gene_expr)>=1]
) %>% unique

ggplot(to.plot2, aes_string(x="motif_enrichment_log_pval", y="diff_gene_expr", size="dot_size")) +
  geom_point(shape=21, fill="gray80") +
  facet_wrap(~celltype) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_vline(xintercept=0, linetype="dashed") +
  labs(x="Motif enrichment", y="TF differential expression") +
  ggrepel::geom_text_repel(data=to.plot.text, aes(x=motif_enrichment_log_pval, y=diff_gene_expr, label=gene), size=3) +
# scale_fill_gradient2(low = "gray50", mid="gray90", high = "red") +
  scale_size_continuous(range = c(0.1,3)) +
  theme_classic() +
  guides(size="none") +
  theme(
    axis.text.x = element_text(color="black", size=rel(0.75)),
    axis.text.y = element_text(color="black")
  )
