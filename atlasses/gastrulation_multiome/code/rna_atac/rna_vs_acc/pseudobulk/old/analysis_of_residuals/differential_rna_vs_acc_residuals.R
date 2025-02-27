library(pheatmap)

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

# I/O
io$basedir <- file.path(io$basedir,"test")
io$rna_diff <- file.path(io$basedir,"results/rna/differential/pseudobulk/celltype/parsed/diff_expr_results.txt.gz")
io$atac_diff <- file.path(io$basedir,"results/atac/archR/differential/pseudobulk/celltype/PeakMatrix/parsed/diff_results.txt.gz")
io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/differential/pseudobulk/residuals")

# Options
opts$celltypes <- c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  # "PGC",
  # "Anterior_Primitive_Streak",
  "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  # "Mixed_mesoderm",
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
  # "Erythroid1",
  # "Erythroid2",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
  # "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm"
  # "Visceral_endoderm",
  # "ExE_endoderm",
  # "ExE_ectoderm",
  # "Parietal_endoderm"
)
# opts$celltypes <- c("Epiblast", "Primitive_Streak", "Caudal_epiblast")

# RNA
opts$min_fold_change_rna <- 2
opts$min_fold_change_atac <- 1.5

##########################################
## Load differential expression results ##
##########################################

rna_diff.dt <- fread(io$rna_diff) %>% 
  .[celltypeA%in%opts$celltypes & celltypeB%in%opts$celltypes] %>%
  .[,.(nhits_positive=sum(logFC>=opts$min_fold_change_rna & padj_fdr<=0.01, na.rm=T), 
       nhits_negative=sum(logFC<=-(opts$min_fold_change_rna) & padj_fdr<=0.01, na.rm=T)),by=c("celltypeA","celltypeB")] %>%
  .[,nhits:=nhits_positive+nhits_negative]

#####################
## Plot DE results ##
#####################

to.plot <- rna_diff.dt
celltype.order <- to.plot %>% .[,.(median(nhits)),by="celltypeA"] %>% setorder(-V1) %>% .$celltypeA
to.plot <- to.plot %>% .[,celltypeA:=factor(celltypeA,levels=celltype.order)]

p <- ggplot(to.plot, aes(x=factor(celltypeA), y=nhits)) +
  # geom_point(aes(fill = celltypeA), shape=21, size=1) +
  geom_boxplot(aes(fill = celltypeA), alpha=0.9, outlier.shape=NA, coef=1) +
  # coord_flip(ylim = c(-5,1500)) +
  geom_hline(yintercept=0, linetype="dashed", size=0.5) +
  scale_fill_manual(values=opts$celltype.colors, drop=F) +
  theme_classic() +
  labs(y="Number of DE genes", x="") +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  )

pdf(sprintf("%s/rna_number_DE_genes_boxplots.pdf",io$outdir), width=5, height=6)
print(p)
dev.off()

#############################################
## Load differential accessibility results ##
#############################################

atac_diff.dt <- fread(io$atac_diff) %>% 
  .[celltypeA%in%opts$celltypes & celltypeB%in%opts$celltypes] %>%
  .[,.(nhits_positive=sum(logFC>=opts$min_fold_change_atac & padj_fdr<=0.01, na.rm=T), 
       nhits_negative=sum(logFC<=-(opts$min_fold_change_atac) & padj_fdr<=0.01, na.rm=T)),by=c("celltypeA","celltypeB")] %>%
  .[,nhits:=nhits_positive+nhits_negative]

#####################
## Plot DA results ##
#####################

to.plot <- atac_diff.dt
celltype.order <- to.plot %>% .[,.(median(nhits)),by="celltypeA"] %>% setorder(-V1) %>% .$celltypeA
to.plot <- to.plot %>% .[,celltypeA:=factor(celltypeA,levels=celltype.order)]

p <- ggplot(to.plot, aes(x=factor(celltypeA), y=nhits)) +
  # geom_point(aes(fill = celltypeA), shape=21, size=1) +
  geom_boxplot(aes(fill = celltypeA), alpha=0.9, outlier.shape=NA, coef=1) +
  # coord_flip(ylim = c(-5,9500)) +
  geom_hline(yintercept=0, linetype="dashed", size=0.5) +
  scale_fill_manual(values=opts$celltype.colors, drop=F) +
  theme_classic() +
  labs(y="Number of DA peaks", x="") +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  )

pdf(sprintf("%s/acc_number_DA_peaks_boxplots.pdf",io$outdir), width=5, height=6)
print(p)
dev.off()

#########################################
## Scatterplot of DE genes vs DA peaks ##
#########################################

to.plot <- merge(rna_diff.dt, atac_diff.dt, by=c("celltypeA","celltypeB"), suffixes=c("_rna","_acc")) %>%
  .[,transition:=sprintf("%s_to_%s",celltypeA,celltypeB)] %>%
  .[,residuals:=lm(formula=nhits_rna~nhits_acc, data=.)[["residuals"]]]

to.plot.text <- to.plot[order(-abs(residuals))] %>% head(n=25)

# to.plot[nhits_acc>7000,nhits_acc:=7000]

p <- ggscatter(to.plot, x="nhits_acc", y="nhits_rna", size=1.5, fill="gray50", shape=21,
          add="reg.line", add.params = list(color="blue", fill="lightgray"), conf.int=TRUE) +
  # scale_fill_manual(values=opts$celltype.colors, drop=F) +
  stat_cor(method = "pearson") +
  # coord_cartesian(ylim=c(0,750)) +
  # ggrepel::geom_text_repel(aes(label=transition), size=3, max.overlaps=Inf, data=to.plot.text) +
  labs(x="Number of differentially expressed genes", y="Number of differentially accessible peaks") +
  theme(
    legend.position = "none",
    axis.text = element_text(size=rel(0.5)),
    axis.title = element_text(size=rel(0.85))
  )

pdf(sprintf("%s/number_DE_genes_vs_number_DA_peaks_scatterplot.pdf",io$outdir), width=6, height=4.5)
print(p)
dev.off()

##########################
## Boxplot of residuals ##
##########################

celltype.order <- to.plot %>% .[,.(mean=mean(residuals)),by="celltypeA"] %>% setorder(-mean) %>% .$celltypeA
to.plot <- to.plot %>% .[,celltypeA:=factor(celltypeA,levels=celltype.order)]

p <- ggplot(to.plot, aes(x=factor(celltypeA), y=residuals)) +
  # geom_point(aes(fill = celltypeA), shape=21, size=1) +
  geom_boxplot(aes(fill = celltypeA), alpha=0.9, outlier.shape=NA, coef=1) +
  # coord_flip(ylim=c(-600,600)) +
  geom_hline(yintercept=0, linetype="dashed", size=0.5) +
  scale_fill_manual(values=opts$celltype.colors, drop=F) +
  theme_classic() +
  labs(y="Residuals", x="") +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black", size=rel(0.75)),
  )

pdf(sprintf("%s/rna_vs_acc_differential_residuals_boxplots.pdf",io$outdir), width=5, height=6)
print(p)
dev.off()
