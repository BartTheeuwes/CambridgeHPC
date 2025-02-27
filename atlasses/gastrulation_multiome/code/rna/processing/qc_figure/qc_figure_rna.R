source(here::here("settings.R"))

#####################
## Define settings ##
#####################

io$outdir <- file.path(io$basedir,"results/rna/qc/fig"); dir.create(io$outdir, showWarnings = F)

opts$samples <- c(
    "E7.5_rep1",
    "E7.5_rep2",
    "E7.75_rep1",
    "E8.0_rep1",
    "E8.0_rep2",
    "E8.5_rep1",
    "E8.5_rep2",
    "E8.75_rep1",
    "E8.75_rep2",
    "E8.5_CRISPR_T_WT",
    "E8.5_CRISPR_T_KO"
)

###################
## Load metadata ##
###################

metadata <- fread(io$metadata) %>% 
    .[pass_rnaQC==TRUE] %>% 
    .[sample%in%opts$samples] %>% .[,sample:=factor(sample,levels=opts$samples)]

###########################
## Boxplot of QC metrics ##
###########################

to.plot <- metadata %>%
    .[nFeature_RNA<=8000 & mitochondrial_percent_RNA<=40 & ribosomal_percent_RNA<=15] %>% # remove outliers for plotting
    melt(id.vars=c("sample","cell","stage"), measure.vars=c("nFeature_RNA","mitochondrial_percent_RNA","ribosomal_percent_RNA"))

facet.labels <- c("nFeature_RNA" = "Num. of genes", "mitochondrial_percent_RNA" = "Mitochondrial %", "ribosomal_percent_RNA" = "Ribosomal %")
    
p <- ggplot(to.plot, aes_string(x="sample", y="value", fill="stage")) +
    geom_boxplot(outlier.shape=NA, coef=1) +
    facet_wrap(~variable, scales="free_y", nrow=1, labeller = as_labeller(facet.labels)) +
    scale_fill_manual(values=opts$stage.colors) +
    guides(x = guide_axis(angle = 90)) +
    theme_classic() +
    theme(
        axis.text.y = element_text(colour="black",size=rel(1)),
        axis.text.x = element_text(colour="black",size=rel(0.8)),
        # axis.text.x = element_text(colour="black",size=rel(0.65), angle=20, hjust=1, vjust=1),
        axis.title.x = element_blank()
    )

pdf(file.path(io$outdir,"qc_metrics_boxplot.pdf"), width=9, height=5)
# pdf(sprintf("%s/qc_metrics_boxplot.pdf",io$outdir))
print(p)
dev.off()


#######################################################
## Plot number of cells that pass QC for each sample ##
#######################################################

to.plot <- metadata %>% .[,.N,by=c("sample","stage")]

p <- ggbarplot(to.plot, x="sample", y="N", fill="stage") +
    scale_fill_manual(values=opts$stage.colors) +
    labs(x="", y="Number of cells that pass QC") +
    theme(
        legend.position = "none",
        axis.text.y = element_text(colour="black",size=rel(0.8)),
        axis.text.x = element_text(colour="black",size=rel(0.65), angle=20, hjust=1, vjust=1),
    )

pdf(file.path(io$outdir,"qc_metrics_barplot.pdf"), width=6, height=5)
print(p)
dev.off()
