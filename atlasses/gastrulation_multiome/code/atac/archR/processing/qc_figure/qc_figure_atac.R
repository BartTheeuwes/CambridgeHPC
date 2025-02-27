source(here::here("settings.R"))

#####################
## Define settings ##
#####################

io$tss_enrichment <- file.path(io$basedir,"results/atac/archR/qc/qc_TSSenrichment.txt.gz")
io$fragment_size <- file.path(io$basedir,"results/atac/archR/qc/qc_FragmentSizeDistribution.txt.gz")

io$outdir <- file.path(io$basedir,"results/atac/archR/qc/fig"); dir.create(io$outdir, showWarnings = F)

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

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE] %>%
  .[sample%in%opts$samples] %>% .[,sample:=factor(sample,levels=opts$samples)]

#########################
## Plot TSS Enrichment ##
#########################

tss_enrichment.dt <- fread(io$tss_enrichment) %>%
  .[,.(value=mean(value)), by = c("sample","x","variable")]

to_plot <- tss_enrichment.dt %>% 
  merge(unique(sample_metadata[,c("sample","stage")])) %>%
  .[,.(value=mean(value)), by = c("stage","x","variable")] %>%
  .[variable=="normValue"]

p <- ggline(to_plot, x="x", y="value", plot_type="l", color="stage") +
  # facet_wrap(~stage, scales="fixed") +
  scale_colour_manual(values=opts$stage.colors) +
  # scale_x_continuous(breaks=seq(-2000,2000,1000)) +
  labs(x="Distance from TSS (bp)", y="TSS enrichment (normalised)") +
  theme(
    axis.text.y = element_text(size=rel(0.85), color="black"),
    axis.text.x = element_text(size=rel(0.85), color="black"),
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank(),
    axis.title = element_text(size=rel(0.75), color="black"),
    legend.position = "right",
    legend.title = element_blank()
  )

pdf(file.path(io$outdir,"qc_TSSenrichment.pdf"), width=5, height=4)
print(p)
dev.off()

#####################################
## Plot Fragment size distribution ##
#####################################

fragment_size.dt <- fread(io$fragment_size) %>%
  .[,.(fragmentPercent=mean(fragmentPercent)), by = c("sample","fragmentSize")]

to_plot <- fragment_size.dt %>% 
  merge(unique(sample_metadata[,c("sample","stage")])) %>%
  .[,.(fragmentPercent=mean(fragmentPercent)), by = c("stage","fragmentSize")] %>%
  .[fragmentSize<=350]

p <- ggline(to_plot, x="fragmentSize", y="fragmentPercent", plot_type="l", color="stage", size=0.75) +
  # facet_wrap(~stage, scales="fixed") +
  # scale_x_continuous(breaks=seq(125,750,125)) +
  scale_colour_manual(values=opts$stage.colors) +
  labs(x="Fragment Size (bp)", y="Percentage of fragments (%)") +
  theme(
    axis.text = element_text(size=rel(0.75)),
    axis.title = element_text(size=rel(0.85)),
    legend.position = "right",
    legend.title = element_blank()
  )

pdf(file.path(io$outdir,"qc_FragmentSizeDistribution.pdf"), width=5, height=4)
print(p)
dev.off()

############################
## Boxplots of QC metrics ##
############################

to.plot <- sample_metadata %>%
  .[nFrags_atac<=150000 & TSSEnrichment_atac<=27] %>% # remove outliers for plotting
  .[,log_nFrags_atac:=log10(nFrags_atac)] %>%
  melt(id.vars=c("sample","cell","sample","stage"), measure.vars=c("TSSEnrichment_atac","log_nFrags_atac"))

facet.labels <- c("log_nFrags_atac" = "Num. of fragments (log10)", "TSSEnrichment_atac" = "TSS enrichment")

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
print(p)
dev.off()


#######################################
## Plot number of cells that pass QC ##
#######################################

to.plot <- sample_metadata %>% .[,.N,by=c("sample","stage")]

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

