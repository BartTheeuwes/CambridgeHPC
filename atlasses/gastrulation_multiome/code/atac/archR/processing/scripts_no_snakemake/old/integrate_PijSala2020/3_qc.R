########################
## Load ArchR project ##
########################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/atac/archR/load_integrated_archR_project.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/atac/archR/load_integrated_archR_project.R")
} else {
  stop("Computer not recognised")
}

#####################
## Define settings ##
#####################

# I/O
# io$metadata <- paste0(io$basedir,"/processed/atac/archR/sample_metadata_after_archR.txt.gz")
io$metadata <- paste0(io$basedir,"/processed/atac/archR_integrated/sample_metadata.txt.gz")
io$outdir <- paste0(io$basedir,"/results/atac/archR_integrated/qc")

# Options
opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E8.5_rep1",
  "E8.5_rep2",
  "E8.25_PijuanSala"
)

# QC thresholds
opts$min.TSSEnrichment <- 9
opts$min.log_nFrags <- 11.25 # 2**11.25 = 2435
# opts$min.log_nFrags <- 11
opts$max.BlacklistRatio <- 0.05

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata)

#############
## Call QC ##
#############

sample_metadata %>%
  .[,pass_atacQC:=TSSEnrichment_atac>=opts$min.TSSEnrichment & log2(nFrags_atac)>=opts$min.log_nFrags & BlacklistRatio_atac<=opts$max.BlacklistRatio] %>%
  .[is.na(pass_atacQC),pass_atacQC:=FALSE]

# print(sample_metadata[,mean(pass_atacQC,na.rm=T),by="sample"])

# Save
outfile <- paste0(io$outdir,"/sample_metadata_after_qc.txt.gz")
fwrite(sample_metadata, outfile, quote=F, na="NA", sep="\t")


##################
## Subset ArchR ##
##################

sample_metadata.filt <- sample_metadata[pass_atacQC==T]
ArchRProject.filt <- ArchRProject[sample_metadata.filt$archR_cell,]

#########################
## Plot TSS Enrichment ##
#########################

to.plot.tss <- plotTSSEnrichment(ArchRProject.filt, groupBy = "Sample", returnDF=T) %>% 
  as.data.table %>% setnames("group","sample") %>% 
  melt(id.vars=c("sample","x")) 

p <- ggline(to.plot.tss[variable=="normValue"], x="x", y="value", plot_type="l") +
  facet_wrap(~sample, scales="fixed") +
  scale_x_discrete(breaks=seq(-2000,2000,1000)) +
  labs(x="Distance from TSS (bp)", y="TSS enrichment (normalised)") +
  theme(
    axis.text = element_text(size=rel(0.75)),
    axis.title = element_text(size=rel(0.75)),
    legend.position = "right",
    legend.title = element_blank()
  )

fwrite(to.plot.tss, sprintf("%s/qc_TSSenrichment.txt.gz",io$outdir))
pdf(sprintf("%s/qc_TSSenrichment.pdf",io$outdir))
print(p)
dev.off()

#####################################
## Plot Fragment size distribution ##
#####################################

to.plot.fragmentsize <- plotFragmentSizes(ArchRProject.filt, groupBy = "Sample", returnDF=T) %>% 
  as.data.table %>% setnames("group","sample") %>% 
  melt(id.vars=c("sample"))

p <- ggline(to.plot.fragmentsize, x="fragmentSize", y="fragmentPercent", plot_type="l") +
  facet_wrap(~sample, scales="fixed") +
  scale_x_continuous(breaks=seq(125,750,125)) +
  labs(x="Fragment Size (bp)", y="Percentage of fragments (%)") +
  theme(
    axis.text = element_text(size=rel(0.55)),
    axis.title = element_text(size=rel(0.75)),
    legend.position = "right",
    legend.title = element_blank()
  )

fwrite(to.plot.fragmentsize, sprintf("%s/qc_FragmentSizeDistribution.txt.gz",io$outdir))
pdf(sprintf("%s/qc_FragmentSizeDistribution.pdf",io$outdir))
print(p)
dev.off()

###########################################
## Plot summary statistics of QC metrics ##
###########################################

to.plot <- sample_metadata %>% copy %>%
    .[,log_nFrags:=log2(nFrags)] %>%
    melt(id.vars=c("sample","cell"), measure.vars=c("TSSEnrichment","log_nFrags","BlacklistRatio"))

# Boxplots
p <- ggboxplot(to.plot, x="sample", y="value") +
    facet_wrap(~variable, scales="free_y") +
    theme(
        axis.text.x = element_text(colour="black",size=rel(0.65), angle=20, hjust=1, vjust=1),  
        axis.title.x = element_blank()
    )
pdf(sprintf("%s/qc_metrics_boxplot.pdf",io$outdir))
print(p)
dev.off()

# Histograms
tmp <- data.table(
  variable = c("TSSEnrichment", "log_nFrags", "BlacklistRatio"),
  value = c(opts$min.TSSEnrichment, opts$min.log_nFrags, opts$max.BlacklistRatio)
)
p <- gghistogram(to.plot, x="value", fill="sample", bins=50) +
  geom_vline(aes(xintercept=value), linetype="dashed", data=tmp) +
  facet_wrap(~variable, scales="free") +
  theme(
    axis.text =  element_text(size=rel(0.8)),
    axis.title.x = element_blank(),
    legend.position = "right",
    legend.text = element_text(size=rel(0.5))
  )
pdf(sprintf("%s/qc_metrics_histogram.pdf",io$outdir))
print(p)
dev.off()


#########################################################
## Plot fraction of cells that pass QC for each sample ##
#########################################################

to.plot <- sample_metadata %>%
  .[,mean(pass_atacQC),by="sample"]

p <- ggbarplot(to.plot, x="sample", y="V1", fill="gray70") +
    labs(x="", y="Fraction of cells that pass QC") +
  coord_cartesian(ylim=c(0,1)) +
    theme(
        axis.text.x = element_text(colour="black",size=rel(0.65), angle=20, hjust=1, vjust=1),  
    )

# pdf(sprintf("%s/qc_metrics_barplot.pdf",io$outdir), width=9, height=7)
pdf(sprintf("%s/qc_metrics_barplot.pdf",io$outdir))
print(p)
dev.off()
