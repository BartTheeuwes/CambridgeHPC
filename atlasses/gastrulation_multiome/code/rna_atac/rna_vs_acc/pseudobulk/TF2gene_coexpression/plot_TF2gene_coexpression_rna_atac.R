
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
io$input.dir <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_Geneexpr")
io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_Geneexpr/pdf"); dir.create(io$outdir, showWarnings = F)


######################################
## Load TF2gene correlation results ##
######################################

opts$min.cor.TF2peak <- 0.25

to.save <- readRDS(sprintf("%s/cor_TFexpr_vs_Geneexpr_SummarizedExperiment_mincor%s.rds",io$outdir,opts$min.cor.TF2peak))
cor.mtx <- dropNA2matrix(assay(to.save,"cor"))
pval.mtx <- dropNA2matrix(assay(to.save,"pvalue"))

tf2gene.dt <- merge(
  x = cor.mtx %>% as.data.table(keep.rownames = T) %>% 
    setnames("rn","TF") %>% melt(id.vars="TF", variable.name="gene", value.name="cor") %>%
    .[!is.na(cor)],
  y = pval.mtx %>% as.data.table(keep.rownames = T) %>% 
    setnames("rn","TF") %>% melt(id.vars="TF", variable.name="gene", value.name="pval") %>%
    .[!is.na(pval)],
  by = c("TF","gene")
)

###########################################################
## Plot number of correlated genes per TF, split by sign ##
###########################################################

to.plot <- tf2gene.dt %>%
  .[,cor_sign:=c("-","+")[(cor>0)+1]] %>%
  .[,.(N=sum(pval<0.20 & abs(cor)>0.25)),by=c("cor_sign","TF")]

TFs.to.plot <- to.plot[,.(N=sum(N)),by="TF"] %>% setorder(-N) %>% head(n=100) %>% .$TF
to.plot2 <- to.plot[TF%in%TFs.to.plot] %>% .[,TF:=factor(TF,levels=TFs.to.plot)]

p <- ggbarplot(to.plot2, x="TF", y="N", fill="cor_sign") +
  coord_flip() +
  labs(x="", y="Number of correlation events with peaks") +
  theme(
    axis.text.y = element_text(size=rel(0.5)),
    axis.text.x = element_text(colour="black",size=rel(0.8)),
    axis.ticks.x = element_blank(),
    legend.position = "top"
  )

pdf(sprintf("%s/barplot_number_correlation_events_perTF_cor%s.pdf",io$outdir,opts$min.cor.TF2peak), width = 6, height = 12)
print(p)
dev.off()

#########################################
## Scatterplot of TF expr vs gene expr ##
#########################################

i <- "RUNX1"
j <- "Trerf1"

cor.mtx[i,j]

to.plot <- data.table(
  TF = logcounts(rna.sce.tf[i,])[1,],
  target_gene = logcounts(rna.sce[j,])[1,],
  celltype = colnames(rna.sce.tf)
)


ggscatter(to.plot, x="TF", y="target_gene", fill="celltype", size=4, shape=21, 
          add="reg.line", add.params = list(color="black", fill="lightgray"), conf.int=TRUE) +
  stat_cor(method = "pearson") +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(x=sprintf("%s expression",i), y=sprintf("%s expression",j)) +
  guides(fill=F) +
  theme(
    plot.title = element_text(hjust = 0.5, size=rel(0.85)),
    axis.text = element_text(size=rel(0.7))
  )

#################
## Exploration ##
#################


######################################
## Load TF2peak correlation results ##
######################################

# io$file <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc/cor_TFexpr_vs_peakAcc_SummarizedExperiment.rds")
# tf2peak_cor.se <- readRDS(io$file)

io$file <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc/cor_TFexpr_vs_peakAcc.txt.gz")
tf2peak_cor.dt <- fread(io$file) %>%
  .[!is.na(cor)] %>%
  .[,cor_sign:=c("negative","positive")[(cor>0)+1]]

#######

foo <- tf2peak_cor.dt[,.N,by=c("cor_sign","TF")] %>%
  dcast(TF~cor_sign, value.var="N", fill=0) %>%
  .[,N:=positive+negative]

bar <- tf2gene.dt %>%
  .[,cor_sign:=c("negative","positive")[(cor>0)+1]] %>%
  .[,.(N=sum(pval<=0.10 & abs(cor)>=0.50)),by=c("cor_sign","TF")] %>%
  dcast(TF~cor_sign, value.var="N", fill=0) %>%
  .[,N:=positive+negative]

foobar <- merge(foo, bar, suffixes=c("_tf2peak","_tf2gene"))
foobar[,fraction_positive_tf2gene:=positive_tf2gene/(positive_tf2gene+negative_tf2gene)]
foobar[,fraction_positive_tf2peak:=positive_tf2peak/(positive_tf2peak+negative_tf2peak)]

p <- ggscatter(foobar, x="N_tf2peak", y="N_tf2gene", color="fraction_positive_tf2peak", add="reg.line", add.params = list(color="blue", fill="lightgray"), conf.int=TRUE) +
  stat_cor(method = "pearson") +
  scale_color_continuous(low="red", high="green") +
  ggrepel::geom_text_repel(data=foobar[N_tf2peak>1000 | N_tf2gene>300], aes(x=N_tf2peak, y=N_tf2gene, label=TF), size=3,  max.overlaps=100) +
  labs(x="Number of correlation events with peaks", y="Number of correlation events with genes") +
  theme_classic() +
  theme(
    # axis.ticks.y = element_blank(),
    # axis.text = element_text(size=rel(0.75), color="black")
  )


pdf(sprintf("%s/scatterplot_ncorrpeaks_vs_ncorrgenes.pdf",io$outdir), width = 6, height = 4)
print(p)
dev.off()




to.plot <- foobar[N_tf2gene>=50] %>% .[,Nlog_tf2gene:=log10(N_tf2gene)]

p <- ggscatter(to.plot, x="fraction_positive_tf2peak", y="fraction_positive_tf2gene", size=3) +
          # add="reg.line", add.params = list(color="blue", fill="lightgray"), conf.int=TRUE) +
  # stat_cor(method = "pearson") +
  scale_color_continuous(low="red", high="green") +
  geom_hline(yintercept=0.5) + geom_vline(xintercept=0.5) +
  ggrepel::geom_text_repel(data=to.plot[fraction_positive_tf2peak<0.40 & fraction_positive_tf2gene>0.60], aes(x=fraction_positive_tf2peak, y=fraction_positive_tf2gene, label=TF), size=3,  max.overlaps=100) +
  ggrepel::geom_text_repel(data=to.plot[fraction_positive_tf2peak<0.40 & fraction_positive_tf2gene<0.40], aes(x=fraction_positive_tf2peak, y=fraction_positive_tf2gene, label=TF), size=3,  max.overlaps=100) +
  ggrepel::geom_text_repel(data=to.plot[fraction_positive_tf2peak>0.80 & fraction_positive_tf2gene>0.80], aes(x=fraction_positive_tf2peak, y=fraction_positive_tf2gene, label=TF), size=3,  max.overlaps=100) +
  ggrepel::geom_text_repel(data=to.plot[fraction_positive_tf2peak>0.75 & fraction_positive_tf2gene<0.50], aes(x=fraction_positive_tf2peak, y=fraction_positive_tf2gene, label=TF), size=3,  max.overlaps=100) +
  labs(x="Fraction of positive correlation events with peaks", y="Fraction of positive correlation events with genes") +
  theme_classic() +
  theme(
    # axis.ticks.y = element_blank(),
    axis.text = element_text(size=rel(0.75), color="black")
  )

pdf(sprintf("%s/scatterplot_fractioncorrpeaks_vs_fractioncorrgenes.pdf",io$outdir), width = 6, height = 4)
print(p)
dev.off()

# tf2gene.dt[TF=="RUNX1" & abs(cor)>0.5] %>% View


