########################
## Load ArchR Project ##
########################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/atac/archR/load_archR_project.R")
} else {
  stop("Computer not recognised")
}

#####################
## Define settings ##
#####################

# I/O
io$outdir <- paste0(io$basedir,"/results/atac/archR/trajectories/blood_trajectory")

# Options

opts$celltypes = c(
  # "Mixed_mesoderm",
  "Haematoendothelial_progenitors",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3"
)

opts$min.expr <- 0.1

#####################
## Update metadata ##
#####################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & pass_rnaQC==TRUE] %>%
  .[celltype.mapped%in%opts$celltypes & sample%in%opts$samples] %>%
  .[,celltype.mapped:=factor(celltype.mapped,levels=opts$celltypes)]

##################
## Subset ArchR ##
##################

ArchRProject.filt <- ArchRProject[sample_metadata$archR_cell,]

###############
## Load data ##
###############

# Load RNA-based trajectory
io$pseudotime <- "/Users/ricard/data/gastrulation_multiome_10x/results/rna/trajectories/blood_trajectory/blood_trajectory.txt.gz"
trajectory.dt <- fread(io$pseudotime) %>% .[,rank_V1:=NULL]

# Load highly variable along the NMP trajectory
io$hvgs <- "/Users/ricard/data/gastrulation_multiome_10x/results/rna/trajectories/blood_trajectory/hvgs.rds"
hvgs <- readRDS(io$hvgs)


##########################
## Load chromVAR scores ##
##########################

atac.deviation.se <- readRDS(io$archR.deviations.se)
atac.deviation.mtx <- atac.deviation.se %>%
  .[,colnames(atac.deviation.se)%in%sample_metadata$archR_cell] %>% 
  assay(.,"z")
dim(atac.deviation.mtx)

######################
## Load peak matrix ##
######################

atac.peak.se <- getMatrixFromProject(ArchRProject.filt, useMatrix = "PeakMatrix", binarize = T)#@assays@data[[1]]

# Define peak names
peak_names <- rowRanges(atac.peak.se) %>% as.data.table %>% .[,id:=sprintf("%s_%s_%s",seqnames,start,end)] %>% .$id
rownames(atac.peak.se) <- peak_names

# Rename cells 
colnames(atac.peak.se) <- sample_metadata %>% 
  .[archR_cell%in%colnames(atac.peak.se)] %>%
  setkey(archR_cell) %>% .[colnames(atac.peak.se)] %>% .$cell

# Prepare data
trajectory.mtx <- trajectory.dt %>%
  # .[,cell:=gsub("rep1_","rep1#",cell)] %>%
  # .[,cell:=gsub("rep2_","rep2#",cell)] %>%
  .[cell%in%colnames(atac.peak.se)] %>%
  .[,c("cell","V1")] %>%
  matrix.please

# Create peak matrix 
peak.mtx <- assay(atac.peak.se)

# Filter cells
peak.mtx <- peak.mtx[,colnames(peak.mtx)%in%rownames(trajectory.mtx)] %>% .[,rownames(trajectory.mtx)]

# Filter peaks
peak.mtx <- peak.mtx[Matrix::rowSums(peak.mtx)>50,]

###############
## Denoising ##
###############


##########################
## Correlation analysis ##
##########################

cor.dt <- psych::corr.test(t(as.matrix(peak.mtx)),trajectory.mtx)[c("r", "p")] %>%
  do.call("cbind",.) %>% as.data.table(keep.rownames = T) %>% setnames(c("peak","r","p")) %>%
  .[,"padj_fdr" := list(p.adjust(p, method="fdr"))] %>%
  .[, sig := padj_fdr <= 0.10] %>% 
  setorder(padj_fdr, na.last = T)

##################
## Volcano plot ##
##################

negative_hits <- cor.dt[sig==TRUE & r<0,peak]
positive_hits <- cor.dt[sig==TRUE & r>0,peak]
all <- nrow(cor.dt)

xlim <- max(abs(cor.dt$r), na.rm=T)
ylim <- max(-log10(cor.dt$padj_fdr+1e-100), na.rm=T)

p <- ggplot(cor.dt, aes(x=r, y=-log10(padj_fdr+1e-100))) +
  # geom_hline(yintercept = -log10(opts$threshold_fdr), color="blue") +
  geom_segment(aes(x=0, xend=0, y=0, yend=ylim-1), color="orange", size=0.5) +
  ggrastr::geom_point_rast(aes(color=sig, size=sig)) +
  scale_color_manual(values=c("black","red")) +
  scale_size_manual(values=c(0.75,1.25)) +
  scale_x_continuous(limits=c(-xlim-0.5,xlim+0.5)) +
  scale_y_continuous(limits=c(0,ylim+6)) +
  annotate("text", x=0, y=ylim+6, size=4, label=sprintf("(%d)", all)) +
  annotate("text", x=-xlim-0.5, y=ylim+6, size=4, label=sprintf("%d (-)",length(negative_hits))) +
  annotate("text", x=xlim+0.5, y=ylim+6, size=4, label=sprintf("%d (+)",length(positive_hits))) +
  labs(x="Pearson correlation", y=expression(paste("-log"[10],"(p.value)"))) +
  theme_classic() +
  theme(
    axis.text = element_text(size=rel(0.75), color='black'),
    axis.title = element_text(size=rel(1.0), color='black'),
    legend.position="none"
  )

# pdf(sprintf("%s/volcano_plots/volcano_pearson_correlation.pdf",io$outdir), width = 9, height = 6)
# png(sprintf("%s/volcano_plots/volcano_pearson_correlation.png",io$outdir), width = 800, height = 500)
print(p)
# dev.off()

# fwrite(cor.dt, sprintf("%s/correlation_results_blood.txt.gz",io$outdir), quote=F, sep="\t", na="NA")

########################################
## Scatterplot of individual examples ##
########################################

peaks.to.plot <-  c(cor.dt %>% .[r>0,peak] %>% head(n=6), cor.dt %>% .[r<0,peak] %>% head(n=3))

to.plot <- peak.mtx[peaks.to.plot,] %>% 
  as.matrix %>% as.data.table(keep.rownames = T) %>%
  setnames("rn","peak") %>% melt(id.vars="peak", variable.name="cell") %>%
  merge(trajectory.dt,by="cell") %>%
  merge(sample_metadata[,c("cell","celltype.mapped")])

ggscatter(to.plot, x="V1", y="value", color="celltype.mapped", size=0.1) +
  facet_wrap(~peak, scales="fixed") +
  stat_smooth(method="lm", color="black", alpha=0.75, span=0.5) +
  scale_color_manual(values=opts$celltype.colors) +
  labs(x="", y="Peak accessibility") +
  theme(
    axis.text = element_text(size=rel(0.75)),
    legend.position = "none"
  )

######################
## Motif enrichment ##
######################

# Load Motif annotations 
peak_annotation.se <- getPeakAnnotation(ArchRProject.filt, name="Motif")[["Matches"]] %>% readRDS
rownames(peak_annotation.se) <- rownames(atac.peak.se)
peak_annotation.se.filt <- peak_annotation.se[rownames(atac.peak.se.filt),]

foreground.peaks <- cor.dt[sig==T & r>0.15,peak]
# background.peaks <- cor.dt[sig==T & r<(-0.15),peak]
background.peaks <- cor.dt[!peak%in%foreground.peaks,peak]

i <- colnames(peak_annotation.se.filt)[1]
i="Klf1_214"

for (i in colnames(peak_annotation.se.filt)) {
  print(i)
  foreground.nmatches <- assay(peak_annotation.se.filt)[,i][foreground.peaks] %>% sum
  background.nmatches <- assay(peak_annotation.se.filt)[,i][background.peaks] %>% sum
  foreground.total = length(foreground.peaks)
  background.total = length(background.peaks)
  
  p <- phyper(
    q = foreground.nmatches, 
    m = length(foreground.peaks), 
    n = length(background.peaks), 
    k = foreground.nmatches+background.nmatches
    )
  print(p)
  
}


peak_annotation.se.filt
pOut <- data.frame(
  feature = rownames(atac.peak.se.filt),
  CompareFrequency = matchCompareTotal,
  nCompare = nrow(matchCompare),
  CompareProportion = matchCompareTotal/nrow(matchCompare),
  BackgroundFrequency = matchBackgroundTotal,
  nBackground = nrow(matchBackground),
  BackgroundProporition = matchBackgroundTotal/nrow(matchBackground)
)

#Enrichment
pOut$Enrichment <- pOut$CompareProportion / pOut$BackgroundProporition

#Get P-Values with Hyper Geometric Test
pOut$mlog10p <- lapply(seq_len(nrow(pOut)), function(x){
  p <- -phyper(pOut$CompareFrequency[x] - 1, # Number of Successes the -1 is due to cdf integration
               pOut$BackgroundFrequency[x], # Number of all successes in background
               pOut$nBackground[x] - pOut$BackgroundFrequency[x], # Number of non successes in background
               pOut$nCompare[x], # Number that were drawn
               lower.tail = FALSE, log.p = TRUE)# P[X > x] Returns LN must convert to log10
  return(p/log(10))
}) %>% unlist %>% round(4)

#Minus Log10 Padj
pOut$mlog10Padj <- pmax(pOut$mlog10p - log10(ncol(pOut)), 0)
pOut <- pOut[order(pOut$mlog10p, decreasing = TRUE), , drop = FALSE]

pOut

###########################################################
## Calculate chromVAR scores using only correlated peaks ##
###########################################################

library(chromVAR)
library(BSgenome.Mmusculus.UCSC.mm10)

# Prepare data for chromVAR
atac.peak.se.filt <- atac.peak.se[rownames(peak.mtx),]
atac.peak.se.filt <- addGCBias(atac.peak.se.filt, genome = BSgenome.Mmusculus.UCSC.mm10)
assayNames(atac.peak.se.filt) <- "counts"
head(rowData(atac.peak.se.filt))

rownames(peak_annotation.se)

# calculate background peaks

# Compute deviations
chromvar.deviations.se <- computeDeviations(
  object = atac.peak.se.filt, 
  annotations = peak_annotation.se.filt
)
head(rowData(chromvar.deviations.se))

# save
# chromvar.deviations.se <- getMatrixFromProject(ArchRProject.filt, "DeviationMatrix_JASPAR")
saveRDS(chromvar.deviations.se, sprintf("%s/deviations_summarized_experiment.rds",io$outdir))
  
###########################
## Plot chromVAR results ##
###########################

# save pseudobulk data.table
chromvar.dt <- assay(chromvar.deviations.se,"z") %>% as.matrix %>% as.data.frame %>%
  as.data.table(keep.rownames = T) %>% setnames("rn","motif") %>%
  .[,motif:=stringr::str_split(motif,"_") %>% map_chr(1) %>% gsub("z:","",.)] %>%
  melt(id.vars="motif", variable.name="cell", value.name="chromvar_zscore") %>%
  merge(sample_metadata[,c("cell","celltype.mapped")], by="cell") %>%
  merge(trajectory.dt,by="cell") %>%
  setnames("celltype.mapped","celltype")

chromvar.dt[,mean(chromvar_zscore),by=c("motif","celltype")] %>% View

# fwrite(dt, sprintf("%s/JASPAR/archr_motif_chromvar_scores_celltype.txt.gz",io$outdir))

motifs.to.plot <- chromvar.dt$motif

for (i in motifs.to.plot) {
  
  to.plot <- chromvar.dt[motif==i]
  
  p1 <- ggplot(to.plot, aes(x=V1, y=chromvar_zscore)) +
    geom_point(aes(fill=celltype), size=1.25, shape=21, stroke=0.1) +
    stat_smooth(method="loess", color="black", alpha=0.75, span=0.5) +
    geom_rug(aes(color=celltype), sides="b") +
    scale_color_manual(values=opts$celltype.colors) +
    scale_fill_manual(values=opts$celltype.colors) +
    guides(fill=F, color=F) +
    labs(x="Pseudotime", y=i) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.title = element_blank(),
      legend.position="top"
    )
  
  p2 <- ggboxplot(to.plot, x="celltype", y="chromvar_zscore", fill="celltype", outlier.shape=NA) +
    scale_fill_manual(values=opts$celltype.colors) +
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
  
  pdf(sprintf("%s/individual_motifs/%s_chromvar_vs_pseudotime_blood.pdf",io$outdir,i), width=8, height=3)
  print(p)
  dev.off()
}

######################
## Motif enrichment ##
######################

##########
## Test ##
##########

# rowData(peak_annotation.se) %>% head

opts$min.motif.score <- 7

motif2peak.positions.se <- getPeakAnnotation(ArchRProject.filt, name="Motif")[["Positions"]] %>% readRDS

motif2peak.positions.dt <- names(motif2peak.positions.se) %>% map(function(i) {
  motif2peak.positions.se[[i]] %>%
    as.data.table() %>%
    setnames("seqnames","chr") %>%
    .[,motif:=factor(i,levels=names(motif2peak.positions.se))] %>%
    return
}) %>% rbindlist# %>%
  # .[,idx:=sprintf("%s:%s-%s",chr,start,end)] %>%
  # .[,motif:=factor(motif,levels=motifs)] %>%
  # .[,chr:=factor(chr,levels=opts$chr)]

motif2peak.positions.dt[,.N,by="motif"] %>% View
foo <- motif2peak.positions.dt[motif=="Klf1_214"] %>% setorder(score) %>% head
foo <- motif2peak.positions.dt[motif=="Klf1_214"] %>% setorder(score) %>% tail

seq <- getSeq(Mmusculus, foo$chr, foo$start, foo$end+1)
seq
