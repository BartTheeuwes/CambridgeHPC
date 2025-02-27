
#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/utils.R")
} else {
  stop("Computer not recognised")
}
io$outdir <- paste0(io$basedir,"/results/atac/archR/PijuanSala2020_comparison")
opts$matrix <- "GeneScoreMatrix"

opts$celltypes = c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  "PGC",
  "Anterior_Primitive_Streak",
  "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  "Mixed_mesoderm",
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
  "Erythroid1",
  "Erythroid2",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
  "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm",
  "Visceral_endoderm",
  "ExE_endoderm",
  "ExE_ectoderm",
  "Parietal_endoderm"
)

opts$match.celltypes <- c(
  "Surface_ectoderm"= "Surface_ectoderm",
  # "Notochord"= "Notochord",
  "Gut"= "Gut",
  "Cardiomyocytes"= "Cardiomyocytes",
  "Mid_Hindbrain"= "Forebrain_Midbrain_Hindbrain",
  "Endothelium"= "Endothelium",
  "Paraxial_mesoderm"= "Paraxial_mesoderm",
  "Spinal_cord"= "Spinal_cord",
  "Somitic_mesoderm"= "Somitic_mesoderm",
  "Erythroid" = "Erythroid3",
  # "Neural_crest"= "Neural_crest",
  # "Mixed_mesoderm"= "Mixed_mesoderm",
  "NMP"= "NMP",
  # "Forebrain"= "Forebrain_Midbrain_Hindbrain",
  "ExE_endoderm"= "ExE_endoderm",
  "Allantois"= "Allantois",
  "Mesenchyme"= "Mesenchyme",
  "Pharyngeal_mesoderm" = "Pharyngeal_mesoderm"
)

opts$samples <- c(
  # "E7.5_rep1",
  # "E7.5_rep2",
  "E8.0_rep1",
  "E8.0_rep2",
  "E8.5_rep1",
  "E8.5_rep2"
)

##########################
## Load sample metadata ##
##########################

# io$metadata <- paste0(io$basedir,"/results/atac/archR/celltype_assignment/sample_metadata_after_archR.txt.gz")

multiome_sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & doublet_call==FALSE & !is.na(celltype.predicted)] %>%
  .[sample%in%opts$samples]

pijuansala_sample_metadata <- fread(paste0(io$pijuansala.basedir,"/cell_metadata.csv.gz")) %>%
  .[,c("cell", "nuclei_type", "celltype")] %>%
  .[,celltype:=stringr::str_replace_all(celltype," ","_")] %>%
  .[,celltype:=stringr::str_replace_all(celltype,"/","_")] %>%
  .[!is.na(celltype)]

##############################
## Define cell types to use ##
##############################

opts$min.cells <- 50

opts$multiome_celltypes <- multiome_sample_metadata[,.N,by=c("celltype.predicted")] %>%
  .[N>opts$min.cells] %>% .[,celltype.predicted]
table(multiome_sample_metadata$celltype.predicted)

opts$pijuansala_celltypes <- pijuansala_sample_metadata[,.N,by=c("celltype")] %>%
  .[N>opts$min.cells] %>% .[,celltype]
table(pijuansala_sample_metadata$celltype)

#######################
## Load marker genes ##
#######################

marker_genes.dt <- fread(io$rna.atlas.marker_genes) %>%
  .[celltype%in%opts$celltypes]

# Filter genes
# marker_genes.dt <- marker_genes.dt[grep("Rik",gene,invert = T)]

#######################################
## Load Multiome pseudobulk matrices ##
#######################################

file <- sprintf("%s/pseudobulk/pseudobulk_%s_summarized_experiment.rds",io$archR.directory,opts$matrix)

multiome_se <- readRDS(file)[,opts$multiome_celltypes]
rownames(multiome_se) <- rowData(multiome_se)$name
multiome_matrix <- assay(multiome_se, opts$matrix)
dim(multiome_matrix)

#########################################
## Load PijuanSala pseudobulk matrices ##
#########################################

file <- sprintf("%s/pseudobulk_matrices/pseudobulk_%s_summarized_experiment.rds",io$pijuansala.archR.directory,opts$matrix)

pijuansala_se <- readRDS(file)[,opts$pijuansala_celltypes]
rownames(pijuansala_se) <- rowData(pijuansala_se)$name
pijuansala_matrix <- assay(pijuansala_se, opts$matrix)
# rownames(pijuansala_matrix) <- rowData(pijuansala_se) %>% as.data.table %>%
#   .[,id:=sprintf("%s_%d_%d",seqnames,start,end)] %>% .[,id]
# colnames(pijuansala_matrix) <- paste0(colnames(pijuansala_matrix),"_PijSal")
dim(pijuansala_matrix)

#############################
## Define feature metadata ##
#############################

feature_metadata <- rowData(multiome_se) %>% as.data.table %>%
  setnames("seqnames","chr") %>% 
  .[,idx:=NULL] %>%
  .[,id:=sprintf("%s_%d_%d",chr,start,end)] %>%
  .[,length:=abs(end-start)]

#######################
## Feature selection ##
#######################

# Remove very long genes
opts$max.length <- 2e6
feature_metadata <- feature_metadata[length<=opts$max.length]

# Remove manually some features
feature_metadata <- feature_metadata[grep("Rik|Mir|Gm",name, invert=T)]

# Filter highly variable features
opts$num_variable_features <- 5000
feature_metadata %>%
  .[,c("multiome_mean","PijSal_mean"):=list(rowMeans(multiome_matrix[feature_metadata$name,]), rowMeans(pijuansala_matrix[feature_metadata$name,]))] %>%
  .[,c("multiome_variance","PijSal_variance"):=list(apply(multiome_matrix[feature_metadata$name,],1,var), apply(pijuansala_matrix[feature_metadata$name,],1,var))]
feature_metadata.filt <- feature_metadata %>% setorder(-multiome_variance) %>% head(opts$num_variable_features)

# Subset matrices
# multiome_matrix_filt <- multiome_matrix[feature_metadata.filt$name,]
# pijuansala_matrix_filt <- pijuansala_matrix[feature_metadata.filt$name,]

genes.to.use <- intersect(
  unique(marker_genes.dt$gene)[unique(marker_genes.dt$gene)%in%rownames(multiome_matrix)],
  unique(marker_genes.dt$gene)[unique(marker_genes.dt$gene)%in%rownames(pijuansala_matrix)]
)
multiome_matrix_filt <- multiome_matrix[genes.to.use,]
pijuansala_matrix_filt <- pijuansala_matrix[genes.to.use,]

# sanity check
cor(feature_metadata$multiome_mean, feature_metadata$PijSal_mean)
cor(feature_metadata$multiome_variance, feature_metadata$PijSal_variance)


##########################
## Scatterplot per gene ##
##########################

stopifnot(names(opts$match.celltypes) %in% colnames(pijuansala_matrix_filt))
stopifnot(opts$match.celltypes %in% colnames(multiome_matrix_filt))

genes.to.plot <- c("Foxa2","Hoxa9","Gbx2")

for (i in genes.to.plot) {
  
  to.plot <- data.table(
    celltype = opts$match.celltypes,
    PijuanSala = pijuansala_matrix_filt[i,names(opts$match.celltypes)] %>% minmax.normalisation,
    Multiome = multiome_matrix_filt[i,opts$match.celltypes] %>% minmax.normalisation
  ) %>% setorder(-Multiome)
  
  p <- ggplot(to.plot, aes(x=PijuanSala, y=Multiome)) +
    geom_point(aes(fill=celltype), color="black", shape=21, size=3, stroke=0.5) +
    stat_cor(method = "pearson") +
    stat_smooth(method="lm", color="black", alpha=0.1, size=0.2) +
    scale_fill_manual(values=opts$celltype.colors) +
    scale_x_continuous(breaks=c(0,0.5,1)) + scale_y_continuous(breaks=c(0,0.5,1)) +
    ggrepel::geom_text_repel(aes(label=celltype), segment.size=0.15, size=2.5, data=head(to.plot,n=5)) +
    labs(x="Gene accessibility (PijuanSala2020)", y="Gene accessibility (This study)", title=i) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(color="black", size=rel(0.7)),
      axis.title = element_text(color="black", size=rel(0.9)),
      legend.position = "none"
    )
  
  pdf(sprintf("%s/scatterplots/per_gene/scatter_%s_%s_PijuanSala2020_comparison.pdf",io$outdir,opts$matrix,i), width=5, height=5)
  print(p)
  dev.off()
  
}

##############################
## Scatterplot per celltype ##
##############################

stopifnot(names(opts$match.celltypes) %in% colnames(pijuansala_matrix_filt))
stopifnot(opts$match.celltypes %in% colnames(multiome_matrix_filt))
# opts$common.celltypes <- intersect(opts$pijuansala_celltypes,opts$multiome_celltypes)

for (i in names(opts$match.celltypes)) {
  
  to.plot <- data.table(
    gene = rownames(pijuansala_matrix_filt),
    PijuanSala = pijuansala_matrix_filt[,i],
    Multiome = multiome_matrix_filt[,opts$match.celltypes[[i]]]
  ) %>% setorder(-Multiome)
  
  p <- ggplot(to.plot, aes(x=PijuanSala, y=Multiome)) +
    geom_point(color="black", fill="gray90", alpha=0.75, shape=21, size=1.5, stroke=0.5) +
      stat_cor(method = "pearson") +
      # geom_abline(slope=1, intercept=0, linetype="dashed", color="black") +
      stat_smooth(method="lm", color="black", alpha=0.15, size=0.5) +
      ggrepel::geom_text_repel(aes(label=gene), segment.size=0.15, size=3, data=head(to.plot,n=10)) +
      labs(x="Gene accessibility (PijuanSala2020)", y="Gene accessibility (This study)", title=i) +
      theme_classic() +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(color="black", size=rel(0.8)),
        legend.position = "none"
      )
  
  pdf(sprintf("%s/scatterplots/scatter_%s_%s_PijuanSala2020_comparison.pdf",io$outdir,opts$matrix,i), width=5, height=5)
  print(p)
  dev.off()
  
}

###################################
## Correlation analysis per gene ##
###################################

foo <- pijuansala_matrix %>% as.data.table(keep.rownames = T) %>%
  melt(id.vars="rn", variable.name="celltype", value.name="accessibility") %>% setnames("rn","gene") %>%
  .[,celltype:=stringr::str_replace_all(celltype,opts$match.celltypes)] %>%
  .[gene%in%unique(marker_genes.dt$gene)]
# %>% .[,dataset:=as.factor("PijuanSala2020")]

bar <- multiome_matrix %>% as.data.table(keep.rownames = T) %>%
  melt(id.vars="rn", variable.name="celltype", value.name="accessibility") %>% setnames("rn","gene") %>%# .[,dataset:=as.factor("This study")]
  .[gene%in%unique(marker_genes.dt$gene)]

atac.dt <- merge(foo,bar,by=c("gene","celltype"), suffixes=c("_PijuanSala2019","_thisStudy"))
opts$pijuansala_celltypes[!opts$pijuansala_celltypes %in% unique(atac.dt$celltype)]
length(unique(atac.dt$gene))
length(unique(atac.dt$celltype))

cor.dt <- atac.dt %>% 
  .[, cor.test(x=accessibility_PijuanSala2019, y=accessibility_thisStudy)[c("estimate","p.value")], by = c("gene")] %>%
  # .[, para := rep(c("r","p"), .N/2)] %>% data.table::dcast(gene ~ para, value.var = "V1") %>%
  .[, padj_fdr := p.adjust(p.value, method="fdr")] %>%
  .[, sig := padj_fdr <= 0.1] %>%
  .[!is.na(estimate)] %>%
  setorder(padj_fdr)

to.plot <- cor.dt

negative_hits <- to.plot[sig==TRUE & estimate<0,gene]
positive_hits <- to.plot[sig==TRUE & estimate>0,gene]
all <- nrow(to.plot)

xlim <- max(abs(to.plot$estimate), na.rm=T)
ylim <- max(-log10(to.plot$p.value), na.rm=T)

p <- ggplot(to.plot, aes(x=estimate, y=-log10(padj_fdr))) +
  labs(title="", x="Pearson correlation (gene accessibility)", y=expression(paste("-log"[10],"(q.value)"))) +
  # geom_hline(yintercept = -log10(opts$threshold_fdr), color="blue") +
  geom_segment(aes(x=0, xend=0, y=0, yend=ylim-1), color="orange") +
  ggrastr::geom_point_rast(aes(color=sig), size=1) +
  scale_color_manual(values=c("black","red")) +
  scale_x_continuous(limits=c(-xlim-0.05,xlim+0.05)) +
  scale_y_continuous(limits=c(0,ylim+1)) +
  annotate("text", x=0, y=ylim+1, size=5, label=sprintf("(%d)", all)) +
  annotate("text", x=-0.75, y=ylim+1, size=5, label=sprintf("%d (-)",length(negative_hits))) +
  annotate("text", x=0.75, y=ylim+1, size=5, label=sprintf("%d (+)",length(positive_hits))) +
  ggrepel::geom_text_repel(data=head(to.plot[sig==T],n=10), aes(x=estimate, y=-log10(padj_fdr), label=gene), max.overlaps=Inf, size=4) +
  theme_classic() +
  theme(
    # plot.title=element_text(size=28, face='bold', margin=margin(0,0,10,0), hjust=0.5),
    axis.text=element_text(size=rel(0.8), color='black'),
    # axis.title=element_text(size=rel(1.95), color='black'),
    legend.position="none"
  )

pdf(sprintf("%s/PijuanSala2020_comparison_volcano_plot_archR.pdf",io$outdir), width = 6, height = 6)
print(p)
dev.off()

################################
## Correlations per cell type ##
###############################

cor_celltype.dt <- atac.dt %>% 
  .[, cor.test(x=accessibility_PijuanSala2019, y=accessibility_thisStudy)[c("estimate","p.value")], by = c("celltype")] %>%
  .[, padj_fdr := p.adjust(p.value, method="fdr")] %>%
  .[, sig := padj_fdr <= 0.1] %>%
  setorder(padj_fdr)

to.plot <- cor_celltype.dt %>% .[,celltype:=factor(celltype,levels=rev(opts$celltypes))]

p <- ggplot(to.plot, aes(x=celltype, y=estimate, fill=celltype)) +
  geom_bar(stat="identity", width=0.25, alpha=0.9) +
  geom_point(shape=21, size=3.5) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme_classic() +
  coord_flip() +
  labs(x="", y="Correlation coefficient") +
  theme(
    axis.text.y = element_text(colour="black",size=rel(1.25)),
    axis.text.x = element_text(colour="black",size=rel(0.85)),
    axis.title.x = element_text(colour="black",size=rel(1.0)),
    legend.position="none"
  )

pdf(sprintf("%s/PijuanSala2020_comparison_correlation_per_celltype.pdf",io$outdir), width = 5, height = 7)
print(p)
dev.off()


#############################
## Heatmap compare samples ##
#############################

r <- cor(multiome_matrix_filt,pijuansala_matrix_filt)

pdf(sprintf("%s/heatmap_%s_PijuanSala_comparison.pdf",io$outdir,opts$matrix), width=9, height=7)
pheatmap::pheatmap(r)
dev.off()


##################
## PCA (IGNORE) ##
##################

# concatenate
concat_matrix <- cbind(multiome_matrix_filt,pijuansala_matrix_filt) %>% t

# PCA
pca <- prcomp(concat_matrix, rank.=5)
pca.var.explained <- 100*(pca$sdev**2 / sum(pca$sdev**2)) %>% round(4)


tail(sort(apply(multiome_matrix_filt,2,sum)))
tail(sort(apply(pijuansala_matrix_filt,2,sum)))

tail(sort(apply(multiome_matrix_filt,1,sum)))
tail(sort(apply(pijuansala_matrix_filt,1,sum)))

# Plot PCA

celltypes <- c(colnames(multiome_matrix),colnames(pijuansala_matrix))
samples <- c(paste0(colnames(multiome_matrix),"_Multiome"),paste0(colnames(pijuansala_matrix),"_PijSal"))
dataset <- c(rep("Multiome",ncol(multiome_matrix)),rep("PijSal",ncol(pijuansala_matrix)))

to.plot <- pca$x %>% as.data.table %>%
  .[,c("celltype","lineage","dataset"):=list(celltypes,lineages,dataset)]

for (i in c("dataset","celltype")) {
  p <- ggplot(to.plot, aes_string(x="PC1", y="PC2", fill=i)) +
    geom_point(size=3, shape=21, stroke=0.1) +
    theme_classic() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
  if (i=="celltype") {
    p <- p + 
      scale_fill_manual(values=opts$celltype.colors) +
      theme(legend.position="none")
  }
  print(p)
}

# pdf(paste0(io$outdir,"/umap_by_celltype.pdf"), width=6, height=4, useDingbats = F)
# print(p)
# dev.off()
