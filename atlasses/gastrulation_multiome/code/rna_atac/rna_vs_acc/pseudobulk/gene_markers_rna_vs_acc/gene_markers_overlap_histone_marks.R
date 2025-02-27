#####################
## Define settings ##
#####################

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))


# I/O
io$tss <- "/Users/argelagr/data/mm10_regulation/genes/TSS.bed"
io$housekeeping.genes <-"/Users/argelagr/data/genesets/manual_genesets/housekeeping/housekeeping.tsv"
io$chip.basedir <- "/Users/argelagr/data/mm10_regulation/Xiang2019/chip/liftover/bedgraph"
io$cpg.islands <- "/Users/argelagr/data/mm10_regulation/genes/old/prom_2000_2000/prom_2000_2000_cgi.bed"
io$cpg_density_promoters <- "/Users/argelagr/data/mm10_regulation/CGI/cpg_density_prom_2000_0.txt.gz"
io$outdir <- file.path(io$basedir,"results_new/rna_atac/rna_vs_acc/pseudobulk/gene_markers_rna_vs_acc/histone_profiles"); dir.create(io$outdir, showWarnings = F)

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
  "Erythroid1",
  "Erythroid2",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
  # "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm"
  # "Visceral_endoderm"
  # "ExE_endoderm",
  # "ExE_ectoderm"
  # "Parietal_endoderm"
)

#######################
## Load gene classes ##
#######################

# marker_genes.dt <- fread(io$rna.atlas.marker_genes) %>%
io$marker_genes_rna <- io$rna.atlas.marker_genes
io$marker_genes_rna <- file.path(io$basedir,"results_new/rna/differential/marker_genes/marker_genes_up.txt.gz")

marker_genes.dt <- fread(io$marker_genes_rna) %>% .[celltype%in%opts$celltypes]

houseekeping.dt <- fread(io$housekeeping.genes, header=F) %>%
  setnames(c("ens_id","gene"))

pluripotency_genes <- c("Dppa2","Dppa3","Dppa4","Morc1","Tex19.1","Zfp981","Dppa5a","Zfp42")

foo <- copy(houseekeping.dt) %>% .[,class:="Housekeping"] %>% .[,class2:=class] %>% .[,c("gene","class")]
bar <- data.table(class = rep("Pluripotency",length(pluripotency_genes)), gene = pluripotency_genes) %>% .[,class2:=class] %>% .[,c("gene","class")]
baz <- marker_genes.dt[,c("celltype","gene")] %>% .[,class:="Marker genes"] %>% setnames("celltype","class2") %>% .[,c("gene","class")]
gene2class <- unique(rbindlist(list(foo,bar,baz))); rm(foo,bar,baz)

########################
## Load gene metadata ##
########################

# Load gene metadata
gene_metadata.dt <- fread(io$gene_metadata) %>% setnames("symbol","gene")

# Load promoter coordinates
promoter.dt <- fread(io$tss) %>% 
  setnames(c("chr","start","end","ens_id","score","strand")) %>% 
  merge(gene_metadata.dt[,c("gene","ens_id")],by="ens_id") %>%
  merge(gene2class,by="gene") %>%
  .[strand=="+",start:=start-2000] %>% .[strand=="-",end:=end+2000] %>%
  .[,c("chr","start","end","gene")]%>% setkey(chr,start,end)

###########################
## Load histone ChIP-seq ##
###########################

H3K27ac.dt <- fread(file.path(io$chip.basedir,"E65Epi_H3K27ac.bed.gz")) %>%
  setnames(c("chr","start","end","value")) %>%
  setkey(chr,start,end)

H3K27me3.dt <- fread(file.path(io$chip.basedir,"E65Epi_H3K27me3.bed.gz")) %>%
  setnames(c("chr","start","end","value")) %>%
  setkey(chr,start,end)

H3K4me3.dt <- fread(file.path(io$chip.basedir,"E65Epi_H3K4me3.bed.gz")) %>%
  setnames(c("chr","start","end","value")) %>%
  setkey(chr,start,end)

###################################################
## Overlap with ChIP-seq data from histone marks ##
###################################################

H3K27ac_quantification.dt <- promoter.dt %>% 
  foverlaps(H3K27ac.dt) %>%
  .[,mean(value,na.rm=T),by="gene"]

H3K27me3_quantification.dt <- promoter.dt %>% 
  foverlaps(H3K27me3.dt) %>%
  .[,mean(value,na.rm=T),by="gene"]

H3K4me3_quantification.dt <- promoter.dt %>% 
  foverlaps(H3K4me3.dt) %>%
  .[,mean(value,na.rm=T),by="gene"]

##########
## Plot ##
##########

to.plot <- rbindlist(
  list(
    H3K27ac_quantification.dt[,mark:=as.factor("H3K27ac")] %>% .[,background:=mean(H3K27ac.dt$value)],
    H3K4me3_quantification.dt[,mark:=as.factor("H3K4me3")] %>% .[,background:=mean(H3K4me3.dt$value)],
    H3K27me3_quantification.dt[,mark:=as.factor("H3K27me3")] %>% .[,background:=mean(H3K27me3.dt$value)]
  )
  ) %>% merge(unique(gene2class[,c("gene","class")]),by="gene") %>%
# ) %>% merge(unique(gene2class[,c("gene","class","class2")]),by="gene", allow.cartesian=T) %>%
  .[,signal:=V1/background]

p <- ggboxplot(to.plot, x="class", y="signal", fill="class", facet="mark", outlier.shape=NA) +
  labs(x="", y="ChIP-seq signal") +
  coord_cartesian(ylim=c(0,9)) +
  scale_fill_brewer(palette="Dark2") +
  theme(
    strip.background = element_rect(colour="black", fill=NA),
    legend.position = "none",
    legend.direction = "vertical",
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=rel(0.75)),
    axis.ticks.x = element_blank()
  )

pdf(file.path(io$outdir,"markers_genes_chipseq.pdf"), width = 6.5, height = 4)
print(p)
dev.off()

# to.plot.scatter <- to.plot %>% dcast(gene+class+class2~mark, value.var="signal")
# ggscatter(to.plot.scatter, x="H3K27ac", y="H3K4me3", fill="class2", shape=21) +
#   facet_wrap(~class, scales="fixed") +
#   # labs(x="", y="ChIP-seq signal") +
#   # coord_cartesian(ylim=c(0,9)) +
#   scale_fill_manual(values=opts$celltype.colors) +
#   theme(
#     legend.position = "none",
#     legend.direction = "vertical",
#     legend.title = element_blank(),
#     axis.text = element_text(size=rel(0.75))
#   )


#################
## CpG density ##
#################

cpg_density_per_gene <- fread(io$cpg_density_promoters) %>%
  merge(gene_metadata.dt[,c("ens_id","gene")], by="ens_id") %>%
  merge(unique(gene2class[,c("gene","class")]), by="gene") %>%
  .[,cpg_density:=100*cpg_density]


p <- ggboxplot(cpg_density_per_gene, x="class", y="cpg_density", fill="class", outlier.shape=NA) +
  scale_fill_brewer(palette="Dark2") +
  labs(x="", y="CpG density (%)") +
  # scale_fill_manual(values=opts$celltype.colors) +
  coord_cartesian(ylim=c(0,6)) +
  stat_summary(fun.data = give.n, geom = "text", size=2.5) +
  theme(
    axis.text.y = element_text(color="black", size=rel(0.75)),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )

pdf(file.path(io$outdir,"markers_genes_cpg_density.pdf"), width = 5, height = 3)
print(p)
dev.off()

##############################
## Overlap with CpG islands ##
##############################

cpg_islands.dt <- fread(io$cpg.islands) %>% setnames(c("chr","start","end","ens_id","score","strand")) %>%
  merge(gene_metadata.dt[,c("ens_id","gene")], by="ens_id")

to.plot <- gene2class[,c("gene","class")] %>% unique %>%
  .[,CGI:=gene%in%cpg_islands.dt$gene] %>%
  .[,mean(CGI),by="class"]

ggbarplot(to.plot, x="class", y="V1", fill="class", position = position_dodge(0.9)) +
  labs(x="", y="Fraction of genes that overlap with a CGI") +
  guides(x = guide_axis(angle = 90)) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    axis.title.y = element_text(size=rel(0.65)),
    axis.text.y = element_text(size=rel(0.75)),
    axis.text.x = element_blank()
  )

pdf(sprintf("%s/markers_genes_cpg_density.pdf",io$outdir), width = 5, height = 3)
print(p)
dev.off()

