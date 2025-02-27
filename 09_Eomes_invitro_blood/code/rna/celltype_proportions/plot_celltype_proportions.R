here::i_am("rna/celltype_proportions/plot_celltype_proportions.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--celltype_label', type="character", help='Cell type label')
p$add_argument('--outdir',          type="character",                               help='Output file')

args <- p$parse_args(commandArgs(TRUE))

args$samples = opts$samples

## START TEST ##
# args$metadata <- file.path(io$basedir,"results/rna/mapping/sample_metadata_after_mapping.txt.gz")
# args$celltype_label <- "celltype"
# args$outdir <- file.path(io$basedir,"results/rna/celltype_proportions")
## END TEST ##

height = 8
if(args$celltype_label == 'celltype_extended.mapped_mnn'){
    opts$celltype.colors  = opts$celltype_extended.colors
    height = 14
}

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE & sample%in%args$samples & !is.na(eval(as.name(args$celltype_label)))] %>%
  setnames(args$celltype_label,"celltype")

################################################
## Calculate cell type proportions per sample ##
################################################

to.plot <- sample_metadata %>%
  .[,N:=.N,by="sample"] %>%
  .[,.(N=.N, celltype_proportion=.N/unique(N)),by=c("sample","day","celltype")] %>%
  setorder(sample)  %>% .[,sample:=factor(sample,levels=args$samples)]


#########################
## Horizontal barplots ##
#########################

# Define colours and cell type order
opts$celltype.colors <- opts$celltype.colors[names(opts$celltype.colors) %in% unique(to.plot$celltype)]
to.plot[,celltype:=factor(celltype, levels=names(opts$celltype.colors))]

for (i in unique(to.plot$day)) {
  p <- ggplot(to.plot[day==i], aes(x=celltype, y=celltype_proportion)) +
    geom_bar(aes(fill=celltype), stat="identity", color="black") +
    scale_fill_manual(values=opts$celltype.colors) +
    facet_wrap(~sample, nrow=1, scales="free_x") +
    coord_flip() +
    labs(y="Proportion of cells") +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(color="black", size=rel(0.5)),
      axis.title.x = element_text(color="black", size=rel(0.9)),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=rel(1), color="black"),
      axis.text.x = element_text(size=rel(1), color="black")
    )
  
  pdf(sprintf("%s/celltype_proportions_%s_horizontal_barplots_per_sample.pdf",args$outdir,i), width=10, height=height)
  print(p)
  dev.off()
}


################################################
## Calculate cell type proportions by genotype ##
################################################

to.plot <- sample_metadata %>%
  .[,N:=.N,by=c("genotype","day")] %>%
  .[,.(N=.N, celltype_proportion=.N/unique(N)),by=c("genotype","day","celltype")] 

# Define colours and cell type order
opts$celltype.colors <- opts$celltype.colors[names(opts$celltype.colors) %in% unique(to.plot$celltype)]
to.plot[,celltype:=factor(celltype, levels=names(opts$celltype.colors))]

for (i in unique(to.plot$day)) {
  p <- ggplot(to.plot[day==i], aes(x=celltype, y=celltype_proportion)) +
    geom_bar(aes(fill=celltype), stat="identity", color="black") +
    scale_fill_manual(values=opts$celltype.colors) +
    facet_wrap(~genotype, nrow=1, scales="free_x") +
    coord_flip() +
    labs(y="Proportion of cells") +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(color="black", size=rel(0.5)),
      axis.title.x = element_text(color="black", size=rel(0.9)),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=rel(1), color="black"),
      axis.text.x = element_text(size=rel(1), color="black")
    ) 

  pdf(sprintf("%s/celltype_proportions_%s_horizontal_barplots_per_genotype.pdf",args$outdir,i), width=10, height=height)
  print(p)
  dev.off()
  
  p <- ggplot(to.plot[day==i], aes(x=celltype, y=celltype_proportion)) +
    geom_bar(aes(fill=genotype), stat="identity", color="black",  position = "dodge") +
    scale_fill_manual(values=opts$genotype.colors) +
    coord_flip() +
    labs(y="Proportion of cells") +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(color="black", size=rel(0.5)),
      axis.title.x = element_text(color="black", size=rel(0.9)),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=rel(1), color="black"),
      axis.text.x = element_text(size=rel(1), color="black")
    ) 
  pdf(sprintf("%s/celltype_proportions_%s_horizontal_barplots_per_genotype_dodged.pdf",args$outdir,i), width=10, height=height)
  print(p)
  dev.off()
}

fwrite(data.table(a=''), sprintf("%s/finished.txt",args$outdir))