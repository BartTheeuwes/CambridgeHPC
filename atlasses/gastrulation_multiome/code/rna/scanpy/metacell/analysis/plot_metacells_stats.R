#####################
## Define settings ##
#####################

# load default setings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
} else {
  stop("Computer not recognised")
}

# Options
opts$number_metacells <- 2500

# I/O
# io$cell2metacell <- paste0(io$basedir,"/results/rna/metacells/cell2metacell_1000metacells.txt.gz")
io$cell2metacell <- sprintf("%s/results/rna_atac/metacells/cell2metacell_%smetacells.txt.gz",io$basedir,opts$number_metacells)
# io$outdir <- paste0(io$basedir,"/results/rna/metacells/pdf"); dir.create(io$outdir, showWarnings = F)


###################
## Load metadata ##
###################

sample_metadata <- fread(io$metadata) %>%
  # .[pass_rnaQC==TRUE & doublet_call==FALSE] %>%
  .[pass_rnaQC==TRUE & pass_atacQC==TRUE & doublet_call==FALSE] %>%
  .[sample%in%opts$samples & celltype.predicted%in%opts$celltypes]

########################
## Load cell2metacell ##
########################

cell2metacell <- fread(io$cell2metacell) %>% 
  merge(sample_metadata[,c("cell","sample","stage","celltype.mapped")] %>% copy %>% setnames("cell","Metacell"))

length(unique(cell2metacell$Metacell))


############################################
## Plot number of metacells per cell type ##
############################################

to.plot <- cell2metacell[,c("Metacell","celltype.mapped")] %>% unique %>%
  .[,.(num_metacells=.N),by=c("celltype.mapped")] %>%
  setorder(-num_metacells) %>% .[,factor:=factor(celltype.mapped,levels=celltype.mapped)]

p1 <- ggbarplot(to.plot, x="celltype.mapped", y="num_metacells", fill="celltype.mapped") +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(x="", y="Number of metacells") +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=rel(0.75)),
    axis.ticks.x = element_blank()
  )

to.plot2 <- to.plot %>% merge(sample_metadata[,.(ncells=.N),by=c("celltype.mapped")])
                             
p2 <- ggscatter(to.plot2, x="ncells", y="num_metacells", fill="celltype.mapped", size=5, shape=21) +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(x="Number of cells per celltype", y="Number of metacells per celltype") +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=rel(0.75)),
    axis.ticks.x = element_blank()
  )

p <- cowplot::plot_grid(p1, p2, ncol = 2, rel_widths = c(1/2,1/2))

pdf(sprintf("%s/ncells_vs_nmetacells.pdf",io$outdir), width=13, height=4)
print(p)
dev.off()