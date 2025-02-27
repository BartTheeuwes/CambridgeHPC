#####################
## Define settings ##
#####################

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

# I/O
# io$virtual_chip.dir <- file.path(io$basedir,"results/rna_atac/virtual_chipseq")
io$outdir <- file.path(io$basedir,"results/rna_atac/virtual_chipseq/stats")

# Options
opts$motif_annotation <- "Motif_cisbp_lenient"
TFs <- list.files(io$virtual_chip.dir, ".bed.gz") %>% str_replace_all(".bed.gz","")
# TFs <- c("T", "ZIC2", "TAL1", "GATA1", "FOXA2", "GATA4", "CDX2")

#######################################
## Load pseudobulk RNA and ATAC data ##
#######################################

source(here::here("rna_atac/load_rna_atac_pseudobulk.R"))

###################################
## Load virtual ChIP-seq library ##
###################################

virtual_chip.dt <- TFs %>% map(function(i) {
  print(i)
  file <- sprintf("%s/%s.bed.gz",io$virtual_chip.dir,i)
  if (file.exists(file)) {
    tmp <- fread(file) 
    if (nrow(tmp)>=10) {
      tmp %>%
        setnames(c("chr","start","end","score")) %>%
        .[,idx:=sprintf("%s:%s-%s",chr,start,end)] %>%
        .[,c("chr","start","end"):=NULL] %>%
        .[,tf:=i] %>%
        return
    }
  }
}) %>% rbindlist


######################################
## Predict TF binding per cell type ##
######################################

rna_tf_pseudobulk.dt[,value:=minmax.normalisation(expr),by="gene"]
atac_chromvar_pseudobulk.dt[,value:=minmax.normalisation(chromvar_zscore),by="gene"]

tmp <- merge(rna_tf_pseudobulk.dt,atac_chromvar_pseudobulk.dt,by=c("celltype","gene")) %>% 
  .[,value:=value.x*value.y] %>%
  setnames("gene","tf")

# i <- "Gut"
celltype_virtual_chip.dt <- opts$celltypes %>% map(function(i) {
  virtual_chip.dt %>% 
    merge(tmp[celltype==i,c("tf","value","celltype")],by="tf") %>%
    setnames(c("score","value"),c("global_score","rna_chromvar_score")) %>%
    .[,score:=round(minmax.normalisation(global_score*rna_chromvar_score),2)] %>%
    .[,c("tf","celltype","idx","score")] %>%
    setorder(-score)
  
  # fwrite(foo[score>=0.25], file.path(io$outdir,sprintf("%s_virtual_chip.txt.gz",i)))
})
names(celltype_virtual_chip.dt) <- opts$celltypes

# fwrite(celltype_virtual_chip.dt, file.path(io$outdir,"celltype_virtual_chip_all.txt.gz"))

celltype_virtual_chip.dt[["Blood_progenitors_2"]] %>% View
##########
## Plot ##
##########

# Score vs number of binding sites per tf

seq.ranges <- seq(0.25,1,by=0.01)

i <- "Blood_progenitors_2"
for (i in opts$celltypes) {
  
  to.plot <- seq.ranges %>% map(function(j) {
    celltype_virtual_chip.dt[[i]][score>=j,log2(.N+1),by="tf"] %>% .[,min_score:=j] %>% return
  }) %>% rbindlist %>% setnames("V1","log2_N")
  
  tfs.to.plot <- to.plot[min_score>=0.40,max(log2_N),by="tf"] %>% .[V1>=6] %>% .$tf
  to.plot2 <- to.plot[tf%in%tfs.to.plot]
  
  p <- ggplot(to.plot2, aes_string(x="min_score", y="log2_N", color="tf")) +
    geom_line(size=1) +
    labs(y="Number of predicted binding sites (log2)", x="Minimum score") +
    # scale_x_continuous(breaks=seq(0,0.75,by=0.10)) +
    # scale_color_brewer(palette="Dark2") +
    theme_classic() +
    theme(
      axis.text = element_text(color="black"),
      legend.position = "top",
      legend.title = element_blank()
    )
  
  pdf(sprintf("%s/%s_virtual_chipseq.pdf",io$outdir,i), width = 6, height = 5)
  print(p)
  dev.off()
}

