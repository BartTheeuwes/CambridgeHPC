
#####################
## Define settings ##
#####################

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

# I/O
io$virtual_chip.dir <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/virtual_chipseq_lenient")
io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/virtual_chipseq_lenient/cobinding")

# Options
opts$TFs <- list.files(io$virtual_chip.dir, "*_filt.bed.gz") %>% str_replace_all("_filt.bed.gz","")

calculate_cobinding_stats <- function(mtx) {
  pval.mtx <- matrix(NA,nrow = nrow(mtx), ncol=ncol(mtx))
  dimnames(pval.mtx) <- dimnames(mtx)
  for (i in 1:nrow(mtx)) {
    m <- sum(mtx[i,])
    for (j in i:ncol(mtx)) {
      if (i!=j) {
        n <- sum(mtx[j,])
        expected = (m*n)/N
        observed = mtx[i,j]
        pval.mtx[i,j] <- pval.mtx[j,i] <-poisson.test(observed, r = expected, alternative = c("greater"))$p.value
      }
    }
  }
  return(pval.mtx)
}

###################################
## Load virtual ChIP-seq library ##
###################################

virtual_chip.dt <- opts$TFs %>% map(function(i) {
  tmp <- fread(sprintf("%s/%s_filt.bed.gz",io$virtual_chip.dir,i)) 
  if (nrow(tmp)>=10) {
    tmp %>% 
      setnames(c("chr","start","end","correlation_score","max_accessibility_score","motif_score","motif_counts","score")) %>%
    .[,idx:=sprintf("%s:%s-%s",chr,start,end)] %>%
    .[,c("chr","start","end"):=NULL] %>%
    .[,tf:=i] %>%
    return
  }
}) %>% rbindlist

opts$TFs <- unique(virtual_chip.dt$tf)

#############################
## Detect cobinding events ##
#############################

motif2peak.mtx <- virtual_chip.dt[,c("tf","idx")] %>% .[,value:=1] %>% dcast(idx~tf, value.var="value", fill = 0) %>%
  matrix.please

# Subset to peaks with N>=2
motif2peak.mtx.filt <- motif2peak.mtx[rowSums(motif2peak.mtx)>=2,] %>% as.matrix

foo <- crossprod(motif2peak.mtx.filt)
foo <- foo[rowSums(foo)>=5,colSums(foo)>=5]

################
## Statistics ##
################

N <- nrow(motif2peak.mtx.filt)

pval.mtx <- matrix(NA,nrow = nrow(foo), ncol=ncol(foo))
dimnames(pval.mtx) <- dimnames(foo)
for (i in 1:nrow(foo)) {
  m <- sum(foo[i,])
  for (j in 1:ncol(foo)) {
    n <- sum(foo[j,])
    expected = (m*n)/N
    observed = foo[i,j]
    pval.mtx[i,j] <- poisson.test(observed, r = expected, alternative = c("greater"))$p.value
  }
}
diag(pval.mtx) <- NA

###################################
## Plot cobinding results per TF ##
###################################

cobinding.dt <- pval.mtx %>% as.data.table(keep.rownames = T) %>%
  melt(id.vars="rn") %>% setnames(c("TF1","TF2","pvalue")) %>%
  .[!is.na(pvalue)] %>%
  .[,log_padj_fdr := -log(p.adjust(pvalue+1e-50, method="fdr")), by="TF1"]

tfs.to.plot <- cobinding.dt %>% .[TF1!=TF2,max(log_padj_fdr),by="TF1"] %>% .[V1>2,TF1]

# i <- "T"
for (i in tfs.to.plot) {
  
  to.plot <- cobinding.dt[TF1==i] %>% 
    setorder(-log_padj_fdr) %>% head(n=15) %>%
    .[,TF2:=factor(TF2,levels=TF2)]
  
  p <- ggplot(to.plot, aes_string(x="TF2", y="log_padj_fdr"), fill="gray70") +
    geom_point(size=2) +
    geom_segment(aes_string(xend="TF2"), size=0.5, yend=0) +
    coord_flip() +
    labs(x="", y="Cobinding score (log10 p-value)") +
    theme_classic() +
    theme(
      axis.ticks.y = element_blank(),
      axis.text = element_text(size=rel(0.75), color="black")
    )
  
  pdf(sprintf("%s/%s_cobinding.pdf",io$outdir,i), width = 5, height = 4)
  print(p)
  dev.off()
}


#############
## Explore ##
#############

names(which(rowSums(motif2peak.mtx.filt[,c("NOTO","FOXA1")])==2))
