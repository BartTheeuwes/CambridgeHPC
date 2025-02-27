#####################
## Define settings ##
#####################

# load default setings
source(here::here("settings.R"))
source(here::here("utils.R"))

# Options
opts$motif_annotation <- "Motif_cisbp_lenient"

# I/O
io$chromvar.diff.pseudobulk <- file.path(io$basedir,"results/atac/archR/chromvar/pseudobulk/differential")
io$outdir <- file.path(io$basedir,"results/atac/archR/chromvar/pseudobulk/differential/pdf")

###########################################
## Load precomputed diff chromVAR scores ##
###########################################

diff.dt <- (1:length(opts$celltypes)) %>% map(function(i) {
  (i:length(opts$celltypes)) %>% map(function(j) {
    if (i!=j) fread(sprintf("%s/%s_vs_%s_chromVAR_pseudobulk.txt.gz", io$chromvar.diff.pseudobulk,opts$celltypes[[i]],opts$celltypes[[j]]))
  }) %>% rbindlist
}) %>% rbindlist

fwrite(diff.dt, file.path(io$outdir,"diff_chromVAR_pseudobulk.txt.gz"))

##########
## Plot ##
##########

i <- "Gut"
j <- "Erythroid3"
# celltypes.to.plot <- c("Gut","Erythroid3")
# genes.to.plot <- c("TAL1")

to.plot <- diff.dt[groupA==i & groupB==j] %>% .[,gene:=factor(gene,levels=rev(gene))]

p <- ggplot(to.plot, aes(x=gene, y=diff)) +
  geom_point(aes(color=abs(diff), alpha=abs(diff))) +
  ggrepel::geom_text_repel(data=head(to.plot[diff>0],n=10), aes(x=gene, y=diff, label=gene), size=10) +
  scale_color_gradient(low = "gray80", high = "red") +
  scale_alpha_continuous(range=c(0.25,1)) +
  theme_classic() +
  labs(y="Differential motif accessibility (chromVAR)", x="") +
  annotate("text", x=35, y=-12, size=4, label=sprintf("(+) %s",i)) +
  annotate("text", x=35, y=12, size=4, label=sprintf("(+) %s",j)) +
  coord_flip() +
  geom_segment(x=0, xend=800, y=0, yend=0, color="black", size=0.25, linetype="dashed") +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size=rel(1.0), color="black")
  )
