
#####################
## Define settings ##
#####################

# Load default settings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
}

io$outdir <- paste0(io$basedir,"/results/rna/soupX")

#######################
## Load soupX output ##
#######################


soup.dt <- opts$samples %>% map(function(i) {
  file <- sprintf("%s/original/%s/filtered_feature_bc_matrix/soupX/soup.tsv.gz",io$basedir,i)
  fread(file) %>% .[,sample:=i]
}) %>% rbindlist

# Save
# fwrite(soup.dt, sprintf("%s/soupX_estimates.txt.gz",io$outdir), sep="\t")


##########
## Plot ##
##########

to.plot <- soup.dt %>% 
  setorder(sample,-est) %>%
  split(.$sample) %>% map(~ head(.,n=35)) %>% rbindlist

gene.order <- to.plot[,mean(est),by="gene"] %>% setorder(-V1) %>% .$gene
to.plot[,gene:=factor(gene,levels=gene.order)]

p <- ggbarplot(to.plot, x="gene", y="est", fill="gray70", ) +
  coord_flip() +
  facet_wrap(~sample, scale="fixed", nrow=1) +
  labs(x="", "Soup relative (?) abundance") +
  theme(
    axis.text = element_text(size=rel(0.74))
  )

pdf(paste0(io$outdir,"/soupX_abundance.pdf"), width=11, height=10)
print(p)
dev.off()