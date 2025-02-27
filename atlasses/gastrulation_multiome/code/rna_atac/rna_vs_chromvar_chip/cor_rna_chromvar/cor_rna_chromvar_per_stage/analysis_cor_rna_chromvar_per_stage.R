
#####################
## Define settings ##
#####################

# Load default settings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
} else {
  stop("Computer not recognised")
}

# I/O
io$input.dir <- paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/all_cells_per_stage")
io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/all_cells_per_stage")

# Options
opts$stages <- c(
  "E7.5",
  "E8.0",
  "E8.5"
)

###############
## Load data ##
###############

# Load correlation results
cor_dt <- opts$stages %>% 
  map(function(x) fread(sprintf(sprintf("%s/%s/cor_rna_vs_chromvar.txt.gz",io$input.dir,x))) %>%
        .[,stage:=x]) %>%
  rbindlist %>% setorder(-stage,padj_fdr)

cor_dt[,stage:=factor(stage,levels=opts$stages)]

###########################################
## Plot number of associations per stage ##
###########################################

to.plot <- cor_dt %>%
  .[,.(N=sum(sig & r>0)), by=c("stage")]

p <- ggbarplot(to.plot, x="stage", y="N", fill="gray70") +
  labs(x="", y="Number of significant (+) associations") +
  theme_classic() +
  theme(
    axis.text.x = element_text(color="black")
  )

pdf(sprintf("%s/number_positive_associations.pdf",io$outdir), width = 5, height=4)
print(p)
dev.off()


###############################################
## Plot genes that change across time points ##
###############################################

# to.plot <- cor_dt %>%
#   .[,N:=sum(sig),by="gene"] %>% .[N<3] %>%
#   dcast(gene~stage, value.var="sig") %>%
#   .[E8.5==TRUE & E7.5==FALSE]

