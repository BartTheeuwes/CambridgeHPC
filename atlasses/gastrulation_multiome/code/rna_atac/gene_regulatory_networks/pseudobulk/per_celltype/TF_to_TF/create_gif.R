library(magick)

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
} else {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
}

# I/O
io$input.dir <- paste0(io$basedir,"/results/rna_atac/gene_regulatory_networks/pseudobulk/per_celltype")
io$output.dir <- paste0(io$basedir,"/results/rna_atac/gene_regulatory_networks/pseudobulk/per_celltype/gif"); dir.create(io$output.dir, showWarnings = F)

# Options

# celltypes to plot
celltypes.to.plot <- c("Gut")


#########
## GIF ##
#########

for (i in celltypes.to.plot) {

  files <- list.files(path=sprintf("%s/%s",io$input.dir,i), pattern="*.png", full.names = T)
  
  files %>%
    map(image_read) %>%
    image_join %>%
    image_animate(fps=1) %>%
    image_write(quality=100, path = sprintf("%s/%s_gene_regulatory_network.gif",io$output.dir,i))
}



