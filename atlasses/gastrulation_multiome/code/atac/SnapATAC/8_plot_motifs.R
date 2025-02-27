library(SnapATAC)
library(purrr)
library(data.table)
library(ggplot2)
library(cowplot)



source(here::here("settings.R"))
source(here::here("atac/SnapATAC/snapatac_settings.R"))










snap <- readRDS(snapio$rds_file)
snap

motifs_all <- colnames(snap@mmat)


length(motifs_all)

plots <- map(seq_along(motifs_all), ~{
  dt <- data.table(x = snap@metaData$cell_type, y = snap@mmat[, .x])
  
  motif_name <- motifs_all[.x]
  print(motif_name)
  outfile <- paste0(snapio$plots_out, "/", motif_name, ".pdf")
  
  
  # violins 
  p <- ggplot(dt, aes(x, y, colour = x, fill = x)) +
    geom_violin(alpha = 0.5) +
    labs(x = "Cell Type", y = "Motif enrichment") +
    ggtitle(motif_name)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5)) +
    guides(fill = FALSE, colour = FALSE)
  
  save_plot(outfile, p)
  
  # umaps
  umap <- plotFeatureSingle(obj=snap,
                            feature.value=snap@mmat[, .x],
                            method="umap", 
                            main=motif_name,
                            point.size=0.2, 
                            point.shape=19, 
                            down.sample=10000,
                            pdf.file.name = gsub(".pdf", "_umap.pdf", outfile),
                            pdf.width = 7,
                            pdf.height = 7,
                            quantiles=c(0.01, 0.99)) # remove outliers)
    
    
  
  list(p, umap)
})








