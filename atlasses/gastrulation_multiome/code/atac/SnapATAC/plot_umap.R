library(SnapATAC)
library(purrr)
library(data.table)
library(ggplot2)
library(cowplot)
library(ggrepel)


source(here::here("settings.R"))
source(here::here("atac/SnapATAC/snapatac_settings.R"))

opts$colour_by <- c("celltype", 
                    "sample", 
                    "atac_fragments", 
                    "gex_umis_count", 
                    "gex_genes_count",
                    "GATA", "FOXA", "SOX", "POU")



umap <- fread(snapio$metadata)

celltypes <- fread(io$metadata)[, .(cell, celltype = celltype.mapped)]
umap <- merge(umap, celltypes, by = "cell", all.x = TRUE)

means <- umap[, .("umap_1" = mean(umap_1), "umap_2" = mean(umap_2)), celltype]

ggplot(umap, aes(umap_1, umap_2)) +
  geom_point(aes(colour = celltype), size = 1) +
  geom_label_repel(data = means, 
                   aes(`umap_1`, `umap_2`, label = celltype), 
                   label.padding = 0.2,
                   size = 2) +
  guides(colour = FALSE) +
  theme_cowplot() 

colour_by <- map(opts$colour_by, ~colnames(umap)[toupper(colnames(umap)) %like% toupper(.x)]) %>% unlist()
colour_by

plots <- map(colour_by, ~{
  outfile <- paste0(snapio$plots_out, "/umap_", .x, ".pdf")
  p <- ggplot(umap, aes(umap_1, umap_2)) +
    geom_point(aes(colour = get(.x)), size = 0.25) +
    geom_label_repel(data = means, 
                     aes(`umap_1`, `umap_2`, label = celltype), 
                     label.padding = 0.2,
                     size = 1.5) +
    guides(colour = FALSE) +
    #labs(colour = .x) +
    ggtitle(.x) +
    theme_cowplot() +
    theme(plot.title = element_text(hjust = 0.5))
  
  if (class(umap[[.x]]) == "numeric") {
    p <- p + scale_colour_viridis_c()
  }
  
  
    
  save_plot(outfile, p, base_height = 5, base_width = 5)
  p
})
plots[[1]]

