library(SnapATAC)
library(purrr)
library(data.table)

library(ggplot2)
library(cowplot)
library(ggrepel)


source(here::here("settings.R"))
source(here::here("public_datasets/Pijuan-Sala_2020/snapatac/snapatac_settings.R"))

snap <- readRDS(snapio$rds_file)
snap

meta <- fread(snapio$metadata) 
meta[, barcode := cell]

umap <- as.data.table(snap@umap)
umap[, cell := snap@barcode]
umap[, cluster := snap@cluster]
umap <- merge(umap, meta, by = "cell")

means <- umap[, .("umap-1" = mean(`umap-1`), "umap-2" = mean(`umap-2`)), celltype]

ggplot(umap, aes(`umap-1`, `umap-2`)) +
  geom_point(aes(colour = celltype), size = 1) +
  geom_label_repel(data = means, aes(label = celltype), label.padding = 0.2) +
  theme_cowplot() 

ggplot(umap, aes(`umap-1`, `umap-2`)) +
  geom_point(aes(colour = cluster), size = 0.5) +
  theme_cowplot()


