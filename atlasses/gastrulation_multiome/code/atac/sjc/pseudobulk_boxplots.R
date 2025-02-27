
library(purrr)
library(data.table)
library(ggplot2)
library(cowplot)


source(here::here("settings.R"))



io$pseudobulk    <- file.path(io$rawdata, "/processed/atac/signac/pseudobulk/NMPs.tsv.gz")

opts$celltypes <- c(
  "Epiblast",
  "Primitive_Streak",
  "NMP",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord"
  
)



dt <- fread(io$pseudobulk)

toplot <- dt[celltype %in% opts$celltypes] %>% 
  .[,celltype := factor(celltype, levels = opts$celltypes)]



ggplot(toplot, aes(celltype, counts, fill = celltype)) +
  geom_boxplot(alpha = 0.5) +
  theme_cowplot() +
  guides(fill = FALSE) +
  facet_wrap(~anno)











