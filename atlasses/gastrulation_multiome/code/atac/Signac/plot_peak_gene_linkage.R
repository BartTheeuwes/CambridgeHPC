library(Seurat)
library(Signac)
library(purrr)
library(data.table)
library(cowplot)
library(ggrepel)




source(here::here("settings.R"))





io$signac          <- file.path(io$basedir, "/processed/atac/signac/signacRNA.rds")
        

io$plots_out       <- "/bi/home/clarks/plots/10X_multiome/Signac/links"

opts$n_genes_plot  <- 50 # top genes to plot by number of sig correlations
opts$lineages      <- c("Gut", "ExE_mesoderm", "Forebrain_Midbrain_Hindbrain", 
                        "Erythroid1")






dir.create(io$plots_out, recursive = TRUE)








signac <- readRDS(io$signac)
signac

idents.plot <- opts$lineages




links <- as.data.table(Links(signac))
links[order(pvalue)]

nhits <- links[, .N, gene] %>%
  setorder(N) %>%
  .[, rank := .I] 

top <- nhits[rank >= (.N - opts$n_genes_plot)]

hits_plot <- ggplot(nhits, aes(rank, N, label = gene)) +
  geom_point() +
  theme_cowplot() +
  geom_text_repel(data = top, size = 3) +
  labs(x = "Ranked Gene List", y = "Number of Correlated Peaks")
hits_plot
save_plot(paste0(io$plots_out, "/number_hits.pdf"), hits_plot)

setorder(links, pvalue)

volc <- ggplot(links, aes(score, -log10(pvalue))) +
  geom_point(alpha = 0.5) +
  theme_cowplot() +
  geom_text_repel(data = links[1:20], size = 3, aes(label = gene))
volc

cov_plots <- map(top$gene, ~{
  p <- CoveragePlot(
    object = signac,
    region = .x,
    features = .x,
    expression.assay = "SCT",
    idents = idents.plot,
    extend.upstream = 1e4,
    extend.downstream = 1e4
  )
  outfile <- paste0(io$plots_out, "/", .x, ".pdf")
  save_plot(outfile, p, base_height = 8, base_width = 8)
  p
})

cov_plots[[1]]
