library(Seurat)
library(Signac)
library(purrr)
library(data.table)




source(here::here("settings.R"))



io$signacRNA           <- file.path(io$basedir, "/processed/atac/signac/signacRNA.rds")
io$signacMotifs        <- file.path(io$basedir, "/processed/atac/signac/signacMotifs_human.rds")

io$plots_out     <- "/bi/home/clarks/plots/10X_multiome/signac/motifs/human"

opts$tfs         <- "All" #c("Sox17")


dir.create(io$plots_out, recursive = TRUE)


signac <- readRDS(io$signacMotifs)
signac

rna <- readRDS(io$signacRNA)
signac[["RNA"]] <- rna[["RNA"]]
signac
rm(rna)


DefaultAssay(signac) <- "peaks"

lineages <- DimPlot(signac, label = TRUE, pt.size = 0.1, repel = TRUE, label.size = 2) + NoLegend()
lineages

DefaultAssay(signac) <- "chromvar"

motif_names <- copy(signac[["chromvar"]][["tf"]]) %>%
  setDT(keep.rownames = "motif")

if (toupper(opts$tfs) == "ALL") {
  opts$tfs <- motif_names$tf
}



walk(opts$tf, ~{
  
  
  .x
  
  motif <- motif_names[tf == .x, motif]
  
  DefaultAssay(signac) <- "RNA"
  
  gene <- paste0(
    toupper(substr(.x, 1, 1)),
    tolower(substr(.x, 2, nchar(.x)))
  )
  
  # check gene is present in RNA matrix
  if (!gene %in% rownames(signac)) return(NULL)
  
  gene_plot <- FeaturePlot(
    object = signac,
    features = gene,
    #min.cutoff = 'q10',
    max.cutoff = 'q50',
    pt.size = 0.1
  ) + ggtitle(paste(.x, "expression"))
  gene_plot
  
  DefaultAssay(signac) <- "chromvar"
  
  motif_plot <- FeaturePlot(
    object = signac,
    features = motif,
    min.cutoff = 'q10',
    max.cutoff = 'q90',
    pt.size = 0.1
  ) + ggtitle(paste(.x, "motif enrichment"))
  
  motif_plot
  
  p <- lineages + motif_plot + gene_plot
  
  outfile <- paste0(io$plots_out, "/", .x, ".pdf")
  cowplot::save_plot(outfile, p, base_height = 5, base_width = 15)
})













