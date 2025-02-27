library(Seurat)
library(Signac)
library(purrr)
library(data.table)
library(future)
library(cowplot)



source(here::here("settings.R"))


io$signac                <- file.path(io$rawdata, "/processed/atac/signac/signac.rds")
io$plots_out             <- "/bi/home/clarks/plots/10X_multiome/signac"

opts$cores               <- 8
opts$mem                 <- 5 # GB
opts$min_frags           <- 2e3
opts$max_frags           <- 2e6
opts$compute_stats       <- FALSE # compute TSS enrichment and nucleosome signal (only needed first time running)




plan("multiprocess", workers = opts$cores)
options(future.globals.maxSize = opts$mem * 1024 ^ 3)
plan()

dir.create(io$plots_out, recursive = TRUE)


signac <- readRDS(io$signac)
signac

signac@meta.data$sample %>% unique()

if (opts$compute_stats){
  print("computing nucleosome signal and TSS enrichment..")
  # compute nucleosome signal score per cell
  signac <- NucleosomeSignal(object = signac)
  
  # compute TSS enrichment score per cell
  signac <- TSSEnrichment(object = signac, fast = FALSE)
  # 
  # saveRDS(signac, io$signac)
}






# signac$high.tss <- ifelse(signac$TSS.enrichment > 2, 'High', 'Low')
# tss <- TSSPlot(signac, group.by = 'high.tss') + NoLegend()
# save_plot(paste0(io$plots_out, "/tss.pdf"), tss)

# signac$nucleosome_group <- ifelse(signac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
# fraghist <- FragmentHistogram(object = signac, group.by = 'nucleosome_group')
# save_plot(paste0(io$plots_out, "/frag_histogram.pdf"), fraghist)



violins <- VlnPlot(
  object = signac,
  features = c("gex_umis_count", "atac_fragments", "TSS.enrichment", "nucleosome_signal"),
  ncol = 2,
  pt.size = 0
)
violins

save_plot(paste0(io$plots_out, "/violins.pdf"), violins, base_height = 10, base_width = 10)

print("filtering cells..")

signac <- subset(
  x = signac,
  subset = atac_fragments <= opts$max_frags &
   # nCount_RNA < 25000 &
    atac_fragments >= opts$min_frags &
    #nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)
signac

print("saving..")
# save 
signac@meta.data$sample %>% unique()
saveRDS(signac, io$signac)

# # plot pass/fails
# qc <- sample_metadata[, .(cell, pass_rnaQC, pass_sigQC = FALSE)]
# qc[cell %in% colnames(signac), pass_sigQC := TRUE]
# qc[, pass_bothQC := as.logical(pass_rnaQC * pass_sigQC)]
# qc <- melt(qc, id.vars = "cell", value.name = "passQC", variable.name = "type")
# 
# p <- ggplot(qc, aes(type, fill = passQC)) + 
#   geom_bar() +
#   theme_cowplot() +
#   labs(x = "", y = "") +
#   #guides(fill = FALSE) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
# p
# 
# save_plot(paste0(io$plots_out, "/passQC.pdf"), p)




