library(Seurat)
library(Signac)
library(purrr)
library(data.table)
library(future)
library(GenomeInfoDb)
library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10)

library(cowplot)

source(here::here("settings.R"))

io$signac           <- file.path(io$rawdata, "/processed/atac/signac/archr_signac.rds")
io$pseudobulk_file <- file.path(io$rawdata, "processed/atac/archR/pseudobulk/pseudobulk_PeakMatrix_summarized_experiment.rds")
io$plots_out     <- "/bi/home/clarks/plots/10X_multiome/default_path"

dir.create(io$plots_out, recursive = TRUE)

# opts$celltypes <- c(
#   "Epiblast",
#   "Primitive_Streak",
#   "Def._endoderm",
#   "Nascent_mesoderm"
# )
# opts$celltype1 <- "Def._endoderm"
# opts$celltype2 <- "Nascent_mesoderm"

# opts$celltypes <- c(
#   "Surface_ectoderm",
#   "Rostral_neurectoderm",
#   "Neural_crest",
#   "Forebrain_Midbrain_Hindbrain" ,
#   "Spinal_cord",
#   "Epiblast",
#   "Caudal_epiblast"
# )
# opts$celltype1 <- "Spinal_cord"
# opts$celltype2 <- "Forebrain_Midbrain_Hindbrain"



# which celltypes to plot
# opts$celltypes <- c(
#   "Epiblast",
#   "Primitive_Streak",
#   "Def._endoderm",
#   "Gut",
#   "Notochord"
# )
# # which celltypes to perform diffacc
# opts$celltype1 <- "Gut"
# opts$celltype2 <- "Notochord"

#which celltypes to plot
opts$celltypes <- c(
  "Epiblast",
  "Primitive_Streak",
  "Rostral_neurectoderm",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Neural_crest",
  "Surface_ectoderm",
  "Nascent_mesoderm"
)
# which celltypes to perform diffacc
opts$celltype1 <- "Epiblast"
opts$celltype2 <- "Forebrain_Midbrain_Hindbrain"


opts$celltypes <- c(
  "Epiblast",
  "Primitive_Streak",
  "Nascent_mesoderm",
  "Haematoendothelial_progenitors",
  "Endothelium",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Blood_progenitors_3",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3"


)
# which celltypes to perform diffacc
opts$celltype1 <- "Blood_progenitors_1"
opts$celltype2 <- "Endothelium"

meta <- fread(io$metadata)
meta[, unique(celltype.predicted)]
meta[, .N, celltype.predicted][order(-rank(N))]



# load signac

signac <- readRDS(io$signac) %>% 
  subset(subset = celltype.mapped %in% c(opts$celltype1, opts$celltype2))
signac


# differential testing


da_peaks <- FindMarkers(
  object = signac,
  ident.1 = opts$celltype1,
  ident.2 = opts$celltype2,
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = "nCount_peaks"
)


da_peaks <- setDT(da_peaks, keep.rownames = "locus")
da_peaks



# filter to keep hyper/hypo-accessible sites


da_peaks <- da_peaks %>% 
  .[p_val_adj < 0.05 & avg_logFC > 0.5, type := paste(opts$celltype1, "hyper-accessible sites")] %>% 
  .[p_val_adj < 0.05 & avg_logFC < -0.5, type := paste(opts$celltype2, "hyper-accessible sites")]

# load pseudobulk

bulk_se <- readRDS(io$pseudobulk_file) 

opts$celltypes_filt <- colnames(bulk_se)[colnames(bulk_se) %in% opts$celltypes]

bulk <- assays(bulk_se)$PeakMatrix %>% 
  cbind(bulk_se@elementMetadata) %>% 
  as.data.table() %>% 
  .[, locus := paste0(seqnames, "-", start, "-", end)] %>% 
  .[, .SD, .SDcol = c("locus", opts$celltypes_filt)] %>% 
  melt(id.vars = "locus", value.name = "acc", variable.name = "celltype")

dt <- merge(bulk, da_peaks[!is.na(type)], by = "locus") %>% 
  .[, celltype := factor(celltype, levels = opts$celltypes)]



hline <- dt[, mean(acc), type]

p <- ggplot(dt, aes(celltype, acc, fill = celltype)) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  guides(fill = FALSE) +
  geom_hline(data = hline, aes(yintercept = V1), linetype = "dashed") +
  ylab("Pseudobulk ATAC") +
  facet_wrap(~type)

p

outfile <- paste0(io$plots_out, "/", opts$celltype1, "_vs_", opts$celltype2, ".pdf")
if (file.exists(outfile)) outfile <- gsub(".pdf", "_new.pdf", outfile)
save_plot(outfile, p, base_width = 8, base_height = 6)

# up <- da_peaks[p_val_adj < 0.05 & avg_logFC > 0.5, locus] %>%
#   signac[., ] %>%
#   #colMeans(.@assays$peaks@data)
#   #colSums(.@assays$peaks@data)
# 
# down  <- da_peaks[p_val_adj < 0.05 & avg_logFC < -0.5, locus] %>%
#   signac[., ] %>%
#   colMeans(.@assays$peaks@data)
#   #colSums(.@assays$peaks@data)

# 
# 
# # check cells match up
# all(names(up) == names(down))
# 
# # plot
# means <- data.table(cell = names(up), up = up, down = down) %>%
#   merge(meta, by = "cell", all.x = TRUE)
# 
# # re-order celltypes
# means[, celltype := factor(celltype.predicted, levels = opts$celltypes)]
# 
# hline <- means[, .(up = mean(up),
#                    down = mean(down),
#                    med_up = median(up),
#                    med_down = median(down))]
# 
# title <- paste("Differentially accesible sites between", opts$celltype1, "and", opts$celltype2)
# 
# p1 <- ggplot(means, aes(celltype, up, fill = celltype)) +
#   #geom_violin(alpha = 0.25) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.5) +
#   theme_cowplot() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   guides(fill = FALSE) +
#   ylab("mean accessibility") +
#   geom_hline(data = hline, aes(yintercept = med_up), linetype = "dashed") +
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank()) +
#   ggtitle(paste(opts$celltype1, "hyper-accessible sites"))
# 
# p2 <- ggplot(means, aes(celltype, down, fill = celltype)) +
#   #geom_violin(alpha = 0.25) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.5) +
#   theme_cowplot() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   guides(fill = FALSE) +
#   geom_hline(data = hline, aes(yintercept = med_down), linetype = "dashed") +
#   ylab("mean accessibility") +
#   ggtitle(paste(opts$celltype2, "hyper-accessible sites"))
# 
# 
# 
# 
# p <- plot_grid(p1, p2, ncol=1, rel_heights = c(2,3))
# 
# outfile <- paste0(io$plots_out, "/", opts$celltype1, "_vs_", opts$celltype2, ".pdf")
# if (file.exists(outfile)) outfile <- gsub(".pdf", "_new.pdf", outfile)
# save_plot(outfile, p, base_height = 12, base_width = 8)
# p
# 
# 
# 
# 
# # up <- da_peaks[p_val_adj < 0.05 & avg_logFC > 1, locus]
# # down  <- da_peaks[p_val_adj < 0.05 & avg_logFC < -1, locus]
# # 
# # sub <- signac[c(up, down), ]
# # sub
# # 
# # gr <- sub@assays$peaks@ranges
# # loci <- paste0(gr@seqnames, "-", gr@ranges)
# # 
# # mat <- copy(da_peaks)[locus %in% c(up, down), .(locus)] %>% 
# #   .[locus %in% up, type := "up"] %>% 
# #   .[locus %in% down, type := "down"] %>% 
# #   .[, v:= 1] %>% 
# #   dcast(locus ~ type, value.var = "v", fill = 0) %>% 
# #   as.matrix(rownames = "locus") %>% 
# #   .[loci, ]
# # 
# # sub <- sub[loci, ]
# # 
# # 
# # 
# # 
# # 
# # motifob <- CreateMotifObject(mat, positions = gr)
# # 
# # chromvar <- CreateChromatinAssay(
# #   counts = sub@assays$peaks@counts,
# #   motifs = motifob,
# #   ranges = sub@assays$peaks@ranges
# #   #fragments = frags,
# #   # annotation = annotation
# # )
# # 
# # chromvar <- RunChromVAR(chromvar, genome = BSgenome.Mmusculus.UCSC.mm10)
# # 
# # toplot <- as.data.table(chromvar@data, keep.rownames = "type") %>% 
# #   melt(id.vars = "type", variable.name = "cell", value.name = "chromvar") %>% 
# #   merge(meta[, .(cell, celltype.predicted)])
# # 
# # ggplot(toplot, aes(celltype.predicted, chromvar, fill = celltype.predicted)) +
# #   geom_boxplot(alpha = 0.75, outlier.shape = NA) +
# #   theme_cowplot() +
# #   guides(fill = FALSE) +
# #   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
# #   facet_wrap(~type)