suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(Matrix)
    library(scran)
    library(Rtsne)
    library(BiocParallel)
    require(irlba)
    library(batchelor)
    library(stringr)
    library(gridExtra)
    library(M3C)
})

out_folder = "/rds/project/bg200/rds-bg200-hphi-gottgens/users/bt392/mouse/Mixl1_KO/data/"
chim_before_umap = as.data.frame(read.csv(paste0(out_folder, 'umap_before_correction_chim.csv'), header=T, row.names=1,sep=","))
chim_after_umap = as.data.frame(read.csv(paste0(out_folder, 'umap_after_correction_chim.csv'), header=T, row.names=1,sep=","))
pca_before_umap = as.data.frame(read.csv(paste0(out_folder, 'umap_before_correction.csv'), header=T, row.names=1,sep=","))
pca_after_umap = as.data.frame(read.csv(paste0(out_folder, 'umap_after_correction.csv'), header=T, row.names=1,sep=","))


#### print plots to png
# plot umaps

# colors
source("/rds/project/bg200/rds-bg200-hphi-gottgens/users/bt392/mouse/Mixl1_KO/plot_colours.R")

# chimaera
options(repr.plot.width = 10, repr.plot.height = 5)
p1 <- ggplot(chim_before_umap, aes(umap1, umap2, color=sample)) +
    geom_point(alpha=0.4, size = 0.4) +
    theme_classic() +
    ggtitle('chimaera only, before correction, color = samples')
p2 <- ggplot(chim_before_umap, aes(umap1, umap2, color=tdTom)) +
    geom_point(alpha=0.4, size = 0.4) +
    theme_classic() +
    ggtitle('chimaera only, before correction, color =  tomato')

p3 <- ggplot(chim_after_umap, aes(umap1, umap2, color=sample)) +
    geom_point(alpha=0.4, size = 0.4) +
    theme_classic() +
    ggtitle('chimaera only, after correction, color = samples')
p4 <- ggplot(chim_after_umap, aes(umap1, umap2, color=tdTom)) +
    geom_point(alpha=0.4, size = 0.4) +
    theme_classic() +
    ggtitle('chimaera only, after correction, color = tomato')

# whole object
p5 <- ggplot(pca_before_umap, aes(umap1, umap2, color=sample)) +
    geom_point(alpha=0.4, size = 0.4) +
    theme_classic() +
    ggtitle('whole object before correction, color = sample')
p6 <- ggplot(pca_after_umap, aes(umap1, umap2, color=sample)) +
    geom_point(alpha=0.4, size = 0.4) +
    theme_classic() +
    ggtitle('whole object after correction, color = sample')

p7 <- ggplot(pca_before_umap, aes(umap1, umap2, color=tdTom)) +
    geom_point(alpha=0.4, size = 0.4) +
    theme_classic() +
    ggtitle('whole object before correction, color = tomato')
p8 <- ggplot(pca_after_umap, aes(umap1, umap2, color=tdTom)) +
    geom_point(alpha=0.4, size = 0.4) +
    theme_classic() +
    ggtitle('whole object after correction, color = tomato')

p9 <- ggplot(pca_before_umap, aes(umap1, umap2, color=stage)) +
    geom_point(alpha=0.4, size = 0.4) +
    theme_classic() +
    scale_colour_manual(values = stage_colours, labels = stage_labels, name = "stage") +
    ggtitle('whole object before correction, color = timepoint')
p10 <- ggplot(pca_after_umap, aes(umap1, umap2, color=stage)) +
    geom_point(alpha=0.4, size = 0.4) +
    theme_classic() +
    scale_colour_manual(values = stage_colours, labels = stage_labels, name = "stage") +
    ggtitle('whole object after correction, color = timepoint')

atlas_before_umap <- pca_before_umap %>% filter(!is.na(celltype))
p11 <- ggplot(atlas_before_umap, aes(umap1, umap2, color=celltype)) +
    geom_point(alpha=0.4, size = 0.4) +
    theme_classic() +
    scale_colour_manual(values = celltype_colours, labels = celltype_labels, name = "celltype") +
    ggtitle('atlas only before correction, color = celltype')
atlas_after_umap <- pca_after_umap %>% filter(!is.na(celltype))
p12 <- ggplot(atlas_after_umap, aes(umap1, umap2, color=celltype)) +
    geom_point(alpha=0.4, size = 0.4) +
    theme_classic() +
    scale_colour_manual(values = celltype_colours, labels = celltype_labels, name = "celltype") +
    ggtitle('atlas only after correction, color = celltype')

pca_before_85 <- pca_before_umap %>% filter(stage == 'E8.5')
p13 <- ggplot(pca_before_85, aes(umap1, umap2, color=tdTom)) +
    geom_point(alpha=0.4, size = 0.4) +
    theme_classic() +
    ggtitle('E8.5 before correction, color = tomato')

pca_after_85 <- pca_after_umap %>% filter(stage == 'E8.5')
p14 <- ggplot(pca_after_85, aes(umap1, umap2, color=tdTom)) +
    geom_point(alpha=0.4, size = 0.4) +
    theme_classic() +
    ggtitle('E8.5 after correction, color = tomato')
    
# plot to png instead of pdf to stop laptop from burning
plot_folder = "/rds/project/bg200/rds-bg200-hphi-gottgens/users/bt392/mouse/Mixl1_KO/plots/"
png(paste0(plot_folder, "umap_batch_correction_comparisons.png"), height=5000, width=3500)
    grid.arrange(p1, p3, p2, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, nrow =7)
dev.off()

png(paste0(plot_folder, "umap_batch_correction_comparisons_chimaera.png"), height=1000, width=1000)
    grid.arrange(p1, p3, p2, p4, nrow =2)
dev.off()

png(paste0(plot_folder, "umap_batch_correction_comparisons_sample.png"), height=1000, width=2500)
    grid.arrange(p5, p6, nrow =1)
dev.off()

png(paste0(plot_folder, "umap_batch_correction_comparisons_tomato.png"), height=1000, width=2500)
    grid.arrange(p7, p8, nrow =1)
dev.off()

png(paste0(plot_folder, "umap_batch_correction_comparisons_timepoint.png"), height=1000, width=2500)
    grid.arrange(p9, p10, nrow =1)
dev.off()

png(paste0(plot_folder, "umap_batch_correction_comparisons_celltype.png"), height=1000, width=2500)
    grid.arrange(p11, p12, nrow =1)
dev.off()

png(paste0(plot_folder, "umap_batch_correction_comparisons_E85.png"), height=1000, width=2500)
    grid.arrange(p13, p14, nrow =1)
dev.off()

png(paste0(plot_folder, "umap_batch_corrected_plots.png"), height=2000, width=3500)
    grid.arrange(p6, p8, p12, p10, nrow =2)
dev.off()