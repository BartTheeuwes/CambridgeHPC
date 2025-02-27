###################
## Load settings ##
###################

source("/Users/ricard/gastrulation_multiome_10x/settings.R")

##########################
## Load sample metadata ##
##########################

# opts$samples <- c(
#   "E7.5_rep1",
#   "E7.5_rep2",
#   "E8.5_rep1",
#   "E8.5_rep2"
# )

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE & doublet_call==FALSE] %>%
  .[,log_nFeature_RNA:=log10(nFeature_RNA)] %>%
  .[,log_nFrags_ATAC:=log10(nFrags_atac)]
  # .[pass_rnaQC==TRUE & doublet_call==FALSE] %>%
  # .[celltype.mapped%in%opts$celltypes & sample%in%opts$samples]

################################################
## Load pre-computed dimensionality reduction ##
################################################

# io$umap <- "/Users/ricard/data/gastrulation_multiome_10x/results/atac/archR/dimensionality_reduction/E8.5/E8.5_rep1-E8.5_rep2_umap_nfeatures50000_ndims30_neigh25_dist0.3.txt.gz"
io$umap <- "/Users/ricard/data/gastrulation_multiome_10x/results/atac/archR/dimensionality_reduction/E7.5/E7.5_rep1-E7.5_rep2_umap_nfeatures50000_ndims30_neigh25_dist0.3.txt.gz"

umap.dt <- fread(io$umap)
to.plot <- umap.dt %>% merge(sample_metadata,by="cell")

to.plot %>% .[,foo:=FALSE] %>% .[(umap1>0 & umap1<4) & (umap2>2 & umap2<6.4),foo:=TRUE] # E7.5
# to.plot %>% .[,foo:=FALSE] %>% .[(umap1>3 & umap1<5) & (umap2>0 & umap2<5),foo:=TRUE] # E7.5
# to.plot %>% .[,foo:=FALSE] %>% .[(umap1>-8 & umap1<(-5)) & (umap2>-1 & umap2<4),foo:=TRUE] # E7.5
# to.plot %>% .[,foo:=FALSE] %>% .[(umap1>-7 & umap1<1.5) & (umap2<10.5 & umap2>2.5),foo:=TRUE] # E8.5

ggplot(to.plot, aes(x=umap1, y=umap2)) +
  geom_point(aes(fill=celltype.mapped), size=1.5, shape=21, color="black", stroke=0.05) +
  scale_fill_manual(values=opts$celltype.colors) +
  theme_classic() +
  theme(legend.position = "none")

ggplot(to.plot, aes_string(x="umap1", y="umap2", color="foo")) +
  geom_point(size=0.5) +
  theme_classic()

ggplot(to.plot, aes_string(x="umap1", y="umap2", color="log_nFrags_ATAC")) +
  geom_point(size=0.5) +
  scale_color_gradient(low = "gray80", high = "red") +
  theme_classic()

mean(to.plot$foo,na.rm=T)

############################
## Update sample metadata ##
############################

sample_metadata <- fread(io$metadata) %>%
  # .[cell%in%to.plot[foo==TRUE,cell],pass_atacQC:=FALSE]
  .[cell%in%to.plot[foo==TRUE,cell],doublet_call:=TRUE]

fwrite(sample_metadata, io$metadata, quote=F, na="NA", sep="\t")

#############################
## Update ArchR's metadata ##
#############################

bar <- foo %>% 
  .[archR_cell%in%rownames(ArchRProject)] %>% setkey(archR_cell) %>% .[rownames(ArchRProject)] %>%
  as.data.frame() %>% tibble::column_to_rownames("archR_cell")

stopifnot(all(bar$TSSEnrichment_atac == getCellColData(ArchRProject, "TSSEnrichment")[[1]]))

for (i in colnames(bar)) {
  ArchRProject <- addCellColData(
    ArchRProject,
    data = bar[[i]], 
    name = i,
    cells = rownames(bar),
    force = TRUE
  )
}

colnames(getCellColData(ArchRProject))

saveArchRProject(ArchRProject)
