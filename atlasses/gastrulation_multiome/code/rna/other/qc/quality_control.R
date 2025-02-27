#####################
## Define settings ##
#####################

source("/Users/ricard/gastrulation_multiome_10x/settings.R")

io$outdir <- paste0(io$basedir,"/results/rna/qc")

# opts$min.nFeature_RNA <- 1000
# opts$max.percent.ribo <- 15
# opts$max.percent.mt <- 40

# opts$min.nFeature_RNA <- 1500
# opts$max.percent.ribo <- 15
# opts$max.percent.mt <- 20

opts$min.nFeature_RNA <- 2000
opts$max.percent.ribo <- 15
opts$max.percent.mt <- 20

###############
## Load data ##
###############

seurat <- readRDS(io$seurat)
metadata <- fread(io$metadata)

# batches <- unique(metadata$batch)

#############################
## Calculate QC statistics ##
#############################

seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "mt-")
ribo.genes <- c(grep(pattern = "^Rpl", x = rownames(seurat), value = TRUE),grep(pattern = "^Rps", x = rownames(seurat), value = TRUE))
seurat[["percent.ribo"]] <- PercentageFeatureSet(seurat, features = ribo.genes)

#####################
## Update metadata ##
#####################

selected_cells <- WhichCells(seurat, expression = nFeature_RNA>=opts$min.nFeature_RNA & percent.ribo<=opts$max.percent.ribo & percent.mt<=opts$max.percent.mt)

seurat[["pass_rnaQC"]] <- colnames(seurat) %in% selected_cells
mean(seurat$pass_rnaQC)

metadata.new <- seurat@meta.data %>% 
  tibble::rownames_to_column("barcode") %>%
  as.data.table
stopifnot(metadata.new$barcode==metadata$barcode)
# fwrite(metadata.new, io$metadata, sep="\t", quote=F)

##########
## Plot ##
##########

to.plot <- seurat@meta.data %>% as.data.table(keep.rownames = T) %>%
  setnames("rn","barcode") %>%
  melt(id.vars=c("barcode"), measure.vars=c("nFeature_RNA","percent.ribo","percent.mt"))

vlines <- data.table(
  nFeature_RNA = opts$min.nFeature_RNA,
  percent.ribo = opts$max.percent.ribo,
  percent.mt = opts$max.percent.mt
) %>% melt(variable.name="variable", value.name="value")

p <- ggplot(to.plot, aes(x=value)) +
  facet_wrap(~variable, nrow=1, scales="free") +
  geom_histogram(position = 'identity', bins=250) +
  labs(x="") +
  geom_vline(aes(xintercept=value), data=vlines, linetype="dashed") +
  theme_classic() +
  theme(
    axis.text = element_text(size = rel(0.75), color="black")
  )
  
pdf(paste0(io$outdir,"/qc_rna.pdf"), width=11, height=5, useDingbats = F)
print(p)
dev.off()