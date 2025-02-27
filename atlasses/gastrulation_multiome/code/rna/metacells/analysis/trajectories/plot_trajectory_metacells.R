# here::i_am("atac/archR/processing/save_archr_matrices.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

## I/O
io$basedir <- file.path(io$basedir,"test")
io$metacell_metadata <- file.path(io$basedir,"results/rna/metacells/trajectories/nmp/metacells_metadata.txt.gz")
io$metacell_sce <- file.path(io$basedir,"results/rna/metacells/trajectories/nmp/SingleCellExperiment_metacells.rds")
io$trajectory <- "nmp"
io$trajectory_file <- file.path(io$basedir,"results/rna/metacells/trajectories/nmp/metacell_trajectory.txt.gz")
io$outdir <- file.path(io$basedir,"results/rna/metacells/trajectories/nmp/pdf"); dir.create(io$outdir, showWarnings = F)

###################
## Load metadata ##
###################

metacell_metadata.dt <- fread(io$metacell_metadata)

###############################
## Load SingleCellExperiment ##
###############################

sce <- readRDS(io$metacell_sce)

#####################
## Load trajectory ##
#####################

# trajectory.dt <- fread(io$atlas_trajectory) %>% setnames(c("cell","V1","V2"))
trajectory.dt <- fread(io$trajectory_file) %>% 
  setnames(c("metacell","V1","V2")) %>% 
  merge(metacell_metadata.dt,by="metacell")

metacells <- intersect(colnames(sce), trajectory.dt$metacell)
trajectory.dt <- trajectory.dt[metacell%in%cells]
sce <- sce[,colnames(sce)%in%metacells]

##########
## Plot ##
##########

to.plot <- trajectory.dt

ggplot(to.plot, aes(x=V1, y=V2)) +
  geom_point(aes(fill=celltype), size=2.5, shape=21, stroke=0.25) +
  # viridis::scale_fill_viridis() +
  scale_fill_manual(values=opts$celltype.colors) +
  labs(x="Force-directed layout (Dim 1)", y="Force-directed layout (Dim 2)") +
  theme_classic() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position="none"
  )

# pdf(sprintf("%s/pseudotime_vs_expression.pdf",io$outdir), width=10, height=8)
print(p)
# dev.off()

##############################
## Plot expression of genes ##
##############################
 
genes.to.plot <- c("Cdx2")

rna.dt <- as.matrix(logcounts(sce[genes.to.plot,])) %>% 
  as.data.table(keep.rownames = "gene") %>% 
  melt(id.vars = "gene", variable.name = "metacell", value.name = "expr")

to.plot <- trajectory.dt %>% merge(rna.dt, by="metacell", allow.cartesian=TRUE)

ggplot(to.plot, aes(x=V1, y=V2)) +
  geom_point(aes(fill=expr), size=2.5, shape=21, stroke=0.25) +
  # facet_wrap(~gene) +
  scale_fill_gradient(low = "gray80", high = "purple") +
  labs(x="Force-directed layout (Dim 1)", y="Force-directed layout (Dim 2)") +
  theme_classic() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position="none"
  )

# pdf(sprintf("%s/pseudotime_vs_expression.pdf",io$outdir), width=10, height=8)
print(p)
# dev.off()


##############################
## Plot interaction of genes #
##############################

geneA <- "T"
geneB <- "Sox2"

tmp <- data.table(
  metacell = colnames(sce),
  expr = logcounts(sce[geneA,])[1,]*logcounts(sce[geneB,])[1,] %>% minmax.normalisation
)

to.plot <- trajectory.dt %>% merge(tmp, by="metacell")

ggplot(to.plot, aes(x=V1, y=V2)) +
  geom_point(aes(fill=expr), size=2.5, shape=21, stroke=0.25) +
  scale_fill_gradient(low = "gray80", high = "purple") +
  labs(x="Force-directed layout (Dim 1)", y="Force-directed layout (Dim 2)") +
  theme_classic() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position="none"
  )

# pdf(sprintf("%s/pseudotime_vs_expression.pdf",io$outdir), width=10, height=8)
print(p)
# dev.off()