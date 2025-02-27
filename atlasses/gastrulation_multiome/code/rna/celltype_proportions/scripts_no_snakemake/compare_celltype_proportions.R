linMap <- function(x, from, to) return( (x - min(x)) / max(x - min(x)) * (to - from) + from )

source("/Users/ricard/gastrulation_multiome_10x/settings.R")

################
## Define I/O ##
################

io$metadata <- paste0(io$basedir,"/results/rna/mapping/sample_metadata_after_mapping.txt.gz")
io$outdir <- paste0(io$basedir,"/results/rna/celltype_proportions")

####################
## Define options ##
####################

opts$query_samples <- c(
  "E7.5_rep1",
  "E7.5_rep2"
  # "E8.5_rep1",
  # "E8.5_rep2"
)

# opts$atlas.stages <- c("E8.25","E8.5")
opts$atlas.stages <- c("E7.5")

opts$to.merge <- c(
  "Erythroid3" = "Erythroid",
  "Erythroid2" = "Erythroid",
  "Erythroid1" = "Erythroid",
  "Blood_progenitors_1" = "Blood_progenitors",
  "Blood_progenitors_2" = "Blood_progenitors",
  "Intermediate_mesoderm" = "Mixed_mesoderm",
  "Caudal_Mesoderm" = "Mixed_mesoderm",
  "Nascent_mesoderm" = "Mixed_mesoderm",
  # "Somitic_mesoderm" = "Mixed_mesoderm",
  # "Paraxial_mesoderm" = "Mixed_mesoderm",
  # "Pharyngeal_mesoderm" = "Mixed_mesoderm",
  "Visceral_endoderm" = "ExE_endoderm"
)

# opts$remove.celltypes <- c("PGC","Caudal_epiblast")

opts$min.ncells <- 25

############################
## Update sample metadata ##
############################

sample_metadata <- fread(io$metadata) %>%
  .[pass_QC==TRUE & sample%in%opts$query_samples & !is.na(celltype.mapped)] %>%
  # .[pass_rnaQC==TRUE & hybrid_call==FALSE & sample%in%opts$samples & !is.na(celltype.mapped)] %>%
  setnames("celltype.mapped","celltype") %>%
  .[,celltype:=stringr::str_replace_all(celltype,opts$to.merge)]

####################################
# Calculate cell type proportions ##
####################################

# Atlas
atlas.proportions <- fread(io$atlas.celltype_proportions) %>%
  .[!celltype%in%opts$remove.celltypes] %>%
  .[,celltype:=stringr::str_replace_all(celltype,opts$to.merge)] %>%
  .[,.(celltype_proportion=mean(celltype_proportion), N=mean(N)), by=c("stage","celltype")] %>% # average across samples
  .[stage%in%opts$atlas.stages] %>%
  .[,foo:=mean(N),by=c("celltype","stage")] %>% .[foo>=opts$min.ncells] %>% .[,foo:=NULL] %>%
  .[,ncells:=sum(N), by="stage"] %>%
  .[,.(proportion=sum(N)/unique(ncells), N=mean(N)),by=c("celltype","stage")]


# Query
query.proportions <- sample_metadata %>%
  .[!celltype%in%opts$remove.celltypes] %>%
  .[,foo:=.N,by=c("celltype")] %>% .[foo>=opts$min.ncells] %>% .[,foo:=NULL] %>%
  .[,ncells:=.N, by="sample"] %>%
  .[,.(proportion=.N/unique(ncells)),by=c("celltype","sample")]

# Sanity check
# unique(query.proportions$celltype)[!unique(query.proportions$celltype) %in% unique(atlas.proportions$celltype)]
# stopifnot(unique(query.proportions$celltype) %in% unique(atlas.proportions$celltype))

# Merge
dt <- merge(
  atlas.proportions, 
  query.proportions, 
 by = c("celltype"), 
 allow.cartesian = T, 
 suffixes = c(".atlas",".query")
)

###############################################
## Plot differences in cell type proportions ##
###############################################

to.plot <- dt %>%
  # .[N.query+N.atlas>25] %>% # only consider cell types with enough observations
  # .[,.(diff_proportion=log2(proportion.query/proportion.atlas), diff_N=N.query-N.atlas), by=c("sample","celltype.mapped")] %>% 
  .[,.(diff_proportion=log2(proportion.query/proportion.atlas)), by=c("sample","celltype","stage")]# %>% 
  # .[,diff_N_norm:=linMap(abs(diff_N), from=0.15, to=2.5)]

to.plot[,stage:=paste(stage,"(atlas)")]

to.plot.atlas_line <- data.table(
  celltype = unique(dt$celltype),
  diff_proportion = log2(1)
  # diff_N = 0 #+  0.01
)

p <- ggplot(to.plot, aes(x=factor(celltype), y=diff_proportion, group=1)) +
  geom_point(aes(color = celltype), size=4, stat = 'identity') + 
  geom_polygon(color="black", fill=NA, alpha=0.5, linetype="dashed", data=to.plot.atlas_line) +
  facet_wrap(~sample, nrow=2) +
  # facet_wrap(~sample+stage, nrow=2) +
  scale_color_manual(values=opts$celltype.colors, drop=F) +
  coord_polar() +
  guides(colour = guide_legend(override.aes = list(size=4), ncol=1)) +
  scale_size(guide = 'none') +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.text = element_text(size=rel(0.75)),
    legend.title = element_blank(),
    axis.title=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.line=element_blank(),
    axis.text.x = element_blank()
    # axis.text.x = element_text(angle= -76 - 360 / length(unique(to.plot.test$celltype)) * seq_along(to.plot.test$celltype))
  )

pdf(paste0(io$outdir,"/celltype_proportions_vs_atlas.pdf"), width=7, height=7)
print(p)
dev.off()

