#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
} else {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
}

io$metadata <- paste0(io$basedir,"/results/rna/mapping/sample_metadata_after_mapping.txt.gz")
io$mapping.results <- paste0(io$basedir,"/results/rna/mapping")
io$outdir <- paste0(io$basedir,"/results/rna/mapping/pdf")


opts$query_samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E8.0_rep1",
  "E8.0_rep2",
  "E8.5_rep1",
  "E8.5_rep2"
)

###############
## Load data ##
###############

# Load sample metadata
dt <- fread(io$metadata) %>% .[pass_rnaQC==TRUE] %>% 
  setnames("celltype.mapped","celltype") %>%
  .[,stage:=substr(sample,1,4)]

# Load mapping results
# mapping.dt <- opts$query_samples %>% 
#   map(~ readRDS(sprintf("%s/mapping_mnn_%s.rds",io$mapping.results,.))[["mapping"]] ) %>%
#   rbindlist %>% merge(sample_metadata,by="cell")

##########
## Plot ##
##########

celltype.colors <- opts$celltype.colors[names(opts$celltype.colors) %in% dt$celltype]
stopifnot(sort(unique(as.character(dt$celltype))) == sort(names(celltype.colors)))

to.plot <- dt %>%
  .[!is.na(celltype) & celltype!="",.N, by=c("stage","celltype","sample")]# %>%
  # .[!celltype%in%c("ExE_ectoderm","ExE_endoderm","Parietal_endoderm","Visceral_endoderm")] %>%
  # .[, celltype:=stringr::str_replace_all( celltype,"_"," ")] %>%
  # .[, celltype:=factor( celltype,levels=names(opts$celltype.colors))] 
  # .[, celltype:=factor(celltype,levels=sort(names(opts$celltype.colors), decreasing = F))]

for (i in unique(to.plot$stage)) {
  p <- ggplot(to.plot[stage==i], aes(x=celltype, y=N)) +
    geom_bar(aes(fill=celltype), stat="identity", color="black") +
    scale_fill_manual(values=opts$celltype.colors, drop=F) + 
    scale_x_discrete(drop=FALSE) +
    facet_wrap(~sample, scales="free_x") +
    coord_flip() +
    labs(y="Number of cells") +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(color="black", size=rel(1.2)),
      axis.title.x = element_text(color="black", size=rel(1.1)),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=rel(1.3), color="black"),
      axis.text.x = element_text(size=rel(1.1), color="black")
    )
  
  # pdf(sprintf("%s/celltype_numbers_%s.pdf",io$outdir,i), width=14, height=14)
  pdf(sprintf("%s/celltype_numbers_%s.pdf",io$outdir,i))
  print(p)
  dev.off()
}
