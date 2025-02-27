
#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/public_datasets/Pijuan-Sala_2020/settings.R")
} else {
  stop("Computer not recognised")
}

# I/O
io$outdir <- paste0(io$basedir,"/results/peaks2genes")

# Options
opts <- list()
opts$gene_window <- 1e6   # window length for the overlap

###############
## Load data ##
###############

# Load gene metadata
gene_metadata <- fread(io$gene_metadata) %>% 
  .[,chr:=as.factor(sub("chr","",chr))] %>%
  setnames(c("ens_id","symbol"),c("id","gene"))

# Load peak metadata
peak_metadata <- fread(io$peak.metadata) %>%
  .[,id:=sprintf("%s:%s-%s",chr,start,end)]

################
## Parse data ##
################

# Prepare metadata for the overlap
gene_metadata <- gene_metadata[, c("chr","start","end","gene","id")] %>%
  # .[,c("start", "end") := list(start-opts$gene_window, end+opts$gene_window)] %>% 
  setkey(chr,start,end)
  
#############
## Overlap ##
#############

ov <- foverlaps(
  peak_metadata %>% setkey(chr,start,end),
  gene_metadata[, c("chr","start","end","gene")],
  nomatch = NA
) %>% 
  setnames(c("i.start","i.end"),c("feature.start","feature.end")) %>%
  setnames(c("start","end"),c("gene.start","gene.end")) %>%
  .[,c("gene.start","gene.end") := list (gene.start+opts$gene_window, gene.end-opts$gene_window)] %>%
  # .[,c("start_dist","end_dist"):=list( abs(gene.end-feature.start), abs(gene.start-feature.end))] %>%
  .[,c("start_dist","end_dist"):=list( gene.end-feature.start, gene.start-feature.end)] %>%
  .[,c("start_dist","end_dist"):=list( ifelse(end_dist<0 & start_dist>0,0,start_dist), ifelse(end_dist<0 & start_dist>0,0,end_dist) )] %>%
  .[,dist:=ifelse(abs(start_dist)<abs(end_dist),abs(start_dist),abs(end_dist))] %>% .[,c("start_dist","end_dist"):=NULL]

# Select nearest gene
ov_nearest <- ov %>%
  .[.[,.I[dist==min(dist)], by=c("id")]$V1] %>%
  .[complete.cases(.)] %>%
  .[!duplicated(id)]

# Sanity check  
ov_nearest$id[(duplicated(ov_nearest$id))]

##########
## Save ##
##########

fwrite(ov, paste0(io$outdir,"/peaks2genes.txt.gz"), sep="\t", na="NA")
fwrite(ov_nearest, paste0(io$outdir,"/peaks2genes_nearest.txt.gz"), sep="\t", na="NA")
