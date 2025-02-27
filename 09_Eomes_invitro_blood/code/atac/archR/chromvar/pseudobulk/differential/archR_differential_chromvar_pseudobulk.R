#####################
## Define settings ##
#####################

# load default setings
source(here::here("settings.R"))
source(here::here("utils.R"))

# Options
opts$motif_annotation <- "Motif_cisbp_lenient"

# I/O
# io$archR.pseudobulk.deviations.se <- sprintf("%s/results/atac/archR/chromvar/pseudobulk/chromVAR_deviations_summarized_experiment_%s_pseudobulk_all_peaks_archr.rds",io$basedir,opts$motif_annotation)
io$outdir <- file.path(io$basedir,"results/atac/archR/chromvar/pseudobulk/differential")

#####################################
## Load pseudobulk chromVAR scores ##
#####################################

# source(here::here("rna_atac/load_rna_atac_pseudobulk.R"))

# Load
atac_pseudobulk_chromvar.se <- readRDS(io$archR.pseudobulk.deviations.se)

# Subset z-scores
if ("z"%in%assayNames(atac_pseudobulk_chromvar.se)) {
  assays(atac_pseudobulk_chromvar.se) <- assays(atac_pseudobulk_chromvar.se)["z"]
} else {
  atac_pseudobulk_chromvar.se <- atac_pseudobulk_chromvar.se[rowData(atac_pseudobulk_chromvar.se)$seqnames=="z",]
}

# Rename motifs
if (any(grepl("_",rownames(atac_pseudobulk_chromvar.se))) | any(grepl("^f",rownames(atac_pseudobulk_chromvar.se)))) {
  rownames(atac_pseudobulk_chromvar.se) <- rowData(atac_pseudobulk_chromvar.se)$name %>% toupper %>% stringr::str_split(.,"_") %>% map_chr(1)
}
rownames(atac_pseudobulk_chromvar.se) <- gsub("TCFAP","TFAP",rownames(atac_pseudobulk_chromvar.se))
rownames(atac_pseudobulk_chromvar.se) <- gsub("NKX2","NKX2-",rownames(atac_pseudobulk_chromvar.se))
rownames(atac_pseudobulk_chromvar.se) <- gsub("NKX3","NKX3-",rownames(atac_pseudobulk_chromvar.se))
rownames(atac_pseudobulk_chromvar.se) <- gsub("NKX6","NKX6-",rownames(atac_pseudobulk_chromvar.se))

# Remove duplicated motifs
atac_pseudobulk_chromvar.se <- atac_pseudobulk_chromvar.se[!duplicated(rownames(atac_pseudobulk_chromvar.se)),]

# Create long data.table
# atac_chromvar_pseudobulk.dt <- assay(atac_pseudobulk_chromvar.se) %>% t %>%
#   as.data.table(keep.rownames = T) %>%
#   setnames("rn","celltype") %>%
#   melt(id.vars=c("celltype"), variable.name="gene", value.name="chromvar_zscore")
  
######################################
## Differential motif accessibility ##
######################################

# i <- 1; j <- 2
for (i in 1:length(opts$celltypes)) {
  for (j in i:length(opts$celltypes)) {
    if (i!=j) {
      foo <- assay(atac_pseudobulk_chromvar.se[,opts$celltypes[[j]]])[,1]
      bar <- assay(atac_pseudobulk_chromvar.se[,opts$celltypes[[i]]])[,1]
      
      chromvar_diff.dt <- data.table(
        gene = names(foo), 
        diff = round(foo-bar,2), 
        groupA = opts$celltypes[[i]], 
        groupB = opts$celltypes[[j]]
      ) %>% sort.abs("diff") 
      # save      
      outfile <- sprintf("%s/%s_vs_%s_chromVAR_pseudobulk.txt.gz", io$outdir,opts$celltypes[[i]],opts$celltypes[[j]])
      fwrite(chromvar_diff.dt, outfile, sep="\t")
    }
  }
}
