source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

opts$matrix <- "PeakMatrix"
opts$group_variable <- "celltype.mapped"

io$diff.pseudobulk <- file.path(io$basedir,sprintf("results/atac/archR/differential/pseudobulk/%s/%s",opts$group_variable,opts$matrix))
io$outdir <- file.path(io$basedir,sprintf("results/atac/archR/differential/pseudobulk/%s/%s",opts$group_variable,opts$matrix)); dir.create(io$outdir, showWarnings = F)

##################
## Load results ##
##################

# source(here::here("atac/archR/differential/pseudobulk/analysis/load_data.R"))

diff.dt <- opts$celltypes %>% map(function(i) { opts$celltypes %>% map(function(j) {
  file <- file.path(io$diff.pseudobulk,sprintf("%s_vs_%s_pseudobulk.txt.gz",i,j))
  if (file.exists(file)) {
    fread(file, select = c(1,2)) %>% 
      .[,c("celltypeA","celltypeB"):=list(i,j)] %>%
      return
  }
}) %>% rbindlist }) %>% rbindlist


##########################
## Create summary table ##
##########################

tmp <- diff.dt %>%
  .[,.(diff=round(mean(diff),2)), by=c("celltypeB","idx")] %>%
  setnames("celltypeB","celltype") %>%
  setorder(celltype,-diff)

##########
## Save ##
##########

# Save marker score for all combination of genes and cell types
length(unique(tmp$gene))
length(unique(tmp$celltype))
fwrite(tmp, file.path(io$outdir,sprintf("differential_atac_%s_pseudobulk_summary.txt.gz",opts$matrix)), sep="\t", na="NA", quote=F)

