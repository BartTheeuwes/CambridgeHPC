#####################
## Define settings ##
#####################

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))


# I/O
io$rna_pseudobulk.sce <- file.path(io$basedir,"results/rna/pseudobulk/celltype_genotype/SingleCellExperiment_pseudobulk.rds")
io$rna_pseudobulk.stats <- file.path(io$basedir,"results/rna/pseudobulk/celltype_genotype/stats.txt")
io$outdir <- file.path(io$basedir,"results/rna/differential/wt_vs_ko/pseudobulk"); dir.create(io$outdir, showWarnings = F)

# Options
opts$min.cells <- 50

##############################
## Load pseudobulk RNA data ##
##############################

# Load SingleCellExperiment
rna_pseudobulk.sce <- readRDS(io$rna_pseudobulk.sce)

# Load stats
rna_pseudobulk.stats <- fread(io$rna_pseudobulk.stats)

#############################
## Differential expression ##
#############################

tmp <- table(strsplit(rna_pseudobulk.stats[N>=opts$min.cells,group], split = "-") %>% map_chr(1))
celltypes.to.use <- tmp[tmp==2] %>% names

# i <- "NMP"
diff.dt <- celltypes.to.use %>% map(function(i) {
  
    foo <- logcounts(rna.sce[,paste0(i,"-WT")])[,1] %>% round(2)
    bar <- logcounts(rna.sce[,paste0(i,"-T_KO")])[,1] %>% round(2)
    tmp <- data.table(
      gene = names(foo), 
      diff = round(foo-bar,2), 
      celltype = i
    ) %>% sort.abs("diff") 
    
    # save      
    # fwrite(tmp, file.path(io$outdir,sprintf("%s_WT_vs_T_KO_pseudobulk.txt.gz",i)), sep="\t")
    
    return(tmp)
}) %>% rbindlist
  

fwrite(diff.dt, file.path(io$outdir,"WT_vs_T_KO_DE_pseudobulk.txt.gz"), sep="\t")
