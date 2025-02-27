here::i_am("atac/archR/chromvar/differential/pseudobulk/differential/wt_vs_ko/run_archR_differential_chromvar_pseudobulk_wt_vs_ko.R")

# load default setings
source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# Options
opts$motif_annotation <- "CISBP"

# I/O
io$archR.chromvar.pseudobulk.se <- sprintf("%s/results/atac/archR/chromvar/pseudobulk/celltype_genotype/chromVAR_deviations_%s_pseudobulk_archr.rds",io$basedir,opts$motif_annotation)
io$atac.pseudobulk.stats <- file.path(io$basedir,"results/atac/archR/pseudobulk/celltype_genotype/stats.txt")
io$outdir <- file.path(io$basedir,"results/atac/archR/chromvar/differential/pseudobulk/wt_vs_ko"); dir.create(io$outdir, showWarnings = F, recursive = T)

###########################
## Load pseudobulk stats ##
###########################

atac_pseudobulk_stats.dt <- fread(io$atac.pseudobulk.stats)
samples.to.use <- atac_pseudobulk_stats.dt[N>=50,group]

#####################################
## Load pseudobulk chromVAR scores ##
#####################################

atac_chromvar_pseudobulk.se <- readRDS(io$archR.chromvar.pseudobulk.se)[,samples.to.use]
assays(atac_chromvar_pseudobulk.se) <- assays(atac_chromvar_pseudobulk.se)["z"]

######################################
## Differential motif accessibility ##
######################################

tmp <- strsplit(samples.to.use,"-") %>% map_chr(1) %>% table
celltypes.to.use <- tmp[tmp==2] %>% names

diff.dt <- list()
# i <- "NMP"
for (i in celltypes.to.use) {
  
  foo <- assay(atac_chromvar_pseudobulk.se[,paste0(i,"-WT")])[,1]
  bar <- assay(atac_chromvar_pseudobulk.se[,paste0(i,"-T_KO")])[,1]
  
  diff.dt[[i]] <- data.table(
    gene = names(foo), 
    diff = round(foo-bar,2), 
    celltype = i
  ) %>% sort.abs("diff") 
  
  # save      
  fwrite(diff.dt[[i]], file.path(io$outdir,sprintf("%s_WT_vs_T_KO_chromVAR_pseudobulk.txt.gz",i)), sep="\t")
}

diff_concat.dt <- rbindlist(diff.dt)
fwrite(diff_concat.dt, file.path(io$outdir,"WT_vs_T_KO_chromVAR_pseudobulk.txt.gz"), sep="\t")
