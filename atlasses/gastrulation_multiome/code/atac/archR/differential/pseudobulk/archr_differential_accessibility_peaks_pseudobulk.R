#####################
## Define settings ##
#####################

# load default setings
source(here::here("settings.R"))
source(here::here("utils.R"))

# I/O
# io$archR.pseudobulk.deviations.se <- sprintf("%s/results/atac/archR/chromvar/pseudobulk/chromVAR_deviations_summarized_experiment_%s_pseudobulk_all_peaks_archr.rds",io$basedir,opts$motif_annotation)
io$pseudobulk_dir <- file.path(io$basedir,"results_new/atac/archR/pseudobulk/celltype.mapped_mnn")
io$outdir <- file.path(io$basedir,"results_new/atac/archR/differential/pseudobulk"); dir.create(io$outdir, showWarnings=F)

# Options
opts$matrix <- "PeakMatrix"

##################################
## Fetch pseudobulk ATAC Matrix ##
##################################

print(sprintf("Fetching pseudobulk ATAC %s matrix...",opts$matrix))

atac_pseudobulk.se <- readRDS(file.path(io$pseudobulk_dir,sprintf("pseudobulk_%s_summarized_experiment.rds",opts$matrix)))#[,opts$celltypes]

# Define feature names (already done in the pseudobulking script)
# if (grepl("peak",tolower(args$matrix),ignore.case=T)) {
#   rownames(atac_pseudobulk.se) <- rowData(atac_pseudobulk.se) %>% as.data.table %>% .[,idx:=sprintf("%s:%s-%s",seqnames,start,end)] %>% .$id
# } else if (grepl("gene",tolower(args$matrix),ignore.case=T)) {
#   rownames(atac_pseudobulk.se) <- rowData(atac_pseudobulk.se)$name
# }

################################
## Differential accessibility ##
################################

stop("REDO DIFF PSEUDOBULK BECAUSE OF SIGN")

# i <- 1; j <- 2
for (i in 1:length(opts$celltypes)) {
  for (j in i:length(opts$celltypes)) {
    if (i!=j) {
      foo <- assay(atac_pseudobulk.se[,opts$celltypes[[j]]])[,1]
      bar <- assay(atac_pseudobulk.se[,opts$celltypes[[i]]])[,1]
      
      atac_diff.dt <- data.table(
        idx = names(foo), 
        diff = round(bar-foo,2) 
        # groupA = opts$celltypes[[i]], 
        # groupB = opts$celltypes[[j]]
      ) %>% sort.abs("diff") 
      
      # save      
      outfile <- sprintf("%s/%s_vs_%s_%s_pseudobulk.txt.gz", io$outdir,opts$celltypes[[i]],opts$celltypes[[j]],opts$matrix)
      fwrite(atac_diff.dt, outfile, sep="\t")
    }
  }
}
