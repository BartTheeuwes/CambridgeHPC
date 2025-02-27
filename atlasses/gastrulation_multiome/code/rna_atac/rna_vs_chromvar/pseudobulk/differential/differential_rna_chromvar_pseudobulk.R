#####################
## Define settings ##
#####################

# Load default settings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/utils.R")
} else {
  stop("Computer not recognised")
}

# Options
opts$motif_annotation <- "Motif_cisbp"

# I/O
io$archR.pseudobulk.deviations.se <- sprintf("%s/results/atac/archR/chromvar/pseudobulk/chromVAR_deviations_summarized_experiment_%s_pseudobulk_correlated_peaks.rds",io$basedir,opts$motif_annotation)
io$outdir <- sprintf("%s/results/rna_atac/rna_vs_chromvar/pseudobulk/differential",io$basedir)

###################################
## Load pseudobulk chromVAR data ##
###################################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/rna_atac/load_rna_atac_pseudobulk.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/rna_atac/load_rna_atac_pseudobulk.R")
} else {
  stop("Computer not recognised")
}

################
## Parse data ##
################

for (i in 1:length(opts$celltypes)) {
  for (j in i:length(opts$celltypes)) {
    if (i!=j) {
      
      chromvar_filt.dt <- chromvar.dt[celltype%in%c(opts$celltypes[[i]],opts$celltypes[[j]])] %>% 
        .[,chromvar_zscore:=round(chromvar_zscore,2)] %>%
        dcast(gene~celltype,value.var="chromvar_zscore") %>%
        setnames(c("gene","groupA","groupB"))
      
      rna_tf_filt.dt <- rna_tf.dt[celltype%in%c(opts$celltypes[[i]],opts$celltypes[[j]])] %>% 
        .[,expr:=round(expr,2)] %>%
        dcast(gene~celltype,value.var="expr") %>%
        setnames(c("gene","groupA","groupB"))
      
      dt <- rbind(
        chromvar_filt.dt[,modality:="chromVAR"],
        rna_tf_filt.dt[,modality:="RNA"]
      )# %>% .[,diff:=round(groupB-groupA,2)] %>% .[,abs_diff:=abs(diff)] %>% setorder(-abs_diff) %>% .[,abs_diff:=NULL]
      
      outfile <- sprintf("%s/rna_chromvar_pseudobulk_%s_vs_%s.txt.gz", io$outdir,opts$celltypes[[i]],opts$celltypes[[j]])
      fwrite(dt, outfile, sep="\t")
    }
  }
}

#######################
## Differential test ##
#######################


##########
## Save ##
##########

