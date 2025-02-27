# i <- opts$celltypes[1]; j <- opts$celltypes[2]
atac_diff_pseudobulk.dt <- opts$celltypes %>% map(function(i) { opts$celltypes %>% map(function(j) {
  file <- file.path(io$diff.pseudobulk,sprintf("%s_vs_%s_pseudobulk.txt.gz",i,j))
  if (file.exists(file)) {
    fread(file, select = c(1,2)) %>% 
      .[,c("celltypeA","celltypeB"):=list(i,j)] %>%
      return
  }
}) %>% rbindlist }) %>% rbindlist
