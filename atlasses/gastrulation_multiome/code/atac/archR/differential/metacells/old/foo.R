################################
## Differential accessibility ##
################################

# if (args$test_mode) {
#   print("Test mode activated, running only a few comparisons...")
#   groups.to.use <- groups.to.use %>% head(n=3)
# }

# i <- 1; j <- 2
# for (i in 1:length(groups.to.use)) {
#   for (j in i:length(groups.to.use)) {
#     if (i!=j) {
#       
#       foo <- assay(atac_metacells.se[,atac_metacells.se$group %in% groups.to.use[[i]]])[,1]
#       bar <- assay(atac_metacells.se[,atac_metacells.se$group %in% groups.to.use[[j]]])[,1]
#       
#       atac_diff.dt <- data.table(
#         idx = names(foo), 
#         diff = round(bar-foo,2) 
#         # groupA = groups.to.use[[i]], 
#         # groupB = groups.to.use[[j]]
#       ) %>% sort.abs("diff") 
#       
#       # save      
#       outfile <- file.path(args$outdir,sprintf("%s_vs_%s_pseudobulk.txt.gz", groups.to.use[[i]],groups.to.use[[j]]))
#       fwrite(atac_diff.dt, outfile, sep="\t")
#     }
#   }
# }




stats.dt <- metacell_metadata.dt[,.N,by="group"] %>% setorder(-N)

groups.to.use <- stats.dt[N>=args$min_metacells,group]
stats.dt[,ignored:=!group%in%groups.to.use] 
cat(sprintf("Groups not considered because having less than %d metacells: %s",args$min_metacells, paste(stats.dt[ignored==T,group],collapse="  ")))

# Save stats
fwrite(stats.dt, file.path(args$outdir,"stats.txt"), quote=F, sep="\t")