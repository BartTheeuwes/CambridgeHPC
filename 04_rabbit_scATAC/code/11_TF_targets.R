here::i_am("01_create_arrow.R")
source(here::here("settings.R"))
source(here::here("load_archr.R"))
library(viridis)
library(igraph)

addArchRThreads(12)
options(repr.plot.width=8, repr.plot.height=8)

matrix = 'GeneScoreMatrix'#'GeneScoreMatrix'
maxDist = 50000

matches = getMatches(ArchRProject)

colnames(matches) = sapply(strsplit(colnames(matches), '_'), '[', 1) # rename columns to actual TF names
colnames(matches@assays@data$matches) = sapply(strsplit(colnames(matches), '_'), '[', 1)

# get gene names
genes = getFeatures(ArchRProj = ArchRProject,
                    useMatrix = matrix)
# All TFs
TFs = sapply(strsplit(colnames(matches), '_'), '[', 1)

# TFs in our data
TFs = intersect(sapply(strsplit(genes, ':'), '[', 2), TFs)

peaks = getPeakSet(ArchRProject)
peaks$idx = 1:length(peaks) # create unique peak IDs

addArchRThreads(1)

genescores = getMatrixFromProject(
  ArchRProj = ArchRProject,
  useMatrix = matrix)
names = rowData(genescores)$name
genescores = assays(genescores)[[matrix]]
rownames(genescores) = names


addArchRThreads(1)

peakacc = getMatrixFromProject(
  ArchRProj = ArchRProject,
  useMatrix = "PeakMatrix")
peakacc = assays(peakacc)$PeakMatrix
rownames(peakacc) = 1:nrow(peakacc)

## Add peak to gene links
# ArchRProject = addPeak2GeneLinks(
#     ArchRProj = ArchRProject,
#     reducedDims = "PeakMatrix",
#     useMatrix = matrix,
#     dimsToUse = 1:30,
#     maxDist = maxDist) # The max distance should probably be greatly reduced

# Retrieve peak to gene links
# peak_gene_links = getPeak2GeneLinks(ArchRProject, returnLoops = FALSE,
#                                 corCutOff = 0.45,
#                                 FDRCutOff = 1e-04, # test with more stringent values
#                                 varCutOffATAC = 0.25,
#                                 varCutOffRNA = 0.25,)

peak_gene_links = readRDS(file.path(io$output.directory, 'Loops/peak_gene_links.rds')) # on GeneScoreMatrix w 50k maxDist

#peak_gene_links = peak_gene_links[peak_gene_links$Correlation > p2g_thr,]
peak_gene_corr = as.data.table(peak_gene_links)[, 1:3] %>% setnames('Correlation', 'peak_gene_cor')

linked_genes = data.table(gene = sapply(strsplit(genes, ':'), '[', 2)[peak_gene_links$idxRNA],
                      idxRNA = peak_gene_links$idxRNA,
                      idxATAC = peak_gene_links$idxATAC)


linked_peaks = as.data.table(peaks[unique(peak_gene_links$idxATAC),])

meta = as.data.table(ArchRProject@cellColData, keep.rownames=TRUE)
pseudobin = data.table('rn'= meta$rn, 'pseudobin' = meta$stage_clusters)

GetTFtargets = function(TF = 'TAL1'){
    if(TF %in% TFs){
        # get peaks that have TF motif
        TF_peak_ids = matches[matches@assays@data$matches[, TF]]@rowRanges$idx

        TF_peaks = linked_peaks[idx %in% TF_peak_ids,]

        # filter peaks that are correlate to expression of TF

        # get imputed TF expression
        TF_expr = genescores[TF,]
        

        # get mean TF expression per pseudotime bin
        pseudobin_expr = merge(pseudobin, as.data.table(TF_expr, keep.rownames=TRUE), 
                        by='rn')

        TF_expr_bin = pseudobin_expr[!is.na(pseudobin),] %>% 
                        .[, sum(TF_expr)/.N, by=pseudobin] %>%
                        setnames('V1', 'TF_expr_bin') %>% 
                        .[order(pseudobin)]

        if(var(TF_expr_bin$TF_expr_bin)>0){ # skip TFs that have no variability in expression
            message(TF)
   
            # get pseudobin accessibility
            peak_subset = as.data.table(t(as.matrix(peakacc[TF_peaks$idx,])), keep.rownames=TRUE)
            pseudobin_peaks = merge(pseudobin, as.data.table(peak_subset, keep.rownames=TRUE), 
                            by='rn')

            peaks_bin = pseudobin_peaks[!is.na(pseudobin),] %>%
                            melt(id.vars = c("rn", "pseudobin")) %>% # wide to long for next calculation
                            .[, sum(value)/.N, by=.(pseudobin, variable)] %>% # calculate mean of accessibility per pseudobin and peak
                            dcast(pseudobin ~ variable, value.var = "V1") %>%  # long to wide
                            .[order(pseudobin)] # order by pseudobin

            # Get correlated peaks, very naive method
            TF_peak_cor = apply(peaks_bin[,2:length(peaks_bin)], 2, function(x){cor(x, TF_expr_bin[,2])})

            TF_peak_cor = as.data.table(TF_peak_cor, keep.rownames=TRUE) %>% 
                            setnames('rn', 'idxATAC') %>% 
                            .[, idxATAC:=as.integer(idxATAC)]

            TF_genes = merge(TF_peak_cor,linked_genes, by = 'idxATAC') %>% 
                            .[order(match(idxATAC,TF_peaks$idx)),] %>%
                            .[,TF:= TF]

            return(TF_genes)
        }
    }
}

addArchRThreads(30)

a = Sys.time()
TF_targets = lapply(TFs, function(x){
    GetTFtargets(x)})
b = Sys.time()
b-a

write.csv(rbindlist(TF_targets), sprintf('/rds/project/rds-SDzz0CATGms/users/bt392/04_Rabbit_ATAC_final/rerun/ArchR/TFNetwork/TF_targets.csv'), row.names=FALSE)