# A set of core functions for the analysis of chimera data sets. 
# Author: Magdalena Strauss

library(scran)
library(irlba)
library(scater)
library(SingleCellExperiment)
library(Matrix)
library(scales)
library(BiocParallel)
library(BiocNeighbors)
library(batchelor)



getHVGs <- function(sce,min_mean = 0.001,FDR = 0.01)
{
  exclude <- c("ENSMUSG00000086503",read.table("data/ygenes.tab", stringsAsFactors = FALSE)[,1],"tomato-td")
  stats <- modelGeneVar(sce,subset.row = setdiff(rownames(sce),exclude),
                        block=sce$sample)
  stats <- stats[stats$mean > min_mean,]
  top_genes <- rownames(stats[stats$FDR < FDR,])
  return(top_genes)
} 

load_atlas_data = function(normalise = TRUE, remove_doublets = FALSE, remove_stripped = FALSE, load_corrected = FALSE){
  
  # This function was taken from 
  # https://github.com/MarioniLab/EmbryoTimecourse2018/blob/master/analysis_scripts/atlas/core_functions.R
  # with minor adaptations
  
  if(load_corrected & (!remove_doublets | !remove_stripped)){
    message("Using corrected PCs, also removing doublets + stripped now.")
    remove_doublets = TRUE
    remove_stripped = TRUE
  }
  
  
  counts = readRDS("/nfs/research/marioni/magda/chimera/fromJonny/atlas_data/raw_counts.rds")
  genes = read.table("/nfs/research/marioni/magda/chimera/fromJonny/atlas_data/genes.tsv", stringsAsFactors = F)
  meta = read.table("/nfs/research/marioni/magda/chimera/data/atlas_meta.tab", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "$")
  
  rownames(counts) = genes[,1] #ensembl
  colnames(counts) = meta$cell
  
  sce = SingleCellExperiment(assays = list("counts" = counts))
  
  if(normalise){
    sfs = read.table("/nfs/research/marioni/magda/chimera/fromJonny/atlas_data/sizefactors.tab", stringsAsFactors = F)[,1]
    sizeFactors(sce) = sfs
    sce = logNormCounts(sce)
  }
  
  if(remove_doublets){
    sce = logNormCounts(sce[,!meta$doublet])
    meta = meta[!meta$doublet,]
  }
  
  if(remove_stripped){
    sce = logNormCounts(sce[,!meta$stripped])
    meta = meta[!meta$stripped, ]
  }
  
  if(load_corrected){
    corrected = readRDS("/nfs/research/marioni/magda/chimera/fromJonny/atlas_data/corrected_pcas.rds")
    assign("corrected", corrected, envir = .GlobalEnv)
    
  }
  
  
  
  
  assign("genes", genes, envir = .GlobalEnv)
  assign("meta", meta, envir = .GlobalEnv)
  assign("sce", sce, envir = .GlobalEnv)
  
  
  invisible(0)
}


getmode <- function(v, dist) {
  # taken from
  # https://github.com/MarioniLab/EmbryoTimecourse2018/blob/master/analysis_scripts/atlas/core_functions.R
  tab = table(v)
  #if tie, break to shortest distance
  if(sum(tab == max(tab)) > 1){
    tied = names(tab)[tab == max(tab)]
    sub = dist[v %in% tied]
    names(sub) = v[v %in% tied]
    return(names(sub)[which.min(sub)])
  } else {
    return(names(tab)[which.max(tab)])
  }
}


mnnMap = function(atlas_pca, atlas_meta, map_pca, meta_chimera, k_map = 10, return.pca = FALSE){
  # taken from (with minor adaptations)
  # https://github.com/MarioniLab/EmbryoTimecourse2018/blob/master/analysis_scripts/atlas/core_functions.R
  
  correct = reducedMNN(atlas_pca, map_pca)$corrected
  atlas = 1:nrow(atlas_pca)
  correct_atlas = correct[atlas,]
  correct_map = correct[-atlas,]
  
  knns = BiocNeighbors::queryKNN(correct_atlas, correct_map, k = k_map, get.index = TRUE, get.distance = FALSE)
  
  #get closest k matching cells
  k.mapped = t(apply(knns$index, 1, function(x) atlas_meta$cell[x]))
  celltypes = t(apply(k.mapped, 1, function(x) atlas_meta$celltype[match(x, atlas_meta$cell)]))
  #celltypes.clustering = t(apply(k.mapped, 1, function(x) atlas_meta$celltype.clustering[match(x, atlas_meta$cell)]))
  stages = t(apply(k.mapped, 1, function(x) atlas_meta$stage[match(x, atlas_meta$cell)]))
  celltype.mapped <- rep("",nrow(celltypes))
  for (j in 1:nrow(celltypes)){
    celltype.mapped[j] <- getmode(celltypes[j,],1:ncol(celltypes))
  }
 # celltype.mapped.clustering <- rep("",nrow(celltypes.clustering))
  # for (j in 1:nrow(celltypes.clustering)){
  #   celltype.mapped.clustering[j] <- getmode(celltypes.clustering[j,],1:ncol(celltypes.clustering))
  # }
  stage.mapped <- rep("",nrow(stages))
  for (j in 1:nrow(stages)){
    stage.mapped[j] <- getmode(stages[j,],1:ncol(stages))
  }

  out = lapply(1:length(celltype.mapped), function(x){
    list(cells.mapped = k.mapped[x,],
         celltype.mapped = celltype.mapped[x],
         stage.mapped = stage.mapped[x])
  })
  
  names(out) = meta_chimera$cell
  
  if(return.pca){
    return(list(out,list("atlas" = correct_atlas,
                         "mapped" = correct_map)))
  }
  
  
  return(out)

}



doBatchCorrect2 = function(counts, 
                           timepoints, 
                           samples, 
                           timepoint_order = c("E6.5", "E6.75", "E7.0", "mixed_gastrulation", "E7.25", "E7.5", "E7.75", "E8.0", "E8.25", "E8.5"), 
                           sample_order=NULL, 
                           npc = 50,
                           pc_override = NULL, 
                           BPPARAM = SerialParam()){
  #taken from 
  # https://github.com/MarioniLab/TChimeras2020/blob/master/core_functions/embryos_core_functions.R
  #default: samples largest to smallest
  if(is.null(sample_order)){
    tab = table(samples)
    sample_order = names(tab)[order(tab), decreasing = TRUE]
  }
  
  #remove timepoints that are not present
  timepoint_order = timepoint_order[timepoint_order %in% timepoints]
  
  if(!is.null(pc_override)){
    pca = pc_override
  } else {
    pca = prcomp_irlba(t(counts), n = npc)$x
    rownames(pca) = colnames(counts)
  }
  
  if(length(unique(samples)) == 1){
    return(pca)
  }
  
  #create nested list
  pc_list = lapply(unique(timepoints), function(tp){
    sub_pc = pca[timepoints == tp, , drop = FALSE]
    sub_samp = samples[timepoints == tp]
    list = lapply(unique(sub_samp), function(samp){
      sub_pc[sub_samp == samp, , drop = FALSE]
    })
    names(list) = unique(sub_samp)
    return(list)
  })
  
  names(pc_list) = unique(timepoints)
  
  #arrange to match timepoint order
  pc_list = pc_list[order(match(names(pc_list), timepoint_order))]
  pc_list = lapply(pc_list, function(x){
    x[order(match(names(x), sample_order))]
  })
  
  #perform corrections within list elements (i.e. within stages)
  correct_list = lapply(pc_list, function(x){
    if(length(x) > 1){
      return(reducedMNN(x, BPPARAM = BPPARAM)$corrected)
    } else {
      return(x[[1]])
    }
  })
  #perform correction over list
  if(length(correct_list)>1){
    correct = reducedMNN(correct_list, BPPARAM = BPPARAM)$corrected
  } else {
    correct = correct_list[[1]]
  }
  
  correct = correct[match(colnames(counts), rownames(correct)),]
  
  return(correct)

}



mapWrap2 <-  function(atlas_sce, atlas_meta, sce_chimera, meta_chimera, return.pca = FALSE, nPC = 50,
                      target_name,latent_time=NULL){
  # taken from the mapWrap function 
  # https://github.com/MarioniLab/TChimeras2020/blob/master/core_functions/embryos_core_functions.R
  
  #prevent duplicate rownames
  colnames(sce_chimera) = paste0("map_", colnames(sce_chimera))
  meta_chimera$cell = paste0("map_", meta_chimera$cell)
  sce_chimera <- logNormCounts(sce_chimera)
  atlas_sce1 <- multiBatchNorm(atlas_sce,batch=atlas_meta$sample)
  big_sce <- multiBatchNorm(cbind(atlas_sce1, sce_chimera),
                            batch=c(rep(1,ncol(atlas_sce1)),rep(2,ncol(sce_chimera))))
  hvgs = getHVGs(big_sce,FDR=0.1)
  big_pca = prcomp_irlba(t(logcounts(big_sce[hvgs,])), n = nPC)$x
  rownames(big_pca) = colnames(big_sce) 
  pca_atlas = big_pca[1:ncol(atlas_sce),]
  map_pca = big_pca[-(1:ncol(atlas_sce)),]
  
  #correct the atlas first
  order_df = atlas_meta[!duplicated(atlas_meta$sample), c("stage", "sample")]
  order_df$ncells = sapply(order_df$sample, function(x) sum(atlas_meta$sample == x))
  order_df$stage = factor(order_df$stage, 
                          levels = rev(c("E8.5", 
                                         "E8.25", 
                                         "E8.0", 
                                         "E7.75", 
                                         "E7.5", 
                                         "E7.25", 
                                         "mixed_gastrulation", 
                                         "E7.0", 
                                         "E6.75", 
                                         "E6.5")))
  order_df = order_df[order(order_df$stage, order_df$ncells, decreasing = TRUE),]
  order_df$stage = as.character(order_df$stage)
  
  set.seed(42)
  atlas_corrected = doBatchCorrect2(counts = logcounts(atlas_sce[hvgs,]), 
                                    timepoints = atlas_meta$stage, 
                                    samples = atlas_meta$sample, 
                                    timepoint_order = order_df$stage, 
                                    sample_order = order_df$sample, 
                                    pc_override = pca_atlas)
  
  
  
  pcas_mapping = mnnMap(atlas_pca = atlas_corrected,
                        atlas_meta = atlas_meta,
                        map_pca = map_pca,
                        meta_chimera = meta_chimera,
                        return.pca = TRUE)
  mapping <- pcas_mapping[[1]]
  pcas = pcas_mapping[[2]]
  save(mapping,pcas,file=paste0(target_name,"_mapping.rda"))
  cn = substr(names(mapping), 5, nchar(names(mapping))) # remove
  ct = sapply(mapping, function(x) x$celltype.mapped)
  #ctc = sapply(mapping, function(x) x$celltype.mapped,clustering)
  st = sapply(mapping, function(x) x$stage.mapped)
  closest = sapply(mapping, function(x) x$cells.mapped[1])
  distance = sapply(1:length(closest),function(n) return(sqrt(sum((pcas$atlas[closest[n],]-pcas$mapped[names(closest)[n],])^2))))
  corr = sapply(1:length(closest),function(n) return(cor(pcas$atlas[closest[n],],
                                                         pcas$mapped[names(closest)[n],])))
  #tm = sapply(mapping, function(x) x$time.mapped)
    out = data.frame(cell = cn, celltype.mapped = ct, stage.mapped = st, closest.cell = closest,distance.to.closest.cell=distance,
                     correlation.closest.cell = corr    )
  return(out)

}







