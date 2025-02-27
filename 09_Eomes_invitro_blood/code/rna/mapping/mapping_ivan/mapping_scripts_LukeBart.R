library(irlba)
library(scran)
library(Matrix)
library(biomaRt)
library(batchelor)
library(BiocParallel)
library(SingleCellExperiment) 

cells  <- readLines("/mnt/b/iimaz/embryos/mesoderm_subset/haem_landscape/haem_full_traj_cells.txt")

sce_atlas  <- readRDS("/mnt/d/iimaz/integrated_sce_sparse_full.rds")[,cells]
meta_atlas <- data.frame(colData(sce_atlas))


setwd("/mnt/d/iimaz/LukesGrantMapping/eomes_multiome")

sce_multiome  <- SingleCellExperiment(list(counts=readRDS("raw_counts.rds")))
meta_multiome <- read.delim("sample_metadata_after_mapping.txt", sep="\t")
meta_multiome <- meta_multiome[match(colnames(sce_multiome), meta_multiome$cell), ]

mouse_ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

gene_map <- getBM(attributes=c("ensembl_gene_id","mgi_symbol"),
filters = "mgi_symbol", values = rownames(sce_multiome), mart = mouse_ensembl)

gene_map_ens <- gene_map[match(rownames(sce_multiome),gene_map$mgi_symbol), ]
rownames(sce_multiome) <- gene_map_ens$ensembl_gene_id

res <- mapWrap2020(atlas_sce = sce_atlas, atlas_meta = meta_atlas,
  map_sce = sce_multiome, map_meta = meta_multiome, k = 15, n_chunks = 20)
saveRDS(res,"/mnt/d/iimaz/LukesGrantMapping/eomes_multiome/multiome2atlas.rds")

setwd("/mnt/d/iimaz/LukesGrantMapping/eomes_public")

files <- list.files()

GEOmtx  <- files[grep("mtx", files)]
GEOcell <- files[grep("barcodes", files)]
GEOgene <- files[grep("features", files)]

for (i in seq(4)){
  counts <- readMM(GEOmtx[i])
  cells  <- readLines(GEOcell[i])
  genes  <- substr(readLines(GEOgene[i]), start = 1, stop = 18)
  colnames(counts) <- cells
  rownames(counts) <- genes
  
  sce_eomes  <- SingleCellExperiment(list(counts = counts))
  meta_eomes <- data.frame(cell = colnames(sce_eomes), stage = rep("Unknown", ncol(sce_eomes)))
  res <- mapWrap2020(atlas_sce = sce_atlas, atlas_meta = meta_atlas,
  map_sce = sce_eomes, map_meta = meta_eomes, k = 15,  n_chunks = 20)
 
  saveRDS(res, file = gsub("_matrix.mtx", ".mapping.rds", GEOmtx[i]))
}

##### Mapping functions

#MAPPING FUNCTIONS

getmode <- function(v, dist) {
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

getcelltypes <- function(v, dist) {
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

getMappingScore <- function(mapping){
  celltypes_accrossK <- matrix(unlist(mapping$celltypes.mapped), 
    nrow=length(mapping$celltypes.mapped[[1]]),
    ncol=length(mapping$celltypes.mapped))
  P <- NULL
  for (i in 1:nrow(celltypes_accrossK)){
    p <- max(table(celltypes_accrossK[i,]))
    index <- which(table(celltypes_accrossK[i,]) == p)
    p <- p/length(mapping$celltypes.mapped)
    P <- c(P,p) 
  }
  return(P)  
}

getHVGs2020 <- function(sce, min.mean = 1e-3, chrY.genes.file = NULL, sparse_matrix = FALSE, p.val_th = 0.05, computer = "ema"){

  if(computer == "ebi"){
    chrY.genes.file <- "/hps/research1/marioni/ivan/EmbryoTimeCourse/atlas/data/ygenes.tab"
  }else if(computer == "ema"){
    chrY.genes.file <- "/mnt/b/iimaz/embryos/ygenes.tab"
  }

  #exclude sex genes and tomato-td
  xist <- "ENSMUSG00000086503"
  if(!is.null(chrY.genes.file)){
    ychr <- read.table(chrY.genes.file, stringsAsFactors = FALSE)[,1]
  }else{
    mouse_ensembl <- biomaRt::useMart("ensembl", dataset="mmusculus_gene_ensembl")
    gene_map      <- biomaRt::getBM(attributes=c("ensembl_gene_id", "chromosome_name"),
      filters = "ensembl_gene_id", values = rownames(decomp), mart = mouse_ensembl)
    ychr <- gene_map[gene_map[,2] == "Y", 1]  
  }  
  other <- c("tomato-td") #for the chimera
  sce   <- sce[!rownames(sce) %in% c(xist, ychr, other),]

  dec <- modelGeneVar(sce)
  out <- getTopHVGs(dec, fdr.threshold = p.val_th)

  message("Number of highly variable genes: ",length(out))

  return(out) 
}  

#ensure counts has columns names for the cells
#match timepoints,samples to the count table
#timepoint_order, sample_order should contain each sample/timepoint ONCE, in correct order for correction
doBatchCorrect = function(counts, timepoints, samples, timepoint_order, sample_order, npc = 50, pc_override = NULL, BPPARAM = SerialParam()){
  require(scran)
  require(irlba)
  require(BiocParallel)
  
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
      return(do.call(fastMNN, c(x, "pc.input" = TRUE, BPPARAM = BPPARAM))$corrected)
    } else {
      return(x[[1]])
    }
  })
  
  #perform correction over list
  if(length(correct_list)>1){
    correct = do.call(fastMNN, c(correct_list, "pc.input" = TRUE, BPPARAM = BPPARAM))$corrected
  } else {
    correct = correct_list[[1]]
  }
  
  correct = correct[match(colnames(counts), rownames(correct)),]
  
  return(correct)
  
}

mnnMap_extended = function(atlas_pca, atlas_meta, map_pca, map_meta, k_map = 10, mixed_gastrulation_map = FALSE){
  require(BiocNeighbors)
  require(scran)
  correct = fastMNN(atlas_pca, map_pca, pc.input = TRUE)$corrected
  atlas = 1:nrow(atlas_pca)
  correct_atlas = correct[atlas,]
  correct_map = correct[-atlas,]
  #browser()
  knns = BiocNeighbors::queryKNN(correct_atlas, correct_map, k = k_map, get.index = TRUE, get.distance = FALSE)  #get closest k matching cells
  k.mapped   = t(apply(knns$index, 1, function(x) atlas_meta$cell[x]))
  celltypes  = t(apply(k.mapped, 1, function(x) atlas_meta$celltype_extended_atlas[match(x, atlas_meta$cell)]))
  stages     = t(apply(k.mapped, 1, function(x) atlas_meta$stage[match(x, atlas_meta$cell)]))
  if(mixed_gastrulation_map){
    stages.ext = t(apply(k.mapped, 1, function(x) atlas_meta$stage.mapped[match(x, atlas_meta$cell)]))
  }
  if (k_map > 1){
    celltype.mapped = apply(celltypes, 1, function(x) getmode(x, 1:length(x)))
    stage.mapped    = apply(stages, 1, function(x) getmode(x, 1:length(x)))
    if(mixed_gastrulation_map){
      stage.mapped.extended = apply(stages.ext, 1, function(x) getmode(x, 1:length(x)))
      out = lapply(1:length(celltype.mapped), function(x){
        list(cells.mapped = k.mapped[x,],
          celltype.mapped = celltype.mapped[x],
          stage.mapped = stage.mapped[x],
          stage.mapped.extended=stage.mapped.extended[x],
          celltypes.mapped= celltypes[x,],
          stages.mapped=stages[x,])
         
      })
    }else{
      out = lapply(1:length(celltype.mapped), function(x){
        list(cells.mapped = k.mapped[x,],
          celltype.mapped = celltype.mapped[x],
          stage.mapped = stage.mapped[x],
          celltypes.mapped= celltypes[x,],
          stages.mapped=stages[x,])

      })
    }  
  }else{
    celltype.mapped = celltypes
    stage.mapped    = stages
    if(mixed_gastrulation_map){
      stage.mapped.extended = stages.ext
      out = lapply(1:length(celltype.mapped), function(x){
        list(cells.mapped = k.mapped[x],
          celltype.mapped = celltype.mapped[x],
          stage.mapped = stage.mapped[x],
          stage.mapped.extended=stage.mapped.extended[x],
          celltypes.mapped= celltypes[x],
          stages.mapped=stages[x])
      })
    }else{
      out = lapply(1:length(celltype.mapped), function(x){
        list(cells.mapped = k.mapped[x],
          celltype.mapped = celltype.mapped[x],
          stage.mapped = stage.mapped[x],
          celltypes.mapped= celltypes[x],
          stages.mapped=stages[x])
      })
    }
}
    
  names(out) = map_meta$cell
  
  return(out)
  
}


mapWrap2020 = function(atlas_sce, atlas_meta, map_sce, map_meta, k=15, nPCs=50, mapstage_x=NULL, map2stage_x=NULL, return.list = FALSE, hvgs.p.val_th = 0.05, normalised = FALSE, mixed_gastrulation_map = FALSE, n_chunks = NULL, computer = "ema"){
 
  message("Matching features... ") 
  genes_shared  <- intersect(rownames(atlas_sce), rownames(map_sce))
  atlas_sce     <- atlas_sce[genes_shared, ]
  map_sce       <- map_sce[genes_shared, ]
  message(nrow(atlas_sce), "\n")
 
  #prevent duplicate rownames
  colnames(map_sce) <- paste0("map_", colnames(map_sce))
  map_meta$cell     <- paste0("map_", map_meta$cell)
  
  if (!normalised){
    message("Normalizing joint dataset...")
    #easier to avoid directly binding sce objects as it is a lot more likely to have issues
    big_sce <- SingleCellExperiment::SingleCellExperiment(
      list(counts=cbind(counts(atlas_sce), counts(map_sce))))
    big_sce <- scater::normalize(big_sce)
    message("Done\n")
  }else{
    message("Joining dataset...")
    big_sce <- SingleCellExperiment::SingleCellExperiment(
      list(logcounts=cbind(logcounts(atlas_sce), logcounts(map_sce))))
    message("Done\n")
  }

  message("Computing highly variable genes...")
  #hvgs <- getHVGs(big_sce, computer = computer, p.val_th = hvgs.p.val_th)
  hvgs <- getHVGs2020(big_sce, computer = computer, p.val_th = hvgs.p.val_th)
  message("Done\n")
  
  message("Performing PCA...")
  big_pca <- irlba::prcomp_irlba(t(logcounts(big_sce[hvgs,])), n = nPCs)$x
  rownames(big_pca) <- colnames(big_sce) 
  atlas_pca <- big_pca[1:ncol(atlas_sce),]
  map_pca   <- big_pca[-(1:ncol(atlas_sce)),] 
  message("Done\n")
    
  message("Atlas batch effect correction...")  
  #correct the atlas first
  order_df        <- atlas_meta[!duplicated(atlas_meta$sample), c("stage", "sample")]
  order_df$ncells <- sapply(order_df$sample, function(x) sum(atlas_meta$sample == x))
  order_df$stage  <- factor(order_df$stage, 
                          levels = rev(c(
                           "E9.5",
                           "E9.25",
                           "E9.0",
                           "E8.75",
                           "E8.5", 
                           "E8.25", 
                           "E8.0", 
                           "E7.75", 
                           "E7.5", 
                           "E7.25", 
                           "mixed_gastrulation", 
                           "E7.0", 
                           "E6.75", 
                           "E6.5")))
  order_df       <- order_df[order(order_df$stage, order_df$ncells, decreasing = TRUE),]
  order_df$stage <- as.character(order_df$stage)
  
  set.seed(42)
  atlas_corrected <- doBatchCorrect(counts         = counts(atlas_sce[hvgs,]), # this should be logcounts but they are override by the pca
                                   timepoints      = atlas_meta$stage, 
                                   samples         = atlas_meta$sample, 
                                   timepoint_order = order_df$stage, 
                                   sample_order    = order_df$sample, 
                                   pc_override     = atlas_pca)
  message("Done\n")
  
  if(is.null(n_chunks)){
    message("MNN mapping...")                        
    # Check if mapping from/to specific embryo stages 
    if(!is.null(map2stage_x)){
      atlas_stage_index <- NULL
      for (i in seq(from = 1, to = length(map2stage_x))){
    	  index1 <- which(atlas_meta$stage == map2stage_x[i])
  	    atlas_stage_index <- c(atlas_stage_index, index1)
      } 
    }
    if(!is.null(mapstage_x)){
      map_stage_index  <- NULL
      for (i in seq(from = 1, to = length(mapstage_x))){  	
  	    index2 <- which(map_meta$stage == mapstage_x[i])
  	    map_stage_index <- c(map_stage_index, index2)
      }
    }
  
    # Select mapping scheme
    if(!is.null(map2stage_x) & !is.null(mapstage_x)){
      mapping <- mnnMap_extended(atlas_pca = atlas_corrected[atlas_stage_index,],
                   atlas_meta = atlas_meta[atlas_stage_index,],
                   map_pca    = map_pca[map_stage_index,],
                   map_meta   = map_meta[map_stage_index,], k_map = k,
                   mixed_gastrulation_map = mixed_gastrulation_map) 
    }else if(!is.null(map2stage_x) & is.null(mapstage_x)){
      mapping <- mnnMap_extended(atlas_pca = atlas_corrected[atlas_stage_index,],
                   atlas_meta = atlas_meta[atlas_stage_index,],
                   map_pca    = map_pca,
                   map_meta   = map_meta, k_map = k,
                   mixed_gastrulation_map = mixed_gastrulation_map)   
    }else if(is.null(map2stage_x) & !is.null(mapstage_x)){
      mapping <- mnnMap_extended(atlas_pca = atlas_corrected,
                   atlas_meta = atlas_meta,
                   map_pca    = map_pca[map_stage_index,],
                   map_meta   = map_meta[map_stage_index,], k_map = k,
                   mixed_gastrulation_map = mixed_gastrulation_map) 
    }else{
      mapping <- mnnMap_extended(atlas_pca = atlas_corrected,
                   atlas_meta = atlas_meta,
                   map_pca    = map_pca,
                   map_meta   = map_meta, k_map = k,
                   mixed_gastrulation_map = mixed_gastrulation_map)
    }
  }else{
    message("Big query MNN mapping...")      
    n     <- nrow(map_pca)
    set.seed(42)
    random_index <- sample(1:n)
    chunk_size     <- ceiling(n/n_chunks)
    chunks_indexes <- split(random_index, ceiling(seq_along(random_index)/chunk_size))
    
    mapping <- NULL
    for (chunk in 1:n_chunks){
       message("Processing chunk ", chunk)
       mapping_tmp <- mnnMap_extended(atlas_pca = atlas_corrected,
                   atlas_meta = atlas_meta,
                   map_pca    = map_pca[chunks_indexes[[chunk]], ],
                   map_meta   = map_meta[chunks_indexes[[chunk]],],
                   k_map = k, mixed_gastrulation_map = mixed_gastrulation_map)                   
      df <- data.frame(
        cell            = names(mapping_tmp), 
        celltype.mapped = sapply(mapping_tmp, function(x) x$celltype.mapped),
        stage.mapped    = sapply(mapping_tmp, function(x) x$stage.mapped),
        closest.cell    = sapply(mapping_tmp, function(x) x$cells.mapped[1])
      )                  
      mapping     <- rbind(mapping, df)
    }
  } 
  message("Done\n")

  if(!is.null(n_chunks)){
    message("Generating output...")
    mapping <- mapping[match(1:n, random_index),]
    message("Done\n")
    return(mapping)
  }
          
  if(return.list){
    message("Generating output...")
    message("Done\n")
    return(mapping)
  }
  
  message("Computing mapping scores...") 
  out <- list()
  for (i in seq(from = 1, to = k)){
   out$closest.cells[[i]]     <- sapply(mapping, function(x) x$cells.mapped[i])
   out$celltypes.mapped[[i]]  <- sapply(mapping, function(x) x$celltypes.mapped[i])
   out$cellstages.mapped[[i]] <- sapply(mapping, function(x) x$stages.mapped[i])
  }  
  celltype.multinomial.prob <- getMappingScore(out)
  message("Done\n")
  
  message("Generating output...") 
  if(mixed_gastrulation_map){
    out$mapping <- data.frame(
      cell            = names(mapping), 
      celltype.mapped = sapply(mapping, function(x) x$celltype.mapped),
      stage.mapped    = sapply(mapping, function(x) x$stage.mapped),
      stage.mapped.extended =  sapply(mapping, function(x) x$stage.mapped.extended),
      closest.cell    = sapply(mapping, function(x) x$cells.mapped[1]))
    out$mapping <- cbind(out$mapping, celltype.multinomial.prob)
  }else{
    out$mapping <- data.frame(
      cell            = names(mapping), 
      celltype.mapped = sapply(mapping, function(x) x$celltype.mapped),
      stage.mapped    = sapply(mapping, function(x) x$stage.mapped),
      closest.cell    = sapply(mapping, function(x) x$cells.mapped[1]))
    out$mapping <- cbind(out$mapping, celltype.multinomial.prob)  
  }
  message("Done\n")
  
  return(out)
  
}
