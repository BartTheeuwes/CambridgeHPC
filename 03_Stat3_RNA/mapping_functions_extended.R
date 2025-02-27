in_folder = "/rds/project/rds-SDzz0CATGms/users/bt392/03_Stat3_RNA/"

library(cowplot)

doBatchCorrect <- function(counts, timepoints, samples, timepoint_order, sample_order, npc = 50, pc_override = NULL, BPPARAM = SerialParam()){
  require(BiocParallel)
  require(batchelor)
  if(!is.null(pc_override)){
    pca = pc_override
  } else {
    pca = irlba::prcomp_irlba(t(counts), n = npc)$x
    rownames(pca) = colnames(counts)
  }
  
  if(length(unique(samples)) == 1){
    return(pca)
  }
  
  #create nested list
  pc_list    <- lapply(unique(timepoints), function(tp){
    sub_pc   <- pca[timepoints == tp, , drop = FALSE]
    sub_samp <- samples[timepoints == tp]
    list     <- lapply(unique(sub_samp), function(samp){
      sub_pc[sub_samp == samp, , drop = FALSE]
    })
    names(list) <- unique(sub_samp)
    return(list)
  })
  
  names(pc_list) <- unique(timepoints)
  
  #arrange to match timepoint order
  pc_list <- pc_list[order(match(names(pc_list), timepoint_order))]
  pc_list <- lapply(pc_list, function(x){
    x[order(match(names(x), sample_order))]
  })
  
  #perform corrections within list elements (i.e. within stages)
  correct_list <- lapply(pc_list, function(x){
    if(length(x) > 1){
        #return(do.call(scran::fastMNN, c(x, "pc.input" = TRUE, BPPARAM = BPPARAM))$corrected)
        return(do.call(reducedMNN, c(x, BPPARAM = BPPARAM))$corrected) # edited 09.02. because of "Error: 'fastMNN' is not an exported object from 'namespace:scran'", 17.02. changed to reducedMNN because otherwise it thinks PCA space is logcounts which would be utter bullcrap
    } else {
      return(x[[1]])
    }
  })
  
  #perform correction over list
  if(length(correct_list)>1){
      #correct <- do.call(scran::fastMNN, c(correct_list, "pc.input" = TRUE, BPPARAM = BPPARAM))$corrected
      correct <- do.call(reducedMNN, c(correct_list, BPPARAM = BPPARAM))$corrected # edited 09.02. because of "Error: 'fastMNN' is not an exported object from 'namespace:scran'", 17.02. changed to reducedMNN because otherwise it thinks PCA space is logcounts which would be utter bullcrap
  } else {
    correct <- correct_list[[1]]
  }
    
  correct <- correct[match(colnames(counts), rownames(correct)),]
  
  return(correct)

}

getHVGs <- function(sce, block, min.mean = 1e-3){
  decomp <- modelGeneVar(sce, block=block)
  decomp <- decomp[decomp$mean > min.mean,]
    
  #exclude sex genes & tomato
  xist = "ENSMUSG00000086503"
  ychr = read.table("/rds/project/rds-SDzz0CATGms/users/bt392/software/ygenes.tab", stringsAsFactors = FALSE)[,1]
  other = c("tomato-td") #for the chimera
  decomp = decomp[!rownames(decomp) %in% c(xist, ychr, other),]  
    
  decomp$FDR <- p.adjust(decomp$p.value, method = "fdr")
  return(rownames(decomp)[decomp$p.value < 0.05])
}

getmode <- function(v, dist) {
  tab <- table(v)
  #if tie, break to shortest distance
  if(sum(tab == max(tab)) > 1){
    tied <- names(tab)[tab == max(tab)]
    sub  <- dist[v %in% tied]
    names(sub) <- v[v %in% tied]
    return(names(sub)[which.min(sub)])
  } else {
    return(names(tab)[which.max(tab)])
  }
}

getcelltypes <- function(v, dist) {
  tab <- table(v)
  #if tie, break to shortest distance
  if(sum(tab == max(tab)) > 1){
    tied <- names(tab)[tab == max(tab)]
    sub  <- dist[v %in% tied]
    names(sub) <- v[v %in% tied]
    return(names(sub)[which.min(sub)])
  } else {
    return(names(tab)[which.max(tab)])
  }
}

getMappingScore <- function(mapping){
    out <- list()
    celltypes_accrossK <- matrix(unlist(mapping$celltypes.mapped),
                                 nrow=length(mapping$celltypes.mapped[[1]]),
                                 ncol=length(mapping$celltypes.mapped))
    
   celltype_originals_accrossK <- matrix(unlist(mapping$celltype_originals.mapped), # ADDED
                                 nrow=length(mapping$celltype_originals.mapped[[1]]),
                                 ncol=length(mapping$celltype_originals.mapped))
    
    cellstages_accrossK <- matrix(unlist(mapping$cellstages.mapped),
                                  nrow=length(mapping$cellstages.mapped[[1]]),
                                  ncol=length(mapping$cellstages.mapped))
    out$celltype.score <- NULL
    for (i in 1:nrow(celltypes_accrossK)){
        p <- max(table(celltypes_accrossK[i,]))
        index <- which(table(celltypes_accrossK[i,]) == p)
        p <- p/length(mapping$celltypes.mapped)
        out$celltype.score <- c(out$celltype.score,p)
    }
    
    out$celltype_original.score <- NULL # ADDED
    for (i in 1:nrow(celltype_originals_accrossK)){
        p <- max(table(celltype_originals_accrossK[i,]))
        index <- which(table(celltype_originals_accrossK[i,]) == p)
        p <- p/length(mapping$celltype_originals.mapped)
        out$celltype_original.score <- c(out$celltype_original.score,p)
    }
    
    out$cellstage.score <- NULL
    for (i in 1:nrow(cellstages_accrossK)){
        p <- max(table(cellstages_accrossK[i,]))
        index <- which(table(cellstages_accrossK[i,]) == p)
        p <- p/length(mapping$cellstages.mapped)
        out$cellstage.score <- c(out$cellstage.score,p)
    }
    return(out)  
}

get_meta <- function(correct_atlas, atlas_meta, correct_map, map_meta, k_map = 10){
  knns <- BiocNeighbors::queryKNN(correct_atlas, correct_map, k = k_map, get.index = TRUE,
    get.distance = FALSE)
  #get closest k matching cells
  k.mapped  <- t(apply(knns$index, 1, function(x) atlas_meta$cell[x]))
  celltypes <- t(apply(k.mapped, 1, function(x) atlas_meta$celltype[match(x, atlas_meta$cell)]))
                       
  celltype_originals <- t(apply(k.mapped, 1, function(x) atlas_meta$celltype_original[match(x, atlas_meta$cell)]))
                       
  stages    <- t(apply(k.mapped, 1, function(x) atlas_meta$stage[match(x, atlas_meta$cell)]))
  celltype.mapped <- apply(celltypes, 1, function(x) getmode(x, 1:length(x)))
                             
  celltype_original.mapped <- apply(celltype_originals, 1, function(x) getmode(x, 1:length(x)))
                           
  stage.mapped    <- apply(stages, 1, function(x) getmode(x, 1:length(x)))
  out <- lapply(1:length(celltype.mapped), function(x){
    list(cells.mapped     = k.mapped[x,],
         celltype.mapped  = celltype.mapped[x],
         
         celltype_original.mapped  = celltype_original.mapped[x],

         stage.mapped     = stage.mapped[x],
         celltypes.mapped = celltypes[x,],
         
         celltype_originals.mapped = celltype_originals[x,],

         stages.mapped    = stages[x,])
  })
  names(out) <- map_meta$cell
  return(out)  
}

mapWrap <- function(atlas_sce, atlas_meta, map_sce, map_meta, order = NULL, k = 30, npcs = 50, genes = NULL, return.list = FALSE){
    
  message("Normalizing joint dataset...")
  
  #easier to avoid directly binding sce objects as it is a lot more likely to have issues
  sce_all <- SingleCellExperiment::SingleCellExperiment(
    list(counts=Matrix::Matrix(cbind(counts(atlas_sce),counts(map_sce)),sparse=TRUE)))
  
  #big_sce <- scater::normalize(sce_all)
  #big_sce <- scater::logNormCounts(sce_all) # edited 09.02. because normalize deprecated in favour of logNormCounts
  big_sce <- multiBatchNorm(sce_all, batch=c(atlas_meta$sample, map_meta$sample)) # edited 17.02. because now multibatchnorm exists
  message("Done\n")
  
  if (is.null(genes)) {
    message("Genes not provided. Computing highly variable genes...")
    hvgs <- getHVGs(big_sce, block=c(atlas_meta$sample, map_meta$sample))
    message("Done\n")
  } else {
    hvgs <- genes
    message(sprintf("%d Genes provided...",length(genes)))
  }
  
  message("Performing PCA...")
  big_pca <- multiBatchPCA(big_sce,
                           batch=c(atlas_meta$sample, map_meta$sample),
                           subset.row = hvgs,
                           d = npcs,
                           preserve.single = TRUE,
                           assay.type = "logcounts")[[1]]
  rownames(big_pca) <- colnames(big_sce) 
  atlas_pca <- big_pca[1:ncol(atlas_sce),]
  map_pca   <- big_pca[-(1:ncol(atlas_sce)),]
  message("Done\n")
  
  message("Batch effect correction for the atlas...")  
  order_df        <- atlas_meta[!duplicated(atlas_meta$sample), c("stage", "sample")]
  order_df$ncells <- sapply(order_df$sample, function(x) sum(atlas_meta$sample == x))
                                         
  order_df$stage  <- factor(order_df$stage, 
                        levels = rev(c("E9.5",
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
  atlas_corrected <- doBatchCorrect(counts         = logcounts(atlas_sce[hvgs,]), 
                                    timepoints      = atlas_meta$stage, 
                                    samples         = atlas_meta$sample, 
                                    timepoint_order = order_df$stage, 
                                    sample_order    = order_df$sample, 
                                    pc_override     = atlas_pca,
                                    npc             = npcs)
  message("Done\n")
  
  message("MNN mapping...")                        
  correct <- reducedMNN(rbind(atlas_corrected, map_pca),
                      batch=c(rep("ATLAS", dim(atlas_meta)[1]), map_meta$sample),
                      merge.order=order)$corrected
  atlas   <- 1:nrow(atlas_pca)
  correct_atlas <- correct[atlas,]
  correct_map   <- correct[-atlas,]

  mapping <- get_meta(correct_atlas = correct_atlas,
                      atlas_meta = atlas_meta,
                      correct_map = correct_map,
                      map_meta = map_meta,
                      k_map = k)
  message("Done\n")
  
  if(return.list){
    return(mapping)
  }
  
  message("Computing mapping scores...") 
  out <- list()
  for (i in seq(from = 1, to = k)) {
    out$closest.cells[[i]]     <- sapply(mapping, function(x) x$cells.mapped[i])
    out$celltypes.mapped[[i]]  <- sapply(mapping, function(x) x$celltypes.mapped[i])
                                         
    out$celltype_originals.mapped[[i]]  <- sapply(mapping, function(x) x$celltype_originals.mapped[i])

    out$cellstages.mapped[[i]] <- sapply(mapping, function(x) x$stages.mapped[i])
  }  
  multinomial.prob <- getMappingScore(out)
  message("Done\n")
  
  message("Writing output...") 
  out$correct_atlas <- correct_atlas
  out$correct_map <- correct_map
  ct <- sapply(mapping, function(x) x$celltype.mapped); is.na(ct) <- lengths(ct) == 0
               
  cto <- sapply(mapping, function(x) x$celltype_original.mapped); is.na(cto) <- lengths(cto) == 0

  st <- sapply(mapping, function(x) x$stage.mapped); is.na(st) <- lengths(st) == 0
  cm <- sapply(mapping, function(x) x$cells.mapped[1]); is.na(cm) <- lengths(cm) == 0
  out$mapping <- data.frame(
      cell            = names(mapping), 
      celltype.mapped = unlist(ct),
      
      celltype_original.mapped = unlist(cto),

      stage.mapped    = unlist(st),
      closest.cell    = unlist(cm))
  
  out$mapping <- cbind(out$mapping,multinomial.prob)
  out$pca <- big_pca
  message("Done\n")
  
  return(out)

}


load_data = function(normalise = TRUE, remove_doublets = FALSE, remove_stripped = FALSE, load_corrected = FALSE){
  
  if(load_corrected & (!remove_doublets | !remove_stripped)){
    message("Using corrected PCs, also removing doublets + stripped now.")
    remove_doublets = TRUE
    remove_stripped = TRUE
  }
  
  require(scran)
  require(scater)
  require(SingleCellExperiment)
  require(Matrix)
  
  counts = readRDS(paste0(in_folder, "03_quality_control/raw_counts.rds"))
  genes = read.table(paste0(in_folder, "02_calledCells/genes.tsv"), stringsAsFactors = F)

if(file.exists(paste0(in_folder, 'meta.tab'))){
    meta = read.csv2(paste0(in_folder, 'meta.tab'), sep='\t')
}else if(file.exists(paste0(in_folder, '05_doublet_detection/meta.tab'))){
        meta = read.csv2(paste0(in_folder, '05_doublet_detection/meta.tab'), sep='\t')
}else{
         meta = read.table(paste0(in_folder, "03_quality_control/meta.tab"), header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "$")
}
  rownames(counts) = genes[,1] #ensembl
  colnames(counts) = meta$cell

  sce = SingleCellExperiment(assays = list("counts" = counts))
  
  if(normalise){
    sfs = read.table(paste0(in_folder, "04_normalisation/sizefactors.tab"), stringsAsFactors = F)[,1]
    sizeFactors(sce) = sfs
    sce = scater::logNormCounts(sce)
  }
  
  if(remove_doublets){
    sce = scater::logNormCounts(sce[,!meta$doublet])
    meta = meta[!meta$doublet,]
  }
  
  if(remove_stripped){
    sce = scater::logNormCounts(sce[,!meta$stripped])
    meta = meta[!meta$stripped, ]
  }
  
  if(load_corrected){
    corrected = readRDS(paste0(in_folder, "corrected_pcas.rds"))
    assign("corrected", corrected, envir = .GlobalEnv)
    
  }
  

  
  
  assign("genes", genes, envir = .GlobalEnv)
  assign("meta", meta, envir = .GlobalEnv)
  assign("sce", sce, envir = .GlobalEnv)
  
  
  invisible(0)
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
  
  atlas_in = '/rds/project/rds-SDzz0CATGms/users/bt392/mouse/Mixl1_KO/atlas/atlas/'
  counts = readMM(paste0(atlas_in, 'raw_counts.mtx')) # readMM for .mtx input 
  genes = read.table(paste0(atlas_in, 'genes.tsv'), stringsAsFactors = F)
  meta = read.table(paste0(atlas_in, 'meta.tab'), header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "$")
  
  rownames(counts) = genes[,1] #ensembl
  colnames(counts) = meta$cell
  
  sce = SingleCellExperiment(assays = list("counts" = counts))
  
  if(normalise){
    sfs = read.table(paste0(atlas_in, 'sizefactors.tab'), stringsAsFactors = F)[,1]
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
    corrected = readRDS(paste0(atlas_in, "corrected_pcas.rds"))
    assign("corrected", corrected, envir = .GlobalEnv)
    
  }
  
  
  
  
  assign("genes", genes, envir = .GlobalEnv)
  assign("meta", meta, envir = .GlobalEnv)
  assign("sce", sce, envir = .GlobalEnv)
  
  
  invisible(0)
}



load_extended_atlas_data = function(normalise = TRUE, remove_doublets = FALSE, remove_stripped = FALSE, load_corrected = FALSE){
  
  # This function was taken from 
  # https://github.com/MarioniLab/EmbryoTimecourse2018/blob/master/analysis_scripts/atlas/core_functions.R
  # with minor adaptations
  
  if(load_corrected & (!remove_doublets | !remove_stripped)){
    message("Using corrected PCs, also removing doublets + stripped now.")
    remove_doublets = TRUE
    remove_stripped = TRUE
  }
  
  atlas_in = '/rds/project/rds-SDzz0CATGms/users/bt392/atlasses/extended/'
  counts = readMM(paste0(atlas_in, 'raw_counts.mtx')) # readMM for .mtx input 
  genes = read.table(paste0(atlas_in, 'genes.tsv'), stringsAsFactors = F)
  meta = read.table(paste0(atlas_in, 'meta.tab'), header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "$")
  
  rownames(counts) = genes[,1] #ensembl
  colnames(counts) = meta$cell
  
  sce = SingleCellExperiment(assays = list("counts" = counts))
  
  if(normalise){
    sfs = read.table(paste0(atlas_in, 'sizefactors.tab'), stringsAsFactors = F)[,1]
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
    corrected = readRDS(paste0(atlas_in, "corrected_pcas.rds"))
    assign("corrected", corrected, envir = .GlobalEnv)
    
  }
   
  
  assign("genes", genes, envir = .GlobalEnv)
  assign("meta", meta, envir = .GlobalEnv)
  assign("sce", sce, envir = .GlobalEnv)
  
  
  invisible(0)
}

               


celltype_colours = c("Epiblast" = "#635547",
                     "Primitive Streak" = "#DABE99",
                     "Caudal epiblast" = "#9e6762",
                     
                     "PGC" = "#FACB12",
                     
                     "Anterior Primitive Streak" = "#c19f70",
                     "Notochord" = "#0F4A9C",
                     "Def. endoderm" = "#F397C0",
                     "Gut" = "#EF5A9D",
                     
                     "Nascent mesoderm" = "#C594BF",
                     "Mixed mesoderm" = "#DFCDE4",
                     "Intermediate mesoderm" = "#139992",
                     "Caudal Mesoderm" = "#3F84AA",
                     "Paraxial mesoderm" = "#8DB5CE",
                     "Somitic mesoderm" = "#005579",
                     "Pharyngeal mesoderm" = "#C9EBFB",
                     "Cardiomyocytes" = "#B51D8D",
                     "Allantois" = "#532C8A",
                     "ExE mesoderm" = "#8870ad",
                     "Mesenchyme" = "#cc7818",
                     
                     "Haematoendothelial progenitors" = "#FBBE92",
                     "Endothelium" = "#ff891c",
                     "Blood progenitors 1" = "#f9decf",
                     "Blood progenitors 2" = "#c9a997",
                     "Erythroid1" = "#C72228",
                     "Erythroid2" = "#f79083",
                     "Erythroid3" = "#EF4E22",
                     
                     "NMP" = "#8EC792",
                     
                     "Rostral neurectoderm" = "#65A83E",
                     "Caudal neurectoderm" = "#354E23",
                     "Neural crest" = "#C3C388",
                     "Forebrain/Midbrain/Hindbrain" = "#647a4f",
                     "Spinal cord" = "#CDE088",
                     
                     "Surface ectoderm" = "#f7f79e",
                     
                     "Visceral endoderm" = "#F6BFCB",
                     "ExE endoderm" = "#7F6874",
                     "ExE ectoderm" = "#989898",
                     "Parietal endoderm" = "#1A1A1A")

# +
celltype_colours_final = c(
"Epiblast" = "#635547",
"Primitive Streak" = "#DABE99",
"Caudal epiblast" = "#9e6762",

"PGC" = "#FACB12",

"Anterior Primitive Streak" = "#c19f70",
"Node"="#153b3d",
"Notochord" = "#0F4A9C",



"Gut tube" = "#EF5A9D",
"Hindgut" = "#F397C0",
"Midgut" = "#ff00b2",
"Foregut" = "#ffb7ff",
"Pharyngeal endoderm"="#95e1ff",
"Thyroid primordium"="#97bad3",

"Nascent mesoderm" = "#C594BF",
"Intermediate mesoderm" = "#139992",
"Caudal mesoderm" = "#3F84AA",
"Lateral plate mesoderm" = "#F9DFE6",
"Limb mesoderm" = "#e35f82",
"Forelimb" = "#d02d75",
"Kidney primordium" = "#e85639",
"Presomitic mesoderm"="#5581ca",#"#0000ff",#blue
"Somitic mesoderm" = "#005579",
"Posterior somitic tissues" = "#5adbe4",#"#40e0d0",#turquoise



"Paraxial mesoderm" = "#8DB5CE",
"Cranial mesoderm" = "#456722",#"#006400",#darkgreen
"Anterior somitic tissues"= "#d5e839",
"Sclerotome" = "#e3cb3a",#"#ffff00",#yellow
"Dermomyotome" = "#00BFC4",#"#a52a2a",#brown



"Pharyngeal mesoderm" = "#C9EBFB",
"Cardiopharyngeal progenitors" = "#556789",
"Anterior cardiopharyngeal progenitors"="#683ed8",



"Allantois" = "#532C8A",
"Mesenchyme" = "#cc7818",
"YS mesothelium" = "#ff7f9c",
"Epicardium"="#f79083",
"Embryo proper mesothelium" = "#ff487d",



"Cardiopharyngeal progenitors FHF"="#d780b0",
"Cardiomyocytes FHF 1"="#a64d7e",
"Cardiomyocytes FHF 2"="#B51D8D",



"Cardiopharyngeal progenitors SHF"="#4b7193",
"Cardiomyocytes SHF 1"="#5d70dc",
"Cardiomyocytes SHF 2"="#332c6c",



"Haematoendothelial progenitors" = "#FBBE92",
"Blood progenitors" = "#6c4b4c",
"Erythroid" = "#C72228",
"Chorioallantoic-derived erythroid progenitors"="#E50000",
"Megakaryocyte progenitors"="#e3cb3a",
"MEP"="#EF4E22",
"EMP"="#7c2a47",



"YS endothelium"="#ff891c",
"YS mesothelium-derived endothelial progenitors"="#AE3F3F",
"Allantois endothelium"="#2f4a60",
"Embryo proper endothelium"="#90e3bf",
"Venous endothelium"="#bd3400",
"Endocardium"="#9d0049",

"NMPs/Mesoderm-biased" = "#89c1f5",
"NMPs" = "#8EC792",

"Ectoderm" = "#ff675c",



"Optic vesicle" = "#bd7300",

"Ventral forebrain progenitors"="#a0b689",
"Early dorsal forebrain progenitors"="#0f8073",
"Late dorsal forebrain progenitors"="#7a9941",
"Midbrain/Hindbrain boundary"="#8ab3b5",
"Midbrain progenitors"="#9bf981",
"Dorsal midbrain neurons"="#12ed4c",
"Ventral hindbrain progenitors"="#7e907a",
"Dorsal hindbrain progenitors"="#2c6521",
"Hindbrain floor plate"="#bf9da8",
"Hindbrain neural progenitors"="#59b545",



"Neural tube"="#233629",



"Migratory neural crest"="#4a6798",
"Branchial arch neural crest"="#bd84b0",
"Frontonasal mesenchyme"="#d3b1b1",



"Spinal cord progenitors"="#6b2035",
"Dorsal spinal cord progenitors"="#e273d6",

"Non-neural ectoderm" = "#f7f79e",
"Surface ectoderm" = "#fcff00",
"Epidermis" = "#fff335",
"Limb ectoderm" = "#ffd731",
"Amniotic ectoderm" = "#dbb400",



"Placodal ectoderm" = "#ff5c00",



"Otic placode"="#f1a262",
"Otic neural progenitors"="#00b000",

"Visceral endoderm" = "#F6BFCB",
"ExE endoderm" = "#7F6874",
"ExE ectoderm" = "#989898",
"Parietal endoderm" = "#1A1A1A"

)
# -

haem_colours = c(
  "Mes1"= "#c4a6b2",
  "Mes2"= "#ca728c",
  
  "Cardiomyocytes" =  "#B51D8D",  
  
  "BP1" = "#6460c5",
  "BP2" = "#96b8e4",
  "Haem3"= "#02f9ff",
  "BP3" = "#07499f",
  "BP4" = "#036ef9",
  
  "Haem1"= "#bb22a7",
  "Haem2" = "#f695e9",
  "Haem4" = "#4c4a81",
  
  "EC1"= "#006737",
  
  "EC2" = "#5eba46",
  "EC3" = "#818068",
  "EC4"="#d6de22",
  "EC5"="#5f7238",
  "EC6"="#2ab572",
  "EC7"="#000000",
  "EC8"="#a0cb3b",
  
  "Ery1"="#f67a58",
  "Ery2" ="#a26724",
  "Ery3"="#cdaf7f",
  "Ery4"= "#625218",
  
  "My" = "#c62127",
  "Mk"= "#f6931d"
)

stage_colours = c("E6.5" = "#D53E4F",
                  "E6.75" = "#F46D43",
                  "E7.0" = "#FDAE61",
                  "E7.25" = "#FEE08B",
                  "E7.5" = "#FFFFBF",
                  "E7.75" = "#E6F598",
                  "E8.0" = "#ABDDA4",
                  "E8.25" = "#66C2A5",
                  "E8.5" = "#3288BD",
                  "mixed_gastrulation" = "#A9A9A9")

stage_colours_final = c(
                  "mixed_gastrulation" = "#A9A9A9",
                  "E6.5" = "#D53E4F",
                  "E6.75" = "#F46D43",
                  "E7.0" = "#FDAE61",
                  "E7.25" = "#FEE08B",
                  "E7.5" = "#FFFFBF",
                  "E7.75" = "#E6F598",
                  "E8.0" = "#ABDDA4",
                  "E8.25" = "#66C2A5",
                  "E8.5" = "#3288BD",
                    "E8.75" = '#3C1ACE',
                    "E9.0" = '#9A28F7' ,
                    "E9.25" =  '#E228F7',
                    "E9.5" = '#F728B8' )

umap_theme =  theme(text = element_text(size = 30),
                axis.text.y   = element_blank(),
                axis.text.x   = element_blank(),
                axis.title.y  = element_blank(),
                axis.title.x  = element_blank(),
                axis.ticks = element_blank(),
                legend.position = "none",
                panel.background = element_blank(),
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.border = element_rect(colour = "black", fill=NA, size=1))

quant_theme = theme(text = element_text(size = 30, color='black'),
                    axis.text.x = element_text(angle = -90, hjust = 0, size = 14, color='black'),
                    legend.position = "none",
                    panel.background = element_blank(),
                    panel.grid.major = element_line(colour='grey93'), 
                    panel.grid.minor = element_line(colour='grey98'),
                    axis.title.x = element_blank(),
                    panel.border = element_rect(colour = "black", fill=NA, size=1),
                   strip.background = element_rect(fill = "white", color = "black", size = 1))


stage_labels = names(stage_colours)
names(stage_labels) = names(stage_colours)
stage_labels[10] = "Mixed"

#removed yellow from position 2 ("#FFFF00")
scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour", "Publication",
                 manual_pal(values = c(
                   "#000000", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
                   "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
                   "#5A0007", "#809693", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
                   "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
                   "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
                   "#372101", "#FFB500", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
                   "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
                   "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
                   "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
                   "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
                   "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
                   "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
                   "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C")), ...)
}

scale_color_Publication = scale_colour_Publication

#removed yellow from position 2 ("#FFFF00")
scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill", "Publication",
                 manual_pal(values = c(
                   "#000000", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
                   "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
                   "#5A0007", "#809693", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
                   "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
                   "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
                   "#372101", "#FFB500", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
                   "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
                   "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
                   "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
                   "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
                   "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
                   "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
                   "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C")), ...)
}


