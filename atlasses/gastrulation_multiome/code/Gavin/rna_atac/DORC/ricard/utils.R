ComputeDORC_ArchR <- function( 
  ArchRProject,
  expression.data,
  regions,
  gene.coords = NULL,
  distance = 5e+03,
  min.cells = 0,
  method = "pearson",
  n_sample = 50,
  pvalue_cutoff = 0.05,
  score_cutoff = 0.05,
  sep = c(":", "-"),
  verbose = TRUE) {
  
  # library(BSgenome.Mmusculus.UCSC.mm10)
  library(GenomicRanges)
  library(future)
  library(pbapply)
  library(future.apply)
  library(Matrix)
  
  
  # Fetch peak matrix
  atac.peaks.se <- getMatrixFromProject(ArchRProject, binarize = TRUE, useMatrix = "PeakMatrix")
  # peak.data <- atac.peaks.se@assays@data$PeakMatrix
  
  # Rename peaks
  rownames(atac.peaks.se) <- paste0(seqnames(rowRanges(atac.peaks.se)),":",ranges(rowRanges(atac.peaks.se)))
  names(regions) <- paste0(seqnames(regions),":",ranges(regions))
  
  common.peaks <- intersect(rownames(atac.peaks.se),names(regions))
  atac.peaks.se <- atac.peaks.se[common.peaks,]
  regions <- regions[common.peaks]
  
  # Add GC content into the SummarizedExperiment (important for calculate_background_peaks_ArchR)
  rowData(atac.peaks.se)$GC <- regions$GC
  
  # Match cells
  cells <- intersect(colnames(expression.data),colnames(atac.peaks.se))
  atac.peaks.se <- atac.peaks.se[,cells]
  expression.data <- expression.data[,cells]
  
  # Match genes
  if (is.null(names(gene.coords))) {
    names(gene.coords) <- gene.coords$gene_name
  }
  genes.to.use <- intersect(rownames(expression.data),names(gene.coords))
  expression.data <- expression.data[genes.to.use,]
  gene.coords <- gene.coords[genes.to.use]
  
  # Filter genes with low counts  
  genecounts <- rowSums(expression.data > 0)
  expression.data <- expression.data[genecounts>min.cells,]
  gene.coords.use <- gene.coords[gene.coords$gene_name %in% rownames(expression.data),]
  
  # Filter peaks with low counts  
  peakcounts <- rowSums(assay(atac.peaks.se) > 0)
  atac.peaks.se <- atac.peaks.se[peakcounts>min.cells, ]
  regions <- regions[peakcounts>min.cells]
  
  # Calculate background peaks
  background.peaks.se <- calculate_background_peaks_ArchR(
    atac_peaks.se = atac.peaks.se,
    method = "ArchR",
    nIterations = 50,
    w = 0.1,
    binSize = 50,
    seed = 1
  )
  # background.peaks.se <- readRDS(io$archR.bgdPeaks)
  # background.peaks.se <- getBgdPeaks(ArchRProject)
  # rownames(background.peaks.se) <- rowRanges(background.peaks.se) %>% as.data.table %>% 
  #   setnames("seqnames","chr") %>% .[,peak:=sprintf("%s:%s-%s",chr,start,end)] %>% .$peak
  # background.peaks.se <- background.peaks.se[common.peaks,]
  
  if (verbose) sprintf("Testing %d genes and %d peaks",nrow(expression.data),length(regions))
  
  # Create peak2gene binary matrix based on genomic distance
  peak2gene_binary.mtx <- DistanceToTSS(
    peaks = regions,
    genes = gene.coords,
    distance = distance
  )
  
  genes.use <- colnames(peak2gene_binary.mtx)
  all.peaks <- rownames(peak2gene_binary.mtx)
  
  peak.data <- t(assay(atac.peaks.se))
  
  coef.vec <- c()
  gene.vec <- c()
  zscore.vec <- c()
  if (nbrOfWorkers() > 1) {
    mylapply <- future_lapply
  } else {
    mylapply <- ifelse(test = verbose, yes = pblapply, no = lapply)
  }
  
  # run in parallel across genes
  res <- mylapply(
    X = seq_along(genes.use),
    FUN = function(i) {
      print(i)
      peak.use <- as.logical(peak2gene_binary.mtx[, genes.use[[i]]])
      gene.expression <- expression.data[genes.use[[i]], ]
      gene.chrom <- as.character(seqnames(gene.coords.use[i]))
      
      if (sum(peak.use) < 2) {
        # no peaks close to gene
        return(list("gene" = NULL, "coef" = NULL, "zscore" = NULL))
      } else {
        coef.result <- cor(
          x = t(as.matrix(assay(atac.peaks.se[peak.use,]))),
          y = as.matrix(gene.expression),
          method = method
        )
        coef.result <- coef.result[coef.result>score_cutoff, , drop = FALSE]
        
        if (nrow(coef.result) == 0) {
          return(list("gene" = NULL, "coef" = NULL, "zscore" = NULL))
        } else {
          
          # select peaks at random with matching GC content and accessibility
          # sample from peaks on a different chromosome to the gene
          peaks.test <- rownames(coef.result)
          trans.peaks <- all.peaks[!grepl(pattern = paste0("^", gene.chrom), x = all.peaks)]
          
          bg.peaks <- lapply(peaks.test,function(x) {
            rownames(background.peaks.se)[assay(background.peaks.se[x,])[1,]]
          }) %>% unlist
          # bg.peaks <- bg.peaks[!is.na(bg.peaks)]
          
          # run background correlations
          bg.coef <- cor(
            x = t(as.matrix(assay(atac.peaks.se[bg.peaks,]))),
            y = as.matrix(gene.expression),
            method = method
          )
          # bg.coef[is.na(bg.coef)] <- 0
          
          zscores <- vector(mode = "numeric", length = length(peaks.test))
          for (j in seq_along(along.with = peaks.test)) {
            coef.use <- bg.coef[(((j - 1) * n_sample) + 1):(j * n_sample), ]
            z <- (coef.result[j] - mean(coef.use)) / sd(coef.use)
            zscores[[j]] <- z
          }
          names(coef.result) <- peaks.test
          names(zscores) <- peaks.test
          zscore.vec <- c(zscore.vec, zscores)
          gene.vec <- c(gene.vec, rep(i, length(coef.result)))
          coef.vec <- c(coef.vec, coef.result)
        }
        gc(verbose = FALSE)
        pval.vec <- pnorm(q = -abs(zscore.vec))
        links.keep <- pval.vec < pvalue_cutoff
        if (sum(links.keep) == 0) {
          return(list("gene" = NULL, "coef" = NULL, "zscore" = NULL))
        } else {
          gene.vec <- gene.vec[links.keep]
          coef.vec <- coef.vec[links.keep]
          zscore.vec <- zscore.vec[links.keep]
          return(list("gene" = gene.vec, "coef" = coef.vec, "zscore" = zscore.vec))
        }
      }
    }
  )
  # combine results
  if(length(res)>1){
    gene.vec <- do.call(what = c, args = sapply(X = res, FUN = `[[`, 1))
    coef.vec <- do.call(what = c, args = sapply(X = res, FUN = `[[`, 2))
    zscore.vec <- do.call(what = c, args = sapply(X = res, FUN = `[[`, 3))
  }else{
    gene.vec=res[[1]][["gene"]]
    coef.vec=res[[1]][["coef"]]
    zscore.vec=res[[1]][["zscore"]]
  }
  
  if (length(coef.vec) == 0) {
    if (verbose) {
      message("No significant links found")
    }
    return(object)
  }
  peak.key <- seq_along(
    along.with = unique(names(coef.vec))
  )
  names(peak.key) <- unique(names(coef.vec))
  coef.matrix <- sparseMatrix(
    i = gene.vec,
    j = peak.key[names(coef.vec)],
    x = coef.vec,
    dims = c(length(genes.use), max(peak.key))
  )
  rownames(coef.matrix) <- genes.use
  colnames(coef.matrix) <- names(peak.key)
  links <- LinksToGRanges(linkmat = coef.matrix, gene.coords = gene.coords.use)
  # add zscores
  z.matrix <- sparseMatrix(
    i = gene.vec,
    j = peak.key[names(zscore.vec)],
    x = zscore.vec,
    dims = c(length(genes.use), max(peak.key))
  )
  rownames(z.matrix) <- genes.use
  colnames(z.matrix) <- names(peak.key)
  z.lnk <- LinksToGRanges(linkmat = z.matrix, gene.coords = gene.coords.use,sep=sep)
  links$zscore <- z.lnk$score
  links$pvalue <- pnorm(q = -abs(links$zscore))
  links <- links[links$pvalue < pvalue_cutoff]
  
  return(links)
}


#' Title
#'
#' @param ranges 
#'
#' @return
#' @export
#'
#' @examples
CollapseToLongestTranscript <- function(ranges) {
  library(data.table)
  range.df <- as.data.table(x = ranges)
  range.df$strand <- as.character(x = range.df$strand)
  range.df$strand <- ifelse(
    test = range.df$strand == "*",
    yes = "+",
    no = range.df$strand
  )
  collapsed <- range.df[
    , .(unique(seqnames),
        min(start),
        max(end),
        strand[[1]],
        gene_biotype[[1]]),
    "gene_name"
  ]
  colnames(x = collapsed) <- c(
    "gene_name", "seqnames", "start", "end", "strand", "gene_biotype"
  )
  gene.ranges <- makeGRangesFromDataFrame(
    df = collapsed,
    keep.extra.columns = TRUE
  )
  return(gene.ranges)
}

#' Title
#'
#' @param peaks 
#' @param genes 
#' @param distance 
#' @param sep 
#'
#' @return
#' @export
#'
#' @examples
DistanceToTSS <- function(peaks, genes, distance = 200000, sep = c("-", "-")) {
  library(GenomicRanges)
  library(BiocGenerics)
  tss <- resize(x = genes, width = 1, fix = 'start')
  genes.extended <- suppressWarnings(
    expr = Extend(
      x = tss, upstream = distance, downstream = distance
    )
  )
  overlaps <- findOverlaps(
    query = peaks,
    subject = genes.extended,
    type = 'any',
    select = 'all'
  )
  hit_matrix <- sparseMatrix(
    i = queryHits(x = overlaps),
    j = subjectHits(x = overlaps),
    x = 1,
    dims = c(length(x = peaks), length(x = genes.extended))
  )
  rownames(x = hit_matrix) <- GRangesToString(grange = peaks, sep = sep)
  colnames(x = hit_matrix) <- genes.extended$gene_name
  return(hit_matrix)
}
library(Matrix)

#' Title
#'
#' @param linkmat 
#' @param gene.coords 
#' @param sep 
#'
#' @return
#' @export
#'
#' @examples
LinksToGRanges <- function(linkmat, gene.coords, sep = c(":", "-")) {
  # get TSS for each gene
  tss <- resize(gene.coords, width = 1, fix = 'start')
  gene.idx <- sapply(
    X = rownames(x = linkmat),
    FUN = function(x) {
      which(x = x == tss$gene_name)[[1]]
    }
  )
  tss <- tss[gene.idx]
  
  # get midpoint of each peak
  peak.ranges <- StringToGRanges(
    regions = colnames(x = linkmat),
    sep = sep
  )
  midpoints <- start(x = peak.ranges) + (width(x = peak.ranges) / 2)
  
  # convert to triplet form
  dgtm <- as(object = linkmat, Class = "dgTMatrix")
  
  # create dataframe
  df <- data.frame(
    chromosome = as.character(x = seqnames(x = peak.ranges)[dgtm@j + 1]),
    tss = start(x = tss)[dgtm@i + 1],
    pk = midpoints[dgtm@j + 1],
    score = dgtm@x,
    gene = rownames(x = linkmat)[dgtm@i + 1],
    peak = colnames(x = linkmat)[dgtm@j + 1]
  )
  
  # work out start and end coords
  df$start <- ifelse(test = df$tss < df$pk, yes = df$tss, no = df$pk)
  df$end <- ifelse(test = df$tss < df$pk, yes = df$pk, no = df$tss)
  df$tss <- NULL
  df$pk <- NULL
  
  # convert to granges
  gr.use <- makeGRangesFromDataFrame(df = df, keep.extra.columns = TRUE)
  return(sort(x = gr.use))
}
#' Title
#'
#' @param ensdb 
#' @param standard.chromosomes 
#' @param biotypes 
#' @param verbose 
#'
#' @return
#' @export
#'
#' @examples
GetGRangesFromEnsDb <- function(
  ensdb,
  standard.chromosomes = TRUE,
  biotypes = c("protein_coding", "lincRNA", "rRNA", "processed_transcript"),
  verbose = TRUE
) {
  library(biovizBase)
  # convert seqinfo to granges
  whole.genome <-  as(object = seqinfo(x = ensdb), Class = "GRanges")
  whole.genome <- keepStandardChromosomes(whole.genome, pruning.mode = "coarse")
  
  # extract genes from each chromosome
  if (verbose) {
    tx <- sapply(X = seq_along(whole.genome), FUN = function(x){
      crunch(
        obj = ensdb,
        which = whole.genome[x],
        columns = c("tx_id", "gene_name", "gene_id", "gene_biotype"))
    })
  } else {
    tx <- sapply(X = seq_along(whole.genome), FUN = function(x){
      suppressMessages(expr = crunch(
        obj = ensdb,
        which = whole.genome[x],
        columns = c("tx_id", "gene_name", "gene_id", "gene_biotype")))
    })
  }
  # combine
  tx <- do.call(what = c, args = tx)
  tx <- tx[tx$gene_biotype %in% biotypes]
  return(tx)}

#' Title
#'
#' @param x 
#' @param upstream 
#' @param downstream 
#' @param from.midpoint 
#'
#' @return
#' @export
#'
#' @examples
Extend <- function(
  x,
  upstream = 0,
  downstream = 0,
  from.midpoint = FALSE
) {
  if (any(strand(x = x) == "*")) {
    warning("'*' ranges were treated as '+'")
  }
  on_plus <- strand(x = x) == "+" | strand(x = x) == "*"
  if (from.midpoint) {
    midpoints <- start(x = x) + (width(x = x) / 2)
    new_start <- midpoints - ifelse(
      test = on_plus, yes = upstream, no = downstream
    )
    new_end <- midpoints + ifelse(
      test = on_plus, yes = downstream, no = upstream
    )
  } else {
    new_start <- start(x = x) - ifelse(
      test = on_plus, yes = upstream, no = downstream
    )
    new_end <- end(x = x) + ifelse(
      test = on_plus, yes = downstream, no = upstream
    )
  }
  ranges(x = x) <- IRanges(start = new_start, end = new_end)
  # x <- trim(x = x)
  return(x)
}
#' Title
#'
#' @param grange 
#' @param sep 
#'
#' @return
#' @export
#'
#' @examples
GRangesToString <- function(grange, sep = c("-", "-")) {
  regions <- paste0(
    as.character(x = seqnames(x = grange)),
    sep[[1]],
    start(x = grange),
    sep[[2]],
    end(x = grange)
  )
  return(regions)
}
#' Title
#'
#' @param regions 
#' @param sep 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
StringToGRanges <- function(regions, sep = c(":", "-"), ...) {
  ranges.df <- data.frame(ranges = regions)
  library(tidyr)
  ranges.df <- separate(
    data = ranges.df,
    col = "ranges",
    sep = paste0(sep[[1]], "|", sep[[2]]),
    into = c("chr", "start", "end")
  )
  granges <- makeGRangesFromDataFrame(df = ranges.df, ...)
  return(granges)
}


calculate_background_peaks_ArchR <- function(
  atac_peaks.se,
  method = "ArchR",
  nIterations = 50,
  w = 0.1,
  binSize = 50,
  seed = 1
) {
  
  # Sanity checks
  stopifnot("GC" %in% colnames(rowData(atac_peaks.se)))
  
  # Create SummarizedExperiment
  se <- SummarizedExperiment(
    assays = SimpleList(counts = as.matrix(rowSums(assay(atac_peaks.se)))),
    rowData = DataFrame(
      bias = rowData(atac_peaks.se)$GC, 
      seqnames = seqnames(rowRanges(atac_peaks.se)), 
      start = start(rowRanges(atac_peaks.se)), 
      end = end(rowRanges(atac_peaks.se))
    )
  )
  
  # Run BGdPeaks algorithm
  if (method=="ArchR") {
    tmp <- .ArchRBdgPeaks(
      object = se,
      bias = rowData(se)$bias, 
      nIterations = nIterations
    )
  } else if (method=="chromVAR") {
    tmp <- chromVAR::getBackgroundPeaks(
      object = se,
      bias = rowData(se)$bias,
      niterations = nIterations,
      w = w,
      bs = binSize
    )
  } else {
    stop("method not recognised")
  }
  
  
  # Prepare output
  bgdPeaks <- SummarizedExperiment(
    assays = SimpleList(bgdPeaks = tmp), 
    rowRanges = GRanges(rowData(se)$seqnames, IRanges(rowData(se)$start,rowData(se)$end), value=assay(se)[,1], GC=rowData(se)$bias)
  )
  rownames(bgdPeaks) <- rownames(se)
  
  biasDF <- data.frame(
    rowSums = Matrix::rowSums(assay(se)),
    bias = rowData(se)$bias,
    length = rowData(se)$end - rowData(se)$start
  )
  
  rowData(bgdPeaks)$bgdSumMean <- round(rowMeans(matrix(biasDF[assay(bgdPeaks),1], nrow = nrow(bgdPeaks))),3)
  rowData(bgdPeaks)$bgdSumSd <- round(matrixStats::rowSds(matrix(biasDF[assay(bgdPeaks),1], nrow = nrow(bgdPeaks))),3)
  
  rowData(bgdPeaks)$bgdGCMean <- round(rowMeans(matrix(biasDF[assay(bgdPeaks),2], nrow = nrow(bgdPeaks))),3)
  rowData(bgdPeaks)$bgdGCSd <- round(matrixStats::rowSds(matrix(biasDF[assay(bgdPeaks),2], nrow = nrow(bgdPeaks))),3)
  
  rowData(bgdPeaks)$bgdLengthMean <- round(rowMeans(matrix(biasDF[assay(bgdPeaks),3], nrow = nrow(bgdPeaks))),3)
  rowData(bgdPeaks)$bgdLengthSd <- round(matrixStats::rowSds(matrix(biasDF[assay(bgdPeaks),3], nrow = nrow(bgdPeaks))),3)
  
  return(bgdPeaks)
}  


.ArchRBdgPeaks <- function(object = NULL, bias = NULL, nIterations = 50){
  
  .cleanSelf <- function(x){
    xn <- matrix(0, nrow = nrow(x), ncol = ncol(x))
    for(i in seq_len(nrow(x))){
      xi <- x[i, ]
      idx <- which(xi != i)
      xn[i, seq_along(idx)] <- xi[idx]
    }
    idx <- which(colSums(xn == 0) > 0)
    if(length(idx) > 0){
      xn <- xn[,-idx]
    }
    xn
  }
  
  #Bias Dataframe
  biasDF <- data.frame(
    rowSums = Matrix::rowSums(assay(object)),
    bias = bias,
    length = rowData(object)$end - rowData(object)$start
  )
  
  #Quantile Normalize
  biasDFN <- apply(biasDF, 2, .getQuantiles)
  
  #Get KNN
  knnObj <- nabor::knn(
    data =  biasDFN,
    k = nIterations + 1
  )[[1]]
  
  #Filter Self
  knnObj <- .cleanSelf(knnObj)
  
  #Shuffle
  idx <- seq_len(ncol(knnObj))
  knnObj2 <- matrix(0, nrow = nrow(knnObj), ncol = ncol(knnObj))
  for(x in seq_len(nrow(knnObj2))){
    knnObj2[x,] <- knnObj[x, sample(idx, length(idx))]
  }
  
  knnObj2
  
}



.getQuantiles <- function(v = NULL, len = length(v)){
  if(length(v) < len){
    v2 <- rep(0, len)
    v2[seq_along(v)] <- v
  }else{
    v2 <- v
  }
  p <- trunc(rank(v2))/length(v2)
  if(length(v) < len){
    p <- p[seq_along(v)]
  }
  return(p)
}
