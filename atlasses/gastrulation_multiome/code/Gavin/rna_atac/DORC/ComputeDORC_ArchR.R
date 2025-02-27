
#' Title
#'
#' @param ArchRProject
#' @param expression.data 
#' @param regions a GRanges regions of peaks
#' @param genome the reference genome
#' @param annotation the genome annotation
#' @param gene.coords 
#' @param distance 
#' @param min.cells 
#' @param method 
#' @param genes.use gene to compute
#' @param n_sample 
#' @param pvalue_cutoff 
#' @param score_cutoff 
#' @param verbose 
#'
#' @return
#' @export
#'
#' @examples
#' annotation <- GetGRangesFromEnsDb(EnsDb.Mmusculus.v79)
#' seqlevelsStyle(annotation) <- "UCSC"
#' ComputeDORC(peak.data = peak.data,expression.data = expression.data,genome = genome,annotation=annotation)


ComputeDORC_ArchR=function( 
  ArchRProject,
  expression.data,
  regions,
  genome,
  annotation,
  gene.coords = NULL,
  distance = 5e+03,
  min.cells = 0,
  method = "pearson",
  genes.use = NULL,
  n_sample = 50,
  pvalue_cutoff =0.05,
  score_cutoff = 0.05,
  sep=c(":", "-"),
  verbose = TRUE){
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(EnsDb.Mmusculus.v79)
  library(future)
  library(pbapply)
  library(GenomicRanges)
  library(future.apply)
  background.peaks.se<-getBgdPeaks(ArchRProject)
  rownames(background.peaks.se) <- rowRanges(background.peaks.se) %>% as.data.table %>% 
    setnames("seqnames","chr") %>% .[,peak:=sprintf("%s:%s-%s",chr,start,end)] %>% .$peak
  atac.peaks.se=getMatrixFromProject(ArchRProject, binarize = TRUE, useMatrix = "PeakMatrix")[,cells.use]
  peak.data<-atac.peaks.se@assays@data$PeakMatrix
  # order the cells
  rowrange=atac.peaks.se@rowRanges
  rownames(peak.data)<-paste0(rowrange@seqnames,':',rowrange@ranges)
  if (!setequal(colnames(expression.data),peak.data@Dimnames[[2]])){stop("the cells are not the same")}
  target=colnames(expression.data)
  peak.data=peak.data[,match(target, peak.data@Dimnames[[2]])]
  if (is.null(x = gene.coords)) {
    gene.coords <- CollapseToLongestTranscript(annotation)
  }
  if(missing(regions)){
    peak_bed = do.call(rbind, strsplit(x = rownames(peak.data), split = '[:-]'))
    peak_bed = as.data.frame(peak_bed)
    colnames(peak_bed)=c('chr','start','end')
    regions=makeGRangesFromDataFrame(peak_bed)
  }
  # meta.features <- ComputeRegionStats(regions=regions,genome = genome)
  # rownames(meta.features)=rownames(peak.data)
  # meta.features<-as.data.frame(meta.features)
  # features.match <- c("GC.percent", "count")
  # if (!("GC.percent" %in% colnames(x = meta.features))) {
  #   stop("GC content per peak has not been computed.\n",
  #        "Run RegionsStats before calling this function.")
  # }
  
  peakcounts <- rowSums(x = peak.data > 0)
  meta.features$count=peakcounts
  genecounts <- rowSums(x = expression.data > 0)
  peaks.keep <- peakcounts > min.cells
  genes.keep <- genecounts > min.cells
  peak.data <- peak.data[peaks.keep, ]
  if(is.null(genes.use)){genes.use=sample(genes.keep[genes.keep],2)}
  genes.keep <- intersect( x = names(x = genes.keep[genes.keep]), y = genes.use)
  if(length(genes.keep)>1){
    expression.data=expression.data[genes.keep,]
  }else{
    cells=colnames(expression.data)
    expression.data=as(t(expression.data[genes.keep, ]),'Matrix') 
    rownames(expression.data)=genes.keep
    colnames(expression.data)=cells
    
  }
  
  
  if (verbose) {
    message(
      "Testing ",
      nrow(x = expression.data),
      " genes and ",
      sum(peaks.keep),
      " peaks"
    )
  }
  genes <- rownames(x = expression.data)
  gene.coords.use <- gene.coords[gene.coords$gene_name %in% genes,]
  peaks <- regions
  peaks <- peaks[peaks.keep]
  library(Matrix)
  peak_distance_matrix <- DistanceToTSS(
    peaks = peaks,
    genes = gene.coords.use,
    distance = distance
  )
  
  genes.use <- colnames(x = peak_distance_matrix)
  all.peaks <- rownames(x = peak.data)
  
  peak.data <- t(x = peak.data)
  
  coef.vec <- c()
  gene.vec <- c()
  zscore.vec <- c()
  library(future)
  library(pbapply)
  if (nbrOfWorkers() > 1) {
    mylapply <- future_lapply
  } else {
    mylapply <- ifelse(test = verbose, yes = pblapply, no = lapply)
  }
  
  # run in parallel across genes
  res <- mylapply(
    X = seq_along(along.with = genes.use),
    FUN = function(i) {
      peak.use <- as.logical(x = peak_distance_matrix[, genes.use[[i]]])
      gene.expression <- expression.data[genes.use[[i]], ]
      gene.chrom <- as.character(x = seqnames(x = gene.coords.use[i]))
      
      if (sum(peak.use) < 2) {
        # no peaks close to gene
        return(list("gene" = NULL, "coef" = NULL, "zscore" = NULL))
      } else {
        peak.access <- peak.data[, peak.use]
        coef.result <- cor(
          x = as.matrix(x = peak.access),
          y = as.matrix(x = gene.expression),
          method = method
        )
        coef.result <- coef.result[x = coef.result > score_cutoff, , drop = FALSE]
        
        if (nrow(x = coef.result) == 0) {
          return(list("gene" = NULL, "coef" = NULL, "zscore" = NULL))
        } else {
          
          # select peaks at random with matching GC content and accessibility
          # sample from peaks on a different chromosome to the gene
          peaks.test <- rownames(x = coef.result)
          trans.peaks <- all.peaks[
            !grepl(pattern = paste0("^", gene.chrom), x = all.peaks)
          ]
          # meta.use <- meta.features[trans.peaks, ]
          # pk.use <- meta.features[peaks.test, ]
          # bg.peaks <- lapply(
          #   X = seq_len(length.out = nrow(x = pk.use)),
          #   FUN = function(x) {
          #     MatchRegionStats(
          #       meta.feature = meta.use,
          #       query.feature = pk.use[x, , drop = FALSE],
          #       features.match = c("GC.percent", "count", "sequence.length"),
          #       n = n_sample,
          #       verbose = FALSE
          #     )
          #   }
          # )
          
          bg.peaks<- lapply(peaks.test,function(x) {rownames(background.peaks.se)[assay(background.peaks.se[x,])[1,]]})
          # run background correlations
          bg.access <- peak.data[, unlist(x = bg.peaks)]
          bg.coef <- cor(
            x = as.matrix(x = bg.access),
            y = as.matrix(x = gene.expression),
            method = method
          )
          zscores <- vector(mode = "numeric", length = length(x = peaks.test))
          for (j in seq_along(along.with = peaks.test)) {
            coef.use <- bg.coef[(((j - 1) * n_sample) + 1):(j * n_sample), ]
            z <- (coef.result[j] - mean(x = coef.use)) / sd(x = coef.use)
            zscores[[j]] <- z
          }
          names(x = coef.result) <- peaks.test
          names(x = zscores) <- peaks.test
          zscore.vec <- c(zscore.vec, zscores)
          gene.vec <- c(gene.vec, rep(i, length(x = coef.result)))
          coef.vec <- c(coef.vec, coef.result)
        }
        gc(verbose = FALSE)
        pval.vec <- pnorm(q = -abs(x = zscore.vec))
        links.keep <- pval.vec < pvalue_cutoff
        if (sum(x = links.keep) == 0) {
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
  
  if (length(x = coef.vec) == 0) {
    if (verbose) {
      message("No significant links found")
    }
    return(object)
  }
  peak.key <- seq_along(
    along.with = unique(x = names(x = coef.vec))
  )
  names(x = peak.key) <- unique(x = names(x = coef.vec))
  coef.matrix <- sparseMatrix(
    i = gene.vec,
    j = peak.key[names(x = coef.vec)],
    x = coef.vec,
    dims = c(length(x = genes.use), max(peak.key))
  )
  rownames(x = coef.matrix) <- genes.use
  colnames(x = coef.matrix) <- names(x = peak.key)
  links <- LinksToGRanges(linkmat = coef.matrix, gene.coords = gene.coords.use)
  # add zscores
  z.matrix <- sparseMatrix(
    i = gene.vec,
    j = peak.key[names(x = zscore.vec)],
    x = zscore.vec,
    dims = c(length(x = genes.use), max(peak.key))
  )
  rownames(x = z.matrix) <- genes.use
  colnames(x = z.matrix) <- names(x = peak.key)
  z.lnk <- LinksToGRanges(linkmat = z.matrix, gene.coords = gene.coords.use,sep=sep)
  links$zscore <- z.lnk$score
  links$pvalue <- pnorm(q = -abs(x = links$zscore))
  links <- links[links$pvalue < pvalue_cutoff]
  return(links)
}
#' Title
#'
#' @param regions 
#' @param genome 
#'
#' @return
#' @export
#'
#' @examples

library(EnsDb.Mmusculus.v79)



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

