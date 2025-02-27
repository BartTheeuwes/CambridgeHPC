#############################
## Commonly-used functions ##
#############################

load_SingleCellExperiment <- function(file, normalise = FALSE, features = NULL, cells = NULL, remove_non_expressed_genes = FALSE) {
  library(SingleCellExperiment); library(scran); library(scater);
  sce <- readRDS(file)
  if (!is.null(cells)) sce <- sce[,cells]
  if (!is.null(features)) sce <- sce[features,]
  if (remove_non_expressed_genes) sce <- sce[which(Matrix::rowSums(counts(sce))>15),]
  if (normalise) sce <- logNormCounts(sce)
  return(sce)
}

load_Seurat <- function(file, assay = "RNA", normalise = FALSE, features = NULL, cells = NULL, remove_non_expressed_genes = FALSE, ...) {
  library(Seurat)
  seurat <- readRDS(file)
  # if (assay%in%Seurat::Assays(seurat)) seurat <- seurat[[assay]]
  if (!is.null(cells)) seurat <- seurat[,cells]
  if (!is.null(features)) seurat <- seurat[features,]
  if (normalise) {
    seurat <- NormalizeData(seurat, normalization.method = "LogNormalize")
    seurat <- ScaleData(seurat, ...)
  }
  if (remove_non_expressed_genes) seurat <- seurat[which(Matrix::rowMeans(seurat@assays[[assay]]@counts)>1e-4),]
  return(seurat)
}

matrix.please<-function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[[1]]
  m
}

#'  method=1: The TF-IDF implementation used by Stuart & Butler et al. 2019. This computes \eqn{\log(TF \times IDF)}.
#'  method=2: The TF-IDF implementation used by Cusanovich & Hill et al. 2018. This computes \eqn{TF \times (\log(IDF))}.
#'  method=3: The log-TF method used by Andrew Hill. This computes \eqn{\log(TF) \times \log(IDF)}.
tfidf <- function(mtx, method = 1, scale.factor = 1e4) {
  npeaks <- colSums(mtx)
  if (any(npeaks == 0)) {
    warning("Some cells contain 0 total counts")
  }

  tf <- Matrix::tcrossprod(mtx, y = Matrix::Diagonal(x=1/npeaks))

  rsums <- rowSums(mtx)
  if (any(rsums == 0)) {
    warning("Some features contain 0 total counts")
  }
  idf <- ncol(mtx) / rsums

  if (method == 2) {
    idf <- log(1 + idf)
  } else if (method == 3) {
    tf <- log1p(tf * scale.factor)
    idf <- log(1 + idf)
  }
  mtx.tfidf <- Matrix::Diagonal(n = length(idf), x = idf) %*% tf

  if (method == 1) {
    mtx.tfidf <- log1p(mtx.tfidf * scale.factor)
  }
  colnames(mtx.tfidf) <- colnames(mtx)
  rownames(mtx.tfidf) <- rownames(mtx)

  # set NA values to 0
  mtx.tfidf[is.na(mtx.tfidf)] <- 0

  return(mtx.tfidf)
}

pdist <- function(tmat){
  # @param tmat A non-negative matrix with samples by features
  # @reference http://r.789695.n4.nabble.com/dist-function-in-R-is-very-slow-td4738317.html
  mtm <- Matrix::tcrossprod(tmat)
  sq <- rowSums(tmat^2)
  out0 <- outer(sq, sq, "+") - 2 * mtm
  out0[out0 < 0] <- 0
  
  sqrt(out0)
}

smoother_aggregate_nearest_nb <- function(mat, D, k){
  # @param mat A matrix in a shape of #genes x #samples.
  # @param D A predefined distance matrix in a shape of #samples x #samples.
  # @param k An integer to choose \code{k} nearest samples (self-inclusive) to
  #  aggregate based on the distance matrix \code{D}.
  denoised_mat <- sapply(seq_len(ncol(mat)), function(cid){
    nb_cid <- head(order(D[cid, ]), k)
    closest_mat <- mat[, nb_cid, drop=FALSE]
    # return(Matrix::rowSums(closest_mat))
    return(Matrix::rowMeans(closest_mat))
  })
  dimnames(denoised_mat) <- dimnames(mat)
  return(denoised_mat)
}


# TO-FINISH.....
smoother_aggregate_nearest_nb_parallel <- function(mat, D, k, cores=1){
  # @param mat A matrix in a shape of #genes x #samples.
  # @param D A predefined distance matrix in a shape of #samples x #samples.
  # @param k An integer to choose \code{k} nearest samples (self-inclusive) to
  #  aggregate based on the distance matrix \code{D}.

  # library(future)
  library(future.apply)
  plan("multiprocess", workers = ncores)


  # sapply(seq_len(ncol(mat)), function(cid){
  future_sapply(seq_len(ncol(mat)), function(cid){
    nb_cid <- head(order(D[cid, ]), k)
    closest_mat <- mat[, nb_cid, drop=FALSE]
    # return(Matrix::rowSums(closest_mat))
    return(Matrix::rowMeans(closest_mat))
  })
}

# regress_covariates <- function(mtx, vars.to.regress) {
#   data <- scale(t(logcounts(sce_filt)), center = T, scale = F)
#   data_regressed <- apply(data, 2, function(x) {
#     lm.out <- lm(formula=expr~covariate, data=data.frame(expr=x, covariate=factor(sce_filt$stage)));
#     residuals <- lm.out[["residuals"]]+lm.out[["coefficients"]][1]
#   })
# }

# Remove unwanted effects from a matrix
#
# @parm mtx An expression matrix to regress the effects of covariates out
# of should be the complete expression matrix in genes x cells
# @param covariates A matrix or data.frame of latent variables, should be cells
# x covariates, the colnames should be the variables to regress
# @param features_idx An integer vector representing the indices of the
# genes to run regression on
# @param model.use Model to use, one of 'linear', 'poisson', or 'negbinom'; pass
# NULL to simply return mtx
# @param verbose Display a progress bar
#' @importFrom stats as.formula lm
#' @importFrom utils txtProgressBar setTxtProgressBar
#
RegressOutMatrix_parallel <- function(mtx, covariates = NULL, features_idx = NULL, split.by = NULL, block.size = 1000, min.cells.to.block = 3000, ncores = 1, verbose = TRUE) {
  
  library(future)
  library(future.apply)
  plan("multiprocess", workers = ncores)
  
  # Check features_idx
  if (is.null(features_idx)) {
    features_idx <- 1:nrow(mtx)
  }
  if (is.character(features_idx)) {
    features_idx <- intersect(features_idx, rownames(mtx))
    if (length(features_idx) == 0) {
      stop("Cannot use features that are beyond the scope of mtx")
    }
  } else if (max(features_idx) > nrow(mtx)) {
    stop("Cannot use features that are beyond the scope of mtx")
  }
  
  # Check dataset dimensions
  if (nrow(covariates) != ncol(mtx)) {
    stop("Uneven number of cells between latent data and expression data")
  }
  
  # Subset
  mtx <- mtx[features_idx,]
  mtx.dimnames <- dimnames(mtx)
  
  # Define chunck points
  chunk.points <- ChunkPoints(dsize = nrow(mtx), csize = block.size)
  
  # Define cell splitting
  split.cells <- split(colnames(mtx), f = split.by %||% TRUE)

   if (nbrOfWorkers() > 1) {

    # Define chuncks
      chunks <- expand.grid(
        names(split.cells),
        1:ncol(chunk.points),
        stringsAsFactors = FALSE
      )

      # Run RegressOutMatrix in parallel
      mtx.resid <- future_lapply(
        X = 1:nrow(chunks),
        FUN = function(i) {
          row <- chunks[i, ]
          group <- row[[1]]
          index <- as.numeric(row[[2]])
          return(RegressOutMatrix(
            mtx = mtx[chunk.points[1, index]:chunk.points[2, index], split.cells[[group]], drop = FALSE],
            covariates = covariates[split.cells[[group]], , drop = FALSE],
            # features_idx = features_idx[chunk.points[1, index]:chunk.points[2, index]],
            verbose = FALSE
          ))
        }
      )

      # Merge splitted cells
      if (length(split.cells) > 1) {
        merge.indices <- lapply(
          X = 1:length(x = split.cells),
          FUN = seq.int,
          to = length(mtx.resid),
          by = length(split.cells)
        )
        mtx.resid <- lapply(
          X = merge.indices,
          FUN = function(x) {
            return(do.call( 'rbind', mtx.resid[x]))
          }
        )
        mtx.resid <- do.call('cbind', mtx.resid)
      } else {
        mtx.resid <- do.call( 'rbind', mtx.resid)
      }
    } else {
      
      mtx.resid <- lapply(
        X = names(split.cells),
        FUN = function(x) {
          if (verbose && length(split.cells) > 1) {
            message("Regressing out variables from split ", x)
          }
          return(RegressOutMatrix(
            mtx = mtx[, split.cells[[x]], drop = FALSE],
            covariates = covariates[split.cells[[x]], , drop = FALSE],
            features_idx = features_idx,
            verbose = verbose
          ))
        }
      )
      mtx.resid <- do.call('cbind', mtx.resid)
    }
    # dimnames(mtx.resid) <- dimnames(mtx)
    return(mtx.resid)
  }

RegressOutMatrix <- function(mtx, covariates = NULL, features_idx = NULL, verbose = TRUE) {
  
  # Check features_idx
  if (is.null(features_idx)) {
    features_idx <- 1:nrow(mtx)
  }
  if (is.character(features_idx)) {
    features_idx <- intersect(features_idx, rownames(mtx))
    if (length(features_idx) == 0) {
      stop("Cannot use features that are beyond the scope of mtx")
    }
  } else if (max(features_idx) > nrow(mtx)) {
    stop("Cannot use features that are beyond the scope of mtx")
  }
  
  # Check dataset dimensions
  if (nrow(covariates) != ncol(mtx)) {
    stop("Uneven number of cells between latent data and expression data")
  }
  
  # Subset
  mtx <- mtx[features_idx,]
  mtx.dimnames <- dimnames(mtx)
  
  # Create formula for regression
  vars.to.regress <- colnames(covariates)
  fmla <- paste('GENE ~', paste(vars.to.regress, collapse = '+')) %>% as.formula

  # In this code, we'll repeatedly regress different Y against the same X
  # (covariates) in order to calculate residuals.  Rather that repeatedly
  # call lm to do this, we'll avoid recalculating the QR decomposition for the
  # covariates matrix each time by reusing it after calculating it once
  regression.mat <- cbind(covariates, mtx[1,])
  colnames(regression.mat) <- c(colnames(covariates), "GENE")
  qr <- lm(fmla, data = regression.mat, qr = TRUE)$qr
  rm(regression.mat)

  # Make results matrix
  data.resid <- matrix(
    nrow = nrow(mtx),
    ncol = ncol(mtx)
  )

  if (verbose) pb <- txtProgressBar(char = '=', style = 3, file = stderr())

  # Extract residuals from each feature by using the pre-computed QR decomposition
  for (i in 1:length(features_idx)) {
    regression.mat <- cbind(covariates, mtx[features_idx[i], ])
    colnames(regression.mat) <- c(vars.to.regress, 'GENE')
    regression.mat <- qr.resid(qr = qr, y = mtx[features_idx[i],])  # The function qr.resid returns the residuals when fitting y to the matrix with QR decomposition.
    data.resid[i, ] <- regression.mat
    if (verbose) {
      setTxtProgressBar(pb = pb, value = i / length(features_idx))
    }
  }

  if (verbose) close(con = pb)

  dimnames(data.resid) <- mtx.dimnames
  
  return(data.resid)
}



# Generate chunk points
#
# @param dsize How big is the data being chunked
# @param csize How big should each chunk be
#
# @return A matrix where each column is a chunk, row 1 is start points, row 2 is end points
#
ChunkPoints <- function(dsize, csize) {
  return(vapply(
    X = 1L:ceiling(dsize / csize),
    FUN = function(i) {
      return(c(
        start = (csize * (i - 1L)) + 1L,
        end = min(csize * i, dsize)
      ))
    },
    FUN.VALUE = numeric(length = 2L)
  ))
}


"%ni%" <- Negate("%in%")

ggplot_theme_NoAxes <- function() {
  theme(
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
}

minmax.normalisation <- function(x)
{
    return((x-min(x,na.rm=T)) /(max(x,na.rm=T)-min(x,na.rm=T)))
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


## sparse -> matrix
dropNA2matrix <- function(x) {
  if(!is(x, "dsparseMatrix")) stop("x needs to be a subclass of dsparseMatrix!")
  
  x <- as(x, "dgCMatrix")
  
  ## remember true NAs
  nas <- Matrix::which(is.na(x), arr.ind=TRUE)
  
  x@x[x@x==0] <- NA
  zeros <- Matrix::which(is.na(x), arr.ind=TRUE)
  x <- as(x, "matrix")
  x[x==0] <- NA
  x[zeros] <- 0
  x[nas] <- NA
  x
}

# dropNA2vector <- function(x) {
#   stop("NEEDS TO BE FIXED")
#   if(!is(x, "dsparseMatrix")) stop("x needs to be a subclass of dsparseMatrix!")
#   x <- as(x, "numeric")

#   # true NAs
#   nas <- which(is.na(x), arr.ind=TRUE)

#   x[x==0] <- NA
#   zeros <- which(is.na(x), arr.ind=TRUE)
#   x[nas] <- NA
#   x[zeros] <- 0
# }

## matrix -> sparse
dropNA <- function(x) {
  if(!is(x, "matrix")) stop("x needs to be a matrix!")
  
  zeros <- which(x==0, arr.ind=TRUE)
  ## keep zeros
  x[is.na(x)] <- 0
  x[zeros] <- NA
  x <- Matrix::drop0(x)
  x[zeros] <- 0
  x
}

# dropNAis.na <- function(x) {
#   if(!is(x, "dsparseMatrix")) stop("x needs to be a subclass of dsparseMatrix!")
#   x <- as(x, "dgCMatrix")
#   
#   ### not represented means NA and 0 means 0
#   ### this coercion keeps 0
#   !as(x, "ngCMatrix")
# }

give.n <- function(x){
  return(c(y = mean(x), label = length(x)))
}

sort.abs <- function(dt, sort.field) dt[order(-abs(dt[[sort.field]]))]


# function to pseudobulk a SingleCellExperiment object
pseudobulk_sce_fn <- function(x, assay = NULL, by, fun = c("sum", "mean", "median"), scale = FALSE) {
  
  # check validity of input arguments
  fun <- match.arg(fun)
  if (is.null(assay))  assay <- assayNames(x)[1] 
  
  # store aggregation parameters &
  # nb. of cells that went into aggregation
  md <- metadata(x)
  md$agg_pars <- list(assay = assay, by = by, fun = fun, scale = scale)
  
  # get aggregation function
  # fun <- switch(fun, sum = "rowSums", mean = "rowMeans", median = "rowMedians")
  
  # drop missing factor levels & tabulate number of cells
  cd <- dplyr::mutate_if(as.data.frame(colData(x)), is.factor, droplevels)
  colData(x) <- DataFrame(cd, row.names = colnames(x),check.names = FALSE)
  md$n_cells <- table(as.data.frame(colData(x)[, by]))
  
  # assure 'by' colData columns are factors so that missing combinations aren't dropped
  for (i in by) 
    if (!is.factor(x[[i]])) 
      x[[i]] <- factor(x[[i]])
  
  # split cells & compute pseudo-bulks
  cs <- .split_cells(x, by)
  # pb <- .pb(x, cs, assay, fun)
  pb <- .pb(x=x, by=by, fun=fun)
  if (scale & length(by) == 2) {
    ls <- lapply(.pb(x, cs, "counts", "rowSums"), colSums)
    pb <- lapply(seq_along(pb), function(i) pb[[i]] / 1e6 * ls[[i]])
    names(pb) <- names(ls)
  }
  
  # construct SCE
  pb <- SingleCellExperiment(pb, metadata = md)
  
  # propagate 'colData' columns that are unique across 2nd 'by'
  if (length(by) == 2) {
    cd <- colData(x)
    ids <- colnames(pb)
    counts <- vapply(ids, function(u) {
      m <- as.logical(match(cd[, by[2]], u, nomatch = 0))
      vapply(cd[m, ], function(u) length(unique(u)), numeric(1))
    }, numeric(ncol(colData(x))))
    cd_keep <- apply(counts, 1, function(u) all(u == 1))
    cd_keep <- setdiff(names(which(cd_keep)), by)
    if (length(cd_keep) != 0) {
      m <- match(ids, cd[, by[2]], nomatch = 0)
      cd <- cd[m, cd_keep, drop = FALSE]
      rownames(cd) <- ids
      colData(pb) <- cd
    }
  }
  return(pb)
}


# split cells by cluster-sample
# auxiliary function to pseudobulk a SingleCellExperiment object
#   - by: character vector specifying colData column(s) to split by. If length(by) == 1, a list of length nlevels(colData$by), else,
#          a nested list with 2nd level of length nlevels(colData$by[2])
.split_cells <- function(x, by) {
  if (is(x, "SingleCellExperiment")) x <- colData(x)
  cd <- data.frame(x[by], check.names = FALSE)
  cd <- data.table(cd, cell = rownames(x)) %>% split(by = by, sorted = TRUE, flatten = FALSE)
  purrr::map_depth(cd, length(by), "cell")
}


# auxiliary function to pseudobulk a SingleCellExperiment object
.pb <- function(x, by, fun) {
  
  # compute pseudobulks
  # y <- scuttle::summarizeAssayByGroup(x, assay.type = assay, ids = (ids <- colData(x)[by]), statistics = fun, BPPARAM = BiocParallel::SerialParam())
  y <- scuttle::summarizeAssayByGroup(x, ids = colData(x)[by], statistics = fun)
  colnames(y) <- y[[by[length(by)]]]
  
  if (length(by) == 1)  return(assay(y))
  
  # reformat into one assay per 'by[1]'
  if (is.factor(ids <- y[[by[1]]]))
    ids <- droplevels(ids)
  is <- split(seq_len(ncol(y)), ids)
  ys <- map(is, ~assay(y)[, .])
  
  # fill in missing combinations
  for (i in seq_along(ys)) {
    fill <- setdiff(unique(y[[by[2]]]), colnames(ys[[i]]))
    if (length(fill != 0)) {
      foo <- matrix(0, nrow(x), length(fill))
      colnames(foo) <- fill
      foo <- cbind(ys[[i]], foo)
      o <- paste(sort(unique(y[[by[2]]])))
      ys[[i]] <- foo[, o]
    }
  }
  return(ys)
}


.summarizeJASPARMotifs <- function(motifs = NULL){

  motifNames <- lapply(seq_along(motifs), function(x){
    namex <- make.names(motifs[[x]]@name)
    if(substr(namex,nchar(namex),nchar(namex))=="."){
      namex <- substr(namex,1,nchar(namex)-1)
    }
    namex <- paste0(namex, "_", x)
    namex
  }) %>% unlist(.)

  motifDF <- lapply(seq_along(motifs), function(x){
    data.frame(
      row.names = motifNames[x],
      name = motifs[[x]]@name[[1]],
      ID = motifs[[x]]@ID,
      strand = motifs[[x]]@strand,
      symbol = ifelse(!is.null(motifs[[x]]@tags$symbol[1]), motifs[[x]]@tags$symbol[1], NA) ,
      family = ifelse(!is.null(motifs[[x]]@tags$family[1]), motifs[[x]]@tags$family[1], NA),
      alias = ifelse(!is.null(motifs[[x]]@tags$alias[1]), motifs[[x]]@tags$alias[1], NA),
      stringsAsFactors = FALSE
    )
  }) %>% Reduce("rbind", .) %>% DataFrame
  
  names(motifs) <- motifNames

  out <- list(motifs = motifs, motifSummary = motifDF)

  return(out)
  
}

.summarizeChromVARMotifs <- function(motifs = NULL){

  motifNames <- lapply(seq_along(motifs), function(x){
    namex <- make.names(motifs[[x]]@name)
    if(grepl("LINE", namex)){
      splitNamex <- stringr::str_split(motifs[[x]]@ID, pattern="\\_", simplify = TRUE)
      namex <- splitNamex[1, grep("LINE",splitNamex[1,]) + 1]
    }
    if(substr(namex,nchar(namex),nchar(namex))=="."){
      namex <- substr(namex,1,nchar(namex)-1)
    }
    namex <- paste0(namex, "_", x)
    namex
  }) %>% unlist(.)

  motifNames2 <- lapply(seq_along(motifs), function(x){
    namex <- make.names(motifs[[x]]@name)
    if(grepl("LINE", namex)){
      splitNamex <- stringr::str_split(motifs[[x]]@ID, pattern="\\_", simplify = TRUE)
      namex <- splitNamex[1, grep("LINE",splitNamex[1,]) + 1]
    }
    if(substr(namex,nchar(namex),nchar(namex))=="."){
      namex <- substr(namex,1,nchar(namex)-1)
    }
    namex
  }) %>% unlist(.)

  motifDF <- lapply(seq_along(motifs), function(x){
    df <- data.frame(
      row.names = motifNames[x],
      name = motifNames2[[x]],
      ID = motifs[[x]]@ID,
      strand = motifs[[x]]@strand,
      stringsAsFactors = FALSE
    )
  }) %>% Reduce("rbind", .) %>% DataFrame

  names(motifs) <- motifNames

  out <- list(motifs = motifs, motifSummary = motifDF)

  return(out)

}

.augment_matrix <-function(mtx, samples) {
  samples <- unique(samples)
  mtx <- t(mtx)
  aug_mtx<-matrix(NA, ncol=ncol(mtx), nrow=length(samples))
  aug_mtx<-mtx[match(samples,rownames(mtx)),,drop=FALSE]
  rownames(aug_mtx)<-samples
  colnames(aug_mtx)<-colnames(mtx)
  return(t(aug_mtx))
}


stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}
                     
                     
# Plotting of browser tracks from bigwigs
plot_browser_track = function(gene = NULL, # Only choose one of gene/input_region, otherwise will use region
                              input_region = NULL, 
                              bigwigs, 
                              plot_expression = F,
                              sce = rna.sce,
                              extend.upstream = 0,
                              extend.downstream = 0, 
                              highlight_regions = NULL,
                              downsample.rate = 1, 
                              no_label = F, 
                              background_color = 'white',
                              highlight_color = 'grey95',
                              bw_height = 2, 
                              gene_height = 1.5, 
                              bw_width = 20,
                              expr_width = 4.5){
    
    
    # x-axis labels function 
    alternate_labels <- function(x) {
      ifelse(seq_along(x) %% 2 == 0, paste("\n", x), x)
    }
    
    if(!is.null(gene)){
        # Get region around gene
        region = ArchRProject@geneAnnotation$genes[grep(paste0('^', gene, '$'), ArchRProject@geneAnnotation$genes$symbol)]
        start(ranges(region)) = start(ranges(region)) - extend.upstream
        end(ranges(region)) = end(ranges(region)) + extend.downstream
        gr <- GRanges(c(seqnames(region)), IRanges(start=start(ranges(region)), end=end(ranges(region))), strand='+')
    }
    if(!is.null(input_region)){
        # Get supplied region as Granges
        gr <- GRanges(c(str_split(input_region, ':') %>% map_chr(1)), 
                      IRanges(start=as.numeric(str_split(str_split(input_region, ':') %>% map_chr(2), '-') %>% map_chr(1)) - extend.upstream, 
                              end=as.numeric(str_split(str_split(input_region, ':') %>% map_chr(2), '-') %>% map_chr(2)) + extend.downstream), 
                        strand='+')
    }
    
    # Get bigwig info
    bw_list = lapply(bigwigs$path, function(x)return(x))
    names(bw_list) = bigwigs$name
    bigwigs.colors = bigwigs$color
    names(bigwigs.colors) = bigwigs$name
                     
    coverages = mclapply(names(bw_list), function(x){
        region_data <- as.data.table(rtracklayer::import(con = bw_list[[x]], 
                    which = gr))
        region_data = rbind(region_data, region_data %>% copy() %>% .[,start:=end]) %>% 
            .[order(start)] %>% 
            .[,c('start', 'score')] %>% 
            .[, bw := names(bw_list[x])]
            # .[, bw := factor(names(bw_list[x]), levels = names(bw_list[x]))] 
        return(region_data)
    }, mc.cores = max(1,detectCores()-2)) %>% rbindlist() %>%
        .[, bw := factor(bw, levels = names(bw_list))] %>% .[order(bw)]
                     
    if(downsample.rate != 1){
        coverages = coverages[start %in% sample(start, nrow(coverages)*downsample.rate)]
    }
                     
    p1 =  ggplot(data = coverages, mapping = aes_string(x = "start", y = "score", fill = "bw")) + 
                     geom_area(color = 'black', linewidth = 0) + 
                     # geom_area() +
                     facet_wrap(facets = ~factor(bw, levels = names(bw_list)), 
                            strip.position = "left", ncol = 1) +
        # ggtitle(label = paste0(seqnames(gr)),
        #     subtitle = paste0(ranges(gr))) + #, ':', ranges(gr))) + 
        geom_hline(yintercept = 0, color='black', linewidth = 0.1) +
        geom_segment(data = as.data.table(gr)[,bw:=bigwigs[1,name]], # Add background segment to make sure x-axis is correct
            mapping = aes(x = start, y = 0, xend = end, yend = 0), color = 'white', alpha=0) +
        scale_x_continuous(expand=c(0,0)) + 
        scale_y_continuous(expand=c(0,0.005)) +
        Signac::theme_browser(axis.text.y = TRUE) + 
        theme(axis.text.y = element_blank(),
                strip.text.y = element_text(hjust=1, color='black', size = 6, vjust=1),
                axis.line = element_blank(),
                axis.line.x = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks = element_blank(), 
                axis.title = element_blank(),
                legend.position = 'none',
                plot.title = element_text(size = 7, hjust=0.5, vjust = -0.5, face='plain'),
                plot.subtitle = element_text(size = 5, hjust=0.5, face='plain'),
                panel.background = element_rect(fill=background_color),
                panel.spacing = unit(0, "lines"),
                plot.margin = margin(0,0,0,0)) 
    
    # Remove y-text if needed
    if(no_label){
        p1 = p1 + theme(strip.text.y = element_text(size = 0), 
                        strip.background = element_blank(),
                        plot.margin = unit(c(0,0,0,0), units = "lines")) 
    }
    
    # Add Highlights
    if(!is.null(highlight_regions)){
        # Keep only highlights in the display range
        highlights = GRanges(c(str_split(highlight_regions, ':') %>% map_chr(1)), 
                      IRanges(start=as.numeric(str_split(str_split(highlight_regions, ':') %>% map_chr(2), '-') %>% map_chr(1)), 
                              end=as.numeric(str_split(str_split(highlight_regions, ':') %>% map_chr(2), '-') %>% map_chr(2))), 
                        strand='+')
        
        highlight_in_region = findOverlaps(highlights, gr, ignore.strand = T, type = 'within')
        highlight_regions = highlight_regions[queryHits(highlight_in_region)]
        
        # Add highlight display
        if(length(highlight_regions) > 0){ # Only if there is any overlap
            highlights.dt = data.table(start = as.numeric(str_split(str_split(highlight_regions, ':') %>% map_chr(2), '-') %>% map_chr(1)),
                                       end = as.numeric(str_split(str_split(highlight_regions, ':') %>% map_chr(2), '-') %>% map_chr(2))) %>% 
                            .[, .(start, end, vector_column = rep(unique(p1$data$bw), each = 1)), by = .(row_index = 1:nrow(.))] %>% 
                            .[, position := start] %>% .[,score := 0]
            
            if(length(highlight_color) == 1){
                highlight_color = rep(highlight_color, length(highlight_regions))
                names(highlight_color) = as.character(unique(highlights.dt$start))
            }
            
            colors = c(bigwigs$color, as.character(highlight_color))
            names(colors) = c(as.character(bigwigs$name), names(highlight_color))

            p1 = p1 + 
                geom_rect(data = highlights.dt, aes(xmin=start, xmax=end, ymin=0, ymax=max(p1$data$score), fill = as.character(start)), 
                                    color=NA, alpha=1) +
                scale_fill_manual(values = colors) 
            # Change layers around so highlight is behind bigwig
            p1$layers = list(p1$layers[[3]], p1$layers[[4]], p1$layers[[1]], p1$layer[[2]])
        }
    }else{ p1 = p1 + scale_fill_manual(values = bigwigs.colors)}

    # Get genome annotation
    annotation.genes = ArchRProject@geneAnnotation$genes 
    annotation.genes = annotation.genes[grep("^Rik|Rik$|^Gm|^Mir", annotation.genes$symbol, invert=T)] # remove uninformative genes from plotting
 
    annotation.exons = ArchRProject@geneAnnotation$exons
    annotation.exons = annotation.exons[grep("^Rik|Rik$|^Gm|^Mir", annotation.exons$symbol, invert=T)]
    
    # Search for genes within region
    tmp1 = findOverlaps(gr, annotation.genes, ignore.strand = T)
    tmp2 = findOverlaps(gr, annotation.exons, ignore.strand = T)
    
    # Add gene track if there are any genes in the window
    if(length(subjectHits(tmp1)) > 0){
        # get overlap among genes
        overlap = as.data.table(findOverlaps(annotation.genes[subjectHits(tmp1)], annotation.genes[subjectHits(tmp1)], type = 'any', ignore.strand = T))

        # check if genes overlap, and if so move different rows
        solver = list()
        for(x in unique(overlap$queryHits)){
            overlaps = overlap[queryHits==x & subjectHits <= x]
            nr_overlaps = nrow(overlaps)
        # No overlap --> Row = 0
            if(nr_overlaps == 1){
                rows = 1
        # 1 overlap --> check which row it was on
            }else{
               rows = setdiff(1:10, rbindlist(solver[overlaps$subjectHits])$row)[1]
            }
        #If overlap is after, do nothing
            solver[[x]] = data.table(queryHits = x,
                                     row = rows)
        }
        solved = rbindlist(solver)

        # keep overlapping genes/exons
        tmp1 = as.data.table(annotation.genes[subjectHits(tmp1)]) %>% 
            .[,row := solved$row] %>% .[width>200] # remove super short genes

        tmp2 = as.data.table(annotation.exons[subjectHits(tmp2)]) %>% 
            merge(., tmp1[,c('symbol', 'row')], by='symbol') # Add rown information

        # Cutoff genes/exon if they fall outside of the range
        tmp1 = tmp1 %>% 
            .[, `:=`(start = ifelse(start <= start(ranges(gr)), start(ranges(gr)), start),
                     end = ifelse(end >= end(ranges(gr)), end(ranges(gr)), end))]

        tmp2 = tmp2 %>% 
            .[, `:=`(start = ifelse(start <= start(ranges(gr)), start(ranges(gr)), start),
                     end = ifelse(end >= end(ranges(gr)), end(ranges(gr)), end))] %>% unique()

        # Add Arrows to indicate orientation
        width = end(ranges(gr)) - start(ranges(gr))

        tmp3 = lapply(1:nrow(tmp1), function(x){
            tmp = tmp1[x,]
            tmpp = data.table(
                          start = seq(tmp$start + 20, tmp$end + 20, by = width/20),
                          end = seq(tmp$start + 40, tmp$end + 40, by = width/20),
                          symbol = tmp$symbol,
                          row = tmp$row)
            return(tmpp)
        }) %>% rbindlist()


        # Plot gene track 
        p = ggplot() + # rbind(tmp1, tmp1 %>% copy() %>% .[,row:=row+1]), aes(start, row)
            geom_segment(data = as.data.table(gr), # Add background segment to make sure x-axis is correct
                         mapping = aes(x = start, y = 1, xend = end, yend = 1), color = 'white') +
            geom_segment(data = tmp1, 
                mapping = aes(x = start, y = row, 
                xend = end, yend = row), color = 'darkblue', show.legend = FALSE, linewidth = unit(0.5, 'mm')) +
            geom_segment(data = tmp2, 
                mapping = aes(x = start, y = row, 
                xend = end, yend = row), color = '#4272f5', show.legend = FALSE, linewidth = unit(2, 'mm'))
        for(i in 1:nrow(tmp1)){
            tmp = tmp1[i,]
            if(as.character(tmp$strand) == '+'){
                p = p + geom_segment(data = tmp3[symbol == tmp$symbol], 
                        mapping = aes(x = start, y = row, 
                xend = end, yend = row), arrow = arrow(ends = "last", 
                        type = "open", angle = 45, length = unit(x = 1, 
                        units = "mm")), show.legend = FALSE, color='black', 
                        linewidth = 0.1)
               }else{
              p = p + geom_segment(data = tmp3[symbol == tmp$symbol], 
                        mapping = aes(x = start, y = row, 
                xend = end, yend = row), arrow = arrow(ends = "first", 
                        type = "open", angle = 45, length = unit(x = 1, 
                        units = "mm")), show.legend = FALSE, color='black', 
                        linewidth = 0.1)  
            }
        }


        p = p + 
            #xlim(start(ranges(gr)), end(ranges(gr))) + 
            ggrepel::geom_text_repel(data = tmp1, aes(label = symbol, x = end - (end - start)/2 , y = row), 
                                     size = 2,
                                     color = 'black', 
                                     direction = "y",
                                     force_pull   = 0, # do not pull toward data points
                                     nudge_y      = 1,
                                     segment.size = 0.2
                                    ) + 
            ylim(0.5, max(c(1, tmp1$row) + 1))

    }else{
        p = ggplot() + # rbind(tmp1, tmp1 %>% copy() %>% .[,row:=row+1]), aes(start, row)
            geom_segment(data = as.data.table(gr), # Add background segment to make sure x-axis is correct
                         mapping = aes(x = start, y = 1, xend = end, yend = 1), color = 'white') +
            ylim(0.5, 2)
        
        
        
    }
    p = p + 
        scale_x_continuous(expand = c(0,0), labels = alternate_labels, breaks = seq(start(ranges(gr)), end(ranges(gr)), 3000)) + # breaks = scales::pretty_breaks(n = 1))
        theme_void() + 
        theme( #axis.text.x = element_text(size=15, color='black', vjust=0.9),
               axis.ticks.x = element_line(color='black', linewidth=0.1), 
               axis.ticks.length.x = unit(0.3, 'mm'),
               axis.line.x = element_line(color='black', linewidth=0.1),
               plot.margin = margin(0,0,0,0) 
    )
    
    if(plot_expression == F){
        patchwork::wrap_plots(plotlist = list(p1, p), ncol = 1, heights = unit(c(nrow(bigwigs)*bw_height, gene_height), c('mm', 'mm')))
    }else{
        stopifnot(!is.null(gene))
        meta.plot = as.data.table(colData(sce)) %>% 
            .[,gene := as.vector(logcounts(sce[gene,]))] %>% 
            merge(., bigwigs, by = c('celltype')) %>% 
            .[,c('genotype', 'celltype', 'name', 'color', 'gene')]

        p2 = ggplot(meta.plot %>% copy(),
                        # .[,celltype_v2 := factor(celltype_v2, levels = rev(WT_plot$name))], 
                    aes(gene, name, fill = name)) + 
            geom_violin(scale = 'width',  size = 0.2) + 
            geom_boxplot(width=0.2, fill='white', linewidth = 0.15, outlier.shape=NA) + 
            scale_fill_manual(values = bigwigs.colors) + 
            scale_x_continuous(expand = c(0, 0)) + 
            scale_y_discrete(expand = c(0, 0), limits = rev(names(bw_list))) + 
            geom_hline(yintercept = seq(0.5, nrow(WT_plot), 1), linewidth = 0.1) + 
            ggtitle(gene) + 
            theme_void() + 
            theme(axis.text.x = element_text(size = 5, vjust = 0),
                  axis.line = element_line(linewidth = 0.1),
                  axis.ticks.x = element_line(linewidth=0.1, color = 'black'),
                  plot.title = element_text(size = 6, hjust=0.5, face='plain'),
                  plot.subtitle = element_text(size = 7, hjust=0.5, face='plain'),
                  panel.background = element_rect(fill = background_color, color = NA),
                  plot.margin = margin(0,0,0,0),
                  legend.position = 'none')

         patchwork::wrap_plots(p1, p2, p, ncol = 2, widths = unit(c(bw_width,expr_width), c('mm','mm')), heights = unit(c(nrow(bigwigs)*bw_height, gene_height), c('mm', 'mm')))
    }
}
                     
do_GO = function(genes, background){
    suppressPackageStartupMessages({
    library(biomaRt)
    library(limma)
    library(org.Mm.eg.db)
    })
    mapping = fread('/rds/project/rds-SDzz0CATGms/users/bt392/atlasses/Mmusculus_genes_BioMart.07052024.txt')
    # Retrieve GO
    go = as.data.table(goana(mapping[match(genes, external_gene_name),entrezgene_id], species="Mm"), keep.rownames=T) %>% setnames('rn', 'GOALL') %>% 
        .[order(P.DE)] %>%
        .[,FDR := p.adjust(P.DE, method = 'fdr')] %>% 
        .[P.DE < 0.1] %>% 
        .[, enrichment := (length(genes)/N) / (N/32285)] %>%  # 32285 = nrow(rna.sce) # Enrichment should be depending on the number of genes in the term
        .[, enrichment := ifelse(enrichment > 100, 100, enrichment)] # set threshold for enrichment at 100
    # Add gene names to annotations
    GO_genes = AnnotationDbi::select(org.Mm.eg.db, keytype="GOALL", keys=go[, GOALL], columns='ENTREZID') %>% 
        as.data.table() %>% 
        .[, EVIDENCEALL := NULL] %>% 
        unique() %>% 
        merge(., mapping[,c('entrezgene_id', 'external_gene_name')] %>% 
                  .[, entrezgene_id := as.character(entrezgene_id)],
              by.x = 'ENTREZID', by.y = 'entrezgene_id')
    # Add genes into GO list
    go = mclapply(1:nrow(go), function(x){ # nrow(go)
        GO_genes_tmp = GO_genes[GOALL == go[x, GOALL], external_gene_name]
        genes_GO = intersect(GO_genes_tmp, genes)
        tmp = go[x,] %>% .[, `:=`(genes = paste0(genes_GO, collapse = ', '),
                                  N = length(intersect(background, GO_genes_tmp)))] %>%  # Ntot is only the number of genes in the GO that are also in the SCE object
        .[,enrichment := (DE/N)/(N/32285)] 
    }, mc.cores = 10) %>% rbindlist()
    
    options(repr.plot.width=10, repr.plot.height=4)
    plot = go %>% copy() %>%
           .[FDR < 0.05] %>% .[Ont == 'BP'] %>% .[N > 20 & N < 1000 & DE > 3 & enrichment > 3] %>% .[order(FDR)] %>% 
           # .[FDR < 0.05] %>% .[Ont == 'BP'] %>% .[N > 50 & N < 1000 & DE > 3 & enrichment > 5] %>% .[order(FDR)] %>% 
            .[,Term := factor(Term, levels = rev(Term))] %>% 
            .[, enrichment := ifelse(enrichment > 100, 100, enrichment)] %>% 
            head(10)
    p = ggplot(plot, 
           aes(enrichment, Term, color = -log10(FDR), fill = -log10(FDR), size = DE)) +  #, 
        geom_bar(stat = 'identity', aes(fill = -log10(FDR)), color = NA,  width = 0.05) +
        geom_point() + 
        geom_text(aes(x = max(plot$enrichment)/100, label = Term), vjust = -0.5, hjust = 0, color = 'black', size = 5 * 0.35) + 
        guides(fill = guide_colourbar(
             barwidth = 0.5, barheight = 2.5,
            frame.colour = "black", 
            ticks.colour = "black"
        )) + 
        scale_fill_gradientn(colors = c('white', 'blue', 'red'), limits = c(0, max(-log10(plot$FDR))), name = '-log10(FDR)') + 
        scale_color_gradientn(colors = c('white', 'blue', 'red'), limits = c(0, max(-log10(plot$FDR))), name = '-log10(FDR)') + 
        scale_x_continuous(limits = c(0, max(plot$enrichment)), expand = c(0,0,0,3)) + 
        theme_bw() + 
        theme(axis.text = element_text(size = 5, color = 'black'),
                axis.text.y = element_blank(), 
                axis.ticks.y = element_blank(),
                axis.title.x = element_text(size = 7),
                axis.title.y = element_blank(),
                panel.border = element_rect(color='black', linewidth=0.5, fill = NA),
                panel.background = element_rect(fill = 'white'),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                plot.margin = margin(0, 0, 0, 0),  # Adjust the margin around the plot
                # legend.margin = margin(0, 0, 0, 0),  # Adjust the margin inside the legend box
                # legend.box.margin = margin(0, 0, 0, 0),  # Adjust the margin outside the legend box
                legend.key.size  = unit(1, "mm"),
                legend.key = element_blank(),
                legend.position = 'right')
    return(list(go, p))
}