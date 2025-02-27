#' PlotDORC
library(RcppRoll)
library(dplyr)
library(scales)
library(fastmatch )
#' @param metadata 
#' @param region A set of genomic coordinates to show. Can be a GRanges object,
#' a string encoding a genomic position, a gene name, or a vector of strings
#' describing the genomic coordinates or gene names to plot. If a gene name is
#' supplied, annotations must be present in the assay.
#' @param annotation 
#' @param features 
#' @param assay 
#' @param show.bulk 
#' @param peak.data 
#' @param expression.data 
#' @param anno 
#' @param Links 
#' @param peaks 
#' @param peaks.group.by 
#' @param ranges 
#' @param ranges.group.by 
#' @param ranges.title 
#' @param links 
#' @param tile 
#' @param tile.size 
#' @param tile.cells 
#' @param group.by 
#' @param window 
#' @param extend.upstream 
#' @param extend.downstream 
#' @param ymax 
#' @param scale.factor 
#' @param cells 
#' @param idents 
#' @param sep 
#' @param heights 
#' @param max.downsample 
#' @param downsample.rate 
#'
#' @return
#' @export
#'
#' @examples
CoveragePlot <- function(
  region,
  annotation,
  Links,
  features = NULL,
  assay = NULL,
  show.bulk = TRUE,
  peak.use,
  bigwig=NULL,
  anno= TRUE,
  peaks = TRUE,
  peaks.group.by = NULL,
  ranges = NULL,
  ranges.group.by = NULL,
  ranges.title = "Ranges",
  links = TRUE,
  tile = FALSE,
  tile.size = 100,
  tile.cells = 100,
  window = 300,
  extend.upstream = 200000,
  extend.downstream = 200000,
  ymax = NULL,
  scale.factor = NULL,
  cells = NULL,
  idents = NULL,
  sep = c("-", "-"),
  heights = NULL,
  bigwig.type='coverage',
  max.downsample = 3000,
  downsample.rate = 0.1,
  add_explot=F
) {
  if( !('type'%in% colnames(annotation@elementMetadata))){
    annotation$type='exon'
  }
  if(!region%in%Links$gene){stop('gene is not in the DORC_output')}
  features=region
  
  # if (!is.null(x = idents)) {
  #   ident.cells <- WhichCells(metadata = meta.data, idents = idents)
  #   cells <- intersect(x = cells, y = ident.cells)
  # }
  region <- FindRegion(
    annotation = annotation,
    region = region,
    sep = sep,
    extend.upstream = extend.upstream,
    extend.downstream = extend.downstream
  )
  # # reads.per.group <- AverageCounts(
  #   peak.data,
  #   meta.data,
  #   verbose = FALSE
  # )
  # # cells.per.group <- CellsPerGroup(
  #   meta.data
  # )
  # group.scale.factors <- suppressWarnings(reads.per.group * cells.per.group)
  #scale.factor <- SetIfNull(
  #   x = scale.factor, y = median(x = group.scale.factors)
  # )
  
  
  
  

  # p <- CoverageTrack(
  #   cutmat = cutmat,
  #   region = region,
  #   group.scale.factors = group.scale.factors,
  #   scale.factor = scale.factor,
  #   window = window,
  #   ymax = ymax,
  #   obj.groups = obj.groups,
  #   downsample.rate = downsample.rate,
  #   max.downsample = max.downsample
  # )
  if(add_explot==TRUE){
    if (!is.null(x = features)) {
      ex.plot <- ExpressionPlot(
        expression.data = expression.data,
        meta.data = meta.data,
        cells = cells,
        features = features,
        idents = idents
      )
      widths <- c(10, length(x = features))
    } else {
      ex.plot <- NULL
      widths <- NULL
    }}
  if (anno) {
    gene.plot <- AnnotationPlot(annotation = annotation, region = region)
  } else {
    gene.plot <- NULL
  }
  if (links) {
    link.plot <- LinkPlot(Links, region = region)
  } else {
    link.plot <- NULL
  }
  if (peaks) {
    peak.plot <- PeakPlot(
      peak.use,
      region = region
    )
  } else {
    peak.plot <- NULL
  }
  if (!is.null(x = bigwig)) {
    if (!inherits(x = bigwig, what = "list")) {
      warning("BigWig should be a list of file paths")
      bigwig <- list("bigWig" = bigwig)
    }
    if (length(x = bigwig.type == 1)) {
      bigwig.type <- rep(x = bigwig.type, length(x = bigwig))
    } else if (length(x = bigwig.type) != length(x = bigwig)) {
      stop("Must supply a bigWig track type for each bigWig file")
    }
    # iterate over list of files
    bigwig.tracks <- list()
    for (i in seq_along(along.with = bigwig)) {
      bigwig.tracks[[i]] <- BigwigTrack(
        region = region,
        bigwig = bigwig[[i]],
        y_label = names(x = bigwig)[[i]],
        type = bigwig.type[[i]]
      )
    }
  } else {
    bigwig.tracks <- NULL
  }
  
  if (!is.null(x = ranges)) {
    range.plot <- PeakPlot(
      object = object,
      assay = assay,
      region = region,
      peaks = ranges,
      group.by = ranges.group.by,
      color = "brown3") +
      ylab(ranges.title)
  } else {
    range.plot <- NULL
  }
  if (tile) {
    # reuse cut matrix
    tile.df <- ComputeTile(
      cutmatrix = cutmat,
      groups = obj.groups,
      window = tile.size,
      n = tile.cells,
      order = "total"
    )
    tile.plot <- CreateTilePlot(
      df = tile.df,
      n = tile.cells
    )
  } else {
    tile.plot <- NULL
  }
  # if (show.bulk) {
  #   metadata$bulk<- "All cells"
  #   bulk=metadata$bulk
  #   names(bulk)=colnames(peak.data)
  #   reads.per.group <- AverageCounts(
  #     peak.data,
  #     meta.data = bulk,
  #     verbose = FALSE
  #   )
  #   cells.per.group <- CellsPerGroup(
  #     bulk
  #   )
  #   bulk.scale.factor <- suppressWarnings(reads.per.group * cells.per.group)
  #   bulk.groups <- rep(x = "All cells", length(x = obj.groups))
  #   names(x = bulk.groups) <- names(x = obj.groups)
  #   
  #   bulk.plot <- CoverageTrack(
  #     cutmat = cutmat,
  #     region = region,
  #     group.scale.factors = bulk.scale.factor,
  #     scale.factor = scale.factor,
  #     window = window,
  #     ymax = ymax,
  #     obj.groups = bulk.groups,
  #     downsample.rate = downsample.rate,
  #     max.downsample = max.downsample
  #   ) +
  #     scale_fill_manual(values = "grey") +
  #     ylab("")
  # } else {
  #   bulk.plot <- NULL
  # }
  nident <- length(x = unique(x = idents))
  bulk.height <- (1 / nident) * 10
  bw.height <- bulk.height * length(x = bigwig.tracks)
  heights <- SetIfNull(
    x = heights, y = c( bw.height, 2, 1, 3)
  )
  # pre-combine bigwig tracks
  if (!is.null(x = bigwig.tracks)) {
    bigwig.tracks <- CombineTracks(
      plotlist = bigwig.tracks,
      heights = rep(x = 1, length(x = bigwig.tracks))
    )
  }
  if(add_explot==TRUE){
    p <- CombineTracks(
      plotlist = list( bigwig.tracks, gene.plot,
                       peak.plot,  link.plot),
      expression.plot = ex.plot,
      heights = heights,
      widths = widths
    ) & theme(
      legend.key.size = unit(x = 1/2, units = "lines"),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8)
    )}else{ 
      p <- CombineTracks(
        plotlist = list( bigwig.tracks, gene.plot, peak.plot,  link.plot),
        # expression.plot = ex.plot,
        heights = heights,
        widths = widths
      ) & theme(legend.key.size = unit(x = 1/2, units = "lines"),
                legend.text = element_text(size = 7),
                legend.title = element_text(size = 8)
      )}
  
  return(p)
}


SetIfNull <- function(x, y) {
  if (is.null(x = x)) {
    return(y)
  } else {
    return(x)
  }
}

WhichCells<- function(
  metadata,
  cells = NULL,
  idents = NULL,
  downsample = Inf,
  seed = 1) 
{
  if (!is.null(x = seed)) {
    set.seed(seed = seed)
  }
  cells <- union(cells,  colnames(expression.data))
  if (is.numeric(x = cells)) {
    cells <- colnames(expression.data)
  }
  cell.order <- cells
  if (!is.null(x = idents)) {
    if (any(!idents %in% levels(x = factor(metadata)))) {
      stop(
        "Cannot find the following identities in the object: ",
        paste(
          idents[!idents %in% levels(x = factor(metadata))],
          sep = ', '
        )
      )
    }
    cells.idents <- unlist(x = lapply(
      X = idents,
      FUN = function(i) {
        cells.use <- which(x = as.vector(x = factor(metadata) == i))
        cells.use <- cells[cells.use]
        return(cells.use)
      }
    ))
    cells <- intersect(x = cells, y = cells.idents)
  }
  
  cells <- CellsByIdentities(metadata = metadata, idents = idents,cells = cells)
  cells <- lapply(
    X = cells,
    FUN = function(x) {
      if (length(x = x) > downsample) {
        x <- sample(x = x, size = downsample, replace = FALSE)
      }
      return(x)
    }
  )
  cells <- as.character(x = na.omit(object = unlist(x = cells, use.names = FALSE)))
  cells <- cells[na.omit(object = match(x = cell.order, table = cells))]
  return(cells)
}

CellsByIdentities <- function(metadata, idents = NULL, cells = NULL) {
  
  
  if (length(x = cells) == 0) {
    stop("Cannot find cells provided")
  }
  
  
  if (length(x = idents) == 0) {
    stop("None of the provided identity class levels were found", call. = FALSE)
  }
  cells.idents <- sapply(
    X = idents,
    FUN = function(i) {
      return(names(metadata)[metadata==i])
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  return(cells.idents)
}

FindRegion <- function(
  annotation,
  region,
  sep = c("-", "-"),
  assay = NULL,
  extend.upstream = 0,
  extend.downstream = 0
) {
  annot.sub<- annotation[annotation$gene_name == region]
  region <- GRanges(seqnames = as.character(x = seqnames(x = annot.sub))[[1]],
                    ranges = IRanges(start = min(start(x = annot.sub)),
                                     end = max(end(x = annot.sub)))) 
  if (is.null(x = region)) {
    stop("Gene not found")
  }
  
  
  region <- suppressWarnings(expr = Extend(
    x = region,
    upstream = extend.upstream,
    downstream = extend.downstream
  )
  )
  return(region)
}

AverageCounts <- function(data, meta.data, verbose = TRUE) {
  # pull nCount_ column
  total.df<-as.data.frame( colSums(data))
  colnames(total.df) <- "readcount"
  total.df$cell <- rownames(x = total.df)
  total.df$group <- meta.data[total.df$cell]
  total.df <- dplyr::group_by(total.df, group)
  group.means <- dplyr::summarize(.data = total.df, mn = mean(x = readcount))
  results <- group.means$mn
  names(x = results) <- group.means$group
  return(results)
}

CellsPerGroup <- function(
  meta.data
) {
  cells.per.group <- table(meta.data, useNA = "always")
  lut <- as.vector(x = cells.per.group)
  names(x = lut) <- names(x = cells.per.group)
  return(lut)
}

CutMatrix <- function(region, peak.use, cells, verbose = TRUE, fragment.path) {
  # run SingleFileCutMatrix for each fragment file and combine results
  cellmap <- peak.use
  names(cellmap)=cellmap
  library(Rsamtools)
  tabix.file <- TabixFile(file = fragment.path)
  open(con = tabix.file)
  # remove regions that aren't in the fragment file
  # seqnames.in.both <-
  
  #   intersect(
  #   x = seqnames(x = region),
  #   y = seqnamesTabix(file = tabix.file)
  # )
  # region <- keepSeqlevels(
  #   x = region,
  #   value = seqnames.in.both,
  #   pruning.mode = "coarse"
  # )
  if (length(x = region) != 0) {
    cm <- SingleFileCutMatrix(
      region = region,
      cellmap = cellmap,
      tabix.file = tabix.file,
      cells = cells,
      verbose = FALSE
    )
    
  }
  close(con = tabix.file)
  return(cm)
}


SingleFileCutMatrix <- function(cellmap, region, cells = NULL, tabix.file, verbose = TRUE) {
  # if multiple regions supplied, must be the same width
  cells <- SetIfNull(x = cells, y = names(x = cellmap))
  if (length(x = region) == 0) {
    return(NULL)
  }
  fragments <- GetReadsInRegion(
    region = region,
    cellmap = cellmap,
    cells = cells,
    tabix.file = tabix.file,
    verbose = verbose
  )
  start.lookup <- start(x = region)
  names(start.lookup) <- seq_along(region)
  # if there are no reads in the region
  # create an empty matrix of the correct dimension
  if (nrow(x = fragments) == 0) {
    cut.matrix <- sparseMatrix(
      i = NULL,
      j = NULL,
      dims = c(length(x = cells), width(x = region)[[1]])
    )
  } else {
    fragstarts <- start.lookup[fragments$ident] + 1
    cut.df <- data.frame(
      position = c(fragments$start, fragments$end) - fragstarts,
      cell = c(fragments$cell, fragments$cell),
      stringsAsFactors = FALSE
    )
    cut.df <- cut.df[
      (cut.df$position > 0) & (cut.df$position <= width(x = region)[[1]]),
      ]
    cell.vector <- seq_along(along.with = cells)
    names(x = cell.vector) <- cells
    cell.matrix.info <- cell.vector[cut.df$cell]
    cut.matrix <- sparseMatrix(
      i = cell.matrix.info,
      j = cut.df$position,
      x = 1,
      dims = c(length(x = cells), width(x = region)[[1]])
    )
  }
  rownames(x = cut.matrix) <- cells
  colnames(x = cut.matrix) <- seq_len(width(x = region)[[1]])
  return(cut.matrix)
}

GetReadsInRegion <- function(
  cellmap,
  region,
  tabix.file,
  cells = NULL,
  verbose = TRUE,
  ...
) {
  file.to.object <- names(x = cellmap)
  names(x = file.to.object) <- cellmap
  
  if (verbose) {
    message("Extracting reads in requested region")
  }
  if (!is(object = region, class2 = "GRanges")) {
    region <- StringToGRanges(regions = region, ...)
  }
  # remove regions that aren't in the fragment file
  common.seqlevels <- intersect(
    x = seqlevels(x = region),
    y = seqnamesTabix(file = tabix.file)
  )
  region <- keepSeqlevels(
    x = region,
    value = common.seqlevels,
    pruning.mode = "coarse"
  )
  reads <- scanTabix(file = tabix.file, param = region)
  reads <- TabixOutputToDataFrame(reads = reads)
  library(fastmatch)
  reads <- reads[
    fmatch(x = reads$cell, table = cellmap, nomatch = 0L) > 0,
    ]
  # convert cell names to match names in object
  reads$cell <- file.to.object[reads$cell]
  if (!is.null(x = cells)) {
    reads <- reads[reads$cell %in% cells, ]
  }
  if (nrow(reads) == 0) {
    return(reads)
  }
  reads$length <- reads$end - reads$start
  return(reads)
}

TabixOutputToDataFrame <- function(reads, record.ident = TRUE) {
  if (record.ident) {
    nrep <- elementNROWS(x = reads)
  }
  reads <- unlist(x = reads, use.names = FALSE)
  df <- fread(
    text = reads,
    sep = "\t",
    header = FALSE,
    fill = TRUE
  )
  if(ncol(df)==4){
    colnames(df)=c("chr", "start", "end", "cell")
    df$count=1
  }else{
    colnames(df)=c("chr", "start", "end", "cell", "count")
  }
  
  if (nrow(x = df) != length(x = reads)) {
    df <- df[!is.na(x = df$start), ]
  }
  if (record.ident) {
    df$ident <- rep(x = seq_along(along.with = nrep), nrep)
  }
  return(df)
}


CoverageTrack <- function(
  cutmat,
  region,
  group.scale.factors,
  scale.factor,
  obj.groups,
  ymax,
  downsample.rate,
  window = 200,
  max.downsample = 3000
) {
  window.size <- width(x = region)
  levels.use <- levels(x = obj.groups)
  coverages <- ApplyMatrixByGroup(
    mat = cutmat,
    fun = colSums,
    groups = obj.groups,
    group.scale.factors = group.scale.factors,
    scale.factor = scale.factor,
    normalize = TRUE
  )
  if (!is.na(x = window)) {
    coverages <- group_by(.data = coverages, group)
    library(RcppRoll )
    coverages <- mutate(.data = coverages, coverage = roll_sum(
      x = norm.value, n = window, fill = NA, align = "center"
    ))
    coverages <- ungroup(x = coverages)
  } else {
    coverages$coverage <- coverages$norm.value
  }
  chromosome <- as.character(x = seqnames(x = region))
  start.pos <- start(x = region)
  end.pos <- end(x = region)
  coverages <- coverages[!is.na(x = coverages$coverage), ]
  coverages <- group_by(.data = coverages, group)
  sampling <- min(max.downsample, window.size * downsample.rate)
  coverages <- slice_sample(.data = coverages, n = sampling)
  
  # restore factor levels
  if (!is.null(x = levels.use)) {
    colors_all <- hue_pal()(length(x = levels.use))
    names(x = colors_all) <- levels.use
    coverages$group <- factor(x = coverages$group, levels = levels.use)
  }
  ymax <- SetIfNull(x = ymax, y = signif(
    x = max(coverages$coverage, na.rm = TRUE), digits = 2)
  )
  ymin <- 0
  
  gr <- GRanges(
    seqnames = chromosome,
    IRanges(start = start.pos, end = end.pos)
  )
  p <- ggplot(
    data = coverages,
    mapping = aes(x = position, y = coverage, fill = group)
  ) +
    geom_area(stat = "identity") +
    geom_hline(yintercept = 0, size = 0.1) +
    facet_wrap(facets = ~group, strip.position = "left", ncol = 1) +
    xlab(label = paste0(chromosome, " position (bp)")) +
    ylab(label = paste0("Normalized accessibility \n(range ",
                        as.character(x = ymin), " - ",
                        as.character(x = ymax), ")")) +
    ylim(c(ymin, ymax)) +
    theme_browser(legend = FALSE)
  if (!is.null(x = levels.use)) {
    p <- p + scale_fill_manual(values = colors_all)
  }
  return(p)
}

ApplyMatrixByGroup <- function(
  mat,
  groups,
  fun,
  normalize = TRUE,
  group.scale.factors = NULL,
  scale.factor = NULL
) {
  if (normalize) {
    if (is.null(x = group.scale.factors) | is.null(x = scale.factor)) {
      stop("If normalizing counts, supply group scale factors")
    }
  }
  all.groups <- as.character(x = unique(x = groups))
  if (any(is.na(x = groups))) {
    all.groups <- c(all.groups, NA)
  }
  ngroup <- length(x = all.groups)
  npos <- ncol(x = mat)
  
  group <- unlist(
    x = lapply(X = all.groups, FUN = function(x) rep(x, npos))
  )
  position <- rep(x = as.numeric(x = colnames(x = mat)), ngroup)
  count <- vector(mode = "numeric", length = npos * ngroup)
  
  for (i in seq_along(along.with = all.groups)) {
    grp <- all.groups[[i]]
    if (is.na(x = grp)) {
      pos.cells <- names(x = groups)[is.na(x = groups)]
    } else {
      pos.cells <- names(x = groups)[groups == all.groups[[i]]]
    }
    if (length(x = pos.cells) > 1) {
      totals <- fun(x = mat[pos.cells, ])
    } else {
      totals <- mat[pos.cells, ]
    }
    count[((i - 1) * npos + 1):((i * npos))] <- totals
  }
  
  # construct dataframe
  coverages <- data.frame(
    "group" = group, "position" = position, "count" = count,
    stringsAsFactors = FALSE
  )
  
  if (normalize) {
    scale.factor <- SetIfNull(
      x = scale.factor, y = median(x = group.scale.factors)
    )
    coverages$norm.value <- coverages$count /
      group.scale.factors[coverages$group] * scale.factor
  } else {
    coverages$norm.value <- coverages$count
  }
  return(coverages)
}
theme_browser <- function(..., legend = TRUE) {
  browser.theme <- theme_classic() +
    theme(
      axis.text.y = element_blank(),
      strip.background = element_blank(),
      strip.text.y.left = element_text(angle = 0)
    )
  if (!legend) {
    browser.theme <- browser.theme +
      theme(
        legend.position = "none"
      )
  }
  return(browser.theme)
}
ExpressionPlot <- function(
  expression.data,
  meta.data,
  cells,
  features,
  idents = NULL
) {
  # get data
  data.plot <- expression.data[features, ]
  obj.groups <- factor(meta.data[cells])
  # if levels set, define colors based on all groups
  levels.use <- levels(x = obj.groups)
  if (!is.null(x = levels.use)) {
    colors_all <- hue_pal()(length(x = levels.use))
    names(x = colors_all) <- levels.use
  }
  if (!is.null(x = idents)) {
    cells.keep <- names(x = obj.groups)[
      fmatch(x = obj.groups, table = idents, nomatch = 0L) > 0
      ]
    if (length(x = features) > 1) {
      data.plot <- data.plot[, cells.keep]
    } else {
      data.plot <- data.plot[cells.keep]
    }
    obj.groups <- obj.groups[cells.keep]
  }
  # construct data frame
  if (length(x = features) == 1) {
    df <- data.frame(
      gene = features,
      expression = data.plot,
      group = obj.groups
    )
  } else {
    df <- data.frame()
    for (i in features) {
      df.1 <- data.frame(
        gene = i,
        expression = data.plot[i, ],
        group = obj.groups
      )
      df <- rbind(df, df.1)
    }
  }
  p.list <- list()
  for (i in seq_along(along.with = features)) {
    df.use <- df[df$gene == features[[i]], ]
    p <- ggplot(data = df.use, aes(x = expression, y = gene, fill = group)) +
      geom_violin(size = 1/4) +
      facet_wrap(~group, ncol = 1, strip.position = "right") +
      theme_classic() +
      scale_y_discrete(position = "top") +
      scale_x_continuous(position = "bottom", limits = c(0, NA)) +
      theme(
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none"
      )
    if (!is.null(x = levels.use)) {
      p <- p + scale_fill_manual(values = colors_all)
    }
    p.list[[i]] <- p
  }
  library(patchwork)
  p <- wrap_plots(p.list, ncol = length(x = p.list))
  return(p)
}

AnnotationPlot <- function(annotation, region) {
  annotation <-annotation
  if (is.null(x = annotation)) {
    return(NULL)
  }
  if (!inherits(x = region, what = "GRanges")) {
    region <- StringToGRanges(regions = region)
  }
  start.pos <- start(x = region)
  end.pos <- end(x = region)
  chromosome <- seqnames(x = region)
  
  # get names of genes that overlap region, then subset to include only those
  # genes. This avoids truncating the gene if it runs outside the region
  annotation.subset <- subsetByOverlaps(x = annotation, ranges = region)
  genes.keep <- unique(x = annotation.subset$gene_name)
  annotation.subset <- annotation[
    fmatch(x = annotation$gene_name, table = genes.keep, nomatch = 0L) > 0L
    ]
  
  if (length(x = annotation.subset) == 0) {
    # make empty plot
    p <- ggplot(data = data.frame())
    
  } else {
    annotation.subset <- split(
      x = annotation.subset,
      f = annotation.subset$gene_name
    )
    annotation.subset=as(annotation.subset,'GRangesList')
    p <- suppressWarnings(expr = suppressMessages(expr = ggbio::autoplot(
      object = annotation.subset,
      GRangesFilter(value = region),
      fill = "darkblue",
      size = 1/2,
      color = "darkblue",
      names.expr = "gene_name"
    )))
    p <- p@ggplot
    # extract y-axis limits and extend slightly so the label isn't covered
    y.limits <- ggplot_build(plot = p)$layout$panel_scales_y[[1]]$range$range
    p <- suppressMessages(p + ylim(y.limits[[1]], y.limits[[2]] + 0.5))
  }
  p=suppressMessages( p +
                        theme_classic() +
                        ylab("Genes") +
                        xlab(label = paste0(chromosome, " position (kb)")) +
                        xlim(start.pos, end.pos) +
                        theme(
                          axis.ticks.y = element_blank(),
                          axis.text.y = element_blank()
                        ))
  return(p)
}
LinkPlot <- function(Links, region, min.cutoff = -1) {
  if (!inherits(x = region, what = "GRanges")) {
    region <- StringToGRanges(regions = region)
  }
  chromosome <- seqnames(x = region)
  
  # extract link information
  links <- Links
  
  # if links not set, return NULL
  if (length(x = links) == 0) {
    return(NULL)
  }
  
  # subset to those in region
  links.keep <- subsetByOverlaps(x = links, ranges = region)
  
  # filter out links below threshold
  link.df <- as.data.frame(x = links.keep)
  link.df <- link.df[abs(x = link.df$score) > min.cutoff, ]
  # link.df=link.df[link.df$start>start(region)&link.df$end>end(region),]
  # remove links outside region
  # link.df <- link.df[link.df$start >= start(x = region) & link.df$end <= end(x = region), ]
  library(ggforce)
  # plot
  if (nrow(x = link.df) > 0) {
    # convert to format for geom_bezier
    link.df$group <- seq_len(length.out = nrow(x = link.df))
    df <- data.frame(
      x = c(link.df$start,
            (link.df$start + link.df$end) / 2,
            link.df$end),
      y = c(rep(x = 0, nrow(x = link.df)),
            log10(link.df$pvalue),
            rep(x = 0, nrow(x = link.df))),
      group = rep(x = link.df$group, 3),
      score = rep(link.df$score, 3),
      correlation=rep(link.df$correlation,3)
    )
    df$correlation<-factor(unfactor(df$correlation),levels = c('positive','negative'))
    p <- ggplot(data = df) +
      geom_bezier2(
        mapping = aes_string(x = "x", y = "y", group = "group", size = "score",color='correlation')
      ) +
      geom_hline(yintercept = 0, color = 'grey') +
      scale_discrete_manual(c('red','blue'))
  } else {
    p <- ggplot(data = link.df)
  }
  library(ggplot2)
  p<- p +
    theme_classic() +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank()) +
    ylab("Links") +
    xlab(label = paste0(chromosome, " position (bp)")) +
    xlim(c(start(x = region), end(x = region)))
  
  return(p)
}
PeakPlot <- function(
  peak.use,
  region,
  peaks = NULL,
  color = "dimgrey"
) {
  
  if (!inherits(x = region, what = "GRanges")) {
    region <- StringToGRanges(regions = region)
  }
  if (is.null(x = peaks)) {
    peak_bed = do.call(rbind, strsplit(x = peak.use, split = '[:-]'))
    peak_bed = as.data.frame(peak_bed)
    colnames(peak_bed)=c('chr','start','end')
    peaks=makeGRangesFromDataFrame(peak_bed)
    # mcols(x = peaks) <- md
  }
  # subset to covered range
  peak.intersect <- subsetByOverlaps(x = peaks, ranges = region)
  peak.df <- as.data.frame(x = peak.intersect)
  start.pos <- start(x = region)
  end.pos <- end(x = region)
  chromosome <- seqnames(x = region)
  
  if (nrow(x = peak.df) > 0) {
    peak.df$start[peak.df$start < start.pos] <- start.pos
    peak.df$end[peak.df$end > end.pos] <- end.pos
    peak.plot <- ggplot(
      data = peak.df
    ) +
      geom_segment(aes(x = start, y = 0, xend = end, yend = 0),
                   size = 2,
                   data = peak.df)
  } else {
    # no peaks present in region, make empty panel
    peak.plot <- ggplot(data = peak.df)
  }
  peak.plot <- peak.plot + theme_classic() +
    ylab(label = "Peaks") +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank()) +
    xlab(label = paste0(chromosome, " position (bp)")) +
    xlim(c(start.pos, end.pos))
  # remove legend, change color
  peak.plot <- peak.plot +
    scale_color_manual(values = color) +
    theme(legend.position = "none")
  
  return(peak.plot)
}

ComputeTile <- function(
  cutmatrix,
  groups,
  window = 200,
  n = 100,
  idents = NULL,
  order = "total"
) {
  # for each group, choose n cells based on total counts
  totals <- rowSums(x = cutmatrix)
  unique.groups <- unique(x = groups)
  cells.use <- vector(mode = "character")
  cell.idx <- vector(mode = "numeric")
  for (i in seq_along(along.with = unique.groups)) {
    tot.use <- totals[names(x = groups[groups == unique.groups[[i]]])]
    cell.keep <- names(x = head(x = sort(x = tot.use, decreasing = TRUE), n))
    cell.idx <- c(cell.idx, seq_along(along.with = cell.keep))
    cells.use <- c(cells.use, cell.keep)
  }
  names(x = cell.idx) <- cells.use
  cutmatrix <- cutmatrix[cells.use, ]
  
  # create sliding window sum of integration sites using RcppRoll
  # note that this coerces to a dense matrix
  smoothed <- apply(
    X = cutmatrix,
    MARGIN = 1,
    FUN = roll_sum,
    n = window,
    by = window
  )
  
  # create dataframe
  smoothed <- as.data.frame(x = smoothed)
  
  # add extra column as bin ID
  smoothed$bin <- seq_len(length.out = nrow(x = smoothed))
  smoothed <- pivot_longer(
    data = smoothed,
    cols = cells.use
  )
  
  smoothed$group <- groups[smoothed$name]
  smoothed$idx <- cell.idx[smoothed$name]
  smoothed$bin <- smoothed$bin + as.numeric(x = colnames(x = cutmatrix)[[1]])
  return(smoothed)
}

CreateTilePlot <- function(df, n, legend = TRUE) {
  # create plot
  p <- ggplot(
    data = df,
    aes_string(x = "bin", y = "idx", fill = "value")) +
    facet_wrap(
      facets = ~group,
      scales = "free_y",
      ncol = 1,
      strip.position = "left"
    ) +
    geom_raster() +
    theme_browser(legend = legend) +
    geom_hline(yintercept = c(0, n), size = 0.1) +
    ylab(paste0("Fragments (", n, " cells)")) +
    scale_fill_gradient(low = "white", high = "darkred") +
    scale_y_reverse() +
    guides(fill = guide_legend(
      title = "Fragment\ncount",
      keywidth = 1/2, keyheight = 1
    )
    ) +
    theme(
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8)
    )
  return(p)
}


CombineTracks <- function(
  plotlist,
  expression.plot = NULL,
  heights = NULL,
  widths = NULL
) {
  # remove any that are NULL
  nullplots <- sapply(X = plotlist, FUN = is.null)
  plotlist <- plotlist[!nullplots]
  heights <- heights[!nullplots]
  
  if (length(x = plotlist) == 1) {
    return(plotlist[[1]])
  }
  
  # remove x-axis from all but last plot
  for (i in 1:(length(x = plotlist) - 1)) {
    plotlist[[i]] <- plotlist[[i]] + theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.line.x.bottom = element_blank(),
      axis.ticks.x.bottom = element_blank()
    )
  }
  
  # combine plots
  if (is.null(x = heights)) {
    # set height of first element to 10x more than other elements
    n.plots <- length(x = plotlist)
    heights <- c(8, rep(1, n.plots - 1))
  } else {
    if (length(x = heights) != length(x = plotlist)) {
      stop("Relative height must be supplied for each plot")
    }
  }
  if (!is.null(x = expression.plot)) {
    # align expression plot with the first element in plot list
    p <- (plotlist[[1]] | expression.plot) +
      plot_layout(widths = widths)
    
    n <- length(x = plotlist)
    heights.2 <- heights[2:n]
    p2 <- wrap_plots(plotlist[2:n], ncol = 1, heights = heights.2)
    
    p <- p + p2 + guide_area() + plot_layout(
      ncol = 2,nrow=2, heights = c(heights[[1]], sum(heights.2)),
      guides = "collect")
  } else {
    p <- wrap_plots(plotlist, ncol = 1, heights = heights)
  }
  return(p)
}


SetIfNull <- function(x, y) {
  if (is.null(x = x)) {
    return(y)
  } else {
    return(x)
  }
}

BigwigTrack <- function(
  region,
  bigwig,
  smooth = 200,
  type = "coverage",
  y_label = "Score",
  max.downsample = 3000,
  downsample.rate = 0.1
) {
  possible_types <- c("line", "heatmap", "coverage")
  if (!(type %in% possible_types)) {
    stop(
      "Invalid type requested. Choose ",
      paste(possible_types, collapse = ", ")
    )
  }
  if (.Platform$OS.type == "windows") {
    message("BigwigTrack not supported on Windows")
    return(NULL)
  }
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    message("Please install rtracklayer. http://www.bioconductor.org/packages/rtracklayer/")
    return(NULL)
  }
  if (!inherits(x = region, what = "GRanges")) {
    stop("region should be a GRanges object")
  }
  region_data <- rtracklayer::import(
    con = bigwig,
    which = region,
    as = "NumericList"
  )[[1]]
  if (!is.null(x = smooth)) {
    region_data <- roll_mean(x = region_data, n = smooth, fill = 0L)
  }
  region_data <- data.frame(
    position = start(x = region):end(x = region),
    score = region_data,
    stringsAsFactors = FALSE
  )
  window.size = width(x = region)
  sampling <- min(max.downsample, window.size * downsample.rate)
  coverages <- slice_sample(.data = region_data, n = sampling)
  p <- ggplot(
    data = coverages,
    mapping = aes_string(x = "position", y = "score")
  )
  if (type == "line") {
    p <- p + geom_line()
  } else if (type == "heatmap") {
    # different downsampling needed for heatmap
    # cut into n bins and average within each bin
    region_data$bin <- floor(x = region_data$position / smooth)
    region_data <- group_by(region_data, bin)
    region_data <- mutate(region_data, score = mean(x = score))
    region_data <- ungroup(region_data)
    region_data <- unique(x = region_data[, c("bin", "score")])
    p <- ggplot(
      data = region_data,
      mapping = aes_string(x = "bin", y = 1, fill = "score")
    ) + geom_tile() + scale_fill_viridis_c()
  } else if (type == "coverage") {
    p <- p + geom_area()
  }
  chromosome <- as.character(x = seqnames(x = region))
  p <- p + theme_browser() +
    xlab(label = paste0(chromosome, " position (bp)")) +
    ylab(label = y_label)
  return(p)
}

