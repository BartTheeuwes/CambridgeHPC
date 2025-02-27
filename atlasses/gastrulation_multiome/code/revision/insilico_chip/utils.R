###################
## Load packages
###################    
suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(SingleCellExperiment)
    library(ggpubr)
    library(stringr)
    library(parallel)
    library(purrr) # might not be needed if I rewrite Ricard's code
})



###################
## General Utils
###################    

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

# Minmax normalisation
minmax.normalisation <- function(x)
{
    return((x-min(x,na.rm=T)) /(max(x,na.rm=T)-min(x,na.rm=T)))
}

# Sort absolute
sort.abs <- function(dt, sort.field) dt[order(-abs(dt[[sort.field]]))]

# Turn data table to matrix
matrix.please<-function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[[1]]
  m
}

# Annotating peaks with genes in gene window
annotate_peaks = function(gene_metadata.dt = gene_metadata.dt,
                          atac.sce = atac.sce,
                          distance = 1.5e5){
    
    # Load gene metadata
    gene_metadata <- copy(gene_metadata.dt) %>% 
      .[,chr:=as.factor(sub("chr","",chr))] %>%
      setnames("symbol","gene") %>%
      .[, c("chr","start","end","gene","ens_id","strand")]

    peakSet.dt <- data.table(peak = rownames(atac.sce)) %>%
        .[,`:=`(chr = as.factor(peak %>% str_split(':') %>% map_chr(1) %>% gsub('chr', '', .)),
                start = as.integer(peak %>% str_split(':') %>% map_chr(2) %>% str_split('-') %>% map_chr(1)),
                end = as.integer(peak %>% str_split(':') %>% map_chr(2) %>% str_split('-') %>% map_chr(2)))] %>%
        .[,c('chr', 'start', 'end', 'peak')] %>% 
      setkey(chr,start,end)

    ## Overlap 
    gene_metadata.ov <- copy(gene_metadata) %>%
      .[strand=="+",c("gene.start","gene.end"):=list(start,end)] %>%
      .[strand=="-",c("gene.start","gene.end"):=list(end,start)] %>%
      .[strand=="+",c("start","end"):=list (gene.start-gene_window, gene.end+gene_window)] %>%
      .[strand=="-",c("end","start"):=list (gene.start+gene_window, gene.end-gene_window)] %>% 
      setkey(chr,start,end)

    stopifnot((gene_metadata.ov$end-gene_metadata.ov$start)>0)

    ov <- foverlaps(
      peakSet.dt,
      gene_metadata.ov,
      nomatch = NA
    ) %>%  .[,c("start","end"):=NULL] %>%
      setnames(c("i.start","i.end"),c("peak.start","peak.end")) %>%
      .[,peak.mean:=(peak.start+peak.end)/2] %>%
      # calculate distance from the peak to the genebody
      .[,dist:=min(abs(gene.end-peak.mean), abs(gene.start-peak.mean)), by=c("gene","ens_id","peak","strand")] %>%
      .[strand=="+" & peak.mean>gene.start & peak.mean<gene.end,dist:=0] %>%
      .[strand=="-" & peak.mean<gene.start & peak.mean>gene.end,dist:=0]
    
    return(ov)   
}

###################
## Correlate TF-expr & Region-accessibility
###################

cor_TF_acc = function(rna.sce, 
                      atac.sce, 
                      TFs_filt = NULL,
                      motifmatcher.se = motifmatcher.se, 
                      motif2gene.dt = motif2gene.dt,
                      correlation_method = "pearson", 
                      remove_motifs = c("T_789"),
                      cores = detectCores()){

    if(is.null(motifmatcher.se)){
        library(motifmatchr)
        stop()
        # to do: 
        # Find motif matches if not supplied
    }

    ###################
    ## Filter TFs 
    ###################

    motifs <- intersect(colnames(motifmatcher.se),motif2gene.dt$motif)
    motifmatcher.se <- motifmatcher.se[,motifs]
    motif2gene.dt <- motif2gene.dt[motif%in%motifs]

    genes <- intersect(toupper(rownames(rna.sce)),motif2gene.dt$gene)
    rna_tf.sce <- rna.sce[str_to_title(genes),]
    rownames(rna_tf.sce) <- toupper(rownames(rna_tf.sce))
    motif2gene.dt <- motif2gene.dt[gene%in%genes]

    # Manually remove some motifs
    if(!is.null(remove_motifs)){
        motif2gene.dt <- motif2gene.dt[!motif %in% remove_motifs]
    }
        
    # Remove duplicated gene-motif pairs
    genes.to.remove <- names(which(table(motif2gene.dt$gene)>1))
    print(sprintf("Removing %d TFs that have duplicated gene-motif pairs:\n%s", length(genes.to.remove), paste(genes.to.remove, collapse=", ")))
    motif2gene.dt <- motif2gene.dt[!gene%in%genes.to.remove]
    rna_tf.sce <- rna_tf.sce[rownames(rna_tf.sce)%in%motif2gene.dt$gene]
    motifmatcher.se <- motifmatcher.se[,colnames(motifmatcher.se)%in%motif2gene.dt$motif]
    stopifnot(table(motif2gene.dt$gene)==1)

    # Sanity checks
    stopifnot(colnames(rna_tf.sce)==colnames(atac.sce))

    TFs <- rownames(rna_tf.sce)

    if(!is.null(TFs_filt)){ 
        TFs = TFs[TFs %in% TFs_filt]
    }
    
    
    ###################
    ## Correlate TF-expr & Region-accessibility
    ###################    
    message(sprintf('Correlate TF-expr & Region-accessibility for %s TFs', length(TFs)))

    # a = Sys.time()
    matrix = mclapply(TFs, function(i){
      cat(i, '\n')
      # Get motif name
      motif_i <- motif2gene.dt[gene==i,motif]

      # Get all peaks with motif
      all_peaks_i <- rownames(motifmatcher.se)[which(assay(motifmatcher.se[,motif_i],"motifMatches")==1)]

      # calculate correlations between TF expression and accessibility of motif containing peaks
      corr_output <- psych::corr.test(
        x = t(logcounts(rna_tf.sce[i,])), 
        y = t(assay(atac.sce[all_peaks_i,],"logcounts")), 
        ci = FALSE,
        method = correlation_method
      )

      # create data.table containing cor & pval
      results.dt = data.table(peak = colnames(corr_output$r),
                              cor = round(corr_output$r[1,],3),
                              pval = round(corr_output$p[1,],10)) %>%
        setnames(c('cor', 'pval'), c(paste0('cor.', i), paste0('pval.', i))) %>% 
            .[match(rownames(atac.sce), peak), ] %>% # order peaks by original order
            .[,peak:=NULL]

        return(results.dt)
    }, mc.cores=cores) %>% dplyr::bind_cols(.) %>% # combine all data.tables
        .[,peak := rownames(atac.sce)] %>%  # add peak names
        as.data.frame() %>% 
        tibble::column_to_rownames('peak') %>% as.matrix() # convert to matrix 

    # Extract correlation matrix
    cor.mtx = matrix[,grep('cor', colnames(matrix))]
    colnames(cor.mtx) = gsub('cor.', '', colnames(cor.mtx))

    # Extract Pvalue matrix
    pvalue.mtx = matrix[,grep('pval', colnames(matrix))]
    colnames(pvalue.mtx) = gsub('pval.', '', colnames(cor.mtx))
    
    output = SummarizedExperiment(
          assays = SimpleList("cor" = dropNA(cor.mtx), "pvalue" = dropNA(pvalue.mtx)),
          rowData = rowData(atac.sce)
        )
    return(output)
}

######################################
## Create virtual chip-seq library 
######################################

silico_chip = function(atac.sce = atac.sce, 
                       tf2peak_cor.se = tf2peak_cor.se,
                       motifmatcher.se = motifmatcher.se,
                       motif2gene.dt = motif2gene.dt, 
                       min_number_peaks = 50,
                       TFs_filt = NULL, 
                       remove_motifs = c("T_789"),
                       cores = detectCores()){
    
    ## Subset peaks 
    peaks <- intersect(rownames(motifmatcher.se),rownames(tf2peak_cor.se))
    print(sprintf("Number of peaks: %s",length(peaks)))

    tf2peak_cor.se <- tf2peak_cor.se[peaks,]
    atac.sce <- atac.sce[peaks,]
    motifmatcher.se <- motifmatcher.se[peaks,]
    
    # Manually remove some motifs
    if(!is.null(remove_motifs)){
        motif2gene.dt <- motif2gene.dt[!motif %in% remove_motifs]
    }

    ## Rename TFs 
    motifs <- intersect(colnames(motifmatcher.se),motif2gene.dt$motif)
    TFs <- intersect(colnames(tf2peak_cor.se),motif2gene.dt$gene)
    if(!is.null(TFs_filt)){ TFs = TFs[TFs %in% TFs_filt]}

    motif2gene_filt.dt <- motif2gene.dt[motif%in%motifs & gene%in%TFs]
    motifs <- motif2gene_filt.dt$motif
    TFs <- motif2gene_filt.dt$gene

    tmp <- TFs; names(tmp) <- motifs

    stopifnot(motif2gene_filt.dt$motif%in%colnames(motifmatcher.se))
    stopifnot(motif2gene_filt.dt$gene%in%colnames(tf2peak_cor.se))

    motifmatcher.se <- motifmatcher.se[,motifs]
    colnames(motifmatcher.se) <- tmp[colnames(motifmatcher.se)]
    tf2peak_cor.se <- tf2peak_cor.se[,TFs]
    stopifnot(colnames(motifmatcher.se)==colnames(tf2peak_cor.se))

    ## Prepare data 
    tf2peak_cor.mtx <- assay(tf2peak_cor.se,"cor")
    motifmatcher.mtx <- assay(motifmatcher.se,"motifScores")
    atac.mtx <- assay(atac.sce[peaks,],"logcounts") %>% round(3)    
    
    ######################################
    ## Create virtual chip-seq library ##
    ######################################

    print("Predicting TF binding sites...")
    stopifnot(!duplicated(TFs))
    print(sprintf("Number of TFs: %s",length(TFs)))

    virtual_chip.dt <- mclapply(TFs, function(i) {

      peaks <- names(which(abs(tf2peak_cor.mtx[,i])>0)) # we only consider chromatin activators # The 'abs' actually makes it so that all non-zero are kept (== all)

      if (length(peaks)>=min_number_peaks) {

        # calculate accessibility score
        max_accessibility_score <- apply(atac.mtx[peaks,],1,max) %>% round(2)

        # calculate correlation score
        correlation_score <- tf2peak_cor.mtx[peaks,i] %>% round(2)
        correlation_score[correlation_score==0] <- NA

        # calculate motif score
        motif_score <- motifmatcher.mtx[peaks,i]
        motif_score <- round(motif_score/max(motif_score),2)

        # calculate motif counts
        predicted_score <- correlation_score * minmax.normalisation(max_accessibility_score * motif_score)

        tmp <- data.table(
          peak = peaks, 
          correlation_score = correlation_score,
          max_accessibility_score = max_accessibility_score,
          motif_score = motif_score,
          score = round(predicted_score,2)
        ) %>% sort.abs("score") %>% 
          .[,c("peak","score","correlation_score","max_accessibility_score","motif_score")]

        to_return.dt <- tmp[!is.na(score),c("peak","score", "correlation_score","max_accessibility_score","motif_score")]  %>% .[,tf:=i] # 
        return(to_return.dt)
      }
    }, mc.cores=cores) %>% rbindlist
    
    output = list()
    output$virtual_chip.dt = virtual_chip.dt
    
    # Create Virtual ChIP-seq matrix
    virtual_chip.mtx <- virtual_chip.dt %>% 
          .[,peak:=factor(peak,levels=rownames(motifmatcher.se))] %>%
          data.table::dcast(peak~tf, value.var="score", fill=0, drop=F) %>%
          matrix.please %>% Matrix::Matrix(.)
    
    ## Update motifmatchr results using the virtual ChIP-seq library
    print("Updating motifmatchr results using the virtual ChIP-seq library...")

    motifmatcher_chip.se <- motifmatcher.se[,colnames(virtual_chip.mtx)]
    # reset motif matches
    stopifnot(rownames(virtual_chip.mtx)==rownames(motifmatcher_chip.se))
    stopifnot(colnames(virtual_chip.mtx)==colnames(motifmatcher_chip.se))
    assay(motifmatcher_chip.se,"VirtualChipScores") <- virtual_chip.mtx
    
    output$motifmatcher_chip.se = motifmatcher_chip.se
    
    return(output)
}

######################################
## ChromVAR - ChIP-seq
######################################

chromVAR_chip = function(atac.sce = atac.sce, 
                           motifmatcher_chip.se = motifmatcher_chip.se,
                           assay = 'VirtualChipScores',
                           background = bgdPeaks.se,
                           genome = BSgenome.Mmusculus.UCSC.mm10, # Only needed if background = NULL
                           positive_only = TRUE, 
                           min_chip_score = 0.15,
                           min_number_peaks = 50,
                           TFs_filt = NULL,
                           test = FALSE,
                           method = 'ArchR', # Option 'ArchR' or 'ChromVAR'
                           cores = detectCores()){
    
    if(method == 'ChromVAR'){
        background = NULL
        cat('Background peaks recalculated for method = ChromVAR \n')
        if(is.null(genome)){
            cat('Please provide genome \n')
        }
    }
    
    if (positive_only){
      print(sprintf("Number of matches before filtering negative TF binding values: %d",sum(assay(motifmatcher_chip.se,"motifMatches"))))
      assay(motifmatcher_chip.se,"motifMatches")[assay(motifmatcher_chip.se, assay)<0] <- F
      print(sprintf("Number of matches after filtering negative TF binding values: %d",sum(assay(motifmatcher_chip.se,"motifMatches"))))
    }

    print(sprintf("Number of matches before filtering based on minimum ChIP-seq score: %d",sum(assay(motifmatcher_chip.se,"motifMatches"))))
    assay(motifmatcher_chip.se,"motifMatches")[abs(assay(motifmatcher_chip.se, assay))<=min_chip_score] <- F
    print(sprintf("Number of matches after filtering based on minimum ChIP-seq score: %d",sum(assay(motifmatcher_chip.se,"motifMatches"))))

    assays(motifmatcher_chip.se) <- assays(motifmatcher_chip.se)["motifMatches"]
    
    
    ## Filter TFs 
    # Filter TFs with too few peaks
    TFs <- which(colSums(assay(motifmatcher_chip.se,"motifMatches")) >= min_number_peaks) %>% names
    TFs.removed <- which(colSums(assay(motifmatcher_chip.se,"motifMatches")) < min_number_peaks) %>% names
    
    # Subset to TFs of interest
    if(!is.null(TFs_filt)){ 
        TFs = TFs[TFs %in% TFs_filt]
    }
    if(length(TFs.removed>0)){
        cat(sprintf("%s TFs removed because they don't have enough binding sites: %s \n", length(TFs.removed), paste(TFs.removed, collapse=" ")))
    }
    
    if(test){
      TFs <- c("FOXA2","MIXL1","GATA1","EOMES","BCL11B","DLX2","FOXC1")
    }
    
    motifmatcher_chip.se <- motifmatcher_chip.se[,TFs]    
    stopifnot(rownames(atac.sce)==rownames(motifmatcher.se))  
    
    if(!is.null(background)){
        ## Load background peaks
        bgdPeaks.se = background
        tmp <- rowRanges(bgdPeaks.se)
        rownames(bgdPeaks.se) <- sprintf("%s:%s-%s",seqnames(tmp), start(tmp), end(tmp))
        bgdPeaks.se <- bgdPeaks.se[rownames(atac.sce),] # subset to same peaks as in filtered atac.sce
        bg = assay(bgdPeaks.se)
    } else{
        ## Determine background peaks
        if(is.null(genome)){
            cat('Please provide genome \n')
            stop()
        }
    
        ## Filter non-accessible peaks
        peaks = rowSums(assay(atac.sce, 'counts'))>0
        atac.sce = atac.sce[peaks,]
        motifmatcher_chip.se = motifmatcher_chip.se[peaks,]   
        
        ## Calculate background peaks
        cat('Calculating background peaks \n')
        peaks = rownames(atac.sce)
        library(GenomicRanges)
        gr <- GRanges(
            seqnames = Rle(peaks %>% str_split(':') %>% map_chr(1)),
            ranges = IRanges(start = as.numeric(peaks %>% str_split(':') %>% map_chr(2) %>% str_split('-') %>% map_chr(1)), 
                             end = as.numeric(peaks %>% str_split(':') %>% map_chr(2) %>% str_split('-') %>% map_chr(2))),
            strand = Rle(rep('*', length(peaks))))
        atac.rse = as(atac.sce, 'RangedSummarizedExperiment')
        atac.rse@rowRanges = gr
        
        ## Adding background peaks
        atac.rse <- addGCBias(atac.rse, genome = genome)
        bg <- getBackgroundPeaks(object = atac.rse)
    }

    
    ## ChromVAR ChIP-seq
    if(method == 'ArchR'){
        cat('Running ChromVAR-ChIP-seq (ArchR implementation) \n')
        
        featureDF <- data.frame(
          rowSums = rowSums(assay(atac.sce))#, 
       #   start = rowData(atac.sce)$start, # -> Check if these can be left out
       #   end = rowData(atac.sce)$end
        )

        # Compute deviations
        chromvar_deviations.se <- ArchR:::.customDeviations(
          countsMatrix = assay(atac.sce),
          annotationsMatrix = as(assay(motifmatcher_chip.se),"dgCMatrix"),
          backgroudPeaks = bg,
          expectation = featureDF$rowSums/sum(featureDF$rowSums),
          prefix = "",
          out = c("deviations", "z"),
          threads = cores,
          verbose = TRUE
        )
    } else if(method == 'ChromVAR'){
        cat('Running ChromVAR-ChIP-seq (ChromVAR implementation) \n')
        library(BiocParallel)
        #stop()
        register(MulticoreParam(cores))
        assayNames(atac.sce) <- "counts"

        # Compute deviations
        chromvar_deviations_chromvar.se <- chromVAR::computeDeviations(
          object = atac.sce,
          annotations = assay(motifmatcher_chip.se),
          background_peaks = bg
        )
        
    } else{
        print('Choose either option "ArchR" or "ChromVAR"')
    }
        
    return(chromvar_deviations.se)
}

######################################
## link TF2genes virtual chip
######################################

TF2Gene = function(rna.sce = rna.sce, 
                   motifmatcher_chip.se = motifmatcher_chip.se,
                   peak2gene.dt = peak2gene.dt,
                   gene_metadata.dt = NULL, # Only needed to find peak-genes if peak2gene.dt not supplied
                   atac.sce = NULL, # Only needed to find peak-genes if peak2gene.dt not supplied
                   filter_genes = TRUE,
                   min_chip_score = 0.15,
                   distance = 5e4,
                   cores = detectCores()){
    
    if(is.null(peak2gene.dt)){
        if(is.null(gene_metadata.dt) | is.null(atac.sce)){
            cat('Please provide gene metadata and atac.sce \n')
            stop()
        }
        cat('Annotating peaks with genes in distance window \n')
        peak2gene.dt = annotate_peaks(gene_metadata.dt = gene_metadata.dt,
                                  atac.sce = atac.sce,
                                  distance = distance)
    }
   peak2gene.dt = peak2gene.dt  %>%
                  .[dist<=distance] %>%
                  .[,peak:=sprintf("chr%s:%s-%s",chr,peak.start,peak.end)]
    
    virtual_chip.mtx = assay(motifmatcher_chip.se, 'VirtualChipScores')

    # Sanity checks
    stopifnot(length(intersect(rownames(virtual_chip.mtx),unique(peak2gene.dt$peak)))>1e5)     

    ## Link TFs to target genes using the virtual ChIP-seq 
    tf2gene_chip.dt <- mclapply(colnames(virtual_chip.mtx), function(i){
    # Select target peaks (note that we only take positive correlations into account)
    target_peaks_i <- names(which(virtual_chip.mtx[,i]>=min_chip_score))

    if (length(target_peaks_i)>=1) {
        tmp <- data.table(
          tf = i,
          peak = target_peaks_i,
          chip_score = virtual_chip.mtx[target_peaks_i,i]
        ) %>% merge(peak2gene.dt[peak %in% target_peaks_i,c("peak","gene","dist")], by="peak")
        return(tmp)
        }
    }, mc.cores=cores) %>% rbindlist        
                       
    ## Filter TFs and genes 

    TFs <- intersect(unique(tf2gene_chip.dt$tf),toupper(rownames(rna.sce)))
    genes <- intersect(unique(tf2gene_chip.dt$gene),rownames(rna.sce))
    
    if(filter_genes){
        # filter out non-informative genes
        genes <- genes[grep("*Rik|^Gm|^mt-|^Rps|^Rpl", genes, invert=T)] 
    }

    tf2gene_chip.dt <- tf2gene_chip.dt[tf%in%TFs & gene%in%genes,]

    # Fetch RNA expression matrices
    rna_tf.mtx <- logcounts(rna.sce)[str_to_title(unique(tf2gene_chip.dt$tf)),]; rownames(rna_tf.mtx) <- toupper(rownames(rna_tf.mtx))
    rna_targets.mtx <- logcounts(rna.sce)[unique(tf2gene_chip.dt$gene),]

    # Filter out lowly variable genes and TFs
    rna_tf.mtx <- rna_tf.mtx[apply(rna_tf.mtx,1,var)>=1,]
    rna_targets.mtx <- rna_targets.mtx[apply(rna_targets.mtx,1,var)>=0.1,]

    TFs <- intersect(unique(tf2gene_chip.dt$tf),rownames(rna_tf.mtx))
    genes <- intersect(unique(tf2gene_chip.dt$gene),rownames(rna_targets.mtx))
    tf2gene_chip.dt <- tf2gene_chip.dt[tf%in%TFs & gene%in%genes,]

    cat(sprintf("Number of TFs: %s \n",length(TFs)))
    cat(sprintf("Number of genes: %s \n",length(genes)))
    
    ## run regression 
    GRN_coef.dt = mclapply(genes, function(i){
        tfs <- unique(tf2gene_chip.dt[gene==i,tf])
        tfs %>% map(function(j) {
          x <- rna_tf.mtx[j,]
          y <- rna_targets.mtx[i,]
          lm.fit <- lm(y~x)
          data.frame(tf=j, gene=i, beta=round(coef(lm.fit)[[2]],3), pvalue=format(summary(lm.fit)$coefficients[2,4], digits=3))
        }) %>% rbindlist
      }, mc.cores=cores) %>% rbindlist
    
    # Add chip score back in so all info is 
    #chip_GRN = merge(tf2gene_chip.dt, GRN_coef.dt, by=c('tf', 'gene'))
    return(GRN_coef.dt)
}