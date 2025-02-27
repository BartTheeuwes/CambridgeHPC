ArchRProject <- loadArchRProject(io$archR.directory)

#######################
# scRNA_INTEGRATION ###
#######################

Mapping_matrix = function(matrix, dims, features){ # MATRIX = Which matrix to use for dim reduction in mapping
    args = list()
    args$Matrix = matrix
    
    # I/O
    io$outdir <- file.path(io$output.directory,"mapping")
    dir.create(file.path(io$outdir), showWarnings = FALSE)

    #### Read RNA ####
    io$RNA_in <- file.path(io$basedir,"RNA/RangedSummarizedExperiment.rds")
    sceRNA = readRDS(io$RNA_in)
    
    #### Mapping function ######
    mapping = function(stage_select){
        # Subset RNA
               
        if(stage_select=='GD9_ExE'){ # Map GD9_ExE to normal GD9
            RNA_stage = sceRNA[,colData(sceRNA)$stage == 'GD9'] # for RangedSummarizedExperiment object
        } else{
            RNA_stage = sceRNA[,colData(sceRNA)$stage == stage_select] # for RangedSummarizedExperiment object
        }
        
        print(stage_select)
            
        
        # Subset ATAC 
        keep = rownames(ArchRProject@cellColData[ArchRProject@cellColData$stage==stage_select,]) # per stage instead of stage_mapping which includes ExE for GD9
        ATAC_stage = ArchRProject[keep]

        # ATAC dim redu
        # Iterative LSI: two iterations
        ATAC_stage <- addIterativeLSI(
          ArchRProj = ATAC_stage,
          useMatrix = args$Matrix, 
          dimsToUse = 1:dims,
          name = "IterativeLSI", 
          firstSelection = "Top",
          depthCol = "nFrags",
          iterations = 2, 
          saveIterations = FALSE,
          varFeatures = features, 
          force = TRUE
        )

        # Add clusters
        ATAC_stage <- addClusters(input = ATAC_stage,  
                                 reducedDims = "IterativeLSI",
                                 resolution = 1.5, 
                                 force=TRUE)
        clusters = paste0(stage_select, '_', getCellColData(ATAC_stage)$Clusters)

        # Add imputation 
        ATAC_stage <- addImputeWeights(ATAC_stage)

        # Integrate
        ATAC_stage <- addGeneIntegrationMatrix(
            ArchRProj = ATAC_stage, 
            useMatrix = "GeneScoreMatrix",
            matrixName = paste0("GeneIntegrationMatrix_", args$Matrix),
            reducedDims = "IterativeLSI",
            seRNA = RNA_stage,
            addToArrow = TRUE,
            force= TRUE,
            groupRNA = "celltype",
            nameCell = paste0("predictedCell_Un_", args$Matrix),
            nameGroup = paste0("predictedGroup_celltype_", args$Matrix),
            nameScore = paste0("predictedScore_celltype_", args$Matrix))

        # Extract mapping
        meta = as.data.frame(ATAC_stage@cellColData[, c(paste0("predictedCell_Un_", args$Matrix), paste0("predictedGroup_celltype_", args$Matrix), paste0("predictedScore_celltype_", args$Matrix))])
        meta$stage_clusters = clusters
        meta$cell = rownames(meta)
        return(meta)
    }

    mapped = suppressWarnings(lapply(unique(ArchRProject@cellColData$stage), mapping)) # for RangedSummarizedExperiment object
    mapped = rbindlist(mapped)

    write.csv(mapped, file.path(io$outdir, paste0('celltype_mapping_', args$Matrix, '.csv')), row.names=FALSE)

    ArchRProject@cellColData = cbind(ArchRProject@cellColData, mapped)
    saveArchRProject(ArchRProj = ArchRProject)
}


###########
#DIM REDU##
###########

dim_redu = function(matrix, nfeatures, batch.variable, ndims, n_neighbors, min_dist, colour_by, final_run){
    outdir <- file.path(io$basedir, 'ArchR/dimensionality_reduction')
    seed <- 42

    dir.create(outdir, showWarnings = FALSE)

    # Options
    opts$lsi.iterations = 2
    opts$lsi.cluster.resolution = 1


    ###########################
    ## Latent Semantic Index ##
    ###########################

    # Iterative LSI: two iterations
    ArchRProject <- addIterativeLSI(
      ArchRProj = ArchRProject,
      useMatrix = matrix, 
      name = matrix, 
      dimsToUse = 1:ndims,
      firstSelection = "Top",
      depthCol = "nFrags",
      iterations = opts$lsi.iterations, 
      saveIterations = FALSE,
      varFeatures = nfeatures, 
      force = TRUE,
      seed = seed
    )

    ############################
    ## LSI + Batch correction ##
    ############################

    if (length(batch.variable)>0) {
      print(sprintf("Applying batch correction for variable: %s",  batch.variable))
      outfile <- sprintf("%s/lsi_%s_nfeatures%d_dims%d_batchcorrection_by_%s.txt.gz",outdir, matrix, nfeatures, ndims,  paste(batch.variable,collapse="-"))

      # Harmony
        ArchRProject <- addHarmony(
          ArchRProj = ArchRProject,
          reducedDims = matrix,
          name = paste0(matrix, "_Harmony"),
          groupBy = "stage",
          force = TRUE)
        lsi.dt <- getReducedDims(ArchRProject, paste0(matrix, "_Harmony")) %>% round(3) %>% 
          as.data.table(keep.rownames = T) %>% setnames("rn","cell")

    } else {
      outfile <- sprintf("%s/lsi_%s_nfeatures%d_ndims%d.txt.gz",outdir, matrix, nfeatures, ndims)
      lsi.dt <- getReducedDims(ArchRProject, matrix) %>% round(3) %>% 
        as.data.table(keep.rownames = T) %>% setnames("rn","cell")
    }

    # Save LSI coordinates on final run
    if(final_run){
        fwrite(lsi.dt, outfile)
    }
    


    ################
    ## clustering ##
    ################
    # Add clusters
    redudim = matrix
    if (length(batch.variable)>0) {
        redudim = paste0(matrix, "_Harmony")
        }

    if(final_run){
 
        ArchRProject <- addClusters(input = ArchRProject, 
                                    name= paste0('Clusters_', matrix),
                                    reducedDims = redudim,
                                    resolution = opts$lsi.cluster.resolution, 
                                    force=TRUE,
                                    seed = seed)
    }
    
    sample_metadata = as.data.table(getCellColData(ArchRProject), keep.rownames=TRUE) %>% setnames('rn', 'cell')


    ##########
    ## UMAP ##
    ##########

    # i <- args$n_neighbors[1]
    # j <- args$min_dist[1]
    for (i in n_neighbors) {
      for (j in min_dist) {


        # Run UMAP
        ArchRProject <- addUMAP(
          ArchRProj = ArchRProject, 
          reducedDims = redudim,
          name = paste0("UMAP_", matrix),
          metric = "cosine",
          nNeighbors = i, 
          minDist = j, 
          seed = seed,
          saveModel = FALSE,
          force = TRUE
        )

        # Fetch UMAP coordinates
        umap.dt <- getEmbedding(ArchRProject,paste0("UMAP_", matrix)) %>%
          round(2) %>%
          as.data.table(keep.rownames = T) %>%
          setnames(c("cell","umap1","umap2"))

        # Save UMAP coordinates
        outfile <- sprintf("%s/umap_%s_nfeatures%d_ndims%d_nn%s_dist%s.txt.gz",outdir, matrix, nfeatures, ndims,i,j)
        if (length(batch.variable)>0) {
              outfile <- sprintf("%s/umap_%s_nfeatures%d_ndims%d_nn%s_dist%s_batchcorrection_by_%s.txt.gz",outdir, matrix, nfeatures,ndims,i,j, paste(batch.variable,collapse="-"))
            }
        
        # save UMAP coordinates only on final run
        if(final_run){
           fwrite(umap.dt, outfile)
        }        
        

        # Plot
        to.plot <- umap.dt %>%
          merge(sample_metadata,by="cell")

        for (k in colour_by) {

          # log10 large numeric values
          if (is.numeric(to.plot[[k]])) {
            if (max(to.plot[[k]],na.rm=T) - min(to.plot[[k]],na.rm=T) > 1000) {
              to.plot[[k]] <- log10(to.plot[[k]]+1)
              to.plot %>% setnames(k,paste0(k,"_log10")); k <- paste0(k,"_log10")
            }
          }

          p <- ggplot(to.plot, aes_string(x="umap1", y="umap2", col=k)) +
            geom_point(size=1, alpha=0.8) +
            theme_void() 

          # Define colormap
          if (is.numeric(to.plot[[k]])) {
            p <- p + scale_color_gradientn(colours = terrain.colors(10))
          }

#           if (grepl("celltype",k)) {
#             p <- p + scale_fill_manual(values=opts$celltype.colors) +
#               theme(
#                 legend.position="none",
#                 legend.title=element_blank()
#               )
#           }

         if (is.character(to.plot[[k]])) {
                p <- p + scale_color_manual(values = as.vector(c(ArchRPalettes[[3]],ArchRPalettes[[5]], ArchRPalettes[[8]],ArchRPalettes[[9]]))[1:length(unique(to.plot[[k]]))])
              }
              
          if (grepl("stage",k)) {
            p <- p + scale_color_manual(values=opts$stage.colors)
          }
           
           

        outfile <- sprintf("%s/umap_%s_nfeatures%d_ndims%d_nn%s_dist%s_%s.pdf",outdir, matrix, nfeatures, ndims,i,j, k)

        if (length(batch.variable)>0) {
              outfile <- sprintf("%s/umap_%s_nfeatures%d_ndims%d_nn%s_dist%s_batchcorrection_by_%s_%s.pdf",outdir, matrix, nfeatures,ndims,i,j, paste(batch.variable,collapse="-"), k)
            }
        
         if (grepl("predictedGroup_celltype",k)){
              pdf(outfile, width=14, height=5)
              print(p)
              dev.off()
          }else{
                pdf(outfile, width=7, height=5)
                print(p)
                dev.off()
             }
        }
        }
      }
    
    
    # save object only on final run
    if(final_run){
        saveArchRProject(ArchRProject) 
        return(ArchRProject)
    }
}