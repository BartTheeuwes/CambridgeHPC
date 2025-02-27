## To do:
# Give input dimensionality reduction instead of recalculating

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(miloR))
suppressPackageStartupMessages(library(patchwork))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--sce',             type="character",                               help='SingleCellExperiment file')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--atlas_metadata',  type="character",                               help='Atlas metadata file')
p$add_argument('--stage',           type="character",                               help='Stages to include')
p$add_argument('--pca_file',        type="character",                               help='Path to PCA file') 
p$add_argument('--n_neighbors',     type="integer",    default=30,                  help='Number of neighbours')
p$add_argument('--tdTom_corr',      type="character",                               help='Keep tdTom+ cells in WT samples (otherwise removed)')
p$add_argument('--prop',            type="double",     default=0.3,                 help='Proportion of cells to sample for MILO')
#p$add_argument('--remove_ExE_cells', action="store_true",                                 help='Remove ExE cells?')
p$add_argument('--outdir',          type="character",                               help='Output file')

args <- p$parse_args(commandArgs(TRUE))



#####################
## Define settings ##
#####################
here::i_am("processing/1_create_seurat_rna.R")
source(here::here("settings.R"))
source(here::here("utils.R"))
source(here::here("mapping/run/mnn/mapping_functions.R"))
test = FALSE

if(test){
## START TEST ##
    args = list()
args$sce <- io$rna.sce
args$metadata <- io$metadata
args$stage <- c('E7.5')#, 'E8.5', 'E9.5')
args$metadata <- paste0(io$basedir,"/results/rna/mapping/sample_metadata_after_mapping.txt.gz")
args$atlas_metadata <- "/rds/project/rds-SDzz0CATGms/users/bt392/atlasses/extended/sample_metadata.txt.gz"
args$pca <- NULL
args$umap <-  NULL
args$n_neighbors = 45 
args$prop = 0.15 # put at 0.15 - 0.2
args$remove_ExE_cells <- FALSE
args$tdTom_corr = NULL
args$outdir <- paste0(io$basedir,"/results/rna/MiloR/test")
## END TEST ##
}

# If passing multiple timepoints split in vector
args$stage = strsplit(args$stage, "_")[[1]] 

dir.create(args$outdir, recursive=TRUE, showWarnings = FALSE)

##########################
## Updata beeswarm code ##
##########################
plotDAbeeswarm = function (da.res, group.by = NULL, alpha = 0.1, subset.nhoods = NULL) 
{
    if (!is.null(group.by)) {
        if (!group.by %in% colnames(da.res)) {
            stop(group.by, " is not a column in da.res. Have you forgot to run annotateNhoods(x, da.res, ", 
                group.by, ")?")
        }
        if (is.numeric(da.res[, group.by])) {
            stop(group.by, " is a numeric variable. Please bin to use for grouping.")
        }
        da.res <- mutate(da.res, group_by = da.res[, group.by])
    }
    else {
        da.res <- mutate(da.res, group_by = "g1")
    }
    if (!is.factor(da.res[, "group_by"])) {
        message("Converting group.by to factor...")
        da.res <- mutate(da.res, factor(group_by, levels = unique(group_by)))
    }
    if (!is.null(subset.nhoods)) {
        da.res <- da.res[subset.nhoods, ]
    }
    da.res %>% mutate(is_signif = ifelse(SpatialFDR < alpha, 
        1, 0)) %>% mutate(logFC_color = ifelse(is_signif == 1, 
        logFC, 0)) %>% arrange(group_by) %>% mutate(Nhood = factor(Nhood, 
        levels = unique(Nhood))) %>% ggplot(aes(group_by, logFC, 
        color = logFC_color)) + scale_color_gradient2(mid='grey90') + guides(color = "none") + 
        xlab(group.by) + ylab("Log Fold Change") + ggbeeswarm::geom_quasirandom(alpha = 1) + 
        coord_flip() + theme_bw(base_size = 22) + theme(strip.text.y = element_text(angle = 0))
}

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
   .[pass_rnaQC==TRUE & doublet_call==FALSE & stage %in% args$stage] %>%
    .[,pool:=stringr::str_replace_all(sample,opts$sample2pool)]


if(args$tdTom_corr!='True'){
    sample_metadata = sample_metadata[tdTom==tdTom_corr]
}

if(test){
    sample_metadata = sample_metadata[sample(1:nrow(sample_metadata), nrow(sample_metadata)/10)]
}

###############
## Load data ##
###############

# Load RNA expression data as SingleCellExperiment object
sce <- load_SingleCellExperiment(args$sce, cells=sample_metadata$cell, normalise = TRUE)

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

# Get PCA
pca = fread(args$pca_file)[cell %in% colnames(sce)]

sce = sce[,pca$cell]
pca = pca[match(colnames(sce), cell)] %>% tibble::column_to_rownames("cell") %>% as.matrix
reducedDim(sce,"PCA") = pca
args$npcs = dim(pca)[2]

# Get UMAP
meta_atlas <- fread(args$atlas_metadata)
umap.dt <- meta_atlas %>%
    .[,c("cell","umapX","umapY")] %>%
    .[match(sce$closest.cell_mnn, cell)] %>%
    .[, cell := colnames(sce)]

umap = umap.dt[match(colnames(sce), cell)] %>% tibble::column_to_rownames("cell") %>% as.matrix
reducedDim(sce,"UMAP") = umap

########################
## Create MILO object ##
########################

# create object
Milo_sce <- Milo(sce)

colData(Milo_sce)$sample_ko = paste0(colData(Milo_sce)$sample, '_',  colData(Milo_sce)$tdTom_corr)

# Build graph
Milo_sce <- buildGraph(Milo_sce, 
                       k = args$n_neighbors, 
                       d = args$npcs, 
                       reduced.dim = "PCA")

# Identify neighbourhoods from NN-cells
# lower prop for larger datasets
Milo_sce <- makeNhoods(Milo_sce, 
                       prop = args$prop, 
                       k = args$n_neighbors, 
                       d = args$npcs, 
                       refined = TRUE, 
                       reduced_dims = "PCA")


# count cells of different samples per neighbourhood
Milo_sce <- countCells(Milo_sce, meta.data = as.data.frame(colData(Milo_sce)), sample="sample_ko")

# Calculate Neighbourhood Distances (most time consuming step)
Milo_sce <- calcNhoodDistance(Milo_sce, d=args$npcs, reduced.dim = "PCA")

# Build neighbourhood graph
Milo_sce <- buildNhoodGraph(Milo_sce)

# Save MILO object
outfile = sprintf("%s/%s_Milo_tdTomcorr%s.rds",args$outdir, paste(args$stage, collapse ='_'), args$tdTom_corr)
saveRDS(Milo_sce, outfile)

# Embryo design table
embryo_design <- data.frame(colData(Milo_sce))[,c("sample_ko", "sample", 'tdTom_corr', 'tdTom', 'stage', 'pool')]

embryo_design <- distinct(embryo_design)
rownames(embryo_design) <- embryo_design$sample_ko

# Test differential abundance per hood
da_results <- testNhoods(Milo_sce, design = ~ pool + tdTom_corr, design.df = embryo_design, reduced.dim="PCA")

# Save MILO results
outfile = sprintf("%s/%s_Milo_tdTomcorr%s_DAresults.csv",args$outdir, paste(args$stage, collapse ='_'), args$tdTom_corr)
write.csv(da_results, outfile, row.names=FALSE)


##################
## Plot results ##
##################

p1 = ggplot(da_results, aes(PValue)) + 
    geom_histogram(bins=50) + 
    theme_bw()

p2 = ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
    geom_point() +
    geom_hline(yintercept = 1) + ## Mark significance threshold (10% FDR)
    theme_bw()


## Plot single-cell UMAP
p3 = plotReducedDim(Milo_sce, dimred = "UMAP", colour_by="tdTom_corr", text_by = "celltype.mapped_mnn", 
                          text_size = 3, point_size=0.3, text_colour = "grey30") +
    scale_color_manual(values=opts$tdTom.color, name = 'tdTom') + 
    theme_void() + 
    theme(legend.position='right') +
    guides(fill="none")

## Plot neighbourhood graph
p4 = plotNhoodGraphDA(Milo_sce, da_results, layout="UMAP",alpha=0.05) 

# plot stats + graphs
outfile = sprintf("%s/%s_Milo_tdTomcorr%s_plots.pdf",args$outdir, paste(args$stage, collapse ='_'), args$tdTom_corr)
pdf(outfile, width=14, height = 7)
    p1 + p2
    p3 + p4 +
      plot_layout(guides="collect")
dev.off()

# If significant results create beeswarm plot
if(TRUE %in% unique(da_results$SpatialFDR < 0.05)){
    # Find Neighbourhood groups
    da_results <- groupNhoods(Milo_sce, da_results, max.lfc.delta = 0.5, overlap = 20)

    # Annotate hoods by celltype
    da_results <- annotateNhoods(Milo_sce, da_results, coldata_col = "celltype.mapped_mnn")
    da_results <- annotateNhoods(Milo_sce, da_results, coldata_col = "celltype_extended.mapped_mnn")

    # Annotate 'mixed' hoods when there is not one main celltype
    da_results$celltype.mapped_mnn <- ifelse(da_results$celltype.mapped_mnn_fraction < 0.4, "Mixed", da_results$celltype.mapped_mnn)
    da_results$celltype_extended.mapped_mnn <- ifelse(da_results$celltype.mapped_mnn_fraction < 0.4, "Mixed", da_results$celltype_extended.mapped_mnn)

    # Save MILO object
    outfile = sprintf("%s/%s_Milo_tdTomcorr%s_DAresults.csv",args$outdir, paste(args$stage, collapse ='_'), args$tdTom_corr)
    write.csv(da_results, outfile, row.names=FALSE)
    
    # plot DA beeswarm
    p5 = plotDAbeeswarm(da_results, group.by = "celltype.mapped_mnn")
    p6 = plotDAbeeswarm(da_results, group.by = "celltype_extended.mapped_mnn")

    # # save beeswarm
    outfile = sprintf("%s/%s_Milo_tdTomcorr%s_beeswarm_celltype.pdf",args$outdir, paste(args$stage, collapse ='_'), args$tdTom_corr)
    pdf(outfile, width=10, height = length(unique(da_results$celltype.mapped_mnn))/2.5)
        print(p5)
    dev.off()
    
        outfile = sprintf("%s/%s_Milo_tdTomcorr%s_beeswarm_celltype_extended.pdf",args$outdir, paste(args$stage, collapse ='_'), args$tdTom_corr)
    pdf(outfile, width=10, height = length(unique(da_results$celltype_extended.mapped_mnn))/2.5)
        print(p6)
    dev.off()
}







