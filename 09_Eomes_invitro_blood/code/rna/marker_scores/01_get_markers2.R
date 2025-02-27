here::i_am("rna/regression/regress_variables.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(edgeR))

BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers = 16

p <- ArgumentParser(description='')
p$add_argument('--celltype_test',    type="character",  help='celltype to test differential genes for')
p$add_argument('--celltypes',    type="character", nargs='+', help='all celltypes')
args <- p$parse_args(commandArgs(TRUE))


# args$celltypes = c("Epiblast",
#         "Primitive Streak",
#         "Caudal epiblast",
#         "PGC",
#         "Anterior Primitive Streak",
#         "Nascent mesoderm",
#         "Intermediate mesoderm",
#         "Caudal mesoderm",
#         "Lateral plate mesoderm",
#         "Limb mesoderm",
#         "Forelimb",
#         "Kidney primordium",
#         "Presomitic mesoderm",
#         "Somitic mesoderm",
#         "Posterior somitic tissues",
#         "Paraxial mesoderm",
#         "Cranial mesoderm",
#         "Anterior somitic tissues",
#         "Sclerotome",
#         "Dermomyotome",
#         "Pharyngeal mesoderm",
#         "Cardiopharyngeal progenitors",
#         "Anterior cardiopharyngeal progenitors",
#         "Allantois",
#         "Mesenchyme",
#         "YS mesothelium",
#         "Epicardium",
#         "Embryo proper mesothelium",
#         "Cardiopharyngeal progenitors FHF",
#         "Cardiomyocytes FHF 1",
#         "Cardiomyocytes FHF 2",
#         "Cardiopharyngeal progenitors SHF",
#         "Cardiomyocytes SHF 1",
#         "Cardiomyocytes SHF 2",
#         "Haematoendothelial progenitors",
#         "Blood progenitors",
#         "Erythroid",
#         "Chorioallantoic-derived erythroid progenitors",
#         "Megakaryocyte progenitors",
#         "MEP",
#         "EMP",
#         "YS endothelium",
#         "YS mesothelium-derived endothelial progenitors",
#         "Allantois endothelium",
#         "Embryo proper endothelium",
#         "Venous endothelium",
#         "Endocardium",
#         "NMPs/Mesoderm-biased",
#         "NMPs"
# )

args$outdir = file.path(io$basedir, 'results/rna/marker_scores/')
dir.create(args$outdir, recursive=TRUE, showWarnings =FALSE)

print(args$celltypes)

print(paste('celltype to test:', args$celltype_test))

extended.sce = readRDS(io$atlas.extended.sce)

colData(extended.sce)$celltype_extended_atlas = gsub(' ', '_', gsub('/', '_', colData(extended.sce)$celltype_extended_atlas)) # remove spaces from cellnames

meta = colData(extended.sce) %>% 
    as.data.table() %>% 
    .[embryo_version=='Original'] %>%
    .[celltype_extended_atlas %in% args$celltypes]

# Subset to original atlas
original.sce = extended.sce[,as.character(meta$cell)]

original.sce

# Remove genes with low total counts
counts = rowSums(assay(original.sce))
original.sce = original.sce[names(counts[counts>=30]),]

# Rename ensemble IDs to gene names in the atlas
gene_metadata <- fread(io$gene_metadata) %>% 
    .[,c("chr","ens_id","symbol")] %>%
     .[symbol!="" & ens_id%in%rownames(original.sce)] %>%
     .[!duplicated(symbol)]

original.sce <- original.sce[rownames(original.sce) %in% gene_metadata$ens_id,]
foo <- gene_metadata$symbol; names(foo) <- gene_metadata$ens_id
rownames(original.sce) <- foo[rownames(original.sce)]

original.sce

# Pseudobulk by celltype.extended + sample
summed <- aggregateAcrossCells(original.sce, 
    id=colData(original.sce)[,c("celltype_extended_atlas", "sample")],
    BPPARAM = BPPARAM)

# Filter pseudobulks
summed.filt <- summed[,summed$ncells >= 15]

# Additionally remove useless genes
summed.filt <- summed.filt[grep("*Rik|^Gm|^Mt-|^Rps|^Rpl|^Olfr",rownames(summed.filt), invert=T),]
summed.filt

# Write function to compare celltype vs all
if(nrow(colData(summed.filt)[colData(summed.filt)$celltype_extended_atlas==args$celltype_test,])>2){ # only if there actually are samples of that celltype
    DEGs_celltype = function(celltype = 'Epiblast'){
            a = Sys.time()
            print(celltype)
            tmp = summed.filt
            tmp$DEG = ifelse(tmp$celltype_extended_atlas == celltype, 'test', 'background')

            between.res <- pseudoBulkDGE(tmp,
                            label=rep("dummy", ncol(tmp)),
                            design=~factor(sample) + DEG,
                            coef="DEGtest")[[1]] %>% 
                                 as.data.table(., keep.rownames=T) %>% 
                                setnames('rn', 'gene') %>%
                                .[order(FDR)] %>% 
                                .[,sig:=ifelse(abs(logFC)>1.25 & FDR<0.05, TRUE, FALSE)] %>%
                                #.[sig==TRUE] %>%
                                .[,celltype:=celltype]
            return(between.res)
            b = Sys.time()
            print('Total time for finding markers:')
            b - a
        }

        # mclapply function for every celltype & combine into one data frame
        out = lapply(args$celltype_test, DEGs_celltype) %>%
            rbindlist()

        add_cdr = function(celltype_i = 'Epiblast', DEGs = out){

            DEGs_tmp = DEGs[celltype == celltype_i]
            tmp = original.sce[DEGs_tmp$gene,original.sce$celltype_extended_atlas == celltype_i]
            cdr = round(rowSums(counts(tmp) > 0) / ncol(tmp), 2)

            tmp.bg = original.sce[DEGs_tmp$gene,original.sce$celltype_extended_atlas != celltype_i]
            cdr.bg = round(rowSums(counts(tmp.bg) > 0) / ncol(tmp.bg), 2)

            DEGs_tmp = DEGs_tmp[,`:=`(cdr=cdr, cdr_bg = cdr.bg)]
            return(DEGs_tmp)
    }

    a = Sys.time()
    out = mclapply(args$celltype_test, add_cdr,  DEGs = out, mc.cores=1) %>% 
        rbindlist()
    b = Sys.time()
    print('Total time for finding markers:')
    b - a
}else{
    out = data.table('gene'=NA,'logFC'=NA,'logCPM'=NA,'F'=NA,'PValue'=NA, 'FDR'=NA,'sig'=NA,'celltype'=NA,'cdr'=NA,'cdr_bg'=NA)
}
# Save all significant hits
fwrite(out, sprintf('%s/celltype_extended_markers_%s.txt.gz', args$outdir, args$celltype_test))

# Try getting a score similar to Ricard
