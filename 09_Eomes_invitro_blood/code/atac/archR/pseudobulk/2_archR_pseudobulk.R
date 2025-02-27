here::i_am("atac/archR/pseudobulk/2_archR_pseudobulk.R")
source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(ArchR))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--archr_directory',    type="character",    help='ArchR directory')
p$add_argument('--metadata',    type="character",    help='metadata file')
p$add_argument('--outdir',     type="character",    help='Output directory')
p$add_argument('--group_by',     type="character",    help='Metadata column to group by')
p$add_argument('--matrices',     type="character",       nargs="+",   help='Matrices to pseudobulk')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$archr_directory <- file.path(io$basedir,"processed/atac/archR")
# args$metadata <- "/bi/group/reik/ricard/data/eomes_10x_multiome/results/atac/archR/qc/sample_metadata_after_qc.txt.gz"
# args$group_by <- "cluster"
# args$matrices <- c("PeakMatrix", "GeneScoreMatrix_distal", "GeneScoreMatrix_TSS")
# args$threads <- 1
# args$outdir <- ""
## END TEST ##

# Parse arguments
dir.create(args$outdir, showWarnings=F, recursive=T)

########################
## Load cell metadata ##
########################

cell_metadata.dt <- fread(args$metadata) %>%
  .[pass_atacQC==TRUE & doublet_call==FALSE | pass_atacQC==TRUE & is.na(doublet_call)]
stopifnot(args$group_by%in%colnames(cell_metadata.dt))
cell_metadata.dt <- cell_metadata.dt[!is.na(cell_metadata.dt[[args$group_by]])]

table(cell_metadata.dt[[args$group_by]])

########################
## Load ArchR project ##
########################

# source(here::here("atac/archR/load_archR_project.R"))

setwd(args$archr_directory)

addArchRGenome("mm10")
addArchRThreads(threads = args$threads)

# Subset
ArchRProject <- loadArchRProject(args$archr_directory)[cell_metadata.dt$cell]

# Sanity checks
stopifnot(args$matrices%in%getAvailableMatrices(ArchRProject))

###########################
## Update ArchR metadata ##
###########################

cell_metadata_to_archr.dt <- cell_metadata.dt %>% 
  .[cell%in%rownames(ArchRProject@cellColData)] %>% setkey(cell)%>% .[rownames(ArchRProject@cellColData)]%>%
  as.data.frame() %>% tibble::column_to_rownames("cell")

stopifnot(all(rownames(cell_metadata_to_archr.dt) == rownames(getCellColData(ArchRProject))))
ArchRProject <- addCellColData(
  ArchRProject,
  data = cell_metadata_to_archr.dt[[args$group_by]],
  name = args$group_by,
  cells = rownames(cell_metadata_to_archr.dt),
  force = TRUE
)

table(getCellColData(ArchRProject,args$group_by)[[1]])

###################################################
## Pseudobulk into a SummarizedExperiment object ##
###################################################

se_list <- list()
for (i in args$matrices) {
  print(sprintf("Calculating pseudobulk matrix for %s",i))

  # summarise
  se_list[[i]] <- getGroupSE(ArchRProject, groupBy = args$group_by, useMatrix = i, divideN = FALSE)
  
  # rename features
  if (grepl("peak",tolower(i),ignore.case=T)) {
    rownames(se_list[[i]]) <- rowData(se_list[[i]]) %>% as.data.table %>% .[,idx:=sprintf("%s:%s-%s",seqnames,start,end)] %>% .$id
  }
  if (grepl("gene",tolower(i),ignore.case=T)) {
    rownames(se_list[[i]]) <- rowData(se_list[[i]])$name
  }

  # save
  saveRDS(se_list[[i]], file.path(args$outdir,sprintf("pseudobulk_%s_summarized_experiment.rds",i)))
}

# Save stats
to_save.dt <- data.table(table(cell_metadata.dt[[args$group_by]])) %>% 
  setnames(c("group","N")) %>% 
  setorder(group)# %>% .[,included:=group%in%groups.to.use]
fwrite(to_save.dt, file.path(args$outdir,"stats.txt"), sep="\t", quote = F)
