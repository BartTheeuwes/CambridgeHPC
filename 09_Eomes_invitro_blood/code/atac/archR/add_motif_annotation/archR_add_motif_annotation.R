suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(argparse))

here::i_am("atac/archR/gene_scores/add_GeneScore_matrices.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',    type="character",    help='metadata file')
p$add_argument('--peakset',    type="character",    help='peakset location file')
# p$add_argument('--outdir',    type="character",    help='Output directory')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$metadata <- "/bi/group/reik/ricard/data/gastrulation_multiome_10x/results/atac/archR/qc/sample_metadata_after_qc.txt.gz"
# args$threads <- 1
## END TEST ##

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

########################
## Load ArchR Project ##
########################

source(here::here("atac/archR/load_archR_project.R"))

addArchRThreads(threads = args$threads)

##########################
## Add motif annotation ##
##########################

# cisbp (stringent threshold)
ArchRProject <- addMotifAnnotations(
  ArchRProject,
  motifSet = "cisbp",
  annoName = "Motif_cisbp",
  cutOff = 5e-05,
  width = 7,
  force = TRUE
)

# # JASPAR2020 human (stringent) - not working for weird reason!
# ArchRProject <- addMotifAnnotations(
#   ArchRProject, 
#   motifSet = "JASPAR2020",      
#   collection = "CORE",  
#   species = "Homo sapiens", # mus musculus only has ~150 motifs in JASPAR2020
#   cutOff = 5e-05,   
#   name = "Motif_JASPAR2020",
#   force = TRUE
# )

# JASPAR2020 human (stringent) - manual fix
args <- list(species = 'Homo sapiens', collection = 'CORE')
motifs <- TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020, args)
obj <- .summarizeJASPARMotifs(motifs)
motifs <- obj$motifs
motifSummary <- obj$motifSummary

ArchRProject = addMotifAnnotations(
  ArchRProj = ArchRProject,
  motifPWMs = TFBSTools::toPWM(motifs),
  motifSet = 'custom',
  annoName = "Motif_JASPAR2020",
  species = NULL,
  collection = "CORE",
  cutOff = 5e-05,
  width = 7,
  force = TRUE
)


################################
## Save peakAnnotation object ##
################################

saveRDS(ArchRProject@peakAnnotation, sprintf("%s/Annotations/peakAnnotation.rds",io$archR.directory))
