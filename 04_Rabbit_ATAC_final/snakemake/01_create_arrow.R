here::i_am("snakemake/01_create_arrow.R")
source(here::here("settings.R"))

# I/O
io$output.directory <- file.path(io$basedir,"ArchR")
dir.create(file.path(io$output.directory), showWarnings = FALSE)

setwd(io$output.directory)

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--sample',        type="character",                               help='input sample')
p$add_argument('--min_fragments',          type="integer",  default=2500,   help='Minimum fragments')
p$add_argument('--min_tss_score',       type="integer",  default=2,  help='Minimum TSS score')
p$add_argument('--outdir',          type="character",                               help='Output directory')

args <- p$parse_args(commandArgs(TRUE))

args$test = FALSE
if(args$test){
args$sample = 'BGRGP1'
args$min_fragments = 2500
args$min_tss_score = 2
args$outdir = io$output.directory
}

# ArchR options
addArchRThreads(threads = 1) 

#genomeAnnotation = readRDS(file.path(io$basedir, 'genomeAnnotation.rds'))
geneAnnotation = readRDS(file.path(io$basedir, 'geneAnnotation_new.rds'))

library(BSgenome.Ocuniculus.NCBI.oryCun2)
genomeAnnotation = createGenomeAnnotation(
  genome = BSgenome.Ocuniculus.NCBI.oryCun2,
  chromSizes = NULL,
  blacklist = NULL,
  filter = FALSE,
  filterChr = c("chrM")
)
# subset annotation bc otherwise ArchR takes ages to run
keep = fread('/rds/project/rds-SDzz0CATGms/users/bt392/04_Rabbit_ATAC_final/data/BGRGP1_chrs.txt')[[1]]
genomeAnnotation$chromSizes = genomeAnnotation$chromSizes[genomeAnnotation$chromSizes@seqnames@values %in% keep]
genomeAnnotation$chromSizes = genomeAnnotation$chromSizes[1:500]
geneAnnotation$genes = geneAnnotation$genes[as.vector(seqnames(geneAnnotation$genes)) %in% genomeAnnotation$chromSizes@seqnames@values]
geneAnnotation$exons = geneAnnotation$exons[as.vector(seqnames(geneAnnotation$exons)) %in% genomeAnnotation$chromSizes@seqnames@values]
geneAnnotation$TSS = geneAnnotation$TSS[as.vector(seqnames(geneAnnotation$TSS)) %in% genomeAnnotation$chromSizes@seqnames@values]


fragment_file_path = paste0(io$basedir, '/data/', args$sample, '_fragments.tsv.gz')

#Create Arrow File
ArrowFiles <- createArrowFiles(
  inputFiles = fragment_file_path,
  sampleNames = paste0('rabbit_', args$sample),
    
  minTSS = args$min_tss_score, #Dont set this too high because you can always increase later
  minFrags = args$min_fragments , 
    
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
    
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation
)


# Calculate doublet scores
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 15, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)
