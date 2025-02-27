here::i_am("snakemake/01_create_arrow.R")
source(here::here("settings.R"))

# I/O
io$output.directory <- file.path(io$basedir,"ArchR")
setwd(io$output.directory)

# Define arguments
p <- ArgumentParser(description='')
p$add_argument('--arrow_files',     type="character",  nargs='+',      help='Arrow files')
args <- p$parse_args(commandArgs(TRUE))

addArchRThreads(threads = 1) 

# Create project
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

ArrowFiles = list.files(io$output.directory, pattern ='arrow')

proj <- ArchRProject(
  ArrowFiles = args$arrow_files, 
  outputDirectory = "Project",
  copyArrows = TRUE, #This is recommened so that you maintain an unaltered copy for later usage.
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation
)

proj <- filterDoublets(ArchRProj = proj)

saveArchRProject(ArchRProj = proj)

write.table('', file=file.path(io$output.directory,'02_completed.txt'))