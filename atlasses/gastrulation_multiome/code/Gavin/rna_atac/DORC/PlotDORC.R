suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(pbapply))
suppressPackageStartupMessages(library(future.apply))
suppressPackageStartupMessages(library(Matrix))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--samples',      type="character",    nargs="+",  help='Samples')
# p$add_argument('--celltypes',  type="character",    nargs="+",  help='Cell type')
p$add_argument('--distance',      type="integer",                  help='Distance')
p$add_argument('--ncores',       type="integer",      default=1,  help='Number of cores')
# p$add_argument('--denoised',     action="store_true",             help='Use denoised ATAC data?')
p$add_argument('--remove_ExE_celltypes', action="store_true",   help='Remove ExE cell types?')
p$add_argument('--test_mode',    action="store_true",             help='Test mode? subset number of cells')
p$add_argument('--outdir',       type="character",                help='Output file')
args <- p$parse_args(commandArgs(TRUE))
