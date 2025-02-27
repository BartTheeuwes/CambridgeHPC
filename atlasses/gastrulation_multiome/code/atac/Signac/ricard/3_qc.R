# NOTE: I HAVE DONE THIS WITH ARCHR, NO NEED TO REPEAT QC, JUST EXTRACT METADATA

#################
## Description ##
#################

# makes Seurat object from accesibility data using merged fragment file
# data is quantified over bins so as to avoid using CellRanger peaks which don't
# match between samples

####################
## Load libraries ##
####################

library(Seurat)
library(Signac)
library(future)

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
} else {
  stop("Computer not recognised")
}


# I/O
io$signac                <- file.path(io$basedir, "/processed/atac/signac/signac.rds")

# Options
opts$cores               <- 8
opts$mem                 <- 25 # GB
opts$max_nCount          <- 1e6
opts$min_nCount          <- 1e3
opts$max_nucleosome      <- 2
opts$min_TSS             <- 1



# Multithreading
plan("multiprocess", workers = opts$cores)
options(future.globals.maxSize = opts$mem * 1024 ^ 3)
plan()

###############
## Load data ##
###############

# Load Signac object
signac <- readRDS(io$signac)
signac


#############################
## Calculate QC statistics ##
#############################

# compute nucleosome signal score per cell
signac <- NucleosomeSignal(object = signac)

# compute TSS enrichment score per cell
signac <- TSSEnrichment(object = signac, fast = FALSE)


signac$high.tss <- ifelse(signac$TSS.enrichment > 2, 'High', 'Low')
tss <- TSSPlot(signac, group.by = 'high.tss') + NoLegend()
save_plot(paste0(io$plots_out, "/tss.pdf"), tss)

# signac$nucleosome_group <- ifelse(signac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
# fraghist <- FragmentHistogram(object = signac, group.by = 'nucleosome_group')
# save_plot(paste0(io$plots_out, "/frag_histogram.pdf"), fraghist)


##########
## Plot ##
##########

violins <- VlnPlot(
  object = signac,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

save_plot(paste0(io$plots_out, "/violins.pdf"), violins, base_height = 5, base_width = 8)

############
## Filter ##
############

# Subset cells from the Signac object
signac <- subset(
  x = signac,
  subset = nCount_ATAC < opts$max_nCount &
    nCount_ATAC > opts$min_nCount &
    nucleosome_signal < opts$max_nucleosome &
    TSS.enrichment > opts$min_TSS
)
signac


# plot pass/fails
# qc <- sample_metadata[, .(cell, pass_rnaQC, pass_sigQC = FALSE)]
# qc[cell %in% colnames(signac), pass_sigQC := TRUE]
# qc[, pass_bothQC := as.logical(pass_rnaQC * pass_sigQC)]
# qc <- melt(qc, id.vars = "cell", value.name = "passQC", variable.name = "type")

# p <- ggplot(qc, aes(type, fill = passQC)) + 
#   geom_bar() +
#   theme_cowplot() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
  
# p
# save_plot(paste0(io$plots_out, "/passQC.pdf"), p)

############################
## Update sample metadata ##
############################

##########
## Save ##
##########

saveRDS(signac, io$signac)

