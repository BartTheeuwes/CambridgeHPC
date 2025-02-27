#########
## I/O ##
#########

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  io$script <- "/Users/ricard/gastrulation_multiome_10x/atac/archR/chromvar/differential/archr_differential_chromvar.R"
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  io$script <- "/homes/ricard/gastrulation_multiome_10x/atac/archR/chromvar/differential/archr_differential_chromvar.R"
  io$tmpdir <- paste0(io$basedir,"/results/atac/archR/chromvar/differential/tmp")
} else {
  stop("Computer not recognised")
}
io$outdir <- paste0(io$basedir,"/results/atac/archR/chromvar/differential")


#############
## Options ##
#############

# Statistical test
opts$statistical.test <- "wilcoxon"

# Define motif annotation
opts$motif_annotation <- c("Motif_JASPAR2020_human", "Motif_cisbp")
# opts$motif_annotation <- "Motif"

# Cell types
opts$celltypes <- c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  # "PGC",
  "Anterior_Primitive_Streak",
  "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  "Mixed_mesoderm",
  "Intermediate_mesoderm",
  "Caudal_Mesoderm",
  "Paraxial_mesoderm",
  "Somitic_mesoderm",
  "Pharyngeal_mesoderm",
  "Cardiomyocytes",
  "Allantois",
  "ExE_mesoderm",
  "Mesenchyme",
  "Haematoendothelial_progenitors",
  "Endothelium",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
  "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm",
  "Visceral_endoderm",
  "ExE_endoderm",
  "ExE_ectoderm",
  "Parietal_endoderm"
)

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_atacQC==TRUE] %>%
  .[sample%in%opts$samples & celltype.predicted%in%opts$celltypes]
table(sample_metadata$celltype.predicted) %>% sort


# subset celltypes with sufficient number of cells
# if (opts$ignore.small.celltypes) {
#   opts$min.cells <- 50
#   sample_metadata <- sample_metadata %>%
#     .[,N:=.N,by=c("celltype.predicted")] %>% .[N>opts$min.cells] %>% .[,N:=NULL]
# }

# Define cell types to use
opts$celltypes <- unique(sample_metadata$celltype.predicted)# %>% head(n=3)

#########
## Run ##
#########

for (i in 1:length(opts$celltypes)) {
  for (j in i:length(opts$celltypes)) {
    for (k in opts$motif_annotation) {
      if (i!=j) {
        groupA <- opts$celltypes[[i]]
        groupB <- opts$celltypes[[j]]
        
        outfile <- sprintf("%s/%s_%s_vs_%s.txt.gz", io$outdir,k,groupA,groupB)
        if (!file.exists(outfile)) {
          # Define LSF command
          if (grepl("ricard",Sys.info()['nodename'])) {
            lsf <- ""
          } else if (grepl("ebi",Sys.info()['nodename'])) {
            lsf <- sprintf("bsub -M 7000 -n 1 -o %s/%s_%s_vs_%s.txt", io$tmpdir,k,groupA,groupB)
          }
          cmd <- sprintf("%s Rscript %s --groupA %s --groupB %s --motif_annotation %s --test %s --outfile %s", lsf, io$script, groupA, groupB, k, opts$statistical.test, outfile)
          
          # Run
          print(cmd)
          system(cmd)
        }
      }
    }
  }
}

