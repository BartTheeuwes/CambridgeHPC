source(here::here("settings.R"))
source(here::here("utils.R"))


io$sce <- "/Users/argelagr/data/gastrulation_multiome_10x/test/processed/rna/SingleCellExperiment.rds"
io$metadata <- "/Users/argelagr/data/gastrulation_multiome_10x/test/results/atac/archR/celltype_assignment/sample_metadata_after_celltype_assignment.txt.gz"

opts$samples <- c(
  "E7.5_rep1",
  "E7.5_rep2",
  "E7.75_rep1",
  "E8.0_rep1",
  "E8.0_rep2",
  "E8.5_rep1",
  "E8.5_rep2",
  "E8.5_CRISPR_T_WT",
  "E8.5_CRISPR_T_KO",
  "E8.75_rep1",
  "E8.75_rep2"
)

opts$celltypes <- c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  "PGC",
  "Anterior_Primitive_Streak",
  "Notochord",
  "Nascent_mesoderm",
  "Mixed_mesoderm",
  "Intermediate_mesoderm",
  "Caudal_Mesoderm",
  "Paraxial_mesoderm",
  "Somitic_mesoderm",
  "Allantois",
  "ExE_mesoderm",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
  "Spinal_cord",
  "Parietal_endoderm"
)

sample_metadata <- fread(io$metadata) %>% 
  .[pass_rnaQC==TRUE & doublet_call==FALSE & sample%in%opts$samples & stage%in%opts$stages & celltype%in%opts$celltypes] %>%
  .[,c("cell", "sample", "nFeature_RNA", "mitochondrial_percent_RNA", "ribosomal_percent_RNA", "celltype","stage", "genotype", "closest.cell")]
  
sce <- readRDS(io$sce)

cells <- intersect(sample_metadata$cell,colnames(sce))
sce <- sce[,cells]
sample_metadata <- sample_metadata[cell%in%cells] %>% setkey(cell) %>% .[cells]

colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

fwrite(sample_metadata, "/Volumes/Untitled/rna/sample_metadata.txt.gz", sep="\t", quote=F, na="NA")
saveRDS(sce, "/Volumes/Untitled/rna/SingleCellExperiment.rds")

