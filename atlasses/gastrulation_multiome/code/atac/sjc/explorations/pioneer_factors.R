library(Seurat)
library(Signac)
library(purrr)
library(data.table)
library(future)
library(chromVAR)

library(motifmatchr)
library(JASPAR)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)

library(ggplot2)
library(cowplot)
library(ggrepel)



source(here::here("settings.R"))




io$signac             <- file.path(io$rawdata, "/processed/atac/signac/signac_peaks.rds")
io$outdir             <- file.path(io$rawdata, "/processed/atac/signac/motifs/")
io$plots_out          <- "/bi/home/clarks/plots/10X_multiome/signac/footprints"


opts$min_cells_exp    <- 0.001 # remove genes (of TFs) which are expressed in few cells

dir.create(io$outdir, recursive = TRUE)


opts$factors     <- c("FOXA1", "FOXA2", "GATA4", "GATA6", "NEUROD1", "ASCL1", 
                      "SP8", "POU5F1", "POU3F1", "TFAP2A", "SOX2", "SOX3", "HNF1B",
                      "GATA1", "ELF3", "SOX17", "TWIST1", "HAND1", "MEIS1", "LEF1",
                      "TEAD4", "ZNF148", "KLF4", "KLF5", "HNF1A", "TEAD2", "KLF16")





plan("multiprocess", workers = 8)
options(future.globals.maxSize = opts$mem * 1024 ^ 3)
plan()

sample_metadata <- fread(io$metadata)

sample_metadata[, .N, celltype.mapped][order(-rank(N))]


# extract position frequency matrices for the motifs
pfm <- getMatrixSet(
  x = JASPAR,
  opts = list(species = 9606, all_versions = FALSE)
)

#pfm <- getMatrixSet(x = JASPAR, opts = list(all_versions = FALSE))
#pfm

pfm_names <- as.data.table(name(pfm), keep.rownames = TRUE) %>%
  setnames(c("motif", "tf")) %>%
  .[, tf := gsub("::|\\(|\\)|\\.", "_", tf) %>% gsub("_$", "", .)] %>%
  .[, rows := motif] 

pfm_select <- pfm_names[tf %in% opts$factors, .(motif, tf)] 

# load RNAseq
# first need to rename rna cells to match atac cells
seurat <- readRDS(io$seurat)

rna_cells <- colnames(seurat) 
atac_cells <- gsub("-1", "", rna_cells) %>% 
  gsub("E7.5_rep1", "rep1_L001_multiome", .) %>% 
  gsub("E7.5_rep2", "rep2_L002_multiome", .) %>% 
  gsub("E8.5_rep1", "multiome1", .) %>% 
  gsub("E8.5_rep2", "multiome2", .)

# quick check
x <- sample(1:length(rna_cells), 10)
data.table(rna_cells[x], atac_cells[x])


rna <- seurat@assays$RNA@data
colnames(rna) <- atac_cells
# subset to only include TFs in our motif database
genes <- toupper(rownames(rna))
tfs <- strsplit(pfm_select$tf, "_") %>% unlist()
sub <- genes %in% tfs
table(sub)
rna <- rna[sub, ]

rna <- as.data.table(rna, keep.rownames = "gene")  %>% 
  melt(id.vars = "gene", value.name = "exp", variable.name = "cell") %>% 
  .[, tf := toupper(gene)] %>%
  .[, cells_exp := sum(exp>0)/.N, tf] %>%
  .[cells_exp > opts$min_cells_exp] %>% 
  dcast(cell ~ gene, value.var = "exp")

tfs <- colnames(rna)[2:ncol(rna)]

# exact matches only here - might want to change to pick up multi-motifs???

pfm_select <- pfm_select[tf %in% toupper(tfs)]



# load acc data

signac <- readRDS(io$signac)
signac










signac <- AddMotifs(
  object = signac,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm[pfm_select$motif]
)
signac



# generate peak x motif matrix (if applicable)
# then run chromVar
opts$find_overlaps <- FALSE
if (opts$find_overlaps) {
  motifmat <- signac@assays$peaks@motifs@data %>% 
    as.data.table(keep.rownames = "peak") %>% 
    setnames(pfm_select$motif, pfm_select$tf) %>% 
    .[, overlap := get(opts$factor1) * get(opts$factor2)] %>% 
    .[overlap == 1, c(opts$factor1, opts$factor2) := .(0, 0)] %>% 
    as.matrix(rownames = "peak") %>% 
    as("dgCMatrix")
  
  
  signac <- RunChromVAR(
    object = signac,
    genome = BSgenome.Mmusculus.UCSC.mm10,
    motif.matrix = motifmat
  )
} else {
  signac <- RunChromVAR(
    object = signac,
    genome = BSgenome.Mmusculus.UCSC.mm10
  )
}


signac

DefaultAssay(signac) <- "chromvar"


saveRDS(signac, file.path(io$outdir, "signac_pioneerTFs.rds"))
signac[["peaks"]] <- NULL
signac
saveRDS(signac, file.path(io$outdir, "signac_pioneerTFs.rds"))








motifs <- as.data.table(signac@assays$chromvar@data, keep.rownames = "motif") %>% 
  melt(id.vars = "motif", value.name = "acc", variable.name = "cell") %>% 
  merge(pfm_select, by = "motif") %>% 
  dcast(cell ~ tf, value.var = "acc")

comb <- merge(rna, motifs, by = c("cell"))# , "tf"))
comb
colnames(comb)

comb[, lm(FOXA1 ~ Foxa1 + Gata4) %>% summary]

ggplot(comb, aes(Foxa1, FOXA1)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_cowplot()

cols <- colnames(comb)[2:ncol(comb)]
genes <- cols[1 : (length(cols)/2)]
motifs <- cols[(length(cols)/2 +1) : length(cols)]

gene_v_motif <- map2(genes, motifs, ~{
  plot <- ggscatter(comb, .x, .y, 
                    add = "reg.line",  # Add regressin line
                    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                    conf.int = TRUE, # Add confidence interval
                    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                    cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
  )
  cors <- comb[, cor.test(get(.x), get(.y))]
  
  list(data.table(R = cors$estimate, p = cors$p.value, tf = .x),
       plot)
})

cors <- map(gene_v_motif, 1) %>% 
  rbindlist()

fwrite_tsv(cors, file.path(io$outdir, "tfs_motif_correlations.tsv.gz"))


fwrite_tsv(comb, file.path(io$outdir, "tfs_motifs.tsv.gz"))
comb=fread(file.path(io$outdir, "tfs_motifs.tsv.gz"))
print("stop")
stop()
cor_test <- function(x, y, ...){
  i <- !is.na(x) & !is.na(y)
  x <- x[i]
  y <- y[i]
  
  
  if (sum(i) < 10) {
    c <- list(p.value = NA, estimate = NA)
  } else {
    c <- cor.test(x, y, ...)
  }
    
  
  list(
    p = c$p.value,
    r = c$estimate,
    mean_x = mean(x),
    sd_x = sd(x),
    mean_y = mean(y),
    sd_y = sd(y),
    N = sum(i)
    
       )
}

cors <- comb[, cor_test(exp, acc), tf]


outfile <- file.path(io$outdir, "tf_motif_cors.tsv.gz")
fwrite_tsv(cors, outfile)

cors <- fread(outfile)

cors[, padj := p.adjust(p, method = "fdr")]

ggplot(cors, aes(r, -log10(padj))) +
  geom_point() +
  geom_text_repel(data = cors[r>0.3], aes(label = tf)) +
  theme_cowplot()



# top hits

cors[order(-rank(r))][!is.na(r)][1:30, tf]

# take a single TF gene and correlate it's expression with all motifs
motif <- comb[, .(cell, tf, acc)] %>% unique()

x <- "GATA4"



multiway <- comb[tf == x, .(cell, exp)] %>% 
  unique() %>% 
  merge(motif, by = "cell", allow.cartesian = TRUE) %>% 
  .[, cor_test(exp, acc, method = "spearman"), tf] %>% 
  .[, padj := p.adjust(p, "fdr")]

ggplot(multiway, aes(r, -log10(padj))) +
  geom_point() +
  geom_text_repel(data = multiway[r>0.15], aes(label = tf), size = 2, max.overlaps = 50) +
  theme_cowplot()

multiway[order(-rank(r))][1:20]
