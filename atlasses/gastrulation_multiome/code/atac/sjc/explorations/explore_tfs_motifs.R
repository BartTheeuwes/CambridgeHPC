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


fwrite_tsv <- partial(fwrite, na = "NA", sep = "\t", quote = FALSE)

io$signac             <- file.path(io$rawdata, "/processed/atac/signac/signacMotifs_human.rds")
io$outdir             <- file.path(io$rawdata, "/processed/atac/signac/motifs/")
io$seurat             <- file.path(io$basedir, "/processed/rna/seurat.rds")
io$plots_out          <- "/bi/home/clarks/plots/10X_multiome/signac/footprints"


opts$min_cells_exp    <- 0.01 # remove genes (of TFs) which are expressed in few cells

dir.create(io$outdir, recursive = TRUE)


seurat <- readRDS(io$seurat)
seurat



signac <- readRDS(io$signac)
signac[["peaks"]] <- NULL
signac

sample_metadata <- fread(io$metadata)

motif_names <- copy(signac[["chromvar"]][["tf"]]) %>%
  setDT(keep.rownames = "motif")

# match up cells 
# first need to rename rna cells to match atac cells


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
sub <- genes %in% motif_names$tf
table(sub)
rna <- rna[sub, ]
rna <- as.data.table(rna, keep.rownames = "gene") %>% 
  melt(id.vars = "gene", value.name = "exp", variable.name = "cell") %>% 
  .[, tf := toupper(gene)] %>% 
  .[, cells_exp := sum(exp>0)/.N, tf] %>% 
  .[cells_exp > opts$min_cells_exp]

motifs <- as.data.table(signac@assays$chromvar@data, keep.rownames = "motif") %>% 
  melt(id.vars = "motif", value.name = "acc", variable.name = "cell") %>% 
  merge(motif_names, by = "motif", all.x = TRUE)

comb <- merge(rna, motifs, by = c("cell", "tf"))
comb=fread(file.path(io$outdir, "tfs_motifs.tsv.gz"))

fwrite_tsv(comb, file.path(io$outdir, "tfs_motifs.tsv.gz"))

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

cors[order(-rank(r))][1:20]

ggplot(cors, aes(r, -log10(padj))) +
  geom_point() +
  geom_text_repel(data = cors[r>0.3], aes(label = tf)) +
  theme_cowplot()



# top hits

cors[order(-rank(r))][!is.na(r)][1:30, tf]

# take a single TF gene and correlate it's expression with all motifs
motif <- comb[, .(cell, tf, acc)] %>% unique()

x <- "HNF4A"



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
