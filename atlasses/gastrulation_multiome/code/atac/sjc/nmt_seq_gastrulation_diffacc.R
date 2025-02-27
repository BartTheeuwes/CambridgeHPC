library(data.table)
library(purrr)
library(furrr)

source(here::here("settings.R"))


io$diffacc     <- file.path(io$rawdata, "results/atac/signac/diffacc_each_vs_all.tsv.gz")
#io$peaks       <- file.path(io$rawdata, "processed/atac/archR/peakSet/peakSet.rds")
io$nmt_meta    <- "/bi/scratch/Stephen_Clark/gastrulation_data/sample_metadata.txt"
io$nmt_acc     <- "/bi/scratch/Stephen_Clark/gastrulation_data/acc/raw/scNMT/"
io$nmt_met     <- "/bi/scratch/Stephen_Clark/gastrulation_data/met/raw/scNMT/"


io$nmt_out     <- file.path(io$rawdata, "public_datasets/Argelaguet_2019/nmt_diffacc_each_vs_all.tsv.gz")

dir.create(dirname(io$nmt_out), recursive = TRUE)

opts$stage_lineages <- c(
  "E4.5_Epiblast", 
  "E5.5_Epiblast",
  "E6.5_Epiblast",
  "E7.5_Epiblast",
  "E7.5_Ectoderm",
  "E7.5_Endoderm",
  "E7.5_Mesoderm"
  )



meta <- fread(io$nmt_meta) %>% 
  .[, stage_lineage := paste0(stage, "_", lineage10x_2)] %>% 
  .[pass_accQC == TRUE & pass_metQC == TRUE & stage_lineage %in% opts$stage_lineages]  
  
cells <- paste0(meta[, id_acc], ".tsv.gz")

files <- list(met = io$nmt_met, acc = io$nmt_acc) %>% 
  map(dir, recursive = TRUE, pattern= ".tsv", full = TRUE) %>% 
  map(~.[basename(.) %in% cells])

map(files, length)
length(cells)

diffacc <- fread(io$diffacc)

# peaks <- readRDS(io$peaks) %>% 
#   as.data.table() %>% 
#   .[, anno := strsplit(GroupReplicate, "\\._.") %>% map_chr(1)] %>% 
#   .[, .(chr = seqnames, start, end, anno, id = paste0(anno, "_", idx))] %>% 
#   setkey(chr, start, end)

# filter to only include open peaks then parse for overlaping
diffacc <- diffacc[p_val_adj < 0.05 & avg_logFC > 0.5] %>%
  .[, c("chr", "start", "end") := strsplit(locus, "-") %>% purrr::transpose() %>% map(unlist)] %>%
  .[, c("start", "end") :=.(as.numeric(start), as.numeric(end))] %>%
  .[, anno := paste0(celltype1)] %>% #, "_vs_", celltype2)] %>%
  split(by = "anno") %>%
  map(~.[, id := paste0(anno, "_", .I)]) %>%
  rbindlist() %>%
  .[, .(chr, start, end, id)] %>%
  setkey(chr, start, end)



#options(future.globals.maxSize = object.size(diffacc) * 2)
cores <- availableCores()
cores
plan(multicore, workers = cores)


dt <- future_map2(files$met, files$acc, ~{
  cellmet <- gsub(".tsv.gz", "", basename(.x))
  
  met <- fread(.x) %>% 
    .[, chr := paste0("chr", chr)] %>% 
    .[, .(chr, start = pos, end = pos, rate)] %>% 
    .[, mean_met := mean(rate)] %>% 
    setkey(chr, start, end) %>% 
    foverlaps(diffacc, nomatch = 0L) %>% 
    .[, .(met = mean(rate), nmet = .N, cell = cellmet), .(chr, start, end, id, mean_met)]
  
  cellacc <- gsub(".tsv.gz", "", basename(.y))
  
  acc <- fread(.y) %>% 
    .[, chr := paste0("chr", chr)] %>% 
    .[, .(chr, start = pos, end = pos, rate)] %>% 
    .[, mean_acc := mean(rate)] %>% 
    setkey(chr, start, end) %>% 
    foverlaps(diffacc, nomatch = 0L) %>% 
    .[, .(acc = mean(rate), nacc = .N, cell = cellacc), .(chr, start, end, id, mean_acc)]
  
  merge(met, acc, by = c("cell", "chr", "start", "end", "id"), all = TRUE)
  
}) %>% 
  rbindlist()

fwrite(dt, io$nmt_out, sep = "\t", na = "NA", quote = FALSE)






