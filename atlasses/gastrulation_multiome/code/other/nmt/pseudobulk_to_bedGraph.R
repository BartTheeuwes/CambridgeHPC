library(data.table)
library(purrr)
library(furrr)

# reads in nmt-seq data (optionally split by celltype) then pseudobulks and saves as a bedGraph file

io <- list()


io$nmt_meta    <- "/bi/scratch/Stephen_Clark/gastrulation_data/sample_metadata.txt"
io$nmt_acc     <- "/bi/scratch/Stephen_Clark/gastrulation_data/acc/raw/scNMT/"
io$nmt_met     <- "/bi/scratch/Stephen_Clark/gastrulation_data/met/raw/scNMT/"


io$outdir      <- "/bi/scratch/Stephen_Clark/gastrulation_data/bedGraph"

dir.create(io$outdir, recursive = TRUE)

opts <- list()
opts$split_by <- "lineage10x_2"

plan(multisession, workers = parallel::detectCores())

meta <- fread(io$nmt_meta) %>% 
  .[, stage_lineage := paste0(stage, "_", get(opts$split_by))] %>% 
  .[pass_accQC == TRUE & pass_metQC == TRUE] %>% 
  split(by = "stage_lineage")



files <- future_map2(meta, names(meta), ~{
  cells <- paste0(.x[, id_acc], ".tsv.gz")
  
  files <- list(met = io$nmt_met, acc = io$nmt_acc) %>% 
    map(dir, recursive = TRUE, pattern= ".tsv", full = TRUE) %>% 
    map(~.[basename(.) %in% cells])
  
  
  met <- map(files$met, fread, select = c(1:2,5)) %>% 
    rbindlist() %>% 
    .[, .(rate = mean(rate)), .(chr, pos)] %>% 
    .[, c("start", "end", "pos") := .(pos, pos, NULL)] %>% 
    setcolorder(c("chr", "start", "end", "rate")) %>% 
    setkey(chr, start, end)
  
  if (met[1, !grepl("chr", chr)]) {
    met[, chr := paste0("chr", chr)]
    
  }
  
  met[chr=="chrMT", chr := "chrM"]
  
  setkey(met, chr, start, end)
  
  met_outfile <- paste0(io$outdir, "/CpG_met_", .y, ".bedGraph")
  print(paste("saving", met_outfile))
  fwrite(met, met_outfile, sep = "\t", col.names = FALSE, quote = FALSE)
  rm(met)
  
  
  acc <- map(files$acc, fread, select = c(1:2,5)) %>% 
    rbindlist() %>% 
    .[, .(rate = mean(rate)), .(chr, pos)] %>% 
    .[, c("start", "end", "pos") := .(pos, pos, NULL)] %>% 
    setcolorder(c("chr", "start", "end", "rate")) %>% 
    setkey(chr, start, end)
  
  if (acc[1, !grepl("chr", chr)]) {
    acc[, chr := paste0("chr", chr)]
    
  }
  
  
  acc[chr=="chrMT", chr := "chrM"]
  
  setkey(acc, chr, start, end)
  
  acc_outfile <- paste0(io$outdir, "/GpC_acc_", .y, ".bedGraph")
  print(paste("saving", acc_outfile))
  fwrite(acc, acc_outfile, sep = "\t", col.names = FALSE, quote = FALSE)
  rm(acc)
  
  c(met_outfile, acc_outfile)
})


