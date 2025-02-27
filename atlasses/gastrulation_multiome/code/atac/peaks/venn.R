library(data.table)
library(purrr)
library(VennDiagram)

source(here::here("settings.R"))


io$peaks_10x     <- c(multiome1 = "/bi/scratch/Stephen_Clark/gastrulation_multiome_10x/multiome1/original/atac_peaks.bed",
                      multiome2 = "/bi/scratch/Stephen_Clark/gastrulation_multiome_10x/multiome2/original/atac_peaks.bed.gz")

io$peaks_test    <- c(snapatac = "/bi/scratch/Stephen_Clark/gastrulation_multiome_10x/processed/atac/snapatac/peaks.bed")

opts$overlap_gap <- 1000L

read_peaks <- function(files, names){
  peaks <- map2(files, names, ~fread(.x)[, c("sample", "i") := .(.y, .I)]) %>%
    rbindlist() %>%
    setnames(c("chr", "start", "end", "sample", "i")) %>%
    .[, chr := gsub("b'|'", "", chr)] %>%
    setkey(chr, start, end)
}

peaks_10x <- read_peaks(io$peaks_10x, names(io$peaks_10x))


peaks_test <- read_peaks(io$peaks_test, names(io$peaks_test))

compare_pairwise <- function(peaks1, peaks2, name1 = "peaks1", name2="peaks2"){
  overlaps <- foverlaps(peaks1, peaks2, nomatch = 0L, maxgap = opts$overlap_gap)
  
  # need to remove multi-overlaps
  overlaps[, n1 := .N, i.i]
  overlaps[, n2 := .N, i]
  overlaps[n1>1 | n2>1]
  
  overlaps[, seq1 := seq(1, .N, 1), by = i.i]
  overlaps[, seq2 := seq(1, .N, 1), by = i]
  
  remove1 <- overlaps[seq1>1, i.i]
  remove2 <- overlaps[seq2>1, i]
  
  peaks1 <- peaks1[!i %in% remove1]
  peaks2 <- peaks2[!i %in% remove2]
  overlaps <- overlaps[!i %in% remove2 & !i.i %in% remove1]
  
  
  missing1 <- peaks1[!i %in% overlaps$i.i] %>% .[, id := paste0("peaks1", i)]
  missing2 <- peaks2[!i %in% overlaps$i] %>% .[, id := paste0("peaks2", i)]
  
  # check number of sites adds up
  nrow(peaks1)
  nrow(overlaps) + nrow(missing1) 
  
  nrow(peaks2)
  nrow(overlaps) + nrow(missing2)
  
  cross <- nrow(overlaps)
  area1 <- nrow(peaks1) + cross
  area2 <- nrow(peaks2) + cross
  
  venn <- draw.pairwise.venn(area1, 
                             area2, 
                             cross, 
                             euler.d = TRUE,
                             fill = c("blue", "red"),
                             category = c(name1, name2))
  list(venn = venn, overlaps = overlaps)
}




compare10x <- compare_pairwise(peaks_10x[sample=="multiome1"],
                               peaks_10x[sample=="multiome2"],
                               "multiome1",
                               "multiome2")

compare_snap <- compare_pairwise(peaks_10x,
                                 peaks_test,
                                 "10X peaks",
                                 "SnapATAC peaks")

compare10x$venn
grid.newpage()
grid.draw(compare_snap$venn)

?draw.pairwise.venn
