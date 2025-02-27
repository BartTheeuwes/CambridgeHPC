library(data.table)
library(purrr)
library(ggplot2)
library(cowplot)


source(here::here("settings.R"))


io$nmt_meta            <- "/bi/scratch/Stephen_Clark/gastrulation_data/sample_metadata.txt"
io$nmt_parsed_data     <- file.path(io$rawdata, "public_datasets/Argelaguet_2019/nmt_diffacc_each_vs_all.tsv.gz")

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

opts$annos  <- c("Forebrain_Midbrain_Hindbrain", "Epiblast", "Surface_ectoderm")

anno_grep <- paste(opts$annos, collapse = "|")

meta <- fread(io$nmt_meta) %>% 
  .[, stage_lineage := paste0(stage, "_", lineage10x_2)] %>% 
  .[pass_accQC == TRUE & stage_lineage %in% opts$stage_lineages] %>% 
  .[, cell := id_acc]

cells <- meta[, cell]


nmt <- fread(io$nmt_parsed_data) %>% 
  setkey(cell) %>% 
  .[cell %in% cells] %>% 
  setkey(id)

toplot <- nmt[id %like% anno_grep] %>% 
  setkey(cell) %>% 
  merge(meta[, .(cell, stage_lineage)], by = "cell") %>% 
  .[, .(acc = mean(acc/mean_acc)), .(id, stage_lineage)] %>% 
  .[, anno := strsplit(id, "_") %>% map_chr(~paste(.[c(1:length(.)-1)], collapse = "_"))] %>% 
  .[, stage_lineage := factor(stage_lineage, levels = opts$stage_lineages)]



ggplot(toplot, aes(stage_lineage, acc, fill = stage_lineage)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  ylab("adjusted accessibility") +
  facet_wrap(~anno)


