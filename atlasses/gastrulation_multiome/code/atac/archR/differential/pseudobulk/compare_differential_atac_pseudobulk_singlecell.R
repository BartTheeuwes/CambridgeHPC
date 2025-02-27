#####################
## Define settings ##
#####################

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

# I/O
io$diff.pseudobulk <- file.path(io$basedir,"results_new/atac/archR/differential/PeakMatrix/pseudobulk")
io$diff.single_cell <- file.path(io$basedir,"results_new/atac/archR/differential/PeakMatrix")
io$outdir <- file.path(io$basedir,"results_new/atac/archR/differential/pseudobulk/comparison"); dir.create(io$outdir, showWarnings = F)

# Options
opts$matrix <- "PeakMatrix"

opts$celltypes <- c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast"
)

#############################
## Load pseudobulk results ##
#############################

# i <- "Epiblast"; j <- "Primitive_Streak"
atac_diff_pseudobulk.dt <- opts$celltypes %>% map(function(i) { opts$celltypes %>% map(function(j) {
  file <- file.path(io$diff.pseudobulk,sprintf("%s_vs_%s_%s_pseudobulk.txt.gz",i,j,opts$matrix))
  if (file.exists(file)) {
    fread(file, select = c(1,2)) %>% 
      .[,c("celltypeA","celltypeB"):=list(i,j)] %>%
      setnames("gene","idx") %>%
      return
  }
}) %>% rbindlist }) %>% rbindlist

##############################
## Load single-cell results ##
##############################

# i <- "Epiblast"; j <- "Primitive_Streak"
atac_diff_cells.dt <- opts$celltypes %>% map(function(i) { opts$celltypes %>% map(function(j) {
  file <- file.path(io$diff.single_cell,sprintf("%s_%s_vs_%s.txt.gz",opts$matrix,i,j))
  if (file.exists(file)) {
    fread(file, select = c(1,2)) %>% 
      .[,c("celltypeA","celltypeB"):=list(i,j)] %>%
      setnames("MeanDiff","diff") %>%
      return
  }
}) %>% rbindlist }) %>% rbindlist


###########
## Merge ##
###########

atac_diff.dt <- merge(
  atac_diff_cells.dt,
  atac_diff_pseudobulk.dt,
  by = c("idx","celltypeA","celltypeB"),
  suffixes = c("_cells","_pseudobulk")
)


##########
## Plot ##
##########

to.plot <- atac_diff.dt[celltypeA=="Primitive_Streak" & celltypeB=="Caudal_epiblast"]
ggscatter(to.plot, x="diff_cells", y="diff_pseudobulk", size=0.5,
          add="reg.line", add.params = list(color="blue", fill="lightgray"), conf.int=TRUE) +
  labs(x="Differential acc. (cells)", y="Differential acc. (pseudobulk)")
