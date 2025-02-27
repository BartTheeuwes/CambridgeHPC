here::i_am("atac/archR/chromvar/differential/pseudobulk/differential/wt_vs_ko/plot_archR_differential_chromvar_pseudobulk_wt_vs_ko.R")

# load default setings
source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# Options
opts$motif_annotation <- "CISBP"

# I/O
io$archR.chromvar.pseudobulk.se <- sprintf("%s/results/atac/archR/chromvar/pseudobulk/celltype_genotype/chromVAR_deviations_%s_pseudobulk_archr.rds",io$basedir,opts$motif_annotation)
io$archR.chromvar.diff <- file.path(io$basedir,"results/atac/archR/chromvar/differential/pseudobulk/wt_vs_ko/WT_vs_T_KO_chromVAR_pseudobulk.txt.gz")
io$atac.pseudobulk.stats <- file.path(io$basedir,"results/atac/archR/pseudobulk/celltype_genotype/stats.txt")
io$outdir <- file.path(io$basedir,"results/atac/archR/chromvar/differential/pseudobulk/wt_vs_ko/pdf"); dir.create(io$outdir, showWarnings = F, recursive = T)

###########################################
## Load precomputed diff chromVAR scores ##
###########################################

chromvar_diff.dt <- fread(io$archR.chromvar.diff)

###########################
## Load pseudobulk stats ##
###########################

atac_pseudobulk_stats.dt <- fread(io$atac.pseudobulk.stats)

#####################################
## Load pseudobulk chromVAR scores ##
#####################################

atac_chromvar_pseudobulk.se <- readRDS(io$archR.chromvar.pseudobulk.se)#[,samples.to.use]
assays(atac_chromvar_pseudobulk.se) <- assays(atac_chromvar_pseudobulk.se)["z"]

##########
## Plot ##
##########

to.plot <- chromvar_diff.dt %>%
  .[,sign:=as.factor(c("Up in T KO","Down in T KO")[(diff>0)+1])] %>%
  .[,rank:=1:.N,by=c("celltype","sign")] %>%
  .[,celltype_sign:=paste(celltype,sign,"_")]

ggplot(to.plot, aes(x=rank, y=abs(diff), color=celltype)) +
  # geom_point() +
  geom_line(size=1) +
  facet_wrap(~sign) +
  scale_color_manual(values=opts$celltype.colors) +
  geom_hline(yintercept = 0, linetype="dashed") +
  labs(y="Differential motif accessibility (chromVAR)", x="Rank") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size=rel(1.0), color="black"),
    axis.text.x = element_text(size=rel(1.0), color="black")
  )
