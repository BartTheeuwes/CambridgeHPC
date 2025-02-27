here::i_am("rna/differential/wt_vs_ko/analysis/barplots.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

# Load utils
source(here::here("rna/differential/wt_vs_ko/analysis/utils.R"))

##############
## Settings ##
##############

# I/O
io$indir <- file.path(io$basedir,"results/rna/differential/wt_vs_ko")
io$outdir <- file.path(io$basedir,"results/rna/differential/wt_vs_ko/pdf/volcano_plots"); dir.create(io$outdir, showWarnings = F, recursive = T)

# Options
opts$min.cells <- 50

###############
## Load data ##
###############

source(here::here("rna/differential/wt_vs_ko/analysis/load_data.R"))

####################
## Filter results ##
####################

# Remove some hits
# diff.dt <- diff.dt[gene!="Xist"]
# diff.dt <- diff.dt[!grepl("mt-",gene)]
# diff.dt <- diff.dt[!grepl("Rps|Rpl",gene)]
# diff.dt <- diff.dt[!grepl("Rik",gene)]
diff.dt <- diff.dt[!grepl("^Hb",gene)]
# diff.dt <- diff.dt[!gene%in%fread(io$gene_metadata)[chr=="chrY",symbol]]

# Filter by minimum number of cells per group
# opts$min.cells <- 30
# diff.dt <- diff.dt[groupA_N>opts$min.cells & groupB_N>opts$min.cells]

# Remove hits that are differentially expressed in all cell type comparisons
# foo <- diff.dt[,mean(sig),by=c("gene")] %>% .[V1>0] %>% setorder(-V1)

# Subset to lineage markers
marker_genes.dt <- fread(io$rna.atlas.marker_genes.up)
diff_markers.dt <- diff.dt[gene%in%unique(marker_genes.dt$gene)]

##########
## Plot ##
##########

# i <- opts$ko.classes[1]; j <- opts$celltypes[1]
celltypes.to.plot <- unique(diff_markers.dt$celltype) %>% as.character
for (i in celltypes.to.plot) {
  to.plot <- diff_markers.dt[celltype==i] %>% .[!is.na(sig)] 
  if (nrow(to.plot)>0) {
    p <- gg_volcano_plot(to.plot, top_genes = 35, groupA = opts$wt.class, groupB = opts$ko.class)
    
    pdf(file.path(io$outdir,sprintf("%s_vs_%s_%s_volcano.pdf",i,opts$wt.class,opts$ko.class)), width=9, height=6)
    print(p)
    dev.off()
  }
}
