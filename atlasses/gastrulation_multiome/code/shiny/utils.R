
#####################
## Plot dimensions ##
#####################

big_plot_width = 9 * 1.5
big_plot_height = 5 * 1.5
narrower_plot_width = 6.5 * 1.5

#####################
## Colour palettes ##
#####################

celltype_colours <- c(
  "Epiblast" = "#635547",
  "Primitive_Streak" = "#DABE99",
  "Caudal_epiblast" = "#9e6762",
  "PGC" = "#FACB12",
  "Anterior_Primitive_Streak" = "#c19f70",
  "Notochord" = "#0F4A9C",
  "Def._endoderm" = "#F397C0",
  "Gut" = "#EF5A9D",
  "Nascent_mesoderm" = "#C594BF",
  "Mixed_mesoderm" = "#DFCDE4",
  "Intermediate_mesoderm" = "#139992",
  "Caudal_Mesoderm" = "#3F84AA",
  "Paraxial_mesoderm" = "#8DB5CE",
  "Somitic_mesoderm" = "#005579",
  "Pharyngeal_mesoderm" = "#C9EBFB",
  "Cardiomyocytes" = "#B51D8D",
  "Allantois" = "#532C8A",
  "ExE_mesoderm" = "#8870ad",
  "Mesenchyme" = "#cc7818",
  "Haematoendothelial_progenitors" = "#FBBE92",
  "Endothelium" = "#ff891c",
  # "Blood_progenitors" = "#c9a997",
  "Blood_progenitors_1" = "#f9decf",
  "Blood_progenitors_2" = "#c9a997",
  # "Erythroid" = "#EF4E22",
  "Erythroid1" = "#C72228",
  "Erythroid2" = "#f79083",
  "Erythroid3" = "#EF4E22",
  "NMP" = "#8EC792",
  # "Neurectoderm" = "#65A83E",
  "Rostral_neurectoderm" = "#65A83E",
  "Caudal_neurectoderm" = "#354E23",
  "Neural_crest" = "#C3C388",
  "Forebrain_Midbrain_Hindbrain" = "#647a4f",
  "Spinal_cord" = "#CDE088",
  "Surface_ectoderm" = "#f7f79e",
  "Visceral_endoderm" = "#F6BFCB",
  "ExE_endoderm" = "#7F6874",
  "ExE_ectoderm" = "#989898",
  "Parietal_endoderm" = "#1A1A1A"
)

celltype_palette = scale_color_manual(values = celltype_colours, name = "", drop=TRUE)
celltype_palette_fill = scale_fill_manual(values = celltype_colours, name = "", drop=TRUE)
sample_palette <- scale_color_brewer(palette="Dark2")
rna_palette <- scale_color_gradient(low = "gray80", high = "red")
atac_palette <- scale_color_gradient(low = "gray80", high = "blue")

stage_colours = c(
  # "E6.5" = "#D53E4F",
  # "E6.75" = "#F46D43",
  # "E7.0" = "#FDAE61",
  # "E7.25" = "#FEE08B",
  "E7.5" = "#FFFFBF",
  # "E7.75" = "#E6F598",
  "E8.0" = "#ABDDA4",
  # "E8.25" = "#66C2A5",
  "E8.5" = "#3288BD"
)

stage_palette = scale_color_manual(values = stage_colours, name = "stage")
stage_palette_fill = scale_fill_manual(values = stage_colours, name = "stage")

###############
## Functions ##
###############

# taken from iSEE
subsetPointsByGrid <- function(X, Y, resolution=200, seed = 42) {
  set.seed(seed)
  # Avoid integer overflow when computing ids.
  resolution <- max(resolution, 1L)
  resolution <- min(resolution, sqrt(.Machine$integer.max))
  resolution <- as.integer(resolution)
  
  # X and Y MUST be numeric.
  rangeX <- range(X)
  rangeY <- range(Y)
  
  binX <- (rangeX[2] - rangeX[1])/resolution
  xid <- (X - rangeX[1])/binX
  xid <- as.integer(xid)
  
  binY <- (rangeY[2] - rangeY[1])/resolution
  yid <- (Y - rangeY[1])/binY
  yid <- as.integer(yid)
  
  # Getting unique IDs, provided resolution^2 < .Machine$integer.max
  # We use fromLast=TRUE as the last points get plotted on top.
  id <- xid + yid * resolution 
  !duplicated(id, fromLast=TRUE)
}


################
## Plot utils ##
################

makeGenePlot = function(gene_name, gene_counts, x_coord, y_coord){
  
  order = order(gene_counts)
  
  p = ggplot(mapping = aes(x = x_coord[order], y = y_coord[order], col = gene_counts[order])) +
    geom_point(size = 1) +
    scale_color_gradient2(name = "Log2\nnormalised\ncounts", mid = "cornflowerblue", low = "gray75", high = "black", midpoint = max(gene_counts)/2) +
    coord_fixed(ratio = 0.8) +
    ggtitle(gene_name) +
    theme_classic() +
    ggplot_theme_NoAxes()
    # theme(axis.title = element_blank(), 
    #       axis.text = element_blank(), 
    #       axis.ticks = element_blank()
    # )
  
  if(max(gene_counts) == 0){
    p = p +
      scale_color_gradient2(name = "Log2\ncounts", mid = "gray75", low = "gray75", high = "gray75", midpoint = max(gene_counts)/2)
  }
  
  return(p)
  
}



minmax.normalisation <- function(x) {
  return((x-min(x,na.rm=T)) /(max(x,na.rm=T)-min(x,na.rm=T)))
}

ggplot_theme_NoAxes <- function() {
  theme(
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
}


matrix.please<-function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[[1]]
  m
}

trajectories <- list(
  "mesoderm" = c("Epiblast", "Primitive_Streak", "Nascent_mesoderm", "Mixed_mesoderm", "Pharyngeal_mesoderm"),
  "endoderm" = c("Epiblast", "Primitive_Streak", "Anterior_Primitive_Streak", "Def._endoderm", "Gut"),
  "ectoderm" = c("Epiblast", "Rostral_neurectoderm","Forebrain_Midbrain_Hindbrain"),
  "blood" = c("Haematoendothelial_progenitors", "Blood_progenitors_1", "Blood_progenitors_2", "Erythroid1", "Erythroid2", "Erythroid3")
)

sort.abs <- function(dt, sort.field) dt[order(-abs(dt[[sort.field]]))]
