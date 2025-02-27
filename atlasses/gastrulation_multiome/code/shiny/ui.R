library(shiny)
library(shinyFiles)
library(ggplot2)
library(DT)
# library(plotly)
library(rintrojs)
library(ggiraph)
library(shinythemes)
require(visNetwork)

#####################
## Define settings ##
#####################

basedir <- "/Users/argelagr/data/gastrulation_multiome_10x/shiny"

source("utils.R")

big_plot_width = "900px"
big_plot_height = "500px"

narrower_plot_width = "650px"

half_plot_width = "450px"
narrower_half_plot_width = "350px"
half_plot_height = "260px"

###############
## Load data ##
###############

TFs <- fread(paste0(basedir,"/TFs.txt"), header=F)[[1]]
celltypes <- fread(paste0(basedir,"/celltypes.txt"), header=F)[[1]]
genes <- fread(paste0(basedir,"/genes.txt"), header=F)[[1]]

##################
## Shiny app UI ##
##################

ui <- shinyUI(fluidPage(
          # sidebarLayout(

            # mainPanel(
            #   id = "main",
            #   width = 10,
            #   titlePanel(
            #     "The ultimate single-cell multi-omics roadmap of mouse gastrulation and early organogenesis"
            #   ),
            navbarPage(
              # tabsetPanel(
                # id = "tabs",
                title = "A single-cell multi-omics roadmap of mouse gastrulation and early organogenesis",
                theme = shinytheme("spacelab"),

                tabPanel(
                  title = "UMAPs", id = "umap",
                  sidebarPanel(width=3,
                    selectInput(inputId = "modality", label = "Modality", choices = c("ATAC" = "atac", "RNA" = "rna"), selected = "atac"),
                    selectInput(inputId = "stage", label = "Cell subset", choices = c("All timepoints"="all", "E7.5"="E7.5", "E7.75"="E7.75", "E8.5"="E8.5"), selected = "all"),
                    selectInput(inputId = "colourby", label = "Plot colour", choices = c("Cell type"="celltype", "Timepoint"="stage", "Sample"="sample", "Gene expression"="gene_expression", "Chromatin accessibility"="chromatin_accessibility"), selected = "celltype"),
                    conditionalPanel(
                      condition = "input.colourby == 'gene_expression'",
                      selectizeInput("gene_umap_rna", "Select gene to show RNA expression", choices = NULL, selected = "T")
                    ),
                    conditionalPanel(
                      condition = "input.colourby == 'chromatin_accessibility'",
                      selectizeInput("gene_umap_atac", "Select gene to show chromatin accessibility", choices = NULL, selected = "T")
                    ),
                    # checkboxInput("numbers", "Annotate clusters in plot"),
                  ),
                  mainPanel(
                    girafeOutput("umap", width = "900px", height = "800px"),
                    # plotOutput("stage_contribution", width = big_plot_width)
                  )
                ),
                
                tabPanel(
                  title = "TF activities", id = "rna_vs_chromvar_per_tf",
                  sidebarPanel(width=3,
                    selectizeInput(inputId = "tf", label = "Select transcription Factor", choices=NULL),
                    plotOutput("rna_vs_chromvar_motif", height="100px", width="300px"),
                    checkboxInput("rna_vs_chromvar_pseudobulk_scatterplot_add_text", "Add cell type labels to the scatterplot"),
                    sliderInput("rna_vs_chromvar_cor_range", label = "RNA vs chromVAR corr. interval", min=-1, max=1, step=0.01, value = c(-1,1))
                  ),
                  mainPanel(
                    girafeOutput("rna_vs_chromvar_per_tf"),
                    verbatimTextOutput('rna_vs_chromvar_print')
                  )
                ),
                
                tabPanel(
                  title = "Celltype TF signatures", id = "tf_activity_per_celltype",
                  sidebarPanel(width=3,
                    selectInput(inputId = "RNA_vs_chromVAR_pseudobulk_per_celltype_colour_by", label = "Colour by",
                                choices = c("Activator/Repressor" = "activator_repressor", "Marker strength" = "marker_strength")),
                    # checkboxInput("ignore_TFs_small_activity", "Ignore TFs with small activity", value=TRUE),
                    selectInput("celltype_tf_signatures", "Select celltype", choices = celltypes, selected = "Gut"),
                    # selectInput("tf_tf_signatures", "Select transcription Factor", choices = TFs, selected = c("FOXA2","SOX17","TAL1"), multiple = TRUE),
                    selectizeInput("tf_tf_signatures", "Select Transcription Factor", choices = NULL, multiple = TRUE),
                    sliderInput("min_rna_threshold_tf_signature", label = "Minimum RNA expression", min = 0, max = 12.5, value = 3),
                    sliderInput("min_chromvar_threshold_tf_signature", label = "Minimum Motif accessibility", min = 0, max = 13, value = 3),
                  ),
                  mainPanel(
                    girafeOutput("rna_vs_chromvar_per_celltype")
                    # verbatimTextOutput('print_celltype_tf_signatures')
                  )
                ),
                
                # tabPanel(
                #   title = "Gene expression vs accessibility", id = "rna_vs_acc",
                #   sidebarPanel(width=3,
                #     selectizeInput("gene_rna_vs_acc", "Select gene", choices = NULL, selected = "T")
                #   ),
                #   mainPanel(
                #     # girafeOutput("plot_expression_vs_accessibility_per_gene", width = big_plot_width, height = big_plot_height),
                #     girafeOutput("plot_expression_vs_accessibility_per_gene")
                #   )
                # ),
                
                tabPanel(
                  title = "In silico ChIP-seq", id = "insilico_chipseq",
                  sidebarPanel(width=3,
                    selectInput("tf_chip", "Transcription Factor", choices = TFs, selected = c("FOXA2","TAL1","ZIC2"), multiple = TRUE),
                    # selectInput("celltype_chip", "Celltype", choices = celltypes, selected = "All celltypes"),
                    sliderInput("min_score_filt", label = "Minimum score threshold", min = 0, max = 1.5, value = 0.50),
                    shinyFilesButton('files', label='File select', title='Please select a file', multiple=T) ,
                    downloadButton("download_chip_bed", "Download BED file"),
                    verbatimTextOutput('download_chip_bed_files_print')
                  ),
                  mainPanel(
                    girafeOutput("chip_stats_plot"),
                  )
                ),
                
                tabPanel(
                  title = "Celltype-specific GRNs", id = "celltype_grn",
                  sidebarPanel(width=3,
                    selectInput("celltype_grn", "Celltype", choices = celltypes, selected = "Gut"),
                    selectInput("celltype_grn_colorby", "Color by", c("Expression"="expression", "Eigenvalue centrality"="eigenvalue_centrality", "Degree centrality"="degree_centrality")),
                    shinyFilesButton('grn_files', label='File select', title='Please select a file', multiple=T) ,
                    downloadButton("download_grn_files", "Download GRN (igraph)"),
                    verbatimTextOutput('download_grn_files_print')
                  ),
                  mainPanel(
                    visNetworkOutput("network", height = "400px"),
                    plotOutput("grn_celltype_stats", height="200px", width="300px")
                  )
                ),
                
                tabPanel(
                  title = "Differential analysis", id = "differential",
                  sidebarPanel(width=3,
                    selectInput(inputId = "celltypeA", label = "Select celltype A", choices = celltypes, selected = "Gut"),
                    selectInput(inputId = "celltypeB", label = "Select celltype B", choices = celltypes, selected = "Erythroid3"),
                    selectInput(inputId = "tf_differential", label = "Select Transcription Factor", choices = NULL),
                    plotOutput("differential_motif", height="100px", width="300px"),
                    # sliderInput("differential_range", label = "Range of differential values", min=-25, max=25, step=0.5, value = c(-25,25)),
                    # sliderInput("highlight_top_n_genes", label = "Highlight top N genes", min=0, max=100, step=1, value = 0)
                  ),
                  mainPanel(
                    girafeOutput("differential_analysis")
                  )
                ),
                tabPanel(
                  title = "Trajectories", id = "trajectories",
                  sidebarPanel(width=3,
                               selectInput(inputId = "trajectory", label = "Select trajectory",
                                           choices = c("Ectoderm"="ectoderm", "Mesoderm"="mesoderm", "Endoderm"="endoderm", "Blood"="blood"), selected = "blood"),
                               selectizeInput(inputId = "tf_trajectory", label = "Select transcription Factor", choices = NULL),
                               # sliderInput("min_logpval_filt", label = "Minimum score threshold", min = 0, max = 1.5, value = 0.50),
                               plotOutput("trajectory_motif", height="100px", width="300px"),
                               sliderInput("trajectory_rna_vs_chromvar_cor_range", label = "RNA vs chromVAR corr. interval", min=-1, max=1, step=0.01, value = c(-1,1))
                  ),
                  mainPanel(
                    girafeOutput("trajectory_rna_vs_chromvar"),
                    verbatimTextOutput('trajectory_rna_vs_chromvar_print')
                  )
                ),
                
                tabPanel(
                  title = "Genome browser", id = "genome_browser",
                  mainPanel(
                    HTML(sprintf("%s: %s you can download a precomputed session to interactively explore the ATAC data on the %s (see snapshot below). Watch %s for instructions", shiny::tags$b("Download IGV session \n"), a("Here", href="igv_session_template.xml"), a("IGV browser", href="https://software.broadinstitute.org/software/igv"), a("this video", href="XXX"))),
                    br(), br(),
                    # shiny::img(src = "igv_session_overview.png", height = 450, width = 900),
                    HTML('<center><img src="igv_session_overview.png", height="450px", width="900px"></center>'),
                    br(), br(),
                    HTML(sprintf("%s: %s you can download a precomputed session to interactively explore the ATAC data on the %s (see snapshot below). Watch %s for instructions", shiny::tags$b("Download UCSC session \n"), a("Here", href="XXX"), a("UCSC browser", href="https://software.broadinstitute.org/software/igv"), a("this video", href="XXX")))
                  )
                )
              )
          )
)