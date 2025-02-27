library(R.utils)
library(shiny)
library(shinyFiles)
library(ggplot2)
library(HDF5Array)
library(data.table)
library(purrr)
library(cowplot)
library(ggrepel)
library(GGally)
library(ggiraph)
library(ggseqlogo)
require(igraph)
require(visNetwork)
require(patchwork)
require(ggpubr) # to remove this dependency?
# library(DT)
# library(plotly)

setwd("/Users/argelagr/gastrulation_multiome_10x/shiny")

basedir <- "/Users/argelagr/data/gastrulation_multiome_10x/shiny"

################
## Load utils ##
################

source("utils.R")

# Updated ggnet2 function
source("/Users/argelagr/gastrulation_multiome_10x/rna_atac/gene_regulatory_networks/pseudobulk/ggnet2.R")

###############
## Load data ##
###############

source("load_data.R")

###############
## Shiny app ##
###############

server <- function(input, output, session) {
  
  
  #######################
  ## Selectize speedup ##
  #######################
  
  updateSelectizeInput(session = session, inputId = 'gene_umap_atac', choices = genes, server = TRUE, selected = "T") 
  updateSelectizeInput(session = session, inputId = 'gene_umap_rna', choices = genes, server = TRUE, selected = "T") 
  updateSelectizeInput(session = session, inputId = 'gene_rna_vs_acc', choices = genes, server = TRUE, selected = "T")
  updateSelectizeInput(session = session, inputId = 'tf', choices = TFs, server = TRUE, selected = "FOXA2")
  updateSelectizeInput(session = session, inputId = 'tf_tf_signatures', choices = TFs, server = TRUE, selected=c("FOXA2"))
  updateSelectizeInput(session = session, inputId = 'tf_trajectory', choices = TFs, server = TRUE, selected = "KLF1")
  updateSelectizeInput(session = session, inputId = 'tf_differential', choices = TFs, server = TRUE, selected = "KLF1")
  
  
  ###########
  ## UMAPs ##
  ###########
  
  plot_UMAP <- reactive({
    
    ## Fetch data ##
    
    if(input$modality == "rna"){
      cells <- sample_metadata[pass_rnaQC==TRUE & doublet_call==FALSE & !is.na(celltype),cell]
      umap.dt <- umap_rna.df[cells,] %>% as.data.table(keep.rownames = T) %>% setnames("rn","cell")
    } else if(input$modality == "atac"){
      cells <- sample_metadata[pass_atacQC==TRUE & doublet_call==FALSE & !is.na(celltype),cell]
      umap.dt <- umap_atac.df[cells,] %>% as.data.table(keep.rownames = T) %>% setnames("rn","cell")
    }
      
    # allowed = get_subset()
    #scramble, and subset if asked
    # color = get_clusters()[allowed]
    
    # color <- sample_metadata[cells][[input$colourby]]
    
    to.plot <- umap.dt %>% 
      merge(sample_metadata[,c("cell","celltype","stage","sample")], by=c("cell"))
    
    if (input$colourby == "gene_expression") {
      tmp <- data.table(
        cell = colnames(link_rna_expr),
        color = as.numeric(link_rna_expr[input$gene_umap_rna,])
      )
      to.plot <- to.plot %>% merge(tmp,by="cell") %>% setorder(color)
    } else if (input$colourby == "chromatin_accessibility") {
        tmp <- data.table(
          cell = colnames(link_gene_acc),
          color = as.numeric(link_gene_acc[input$gene_umap_atac,])
        )
        to.plot <- to.plot %>% merge(tmp,by="cell") %>% setorder(color)
    } else {
      to.plot$color <- to.plot[[input$colourby]]
    }

    ## Plot PAGA ##
    
    celltypes <- sapply(paga$val,"[[","vertex.names")
    alphas <- rep(1.0,length(celltypes)); names(alphas) <- celltypes
    sizes <- rep(8,length(celltypes)); names(sizes) <- celltypes
    
    # p1 <- ggnet2(
    p.paga <- ggnet2_interactive(
      net = paga,
      mode = c("x", "y"),
      color = celltype_colours[celltypes],
      node.alpha = alphas,
      node.size = sizes,    
      edge.size = 0.15,
      edge.color = "grey",
      label = TRUE,
      label.size = 3
    )
    
    ## Plot UMAP ##
    
    to.plot <- to.plot[sample(1:.N, size = 10000)]
    
    p.umap <- ggplot(to.plot, aes(x = X, y = Y, color = color)) +
      # geom_point(size = 1, alpha = 0.9) +
      geom_point_interactive(aes(tooltip = celltype, data_id = celltype), size=1, alpha=0.9) +
      coord_fixed(ratio = 0.8) +
      theme_classic() +
      ggplot_theme_NoAxes() +
      theme(
        legend.position = "none"
      )

    # Modify legends    
    if (input$colourby%in%c("celltype","stage","sample")) {
      p.umap <- p.umap + 
        guides(colour = guide_legend(override.aes = list(size=5, alpha = 1)))
    } else {
      p.umap <- p.umap +
        theme(
          legend.title = element_blank()
        )
    }
    
    # Define palette
    palette <- switch(input$colourby, 
      "celltype" = celltype_palette, 
      "stage" = stage_palette, 
      "sample" = sample_palette, 
      "gene_expression" = rna_palette, 
      "chromatin_accessibility" = atac_palette
    )
    p.umap <- p.umap + palette
      
      # if (input$numbers) {
      #   centroids = get_cluster_centroids()
      #   centroids$num = gsub(" ", "\n", centroids$num)
      #   p = p + geom_label_repel(data = centroids, mapping = aes(x = X, y = Y, label = num), col = "black", alpha = 0.8, size = 4)
      # }
      
    ## Barplot/Violin plot with statistics ##
    
    if (input$colourby=="celltype") {
      tab <- table(to.plot$celltype, droplevels(to.plot$stage))
      fractions <- sweep(tab, 1, rowSums(tab), "/")
      means <- apply(fractions, 1, function(x) sum(x * 1:length(x)))
      to.plot2 <- fractions %>% as.data.table(keep.rownames = T) %>% setnames(c("celltype","stage","value"))
      
      p2 <- ggplot(to.plot2, aes(x = factor(celltype, levels = names(means)[order(means)]), y = value, fill = stage)) +
        # geom_bar(stat = "identity") +
        geom_bar_interactive(aes(tooltip=celltype, data_id=celltype), stat = "identity") +
        labs(y = "Fraction of cells") +
        stage_palette_fill +
        theme_classic() +
        theme(
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black", size=rel(0.75)),
          axis.text.y = element_text(color="black", size=rel(0.8)),
          axis.title.x = element_blank()
        )
    } else if (input$colourby%in%c("stage","sample")) {
      
      to.plot2 <- to.plot[,.N,by=c(input$colourby,"celltype")]
      
      p2 <- ggplot(to.plot2, aes_string(x = input$colourby, y = "N", fill = "celltype")) +
        # geom_bar(stat = "identity") +
        geom_bar_interactive(aes(tooltip=celltype, data_id=celltype), stat = "identity") +
        labs(y = "Fraction of cells") +
        celltype_palette_fill +
        theme_classic() +
        theme(
          legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black"),
          axis.text.y = element_text(color="black", size=rel(0.8)),
          axis.title.x = element_blank()
        )
      
    } else if (input$colourby%in%c("gene_expression","chromatin_accessibility")) {
      
      clust.sizes <- table(to.plot$celltype)
      
      p2 <- ggplot(to.plot, aes(x = celltype, y = color, fill = celltype)) +
        geom_violin(scale = "width", alpha=0.8) +
        # geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.8) +
        geom_boxplot_interactive(aes(tooltip=celltype, data_id=celltype), width=0.5, outlier.shape=NA, alpha=0.8) +
        labs(y = "Log2 normalised counts") +
        annotate(
          geom = "text",
          # x = factor(names(clust.sizes)),
          x = names(clust.sizes),
          y = rep_len(c(max(to.plot$color)*1.1, max(to.plot$color) * 1.2), length.out = length(clust.sizes)),
          label = as.vector(clust.sizes),
          size = 3
        ) +
        celltype_palette_fill +
        theme_classic() +
        theme(
          axis.title = element_text(size = 11, color="black"),
          axis.text.y = element_text(size = 9, color="black"),
          axis.text.x = element_text(size = 12, color="black", angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "none",
          axis.title.x = element_blank()
        )
    }

  # p <- plot_grid(plotlist = list(p.umap, p2), nrow=2, scale=1, rel_heights = c(2/3,1/3))
    
    layout <- "
    ABB
    ABB
    ABB
    ABB
    CCC
    "
    girafe(
      # code = print((p.paga+p.umap)/p2 + plot_layout(widths = c(1,10), heights=c(4,1))),
      code = print(p.paga+p.umap+p2 + plot_layout(design = layout)),
      width_svg = 13, height_svg = 9,
      options = list( 
        opts_sizing(rescale = FALSE),
        # opts_selection(type = "single", css = "cursor:pointer;fill:magenta;color:magenta"),
        opts_selection(type = "single", css = ""),
        # opts_hover_inv(css = "opacity:0.45;"),
        opts_hover(css = "cursor:pointer;fill:magenta;color:magenta")
      )
    ) %>% return(.)
  
  })
  
  output$umap = renderGirafe({
    shiny::validate(need(input$gene_umap_rna%in%genes, "" ))
    shiny::validate(need(input$gene_umap_atac%in%genes, "" ))
    plot_UMAP()
  })
  
  ##############################################################
  ## RNA expression versus chromatin accessibility (per gene) ##
  ##############################################################
  
  # plot_expression_vs_accessibility_per_gene = reactive({
  # 
  #   # input$gene_rna_vs_acc <- "T"
  #   
  #   cells <- cells_rna
  #   
  #   dat <- umap_rna.df[cells,]
  # 
  #   ## PAGA ##
  #   
  #   celltypes <- sapply(paga$val,"[[","vertex.names")
  #   colors <- celltype_colours[celltypes]
  #   alphas <- rep(1.0,length(celltypes)); names(alphas) <- celltypes
  #   sizes <- rep(11,length(celltypes)); names(sizes) <- celltypes
  #   
  #   p.paga <- ggnet2_interactive(
  #     net = paga,
  #     mode = c("x", "y"),
  #     color = colors,
  #     node.alpha = alphas,
  #     node.size = sizes,    
  #     edge.size = 0.15,
  #     edge.color = "grey",
  #     label = TRUE,
  #     label.size = 3
  #   )
  #   
  #   ## RNA vs ACC scatterplot
  #   to.plot <- gene_rna_vs_atac_pseudobulk.dt[gene==input$gene_rna_vs_acc]
  #   
  #   to.plot$tooltip <- to.plot$celltype
  #   
  #   p.scatter <- ggplot(to.plot, aes(x=acc, y=expr, fill=celltype)) +
  #     geom_point_interactive(aes(tooltip = tooltip, data_id = tooltip), size=6, shape=21) +
  #     scale_fill_manual(values=celltype_colours) +
  #     # geom_text(aes(label=gene), size=3, data=to.plot.text) +
  #     labs(x="Gene accessibility (promoter region)", y="Gene expression") +
  #     theme_classic() +
  #     theme(
  #       legend.position = "none",
  #       axis.text = element_text(size=rel(0.7))
  #     )
  #   
  #   # if (input$rna_vs_chromvar_pseudobulk_scatterplot_add_text) {
  #   # p.scatter <- p.scatter + ggrepel::geom_text_repel(aes(label=celltype), size=3, max.overlaps=Inf)
  #   # }
  #   
  #   # Concatenate plots
  #   # p1 <- plot_grid(plotlist=list(p.paga,p.scatter,p.umap), scale=0.95, nrow=2, rel_heights = c(1/2,1/2))
  #   
  #   girafe(
  #     # ggobj = plot_grid(p1, p.umap, nrow=2, scale=0.95, rel_widths = c(1/2,1/2)), 
  #     ggobj = plot_grid(plotlist=list(p.paga,p.scatter), scale=0.95, ncol=2), 
  #     width_svg = 14, height_svg = 5,
  #     options = list(
  #       opts_sizing(rescale = FALSE),
  #       # opts_selection(type = "single", css = "cursor:pointer;r:25px"),
  #       opts_selection(type = "single", css = "cursor:pointer;r:2%"),
  #       opts_hover_inv(css = "opacity:0.65;cursor:pointer;r:.9%"),
  #       opts_hover(css = "cursor:pointer;r:2%")
  #     )
  #   ) %>% return(.)
  #   
  # })
  # 
  # output$plot_expression_vs_accessibility_per_gene = renderGirafe({
  #   shiny::validate(need(input$gene_rna_vs_acc%in%genes, "Please select a gene; if you have already selected one, this gene is not in our annotation." ))
  #   return(plot_expression_vs_accessibility_per_gene())
  # })



  #############################
  ## RNA vs chromVAR+ per TF ##
  #############################
  
  selected_gene_RNA_vs_chromVAR <- reactive({
    input$rna_vs_chromvar_per_tf_selected
  })
  
  plot_RNA_vs_chromVAR_pseudobulk_per_TF = reactive({
  
    ## Volcano plot ##
    # input$rna_vs_chromvar_cor_range <- c(-1,-0.3)
    xlim_min <- input$rna_vs_chromvar_cor_range[1]
    xlim_max <- input$rna_vs_chromvar_cor_range[2]
    
    to.plot <- cor_rna_vs_chromvar_per_gene.dt[r>=input$rna_vs_chromvar_cor_range[1] & r<=input$rna_vs_chromvar_cor_range[2]]
    ylim_min <- min(to.plot$log_pval,na.rm=T)
    ylim_max <- max(to.plot$log_pval,na.rm=T)
    
    negative_hits <- to.plot[sig==TRUE & r<0 & r>input$rna_vs_chromvar_cor_range[1],gene]
    positive_hits <- to.plot[sig==TRUE & r>0 & r<input$rna_vs_chromvar_cor_range[2],gene]
    all <- nrow(to.plot)
    
    # Select dot size
    if (all>=500) { dot_size <- 2.5 } else if (all>100 & all<500) { dot_size <- 5 } else if (all<=100) { dot_size <- 8 }
    
    p.volcano <- ggplot(to.plot, aes(x=r, y=log_pval)) +
      # geom_segment(aes(x=0, xend=0, y=0, yend=ylim_max-1), color="orange", size=0.25) +
      geom_jitter_interactive(aes(fill=log_pval, alpha=log_pval, tooltip=gene, data_id=gene, onclick=gene), width=0.03, height=0.03, shape=21, size=dot_size) + 
      scale_fill_gradient(low = "gray80", high = "red") +
      scale_alpha_continuous(range=c(0.25,1)) +
      scale_x_continuous(limits=c(xlim_min-0.15,xlim_max+0.15)) +
      scale_y_continuous(limits=c(ylim_min,ylim_max+3)) +
      annotate("text", x=median(to.plot$r,na.rm=T), y=ylim_max+3, size=4, label=sprintf("(%d)", all)) +
      annotate("text", x=xlim_min-0.05, y=ylim_max+2, size=4, label=sprintf("%d (-)",length(negative_hits))) +
      annotate("text", x=xlim_max+0.05, y=ylim_max+2, size=4, label=sprintf("%d (+)",length(positive_hits))) +
      labs(x="Pearson correlation\n(RNA expression vs Motif accessibility)", y=expression(paste("-log"[10],"(p.value)"))) +
      theme_classic() +
      theme(
        plot.margin = margin(t = 0, r = 25, b = 0, l = 0, unit = "pt"),
        axis.text = element_text(size=rel(1.1), color='black'),
        axis.title = element_text(size=rel(1.25), color='black'),
        legend.position="none"
      )
    
    if (cor_rna_vs_chromvar_per_gene.dt[gene==input$tf,r]>0) { nudge_x = -0.55 } else { nudge_x = 0.55  }
    p.volcano <- p.volcano + 
      ggrepel::geom_text_repel(data=cor_rna_vs_chromvar_per_gene.dt[gene==input$tf], aes(x=r, y=log_pval, label=gene), size=6, 
        nudge_x = nudge_x,
        box.padding = 0.5,
        nudge_y = 1,
        segment.curvature = -0.1,
        segment.ncp = 3,
        segment.angle = 20,
        arrow = arrow(length = unit(0.02, "npc"))
      )

    ## Scatterplot ##
    
    # input <- list(); input$tf <- "FOXA2"; input$rna_vs_chromvar_pseudobulk_scatterplot_add_text <- FALSE
    to.plot <- rna_vs_chromvar_pseudobulk.dt[gene==input$tf]
    
    to.plot[,tooltip:=sprintf("%s\nRNA expression = %.2f\nMotif accessibility = %.2f",celltype,expr,chromvar_zscore)]
    
    p.scatter <- ggplot(to.plot, aes(x=chromvar_zscore, y=expr, fill=celltype)) +
      geom_point_interactive(aes(tooltip=tooltip), size=6, shape=21) +
      celltype_palette_fill +
      # geom_text(aes(label=gene), size=3, data=to.plot.text) +
      labs(x="Motif accessibility (z-score)", y="Gene expression", title=input$tf) +
      theme_classic() +
      theme(
        plot.margin = margin(t = 0, r = 0, b = 25, l = 0, unit = "pt"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size=rel(1.5), color="black"),
        axis.title = element_text(size=rel(1.10), color="black"),
        axis.text = element_text(size=rel(1.10), color="black")
      )
    
    if (input$rna_vs_chromvar_pseudobulk_scatterplot_add_text) {
      p.scatter <- p.scatter + ggrepel::geom_text_repel(aes(label=celltype), size=3, max.overlaps=Inf)
    }

    ## motif seqLogo ##
    # m <- 0.25*exp(as.matrix(motifs[[input$tf]]))
    # p.motif <- ggseqlogo(m)
    
    ## PAGA ##
    
    to.plot2 <- to.plot %>% copy %>%
      .[,expr:=minmax.normalisation(expr)] %>%
      .[,chromvar_zscore:=minmax.normalisation(chromvar_zscore)]
      
    # Define colors
    rna.col.seq <- chromvar.col.seq <- round(seq(0,1,0.1), 2)
    rna.colors <- colorRampPalette(c("gray92", "darkgreen"))(length(rna.col.seq))
    chromvar.colors <- colorRampPalette(c("gray92", "purple"))(length(chromvar.col.seq)) 
    
    celltypes <- sapply(paga$val,"[[","vertex.names")
    colors <- celltype_colours[celltypes]
    alphas <- rep(1.0,length(celltypes)); names(alphas) <- celltypes
    sizes <- rep(6.5,length(celltypes)); names(sizes) <- celltypes
    
    p.paga <- ggnet2(
    # p.paga <- ggnet2_interactive(
      net = paga,
      mode = c("x", "y"),
      color = colors,
      node.alpha = alphas,
      node.size = sizes,    
      edge.size = 0.15,
      edge.color = "grey",
      label = TRUE,
      label.size = 2
    )
    
    # Interactive mode doesn't work well    
    p.base <- ggnet2(
      net = paga,
      mode = c("x", "y"),
      node.size = 0,
      edge.size = 0.15,
      edge.color = "grey",
      label = FALSE
    )
    
    expr.values <- to.plot2[,c("celltype","expr")] %>% matrix.please %>% .[celltypes,]
    expr.colors <- round(expr.values,1) %>% map(~ rna.colors[which(rna.col.seq == .)]) %>% unlist
    
    p.rna <- p.base + geom_text(label = "\u25D0", aes(x=x, y=y), color=expr.colors, size=12, family = "Arial Unicode MS",
                        data = p.base$data[,c("x","y")] %>% dplyr::mutate(expr=expr.colors)) +
      # scale_colour_manual(values=expr.colors) + 
      labs(title="RNA expression") +
      theme(
        plot.title = element_text(hjust = 0.5)
      )
    
    acc.values <- to.plot2[,c("celltype","chromvar_zscore")] %>% matrix.please %>% .[celltypes,]# %>% .[celltypes,]
    acc.colors <- round(acc.values,1) %>% map(~ chromvar.colors[which(chromvar.col.seq == .)]) %>% unlist
    
    p.acc <- p.base + geom_text(label = "\u25D1", aes(x=x, y=y), color=acc.colors, size=12, family = "Arial Unicode MS",
                        data = p.base$data[,c("x","y")] %>% dplyr::mutate(acc=acc.colors)) +
      # scale_fill_manual(values=acc.colors) + 
      labs(title="Motif accessibility") +
      theme(
        plot.title = element_text(hjust = 0.5)
      )
    
    # Concatenate plots
    # p1 <- plot_grid(plotlist=list(p.scatter,p.motif), nrow=1, rel_widths = c(3/5,2/5))
    # p2 <- plot_grid(plotlist=list(p.paga,p.rna,p.acc), nrow=1, scale = 0.95)
    
    # htmlwidget call
    x <- girafe(
      # ggobj = plot_grid(p1, p2, nrow=2, scale=0.95, rel_heights = c(2/5,3/5)), 
      code = print(((p.volcano)|(p.scatter/(p.paga+p.rna+p.acc))) + plot_layout(width=c(1,2))),
      width_svg = 16, height_svg = 10,
      options = list(
        opts_zoom(min = 1, max = 5),
        opts_sizing(rescale = TRUE),
        # opts_selection(type = "single", css = "cursor:pointer;r:1.5%"),
        # opts_hover_inv(css = "opacity:0.65;cursor:pointer;r:.7%"),
        # opts_hover(css = "cursor:pointer;r:1.5%")
        opts_selection(type = "single", css = ""),
        opts_hover(css = "cursor:pointer;fill:magenta;color:magenta")
      )
    ) %>% return(.)
  })
  

  output$rna_vs_chromvar_per_tf = renderGirafe({
    shiny::validate(need(input$tf%in%TFs, "Please select a TF from our annotation" ))
    return(plot_RNA_vs_chromVAR_pseudobulk_per_TF())
  })
  
  observeEvent(selected_gene_RNA_vs_chromVAR(), {
    updateTextInput(session = session, "tf", value = selected_gene_RNA_vs_chromVAR() )
  })
  
  output$rna_vs_chromvar_print <- renderPrint({
    if (length(input$tf)>0) {
      tmp <- cor_rna_vs_chromvar_per_gene.dt %>%
        .[r>=input$rna_vs_chromvar_cor_range[1] & r<=input$rna_vs_chromvar_cor_range[2]] %>%
        .[,abs_r:=abs(r)] %>% setorder(-abs_r) %>% head(n=15)
      cat("Top 15 TFs with largest absolute RNA vs chromVAR+ correlation (within the selected range):\n")
      cat(sprintf("%s: r=%.2f\t", tmp$gene, round(tmp$r,2)))
    }
  })
  
  output$rna_vs_chromvar_motif <- renderPlot({
    shiny::validate(need(input$tf%in%TFs, "Please select a TF from our annotation" ))
    m <- 0.25*exp(as.matrix(motifs[[input$tf]]))
    ggseqlogo(m) +
      theme(
        panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank()
      )
  }, bg="transparent", execOnResize = TRUE)
  

  
  ############################################
  ## RNA vs chromVAR+ per TF (per celltype) ##
  ############################################

  plot_RNA_vs_chromVAR_pseudobulk_per_celltype = reactive({
    
    # input <- list(); input$celltype_tf_signatures <- "Gut"; input$tf_tf_signatures <- c("FOXA2","SOX17","TAL1"); 
    # input$ignore_TFs_small_activity <- TRUE; input$RNA_vs_chromVAR_pseudobulk_per_celltype_colour_by <- "activator_repressor"
    # input$min_chromvar_threshold_tf_signature <- 3; input$min_rna_threshold_tf_signature <- 3
    
    celltype_virtual_chip.dt <- fread(sprintf("%s/insilico_chipseq/celltype/%s_virtual_chip.txt.gz",basedir,input$celltype_tf_signatures))
    
    ## Lineplot for TF binding sites ##
    
    seq.ranges <- seq(0.25,1,by=0.01)
    celltype_virtual_chip_to_plot.dt <- seq.ranges %>% map(function(j) {
      celltype_virtual_chip.dt[score>=j,log2(.N+1),by="tf"] %>% .[,min_score:=j] %>% return
    }) %>% rbindlist %>% setnames("V1","log2_N")
    
    active.tfs <- celltype_virtual_chip_to_plot.dt[min_score>=0.40,max(log2_N),by="tf"] %>% .[V1>=6] %>% .$tf
    updateTextInput(session = session, "tf_tf_signatures", value = active.tfs)
    celltype_virtual_chip_to_plot.dt <- celltype_virtual_chip_to_plot.dt[tf%in%active.tfs]
    
    p.lineplot <- ggplot(celltype_virtual_chip_to_plot.dt, aes_string(x="min_score", y="log2_N", color="tf")) +
      geom_line_interactive(size=1, aes(tooltip=tf, data_id=tf)) +
      labs(y="Number of predicted binding sites (log2)", x="Minimum predicted binding score") +
      theme_classic() +
      theme(
        axis.text = element_text(color="black"),
        # legend.position = "right",
        legend.position = c(.9,.65),
        legend.title = element_blank()
      )
    
    
    ## Scatterplot ##
    
    max.chromvar <- 13
    max.expr <- 12
    
    to.plot <- rna_vs_chromvar_pseudobulk.dt %>%
      .[expr>=input$min_rna_threshold_tf_signature & chromvar_zscore>=input$min_chromvar_threshold_tf_signature] %>%
      .[celltype==input$celltype_tf_signatures] %>%
      merge(cor_rna_vs_chromvar_per_gene.dt[,c("gene","cor_sign")], by="gene") %>%
      .[chromvar_zscore<0,chromvar_zscore:=0] %>%
      .[chromvar_zscore>=max.chromvar,chromvar_zscore:=max.chromvar] %>%
      .[expr>=max.expr,expr:=max.expr] %>% 
      .[,dot_size:=minmax.normalisation(expr)*minmax.normalisation(chromvar_zscore)]
    
    # if (input$ignore_TFs_small_activity) {
    #   to.plot <- to.plot[expr>=2 & chromvar_zscore>=2] 
    # }
    
    to.plot.text <- to.plot[expr>=input$min_rna_threshold_tf_signature & chromvar_zscore>=input$min_chromvar_threshold_tf_signature] 
    to.plot.dots <- to.plot[!gene%in%to.plot.text$gene] 
    
    if (input$RNA_vs_chromVAR_pseudobulk_per_celltype_colour_by=="activator_repressor") {
      to.plot$colour_by <- to.plot$cor_sign
    } else if (input$RNA_vs_chromVAR_pseudobulk_per_celltype_colour_by=="marker_strength") {
      # to.plot$colour_by <- to.plot$cor_sign
      stop("Not implemented")
    }
    
    p.scatter <- ggplot(to.plot, aes(x=chromvar_zscore, y=expr)) +
      ggrepel::geom_text_repel(data=to.plot.text, aes(label=gene), size=3, max.overlaps=Inf) +
      # scale_size_continuous(range = c(0.25,6)) +
      coord_cartesian(
        xlim = c(min(to.plot$chromvar_zscore)-0.1,max.chromvar+0.1),
        ylim = c(min(to.plot$expr)-0.1,max.expr+0.1)
      ) +
      labs(x="Motif accessibility (chromVAR+, z-score)", y="Gene expression") +
      guides(size="none", fill = guide_legend(override.aes = list(size=2))) +
      # scale_fill_brewer(palette="Dark2") +
      theme_classic() +
      theme(
        legend.position = "top",
        legend.title = element_blank(),
        axis.text = element_text(size=rel(0.75), color="black")
      )
    
    # Add dots
    if (length(input$tf_tf_signatures)>0) {
      to.plot[,selected:=gene%in%input$tf_tf_signatures]
      p.scatter <- p.scatter +
        geom_jitter_interactive(data=to.plot, aes(tooltip=gene, data_id=gene, alpha=selected, fill=colour_by), size=5, shape=21, width=0.35, height=0.35) +
        scale_alpha_manual(values=c(0.30,1)) + guides(alpha="none")
    } else {
      p.scatter <- p.scatter + 
        geom_jitter_interactive(aes(tooltip=gene, data_id=gene, fill=colour_by), size=5, shape=21, width=0.35, height=0.35)
        # geom_jitter_interactive(aes(tooltip=gene, data_id=gene, size=dot_size, fill=colour_by), shape=21, width=0.35, height=0.35) +
    }
    
    ## Plot stacked barplots ##
    if (length(input$tf_tf_signatures)>0) {
      
      to.plot <- rna_vs_chromvar_pseudobulk.dt[gene%in%input$tf_tf_signatures] %>%
        .[chromvar_zscore<0,chromvar_zscore:=0] %>%
        .[chromvar_zscore>=max.chromvar,chromvar_zscore:=max.chromvar] %>%
        .[expr>=max.expr,expr:=max.expr] %>%
        .[,chromvar_zscore:=chromvar_zscore/sum(chromvar_zscore), by="gene"] %>%
        .[,expr:=exp(expr)/sum(exp(expr)), by="gene"] %>%
        melt(id.vars=c("celltype","gene"), variable.name="modality")
        
      to.plot <- to.plot[value>0.05]
      facet.labels <- c(expr = "RNA expression", chromvar_zscore = "Motif accessibility (z-score)")
        
      to.plot[,tooltip:=sprintf("%s:%s",gene,celltype)]
      
      p.barplot <- ggplot(to.plot, aes(x=gene, y=value, fill=celltype)) +
        geom_bar_interactive(aes(tooltip=tooltip, data_id=gene), stat="identity", color="black", position="fill") +
        celltype_palette_fill +
        facet_wrap(~modality, ncol=2, labeller = as_labeller(facet.labels)) +
        # scale_fill_manual(values=opts$celltype.colors[names(opts$celltype.colors)%in%unique(to.plot$celltype)]) +
        theme_classic() +
        labs(x="", y="") +
        theme(
          legend.position = "none",
          legend.title = element_blank(),
          # axis.text.x = element_text(color="black", size=rel(0.75)),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank()
        )
    }
    
    layout <-
    "
    AAB
    AAB
    AAB
    CCC
    CCC
    "
    
    girafe(
      # code = ifelse(length(input$tf_tf_signatures)>0, print(p.scatter+p.barplot + plot_layout(nrow=1, widths=c(2.5,1))), print(p.scatter)),
      code = print(p.scatter+p.barplot+p.lineplot + plot_layout(design = layout)),
      width_svg = 12, height_svg = 10,
      options = list(
        opts_zoom(min = 1, max = 3),
        opts_sizing(rescale = FALSE),
        opts_selection(type = "single", css = ""),
        opts_hover_inv(css = "opacity:0.25"),
        opts_hover(css = "")
      )
    ) %>% return(.)
    
  })
  
  output$rna_vs_chromvar_per_celltype <- renderGirafe({
    shiny::validate(need(input$celltype_tf_signatures%in%celltypes, "Please select one or multiple celltypes from our annotation." ))
    return(plot_RNA_vs_chromVAR_pseudobulk_per_celltype())
  })
  
  # output$print_celltype_tf_signatures <- renderPrint({
  #   rna_vs_chromvar_pseudobulk.dt %>%
  #     .[expr>=input$min_rna_threshold_tf_signature & chromvar_zscore>=input$min_chromvar_threshold_tf_signature] %>%
  #     .[celltype==input$celltype_tf_signatures] %>% 
  #     .[,tmp:=expr*chromvar_zscore] %>% setorder(-tmp) %>%
  #     .$gene %>% unique
  # })
    
  ########################
  ## In silico ChIP-seq ##
  ########################
  
  plot_InSilicoChIP = reactive({
    
    # input$tf_chip <- c("FOXA2","TAL1")
    # input$min_score_filt <- 0.25
    virtual_chip.dt <- input$tf_chip %>% map(function(i) {
      fread(sprintf("%s/insilico_chipseq/%s.txt.gz",basedir,i)) %>%
        setnames(c("chr","start","end","score","correlation_score","max_accessibility_score","motif_score","motif_counts")) %>%
        .[!is.na(score)] %>%
        .[,idx:=sprintf("%s:%s-%s",chr,start,end)] %>%
        .[,c("chr","start","end"):=NULL] %>%
        .[,tf:=i] %>%
        return
    }) %>% rbindlist
    
    ## Plot minimum score vs number of peaks ##
    
    to.plot <- insilico_chip_stats[tf%in%input$tf_chip]
    
    p1 <- ggplot(to.plot, aes_string(x="min_score", y="log2_N", color="tf")) +
      geom_line_interactive(size=1.5, aes(data_id=tf)) +
      labs(y="Number of predicted binding sites (log2)", x="Minimum score") +
      geom_vline(xintercept=input$min_score_filt, linetype="dashed") +
      scale_x_continuous(breaks=seq(0,0.75,by=0.10)) +
      theme_classic() +
      theme(
        axis.text = element_text(color="black"),
        legend.position = "top",
        legend.title = element_blank()
      )
    
    ## Plot fraction of positively correlated peaks ##

    to.plot <- virtual_chip.dt[abs(score)>=input$min_score_filt] %>% 
      .[,sign:=factor(c("-","+"))[(sign(score)>0)+1]] %>%
      .[,.N,by=c("tf","sign")]
    
    p2 <- ggplot(to.plot, aes(x=tf, y=N, fill=sign)) +
      geom_bar_interactive(stat="identity", aes(data_id=tf), color="black") +
      labs(x="", y="Number of putative binding sites") +
      scale_fill_brewer(palette="Dark2") +
      theme_classic() +
      theme(
        legend.title = element_blank(),
        axis.text.y = element_text(size=rel(0.75), color="black"),
        axis.text.x = element_text(size=rel(1), color="black")
      )
    
    # plot_grid(plotlist=list(p1,p2), scale=0.95, ncol=2, rel_widths = c(1/2,1/2)) %>% return(.)
    
    girafe(
      code = print(p1+p2),
      width_svg = 12, height_svg = 8,
      options = list( 
        opts_sizing(rescale = FALSE),
        opts_hover_inv(css = "opacity:0.2;"),
        opts_hover(css = "")
      )
    ) %>% return(.)
  })
  
  output$chip_stats_plot = renderGirafe({
    shiny::validate(need(input$tf_chip%in%TFs, "Please select one or multiple TFs from our annotation." ))
    return(plot_InSilicoChIP())
  })
  
  ## Download ChIP-seq data ##
  roots1 <- c(wd = file.path(basedir,'insilico_chipseq/bed'))
  shinyFileChoose(input, 'files', roots = roots1, filetypes=c('gz'))
  
  output$download_chip_bed_files_print <- renderPrint({
    if (length(input$files)>1) {
      tmp <- sapply(input$files[["files"]], function(x) x[[2]]) %>% unname
      cat(paste(tmp,collapse="\n"))
    }
  })
                
  output$download_chip_bed <- downloadHandler(
    filename = function() {
      files.to.download <- as.character(parseFilePaths(roots1, input$files)$datapath)
      if (length(files.to.download)>1) {
        "insilico_chip_seq.zip"
      } else {
        input$files[[1]][[1]][[2]]
      }
    },
    content = function(file) {
      files.to.download <- as.character(parseFilePaths(roots1, input$files)$datapath)
      if (length(files.to.download)>1) {
        utils::zip(file, files=files.to.download, flags = "-jr9X")
        # tar(file, files.to.download)
      } else {
        file.copy(from=files.to.download, to=file)
      }
    },
    contentType = "application/zip"
  )
  
  ############################
  ## Celltype-specific GRNs ##
  ############################
  
  plot_GRN_celltype <- reactive({
    
    # to.plot <- insilico_chip_stats[tf%in%input$tf_chip]
    file <- sprintf("%s/gene_regulatory_networks/per_celltype/%s_network.rds",basedir,input$celltype_grn)
    net <- readRDS(file)
    
    TFs_grn <- names(V(net))[V(net)$class=="TF"]
    genes_grn <- names(V(net))[V(net)$class=="gene"]
    
    pal.TFs <- grDevices::colorRamp(c("gray62", "purple"))( (1:100)/100 )
    pal.TFs <- cbind(pal.TFs, seq(100, 255, length.out = 100))
    pal.genes <- grDevices::colorRamp(c("gray60", "darkgreen"))( (1:100)/100 )
    pal.genes <- cbind(pal.genes, seq(100, 255, length.out = 100))
    
    # input$celltype_grn
    if (input$celltype_grn_colorby=="eigenvalue_centrality") {
      eigenvalue_centrality_scores <- igraph::eigen_centrality(net)$vector[TFs_grn]
      color_TF_nodes <- colourvalues::colour_values(eigenvalue_centrality_scores, palette = pal.TFs)
      names(color_TF_nodes) <- TFs_grn
      color_genes_nodes <- rep("#FFFFFF",length(genes_grn))
      names(color_genes_nodes) <- genes_grn
      V(net)$color <- c(color_TF_nodes, color_genes_nodes)[names(V(net))]
    } else if (input$celltype_grn_colorby=="degree_centrality") {
      degree_centrality_scores <- igraph::degree(net)[TFs_grn]
      color_TF_nodes <- colourvalues::colour_values(degree_centrality_scores, palette = pal.TFs)
      names(color_TF_nodes) <- TFs_grn
      color_genes_nodes <- rep("#FFFFFF",length(genes_grn))
      names(color_genes_nodes) <- genes_grn
      V(net)$color <- c(color_TF_nodes, color_genes_nodes)[names(V(net))]
    } else if (input$celltype_grn_colorby=="expression") {
      expression_values <- V(net)$expr; names(expression_values) <- names(V(net))
      expression_values.tfs <- expression_values[TFs_grn]
      expression_values.genes <- expression_values[genes_grn]
      colors.tfs <- colourvalues::colour_values(minmax.normalisation(expression_values.tfs), palette = pal.TFs)
      names(colors.tfs) <- names(expression_values.tfs)
      colors.genes <- colourvalues::colour_values(minmax.normalisation(expression_values.genes), palette = pal.genes)
      names(colors.genes) <- names(expression_values.genes)
      V(net)$color <- c(colors.tfs,colors.genes)[names(V(net))]
    }
    
    # Modify vertex attributes
    V(net)$group <- V(net)$class
    V(net)$shape <- stringr::str_replace_all(V(net)$class,c("gene"="triangle","TF"="circle"))
    
    visIgraph(net, randomSeed=42) %>%
      visIgraphLayout(randomSeed=42, physics = TRUE) %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE, selectedBy = "group") %>%
      visNodes(shadow = TRUE) %>%
      visEdges(width = 1, color = "black", smooth = FALSE) %>%
      visGroups(groupname = "TF", shape = "circle", size=35, font = list(color="black", size=30), color = list(background="#63B8FF", hover="#4876FF", border="#63B8FF")) %>%
      # visGroups(groupname = "TF", shape = "circle", size=35, font = list(color="black", size=30), color = list(background=color_TF_nodes)) %>%
      visGroups(groupname = "gene", shape = "triangle", size=20, font = list(color="black", size=35), color = list(background="#EE6363", hover="red", border="#EE6363")) %>%
      visPhysics(
        solver = "forceAtlas2Based",
        minVelocity = 0.50,
        forceAtlas2Based = list(gravitationalConstant = -250),
        stabilization = TRUE # By default, vis.js computes coordinates dynamically and waits for stabilization before rendering
      ) %>%
      visInteraction(
        hover = TRUE,
        dragNodes = TRUE,
        dragView = TRUE,
        zoomView = TRUE
      )
    })
  
  output$network <- renderVisNetwork({
    # shiny::validate(need(input$tf_chip%in%TFs, "Please select one or multiple TFs from our annotation." ))
    return(plot_GRN_celltype())
  })
  
  
  plot_GRN_celltype_stats <- reactive({
    
    # input$celltype_grn <- "Gut"
    file <- sprintf("%s/gene_regulatory_networks/per_celltype/%s_network.rds",basedir,input$celltype_grn)
    net <- readRDS(file)
    
    TFs_grn <- names(V(net))[V(net)$class=="TF"]
    
    if (input$celltype_grn_colorby=="eigenvalue_centrality") {
      eigenvalue_centrality_scores <- igraph::eigen_centrality(net)$vector[TFs_grn]
      to.plot <- data.frame(value=eigenvalue_centrality_scores) %>% as.data.table(keep.rownames = T) %>% setnames("rn","tf") %>%
        setorder(-value) %>% .[,tf:=factor(tf, levels=tf)]
      ylabel <- "Eigenvalue centrality"
    } else if (input$celltype_grn_colorby=="degree_centrality") {
      degree_centrality_scores <- igraph::degree(net)[TFs_grn]
      to.plot <- data.frame(value=degree_centrality_scores) %>% as.data.table(keep.rownames = T) %>% setnames("rn","tf") %>%
        setorder(-value) %>% .[,tf:=factor(tf, levels=TFs.to.plot)]
      ylabel <- "Degree centrality"
    } else if (input$celltype_grn_colorby=="expression") {
      expression_values <- V(net)$expr; names(expression_values) <- names(V(net))
      to.plot <- data.frame(value=expression_values[TFs_grn]) %>% as.data.table(keep.rownames = T) %>% setnames("rn","tf") %>%
        setorder(-value) %>% .[,tf:=factor(tf, levels=tf)]
      ylabel <- "RNA expression (scaled)"
    }
    
    ggbarplot(to.plot, x="tf", y="value", fill="value", width=0.55) +
      coord_flip() +
      scale_fill_gradient(low = "gray62", high = "purple") +
      # scale_fill_manual(values=c("-"="blue", "+"="red"), drop=F) +
      labs(x="", y=ylabel) +
      theme(
        axis.text.y = element_text(colour="black", size=rel(1)),
        axis.text.x = element_text(colour="black", size=rel(0.8)),
        axis.ticks.x = element_line(size=rel(0.75)),
        legend.position = "none"
      )
  })
  
  output$grn_celltype_stats = renderPlot({
    # shiny::validate(need(input$tf_chip%in%TFs, "Please select one or multiple TFs from our annotation." ))
    return(plot_GRN_celltype_stats())
  })
  
  ## Download GRNs ##
  roots <- c(wd = file.path(basedir,'gene_regulatory_networks/per_celltype'))
  shinyFileChoose(input, 'grn_files', roots =  roots, filetypes=c('rds'))
  
  output$download_grn_files_print <- renderPrint({
    if (length(input$grn_files)>1) {
      tmp <- sapply(input$grn_files[["files"]], function(x) x[[2]]) %>% unname %>% stringr::str_replace(".rds","")
      cat(paste(tmp,collapse=".rds\n"))
    }
  })
  
  output$download_grn_files <- downloadHandler(
    filename = function() {
      files.to.download <- as.character(parseFilePaths(roots, input$grn_files)$datapath)
      if (length(files.to.download)>1) {
        "gene_regulatory_networks_per_celltype.zip"
      } else {
        input$grn_files[[1]][[1]][[2]]
      }
    },
    content = function(file) {
      files.to.download <- as.character(parseFilePaths(roots, input$grn_files)$datapath)
      if (length(files.to.download)>1) {
        utils::zip(file, files=files.to.download, flags = "-jr9X")
      } else {
        file.copy(from=files.to.download, to=file)
      }
    },
    contentType = "application/zip"
  )

  ###########################
  ## Differential analysis ##
  ###########################
  
  selected_gene_differential <- reactive({
    input$differential_analysis_selected
  })
  
  plot_differential <- reactive({
    
    # Define settings for testing
    # input$celltypeA <- "Gut"; input$celltypeB <- "Neural_crest"; input$differential_range <- c(-25,25); input$modality_differential <- "atac"
    
    # Load precomputed differential analysis results
    diff_rna.dt <- fread(file.path(basedir,"diff_rna_pseudobulk.txt.gz"))
    diff_chromvar.dt <- fread(file.path(basedir,"diff_chromVAR_pseudobulk.txt.gz"))
    
    
    # Select cell types of interest
    diff_rna_filt.dt <- diff_rna.dt[groupA==input$celltypeA & groupB==input$celltypeB] 
    if (nrow(diff_rna_filt.dt)==0) {
      diff_rna_filt.dt <- diff_rna.dt[groupA==input$celltypeB & groupB==input$celltypeA]  %>% 
        setnames(c("gene","diff","groupB","groupA")) %>% .[,diff:=-diff]
    }
    
    diff_chromvar_filt.dt <- diff_chromvar.dt[groupA==input$celltypeA & groupB==input$celltypeB] 
    if (nrow(diff_chromvar_filt.dt)==0) {
      diff_chromvar_filt.dt <- diff_chromvar.dt[groupA==input$celltypeB & groupB==input$celltypeA]  %>% 
        setnames(c("gene","diff","groupB","groupA")) %>% .[,diff:=-diff]
    }

  to.plot <- merge(diff_rna_filt.dt, diff_chromvar_filt.dt, by=c("gene","groupA","groupB"), suffixes=c("_rna","_chromvar"))
    
  # if (input$highlight_top_n_genes>0) {
  #   to.plot2 <- to.plot %>% head(n=input$highlight_top_n_genes)
  # } else {
  #   to.plot2 <- to.plot
  # }
  
  to.plot <- to.plot %>% 
    .[,dot_size:=minmax.normalisation(abs(diff_rna)*abs(diff_chromvar))] %>% 
    sort.abs("dot_size")
  
  to.plot[,tooltip:=sprintf("%s\ndiff_rna = %.2f\ndiff_chromvar = %.2f|\n Higher expression in %s\nHigher motif acc in %s",gene,diff_rna,diff_chromvar,ifelse(diff_rna>0,groupB,groupA),ifelse(diff_chromvar>0,groupB,groupA))]
  
  p.scatter <- ggplot(to.plot, aes(x=diff_rna, y=diff_chromvar)) +
    geom_jitter_interactive(aes(size=dot_size, data_id=gene, tooltip=tooltip, onclick=gene), width=0.15, height=0.5) +
    scale_size_continuous(range = c(0.25,4)) +
    # geom_point(size=0.5, data=to.plot.i[sig==F], alpha=0.2, color="grey") +
    geom_segment(x=min(to.plot$diff_rna,na.rm=T), xend=max(to.plot$diff_rna,na.rm=T), y=0, yend=0, size=0.35, color="orange") +
    geom_segment(x=0, xend=0, y=min(to.plot$diff_chromvar,na.rm=T), yend=max(to.plot$diff_chromvar,na.rm=T), size=0.35, color="orange") +
    ggrepel::geom_text_repel(data=head(to.plot,n=30), aes(x=diff_rna, y=diff_chromvar, label=gene), size=3, max.overlaps = Inf) +
    labs(x="Differential RNA expression", y="Differential motif acc. (chromVAR+ z-score)") +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text = element_text(size=rel(0.8), color="black"),
      axis.title = element_text(size=rel(1), color="black")
    )
  
  ## PAGA ##
  colors <- celltype_colours[celltypes]
  colors[!names(colors)%in%c(input$celltypeA,input$celltypeB)] <- "gray90"

  alpha <- rep(1,length(celltypes))
  names(alpha) <- celltypes
  alpha[!names(alpha)%in%c(input$celltypeA,input$celltypeB)] <- 0.35

  size <- rep(6,length(celltypes))
  names(size) <- celltypes
  size[!names(size)%in%c(input$celltypeA,input$celltypeB)] <- 1.5

  text_size <- rep(3,length(celltypes))
  names(text_size) <- celltypes
  text_size[!names(text_size)%in%c(input$celltypeA,input$celltypeB)] <- 0

  p.paga <- ggnet2(
    net = paga,
    mode = c("x", "y"),
    color = colors,
    # color = celltype.colors[celltypes],
    node.size = size,
    edge.size = 0.15,
    edge.color = "grey",
    alpha = alpha,
    label = TRUE,
    label.size = text_size
  ) + guides(size="none")
  
  
  p.paga_base <- ggnet2(
    net = paga,
    mode = c("x", "y"),
    node.size = 0,
    edge.size = 0.15,
    edge.color = "grey",
    label = FALSE
  )
  
  tmp <- rna_vs_chromvar_pseudobulk.dt[gene==input$tf_differential]
  
  # Define colors  
  rna.col.seq <- chromvar.col.seq <- round(seq(0,1,0.1), 2)
  rna.colors <- colorRampPalette(c("gray92", "darkgreen"))(length(rna.col.seq))
  chromvar.colors <- colorRampPalette(c("gray92", "purple"))(length(chromvar.col.seq)) 
  values.rna <- minmax.normalisation(tmp$chromvar_zscore); names(values.rna) <- tmp$celltype; values.rna <- values.rna[celltypes]
  colors.rna <- round(values.rna,1) %>% map(~ chromvar.colors[which(chromvar.col.seq == .)]) %>% unlist
  values.chromvar <- minmax.normalisation(tmp$expr); names(values.chromvar) <- tmp$celltype; values.chromvar <- values.chromvar[celltypes]
  colors.chromvar <- round(values.chromvar,1) %>% map(~ rna.colors[which(rna.col.seq == .)]) %>% unlist
  
  p.paga_rna <- p.paga_base + 
    geom_text(label = "\u25D0", aes(x=x, y=y), color=colors.rna, size=12, family = "Arial Unicode MS", data = p.paga_base$data[,c("x","y")])
  
  p.paga_chromvar <- p.paga_base + 
    geom_text(label = "\u25D0", aes(x=x, y=y), color=colors.chromvar, size=12, family = "Arial Unicode MS", data = p.paga_base$data[,c("x","y")])
  
  
  girafe(
    code = print(((p.scatter)|(p.paga/p.paga_rna/p.paga_chromvar)) + plot_layout(width=c(4,1))),
    width_svg = 12, height_svg = 8,
    options = list( 
      opts_zoom(min = 1, max = 4),
      opts_sizing(rescale = FALSE),
      opts_selection(type = "single", css = ""),
      opts_hover(css = "cursor:pointer;fill:magenta;color:magenta")
    )
  ) %>% return(.)
    
  })
  
  output$differential_motif <- renderPlot({
    shiny::validate(need(input$tf_differential%in%TFs, "Please select one or multiple TFs from our annotation." ))
    m <- 0.25*exp(as.matrix(motifs[[input$tf_differential]]))
    ggseqlogo(m) +
      theme(
        panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank()
      )
  }, bg="transparent", execOnResize = TRUE)
  
  output$differential_analysis <- renderGirafe({
    shiny::validate(need(input$tf_differential%in%TFs, "Please select one or multiple transcription factors from our annotation." ))
    shiny::validate(need(input$celltypeA%in%celltypes, "Please select one or multiple celltypes from our annotation." ))
    shiny::validate(need(input$celltypeB%in%celltypes, "Please select one or multiple celltypes from our annotation." ))
    return(plot_differential())
  })
  
  observeEvent(selected_gene_differential(), {
    updateTextInput(session = session, "tf_differential", value = selected_gene_differential() )
  })
  
  ##################
  ## trajectories ##
  ##################
  
  
  selected_gene_trajectories <- reactive({
    input$trajectory_rna_vs_chromvar_selected
  })
  
  plot_trajectories <- reactive({
      
    # Define settings for testing
    # input$trajectory <- "blood"; input$tf_trajectory <- "GATA1"
    i <- input$trajectory
    
    # Load data
    trajectory.dt <- fread(sprintf("%s/trajectories/%s_trajectory_knn50/trajectory.txt.gz",basedir,i))
    chromvar_rna.dt <- fread(sprintf("%s/trajectories/%s_trajectory_knn50/chromvar_rna.txt.gz",basedir,i)) %>%
      setnames("celltype.predicted","celltype")
    cor.dt <- fread(sprintf("%s/trajectories/%s_trajectory_knn50/correlation_results.txt.gz",basedir,i), select=c(1,3,4,5))
    
    ## PAGA ##
    colors <- celltype_colours[celltypes]
    colors[!names(colors)%in%trajectories[[i]]] <- "gray70"
    
    alpha <- rep(1,length(celltypes))
    names(alpha) <- celltypes
    alpha[!names(alpha)%in%trajectories[[i]]] <- 0.4
    
    size = rep(6,length(celltypes))
    names(size) <- celltypes
    size[!names(size)%in%trajectories[[i]]] <- 2.5
    
    text_size = rep(2.3,length(celltypes))
    names(text_size) <- celltypes
    text_size[!names(text_size)%in%trajectories[[i]]] <- 0
    
    p.paga <- ggnet2_interactive(
      net = paga,
      mode = c("x", "y"),
      color = colors,
      # color = celltype.colors[celltypes],
      node.size = size,
      edge.size = 0.15,
      edge.color = "grey",
      alpha = alpha,
      label = TRUE,
      label.size = text_size
    ) + guides(size="none")
    
    ## RNA vs chromVAR pseudotime plots ##
    
    to.plot <- chromvar_rna.dt[gene==input$tf_trajectory] %>%
      melt(id.vars=c("cell","PC1","celltype"), measure.vars=c("chromvar_zscore","expr"), variable.name="modality")
    
    facet.labels <- c(expr = "RNA expression", chromvar_zscore = "Motif accessibility (chromVAR+, z-score)")
    
    p.pseudotime <- ggplot(to.plot, aes(x=PC1, y=value)) +
      geom_point_interactive(aes(fill=celltype, tooltip=celltype, data_id=celltype), size=1.75, shape=21, stroke=0.1) +
      stat_smooth(method="loess", color="black", alpha=0.75, span=0.5) +
      geom_rug(aes(color=celltype), sides="b") +
      facet_wrap(~modality, nrow=2, scales="free_y", labeller = as_labeller(facet.labels)) +
      celltype_palette + celltype_palette_fill + 
      guides(fill="none", color="none") +
      labs(x="Pseudotime", y="") +
      theme_classic() +
      theme(
        axis.text.x = element_blank(),
        axis.text.y = element_text(color="black", size=rel(1.0)),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "none"
      )
    
    ## Volcano plot ##
    
    cor.dt[,log_pval:=-log10(padj_fdr+1e-100)]
      
    # xlim <- max(abs(cor.dt$r), na.rm=T)
    # ylim <- max(-log10(cor.dt$padj_fdr+1e-100), na.rm=T)
    xlim_min <- input$rna_vs_chromvar_cor_range[1]
    xlim_max <- input$rna_vs_chromvar_cor_range[2]

    to.plot <- cor.dt[r>=input$trajectory_rna_vs_chromvar_cor_range[1] & r<=input$trajectory_rna_vs_chromvar_cor_range[2]]
    ylim_min <- min(to.plot$log_pval,na.rm=T)
    ylim_max <- max(to.plot$log_pval,na.rm=T)
    
    negative_hits <- to.plot[sig==TRUE & r<0 & r>input$trajectory_rna_vs_chromvar_cor_range[1],gene]
    positive_hits <- to.plot[sig==TRUE & r>0 & r<input$trajectory_rna_vs_chromvar_cor_range[2],gene]
    all <- nrow(to.plot)
    
  
    if (all>=500) { dot_size <- 1.5 } else if (all>=100 & all<=500) { dot_size <- 2.5 } else if (all<100) { dot_size <- 5}
    
    p.volcano <- ggplot(to.plot, aes(x=r, y=log_pval)) +
      geom_segment(x=0, xend=0, y=0, yend=ylim_max-1, color="orange", size=0.25) +
      geom_jitter_interactive(aes(fill=log_pval, alpha=log_pval, tooltip=gene, data_id=gene, onclick=gene), width=0.05, height=3, shape=21, size=dot_size) + 
      scale_fill_gradient(low = "gray80", high = "red") +
      scale_alpha_continuous(range=c(0.25,1)) +
      scale_x_continuous(limits=c(xlim_min-0.35,xlim_max+0.35)) +
      scale_y_continuous(limits=c(ylim_min,ylim_max+5)) +
      annotate("text", x=median(to.plot$r,na.rm=T), y=ylim_max+5, size=4, label=sprintf("(%d)", all)) +
      annotate("text", x=xlim_min-0.20, y=ylim_max+4, size=4, label=sprintf("%d (-)",length(negative_hits))) +
      annotate("text", x=xlim_max+0.20, y=ylim_max+4, size=4, label=sprintf("%d (+)",length(positive_hits))) +
      labs(x="Pearson correlation\n(RNA expression vs Motif accessibility)", y=expression(paste("-log"[10],"(p.value)"))) +
      theme_classic() +
      theme(
        axis.text = element_text(size=rel(0.75), color='black'),
        axis.title = element_text(size=rel(1.0), color='black'),
        legend.position="none"
      )
      
    p.volcano <- p.volcano + ggrepel::geom_text_repel(data=to.plot[gene==input$tf_trajectory], aes(x=r, y=log_pval, label=gene), size=3)

    # add motif to the volcano plot
    # p.volcano <- p.volcano + inset_element(p.motif, left=0.7, top=0.20, right = 1, bottom = -0.05)
    
    girafe(
      code = print((p.paga/p.volcano)|p.pseudotime),
      width_svg = 11, height_svg = 8,
      options = list( 
        opts_zoom(min = 1, max = 3),
        opts_sizing(rescale = FALSE),
        # opts_selection(type = "single", css = "cursor:pointer;fill:magenta;color:magenta"),
        opts_selection(type = "single", css = ""),
        # opts_hover_inv(css = "opacity:0.45;"),
        opts_hover(css = "cursor:pointer;fill:magenta;color:magenta")
      )
    ) %>% return(.)
    
  })
  
  
  observeEvent(selected_gene_trajectories(), {
    updateTextInput(session = session, "tf_trajectory", value = selected_gene_trajectories() )
  })
  
  
  output$trajectory_motif <- renderPlot({
    shiny::validate(need(input$tf_trajectory%in%TFs, "Please select one or multiple TFs from our annotation." ))
    m <- 0.25*exp(as.matrix(motifs[[input$tf_trajectory]]))
    ggseqlogo(m) +
      theme(
        panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank()
      )
  }, bg="transparent", execOnResize = TRUE)
  
  output$trajectory_rna_vs_chromvar <- renderGirafe({
    shiny::validate(need(input$trajectory%in%names(trajectories), "Please select one or multiple trajectories from our annotation." ))
    shiny::validate(need(input$tf_trajectory%in%TFs, "Please select one or multiple TFs from our annotation." ))
    return(plot_trajectories())
  })
  
  output$trajectory_rna_vs_chromvar_print <- renderPrint({
    if (length(input$tf_trajectory)>0) {
      tmp <- fread(sprintf("%s/trajectories/%s_trajectory_knn50/correlation_results.txt.gz",basedir,input$trajectory), select=c(1,3,4,5)) %>%
        .[r>=input$trajectory_rna_vs_chromvar_cor_range[1] & r<=input$trajectory_rna_vs_chromvar_cor_range[2]] %>%
        .[,abs_r:=abs(r)] %>% setorder(-abs_r) %>% head(n=15)
      cat("Top 15 TFs with largest absolute RNA vs chromVAR+ correlation (among the selected range):\n")
      cat(sprintf("%s: r=%.2f\t", tmp$gene, round(tmp$r,2)))
    }
  })
  
}
