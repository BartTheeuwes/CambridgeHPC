
#####################
## Define settings ##
#####################

# Load default settings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/utils.R")
} else {
  stop("Computer not recognised")
}

# Options
opts$celltypes = c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  "PGC",
  "Anterior_Primitive_Streak",
  "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  "Mixed_mesoderm",
  "Intermediate_mesoderm",
  "Caudal_Mesoderm",
  "Paraxial_mesoderm",
  "Somitic_mesoderm",
  "Pharyngeal_mesoderm",
  "Cardiomyocytes",
  "Allantois",
  "ExE_mesoderm",
  "Mesenchyme",
  "Haematoendothelial_progenitors",
  "Endothelium",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
  "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm",
  "Visceral_endoderm",
  "ExE_endoderm",
  "ExE_ectoderm",
  "Parietal_endoderm"
)

# opts$motif_annotation <- "Motif_cisbp"
opts$motif_annotation <- "Motif_cisbp_lenient"

# Scatterplots of TF expression vs peak accessibility?
opts$scatterplots <- FALSE

# I/O
io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_acc/pseudobulk/TFexpr_vs_peakAcc"); dir.create(io$outdir, showWarnings = F)
# io$motifmatcher.se <- sprintf("%s/Annotations/%s-Matches-In-Peaks.rds",io$archR.directory,opts$motif_annotation)
if(opts$scatterplots) dir.create(paste0(io$outdir,"/scatterplots"), showWarnings = F)

##################################
## Load pseudobulk RNA and ATAC ##
##################################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/rna_atac/load_rna_atac_pseudobulk.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/rna_atac/load_rna_atac_pseudobulk.R")
} else {
  stop("Computer not recognised")
}

###############################
## Load motifmatcher results ##
###############################

stop("Make sure NKX is there")

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/load_motifmatchR.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/load_motifmatchR.R")
} else {
  stop("Computer not recognised")
}

# Subset peaks
# motifmatcher.se <- motifmatcher.se[rownames(atac.peakMatrix.se),]

#########################################################
## Correlate peak accessibility with TF RNA expression ##
#########################################################

stopifnot(colnames(rna.sce.tf)==colnames(atac.peakMatrix.se))

TFs <- intersect(colnames(motifmatcher.se),rownames(rna.sce.tf))# %>% head(n=5)
# TFs <- c("GATA1","TAL1")

# Prepare output data objects
cor.dt <- list()
cor.mtx <- matrix(as.numeric(NA), nrow=nrow(atac.peakMatrix.se), ncol=length(TFs))
pvalue.mtx <- matrix(as.numeric(NA), nrow=nrow(atac.peakMatrix.se), ncol=length(TFs))
rownames(cor.mtx) <- rownames(atac.peakMatrix.se); colnames(cor.mtx) <- TFs
dimnames(pvalue.mtx) <- dimnames(cor.mtx)


for (i in TFs) {
  print(i)
  all_peaks_i <- rownames(motifmatcher.se)[which(assay(motifmatcher.se[,i],"motifMatches")==1)]
  # TF_peak.cor <- cor(assay(rna.sce.tf[i,])[1,], t(assay(atac.peakMatrix.se[all_peaks_i,])))[1,]
  
  # calculate correlations
  corr_output <- psych::corr.test(t(logcounts(rna.sce.tf[i,])), t(assay(atac.peakMatrix.se[all_peaks_i,])), ci=FALSE)
  
  # Fill matrices
  cor.mtx[all_peaks_i,i] <- round(corr_output$r[1,],3)
  pvalue.mtx[all_peaks_i,i] <- round(corr_output$p[1,],5)
  
  # Define significant peaks
  # sig_peaks <- which(corr_output$p[1,]<0.10 & abs(corr_output$r[1,])>0.25)
  sig_peaks <- which(corr_output$p[1,]<0.10)
  
  if (length(sig_peaks)>0) {
    
    # Plot
    if (opts$scatterplots) {
      for (j in names(sig_peaks)) {

        to.plot <- data.table(
          rna = logcounts(rna.sce.tf[i,])[1,], 
          atac = assay(atac.peakMatrix.se[j,])[1,],
          celltype = colnames(rna.sce.tf)
        )
        
        p <- ggscatter(to.plot, x="rna", y="atac", fill="celltype", size=4, shape=21, 
                        add="reg.line", add.params = list(color="black", fill="lightgray"), conf.int=TRUE) +
          stat_cor(method = "pearson") +
          scale_fill_manual(values=opts$celltype.colors) +
          labs(x=sprintf("%s expression",i), y="chromatin accessibility", title=j) +
          guides(fill=F) +
          theme(
            plot.title = element_text(hjust = 0.5, size=rel(0.85)),
            axis.text = element_text(size=rel(0.7))
          )
        
        pdf(sprintf("%s/scatterplots/%s_expr_vs_%s_accessibility_pseudobulk.pdf",io$outdir,i,gsub(":","-",j)), width = 6, height = 5)
        print(p)
        dev.off()
        
      }
    }
    
    # Return data.table
    cor.dt[[i]] <- data.table(
      TF = i, 
      peak = names(sig_peaks), 
      cor = round(corr_output$r[1,][sig_peaks],2),  
      p = round(corr_output$p[1,][sig_peaks],5)
    ) %>% setorder(p)
  }
}

# Save data.table with significant correlations
cor.dt <- rbindlist(cor.dt)
fwrite(cor.dt, paste0(io$outdir,"/cor_TFexpr_vs_peakAcc_lenient.txt.gz"), sep="\t", quote=F)

# Save SummarizedExperiment object
to.save <- SummarizedExperiment(
  assays = SimpleList("cor" = dropNA(cor.mtx), "pvalue" = dropNA(pvalue.mtx)),
  rowData = rowData(atac.peakMatrix.se)
)
saveRDS(to.save, paste0(io$outdir,"/cor_TFexpr_vs_peakAcc_SummarizedExperiment_lenient_v2.rds"))

##########
## TEST ##
##########

# target_peaks_i <- which(!is.na(dropNA2matrix(assay(to.save[,i],"pvalue"))[,1]))

