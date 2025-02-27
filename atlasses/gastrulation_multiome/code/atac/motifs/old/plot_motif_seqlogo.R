library(ggseqlogo)

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

io$outdir <- paste0(io$basedir,"/results/atac/archR/motif_seqlogo")

########################
## Load ArchR Project ##
########################

source(here::here("atac/archR/load_archR_project.R"))

#########################
## Load sequence Logos ##
#########################

motifs <- getPeakAnnotation(ArchRProject)[["motifs"]]

# Rename motifs
names(motifs) <- names(motifs) %>% toupper %>% stringr::str_split(.,"_") %>% map_chr(1)
names(motifs) <- gsub("TCFAP","TFAP",names(motifs))
names(motifs) <- gsub("NKX2","NKX2-",names(motifs))
names(motifs) <- gsub("NKX3","NKX3-",names(motifs))
names(motifs) <- gsub("NKX6","NKX6-",names(motifs))

# Remove duplicated motifs
motifs <- motifs[!duplicated(names(motifs))]
saveRDS(motifs, paste0(io$outdir,"/PWMatrixList.rds"))

# Find motifs
grep("KLF",names(motifs), value=T)

#########################
## Plot Sequence Logos ##
#########################

# Load precomputed PWMs
motifs <- readRDS(paste0(io$outdir,"/PWMatrixList.rds"))

# i <- "FOXA2"
for (i in names(motifs)) {
  m <- 0.25*exp(as.matrix(motifs[[i]]))
  # seqLogo::seqLogo(m)
  p <- ggseqlogo(m) + 
    theme(
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.y = element_blank()
    )
  pdf(sprintf("%s/seqlogo_%s.pdf",io$outdir,i), width=4, height=2)
  print(p)
  dev.off()
}

#######################
## Use JASPAR motifs ##
#######################

library(JASPAR2020)
library(TFBSTools)


# .summarizeJASPARMotifs is a function inside ArchR::AnnotationPeaks.R
jaspar_motifs <- .summarizeJASPARMotifs(getMatrixSet(JASPAR2020, list(species="Homo sapiens", collection="CORE", matrixtype="PWM")))$motifs

names(jaspar_motifs) <- names(jaspar_motifs) %>% toupper %>% stringr::str_split(.,"_") %>% map_chr(1)

grep("DUX",names(jaspar_motifs),value=T)

motifs.to.plot <- c("DUX4","KLF17")

for (i in motifs.to.plot) {
  m <- 0.25*exp(as.matrix(jaspar_motifs[[i]]))
  # seqLogo::seqLogo(m)
  p <- ggseqlogo(m) + 
    theme(
      axis.line = element_line(size=rel(0.5), color="black"),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size=rel(0.75)),
      axis.title.y = element_text(size=rel(0.75)),
      # axis.title.y = element_blank()
    )
  pdf(sprintf("%s/seqlogo_JASPAR_%s.pdf",io$outdir,i), width=5, height=2.2)
  print(p)
  dev.off()
}

########################
## Use chromVARmotifs ##
########################

# library(chromVARmotifs)
# 
# data("human_pwms_v2")
# grep("DUX",names(human_pwms_v2),value=T)
# 
# data("mouse_pwms_v1")
# data("mouse_pwms_v2")
# data("encode_pwms")
# data("homer_pwms")
# names(encode_pwms)
# grep("KLF1",names(encode_pwms),value=T)
