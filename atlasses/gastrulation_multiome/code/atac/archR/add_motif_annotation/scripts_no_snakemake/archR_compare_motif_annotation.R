#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
} else {
  stop("Computer not recognised")
}

###########################
## Load motif annotation ##
###########################

peakAnnotation.files <- readRDS(sprintf("%s/Annotations/peakAnnotation.rds",io$archR.directory))

motifmatcher1.se <- readRDS(peakAnnotation.files$Motif_cisbp$Matches)
motifmatcher2.se <- readRDS(peakAnnotation.files$Motif_cisbp_lenient$Matches)

motifpositions1.se <- readRDS(peakAnnotation.files$Motif_cisbp$Positions)
motifpositions2.se <- readRDS(peakAnnotation.files$Motif_cisbp_lenient$Positions)

# Rename TFs
colnames(motifmatcher1.se) <- colnames(motifmatcher1.se) %>% toupper %>% stringr::str_split(.,"_") %>% map_chr(1)
colnames(motifmatcher2.se) <- colnames(motifmatcher2.se) %>% toupper %>% stringr::str_split(.,"_") %>% map_chr(1)
names(motifpositions1.se) <- names(motifpositions1.se) %>% toupper %>% stringr::str_split(.,"_") %>% map_chr(1)
names(motifpositions2.se) <- names(motifpositions2.se) %>% toupper %>% stringr::str_split(.,"_") %>% map_chr(1)

# Rename peaks
tmp <- rowRanges(motifmatcher1.se)
rownames(motifmatcher1.se) <- sprintf("%s:%s-%s",seqnames(tmp), start(tmp), end(tmp))

tmp <- rowRanges(motifmatcher2.se)
rownames(motifmatcher2.se) <- sprintf("%s:%s-%s",seqnames(tmp), start(tmp), end(tmp))

mean(assay(motifmatcher1.se)==1)
mean(assay(motifmatcher2.se)==1)

#############
## Explore ##
#############

foo <- rowSums(assay(motifmatcher1.se))
bar <- rowSums(assay(motifmatcher2.se))

to.plot <- data.table(motifmatchr1=foo, motifmatcher2=bar) %>%
  melt()

gghistogram(to.plot, x="value", fill="variable", bins=60)



foo <- colSums(assay(motifmatcher1.se))
bar <- colSums(assay(motifmatcher2.se))

to.plot <- data.table(motifmatchr1=foo, motifmatcher2=bar) %>% melt()

gghistogram(to.plot, x="value", fill="variable", bins=60)



assay(motifmatcher2.se)[,"TAL1"][,1]

motifpositions2.se[["TAL1"]]


assay(motifmatcher1.se)["chr8:122719672-122720272","TAL1"]
assay(motifmatcher2.se)["chr8:122719672-122720272","TAL1"]
