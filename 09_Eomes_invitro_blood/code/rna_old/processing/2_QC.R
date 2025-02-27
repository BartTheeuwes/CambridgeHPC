here::i_am("rna/processing/2_QC.R")

source(here::here("settings.R"))

#####################
## Define arguments ##
#####################

p <- ArgumentParser(description='')
p$add_argument('--metadata',       type="character",                    help='Metadata')
p$add_argument('--outputdir',       type="character",                    help='Output directory')
p$add_argument('--min_nFeature_RNA',       type="integer",                    help='Minimum number of expressed genes')
p$add_argument('--max_nFeature_RNA',       type="integer",                    help='Maximum number of expressed genes')
p$add_argument('--mitochondrial_percent_RNA',       type="integer",                    help='Maximum percentage of mitochondrial reads')
p$add_argument('--ribosomal_percent_RNA',       type="integer",                    help='Maximum percentage of ribosomal reads')
args <- p$parse_args(commandArgs(TRUE))


#####################
## Define settings ##
#####################

## START TEST ##
# args <- list()
# args$outputdir <- paste0(io$basedir,"/results/rna_new/qc")
# args$min_nFeature_RNA <- 500
# args$max_nFeature_RNA <- 8000
# args$mitochondrial_percent_RNA <- 10
# args$ribosomal_percent_RNA <- 10
# args$metadata <- paste0(io$basedir,"/processed/rna_new/metadata.txt.gz")
## END TEST ##

dir.create(args$outputdir, showWarnings=F)

###############
## Load data ##
###############

metadata <- fread(args$metadata) %>% 
    .[,pass_rnaQC:=nFeature_RNA<=args$max_nFeature_RNA & nFeature_RNA>=args$min_nFeature_RNA & mitochondrial_percent_RNA<args$mitochondrial_percent_RNA & ribosomal_percent_RNA<args$ribosomal_percent_RNA]

table(metadata$pass_rnaQC)

#####################
## Plot QC metrics ##
#####################

to.plot <- metadata %>% copy %>%
    # .[,log_nCount_RNA:=log2(nCount_RNA)] %>%
    melt(id.vars=c("sample","cell"), measure.vars=c("nFeature_RNA","mitochondrial_percent_RNA","ribosomal_percent_RNA"))

## Box plot 

# p <- ggboxplot(to.plot, x="sample", y="value") +
#     facet_wrap(~variable, scales="free_y", nrow=1) +
#     theme(
#         axis.text.x = element_text(colour="black",size=rel(0.45), angle=20, hjust=1, vjust=1),  
#         axis.title.x = element_blank()
#     )
# 
# pdf(sprintf("%s/qc_metrics_boxplot.pdf",args$outputdir), width=12, height=6)
# # pdf(sprintf("%s/qc_metrics_boxplot.pdf",args$outputdir))
# print(p)
# dev.off()

## histogram 

tmp <- data.table(
    variable = c("nFeature_RNA", "mitochondrial_percent_RNA", "ribosomal_percent_RNA"),
    value = c(args$min_nFeature_RNA, args$mitochondrial_percent_RNA, args$ribosomal_percent_RNA)
)

p <- gghistogram(to.plot, x="value", fill="sample", bins=50) +
    geom_vline(aes(xintercept=value), linetype="dashed", data=tmp) +
    facet_wrap(~variable, scales="free", nrow=1) +
    theme(
        axis.text =  element_text(size=rel(0.8)),
        axis.title.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(size=rel(0.5))
    )
    
pdf(sprintf("%s/qc_metrics_histogram.pdf",args$outputdir), width=13, height=6)
# pdf(sprintf("%s/qc_metrics_histogram.pdf",args$outputdir))
print(p)
dev.off()


########################################################
## Plot number of cells that pass QC for each sample ##
########################################################

to.plot <- metadata %>%
    .[,sum(pass_rnaQC),by="sample"]

p <- ggbarplot(to.plot, x="sample", y="V1", fill="gray70") +
    labs(x="", y="Number of cells that pass QC") +
    # facet_wrap(~stage)
    theme(
        legend.position = "none",
        axis.text.y = element_text(colour="black",size=rel(0.8)),
        axis.text.x = element_text(colour="black",size=rel(0.65), angle=20, hjust=1, vjust=1),  
    )

pdf(sprintf("%s/qc_metrics_barplot.pdf",args$outputdir))
print(p)
dev.off()

##########
## Save ##
##########

fwrite(metadata, paste0(args$outputdir,"/sample_metadata_after_qc.txt.gz"), quote=F, na="NA", sep="\t")

