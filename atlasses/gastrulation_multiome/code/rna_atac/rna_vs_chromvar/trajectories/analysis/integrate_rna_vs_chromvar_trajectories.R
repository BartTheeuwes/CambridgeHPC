library(pheatmap)
library(RColorBrewer)

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

# I/O
io$outdir <- paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/trajectories/integrated"); dir.create(io$outdir, showWarnings = F)
io$trajectories.inputdir <- c(
  "ectoderm" = paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/trajectories/ectoderm_trajectory_knn50"),
  "endoderm" = paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/trajectories/endoderm_trajectory_knn50"),
  "mesoderm" = paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/trajectories/mesoderm_trajectory_knn50"),
  "blood" = paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/trajectories/blood_trajectory_knn50")
)

# Options

###############
## Load data ##
###############

# Load data.tables with RNA + chromVAR estimates
# chromvar_rna_dt <- names(io$trajectories.inputdir) %>% 
#   map(function(x) fread(sprintf(sprintf("%s/chromvar_rna_%s.txt.gz",io$trajectories.inputdir[[x]],x))) %>%
#         .[,trajectory:=x]) %>%
#   rbindlist

# Load correlation results
cor_dt <- names(io$trajectories.inputdir) %>% 
  map(function(x) fread(sprintf(sprintf("%s/correlation_results.txt.gz",io$trajectories.inputdir[[x]]))) %>%
        .[,trajectory:=x]) %>%
  rbindlist

################
## Parse data ##
################

opts$min.cor <- 0.35
opts$max.pval <- 0.01

# Filter correlation results by statistical significance
cor_dt.filt <- cor_dt %>% copy %>%
  .[,sig:=padj_fdr<opts$max.pval & abs(r)>=opts$min.cor] %>%
  .[,foo:=sum(sig==TRUE),by=c("gene")] %>% .[foo>=1] %>% .[,foo:=NULL] %>%
  .[,class:=c("Repressor","Activator")[as.numeric(r>0)+1]]
  
  
# Create matrix
# cor.mtx <- cor_dt.filt %>% 
#   dcast(gene~trajectory, value.var = "r", fill = 0) %>%
#   matrix.please

#####################################
## Identify upregulated activators ##
#####################################

cor_upregulated_activators.dt <- cor_dt.filt %>%
  .[class=="Activator" & rna_sign=="Up"] %>%
  .[,foo:=sum(sig==TRUE),by=c("gene")] %>% .[foo>=1] %>% .[,foo:=NULL] %>%
  .[,foo:=mean(r>=0),by=c("gene")] %>% .[foo==1] %>% .[,foo:=NULL]

to.plot <- cor_upregulated_activators.dt %>% 
  dcast(gene~trajectory, value.var = "r", fill = 0) %>%
  matrix.please

to.plot[to.plot<0.05] <- 0

pheatmap(to.plot,
# pheatmap(to.plot[sample(1:nrow(to.plot), size = 50),], 
         cluster_cols = F, 
         color = colorRampPalette(c("#F2F2F2", "#CD0000"))(100),
         fontsize_row = 10,
         fontsize_col = 11,
         angle_col = 0,
         treeheight_row = 0,
         filename = sprintf("%s/upregulated_activators_cor%s.pdf",io$outdir,opts$min.cor),
         width = 5, height = 11.5
)


#######################################
## Identify downregulated activators ##
#######################################

cor_downregulated_activators.dt <- cor_dt.filt %>%
  .[class=="Activator" & rna_sign=="Down"] %>%
  .[,foo:=sum(sig==TRUE),by=c("gene")] %>% .[foo>=1] %>% .[,foo:=NULL] %>%
  .[,foo:=mean(r>=0),by=c("gene")] %>% .[foo==1] %>% .[,foo:=NULL]

to.plot <- cor_downregulated_activators.dt %>% 
  dcast(gene~trajectory, value.var = "r", fill = 0) %>%
  matrix.please

to.plot[to.plot<0.05] <- 0

pheatmap(to.plot,
# pheatmap(to.plot[sample(1:nrow(to.plot), size = 50),], 
         cluster_cols = F, 
         color = colorRampPalette(c("#F2F2F2", "#CD0000"))(100),
         fontsize_row = 10,
         fontsize_col = 11,
         angle_col = 0,
         treeheight_row = 0,
         filename = sprintf("%s/downregulated_activators_cor%s.pdf",io$outdir,opts$min.cor),
         width = 5, height = 11.5
)


#####################################
## Identify upregulated repressors ##
#####################################

cor_upregulated_repressors.dt <- cor_dt.filt %>%
  .[class=="Repressor" & rna_sign=="Up"] %>%
  .[,foo:=sum(sig==TRUE),by=c("gene")] %>% .[foo>=1] %>% .[,foo:=NULL] %>%
  .[,foo:=mean(r<=0),by=c("gene")] %>% .[foo==1] %>% .[,foo:=NULL]

to.plot <- cor_upregulated_repressors.dt %>% 
  dcast(gene~trajectory, value.var = "r", fill = 0) %>%
  matrix.please

to.plot[abs(to.plot)<0.05] <- 0

pheatmap(to.plot, 
         cluster_cols = F, 
         color = colorRampPalette(c("#104E8B", "#F2F2F2"))(100),
         fontsize_row = 10,
         fontsize_col = 11,
         angle_col = 0,
         treeheight_row = 0,
         filename = sprintf("%s/upregulated_repressors_cor%s.pdf",io$outdir,opts$min.cor),
         width = 5, height = 11.5
)

#####################################
## Identify downregulated repressors ##
#####################################

cor_downregulated_repressors.dt <- cor_dt.filt %>%
  .[class=="Repressor" & rna_sign=="Down"] %>%
  .[,foo:=sum(sig==TRUE),by=c("gene")] %>% .[foo>=1] %>% .[,foo:=NULL] %>%
  .[,foo:=mean(r<=0),by=c("gene")] %>% .[foo==1] %>% .[,foo:=NULL]

to.plot <- cor_downregulated_repressors.dt %>% 
  dcast(gene~trajectory, value.var = "r", fill = 0) %>%
  matrix.please

to.plot[abs(to.plot)<0.05] <- 0

pheatmap(to.plot, 
         cluster_cols = F, 
         color = colorRampPalette(c("#104E8B", "#F2F2F2"))(100),
         fontsize_row = 10,
         fontsize_col = 11,
         angle_col = 0,
         treeheight_row = 0,
         filename = sprintf("%s/downregulated_repressors_cor%s.pdf",io$outdir,opts$min.cor),
         width = 5, height = 11.5
)

####################
## Identify mixed ##
####################


mixed.TFs <- cor_dt.filt %>% .[sig==T] %>% .[,.(N=length(unique(class))),by="gene"] %>% .[N>1,gene]
cor_dt.mixed <- cor_dt.filt %>% .[gene%in%mixed.TFs]

to.plot <- cor_dt.mixed %>% 
  dcast(gene~trajectory, value.var = "r", fill = 0) %>%
  matrix.please

# to.plot[abs(to.plot)<0.25] <- 0

pheatmap(to.plot, 
         cluster_cols = F, 
         color = colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100),
         fontsize_row = 10,
         fontsize_col = 11,
         angle_col = 0,
         treeheight_row = 0,
         filename = sprintf("%s/chromatin_mixed_cor%s.pdf",io$outdir,opts$min.cor),
         width = 5, height = 11.5
)

###########################################################
## Scatterplots of RNA difference vs chromVAR difference ##
###########################################################

trajectories.to.plot <- c("endoderm","mesoderm","blood")
to.plot <- cor_dt %>% 
  .[trajectory%in%trajectories.to.plot]

to.plot[,dot_size:=0.5*as.numeric(cut(abs(r), breaks=seq(0,1,by=0.05)))]
to.plot.text <- to.plot[abs(r)>0.75]

p <- ggplot(to.plot, aes_string(x="rna_diff", y="chromvar_diff")) +
  geom_point(aes_string(size="dot_size", fill="r", alpha="abs_r"), shape=21) +
  facet_wrap(~trajectory, nrow=3, scales="fixed") +
  scale_size_continuous(range = c(0.25,3.5)) +
  geom_vline(xintercept=0, linetype="dashed", size=0.5) +
  geom_hline(yintercept=0, linetype="dashed", size=0.5) +
  scale_fill_gradient2(low = "yellow", mid="gray90", high = "purple") +
  ggrepel::geom_text_repel(aes(label=gene), size=3.5, max.overlaps=Inf, data=to.plot.text) +
  # coord_cartesian(ylim=c(0,750)) +
  labs(x="Differential RNA expression", y="Differential motif accessibility") +
  guides(size=F, alpha=F) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size=rel(0.5)),
    axis.title = element_text(size=rel(0.85))
  )

pdf(sprintf("%s/diffrna_vs_diffacc_scatterplot.pdf",io$outdir), width=4.5, height=11)
print(p)
dev.off()

##########################################
## Plot number of events per trajectory ##
##########################################

to.plot <- cor_dt %>%
  .[,.(N=sum(sig)),by=c("trajectory","rna_sign","chromvar_sign")]

p <- ggplot(to.plot, aes_string(x="rna_sign",y="chromvar_sign")) +
  geom_tile(aes(fill=N), color="black") +
  scale_fill_gradient(low = "gray80", high = "purple") +
  facet_wrap(~trajectory, nrow=1) +
  labs(x="RNA expression directionality", y="Motif accessibility directionality") +
  theme_classic() +
  theme(
    axis.title = element_text(color="black"),
    axis.text = element_text(color="black")
  )


pdf(sprintf("%s/diffrna_vs_diffacc_number_associations_tile.pdf",io$outdir), width=8, height=3)
print(p)
dev.off()


to.plot <- cor_dt.filt %>%
  .[,.(N=sum(sig)),by=c("trajectory","rna_sign","class")]

ggplot(to.plot, aes_string(x="rna_sign", y="N", fill="class")) +
  geom_bar(stat="identity", color="black") +
  # scale_fill_gradient(low = "gray80", high = "purple") +
  facet_wrap(~trajectory, nrow=1) +
  # labs(x="RNA expression directionality", y="Motif accessibility directionality") +
  theme_classic() +
  theme(
    axis.title = element_text(color="black"),
    axis.text = element_text(color="black")
  )



#############################################################################
## Compare trajectory-specific correlation vs global correlation estimates ##
#############################################################################

# Load global correlation estimates
io$rna_vs_chromar_per_gene <- paste0(io$basedir,"/results/rna_atac/rna_vs_chromvar/pseudobulk/per_gene/cor_rna_vs_chromvar_correlated_peaks_pseudobulk.txt.gz")

rna_vs_chromvar_global.dt <- fread(io$rna_vs_chromar_per_gene) 

to.plot <- merge(rna_vs_chromvar_global.dt[,c("gene","r")], cor_dt.filt[,c("gene","r","trajectory")], by=c("gene"), suffixes = c(".global", ".local"))

to.plot[,dot_size:=as.numeric(cut(abs(r.local), breaks=seq(0,1,by=0.05)))]

# to.plot[,residual:=lm(formula=r.global~r.local)[["residuals"]], by=c("trajectory")] %>%
#   .[,abs_residual:=abs(residual)] %>%
#   .[,scale_residual:=minmax.normalisation(abs_residual)]

to.plot %>%
  .[,diff_cor:=r.local-r.global] %>% 
  .[,abs_diff_cor:=abs(diff_cor)] %>% 
  .[,scale_diff:=minmax.normalisation(diff_cor)]

p <- ggscatter(to.plot, x="r.local", y="r.global", fill="diff_cor", alpha="scale_diff", size="abs_diff_cor", shape=21) +
          # add="reg.line", add.params = list(color="black", fill="lightgray"), conf.int=TRUE) +
  facet_wrap(~trajectory, scales="fixed", nrow=1) +
  scale_size_continuous(range = c(0.25,3.5)) +
  geom_vline(xintercept=0, linetype="dashed", size=0.5) +
  geom_hline(yintercept=0, linetype="dashed", size=0.5) +
  # geom_abline(slope=1, intercept=0, linetype="dashed") +
  scale_fill_gradient2(low = "yellow", mid="gray90", high = "purple") +
  # ggrepel::geom_text_repel(aes(label=gene), size=3.5, max.overlaps=Inf, data=to.plot[abs_residual>=0.70]) +
  ggrepel::geom_text_repel(aes(label=gene), size=3.5, max.overlaps=Inf, data=to.plot[diff_cor>=0.6]) +
  ggrepel::geom_text_repel(aes(label=gene), size=3.5, max.overlaps=Inf, data=to.plot[diff_cor<=(-0.9)]) +
  # ggrepel::geom_text_repel(aes(label=gene), size=3.5, max.overlaps=Inf, data=to.plot[diff_cor<(-0.5)] +
  # ggrepel::geom_text_repel(aes(label=gene), size=3.5, max.overlaps=Inf, data=to.plot[diff_cor<(-0.8)]) +
  coord_cartesian(ylim=c(-1,1)) +
  labs(x="Trajectory-specific correlation", y="Global correlation") +
  guides(size=F, alpha=F) +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text = element_text(size=rel(1), color="black"),
    axis.title = element_text(size=rel(1.25), color="black")
  )

pdf(sprintf("%s/global_vs_local_correlation_scatterplot_residuals.pdf",io$outdir), width=18, height=5)
print(p)
dev.off()
