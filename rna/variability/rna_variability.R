library(scater)
library(RColorBrewer)
library(scran)

#####################
## Define settings ##
#####################

source("/Users/ricard/scnmt_gastrulation/settings.R")
io$outdir <- paste0(io$basedir,"/rna/results/variability")

## Define options ##
# Define stage and lineage
opts$stage_lineage <- c(
  "E4.5_Epiblast",
  "E5.5_Epiblast",
  "E6.5_Epiblast",
  "E6.5_Primitive_Streak",
  # "E7.5_Epiblast",
  # "E7.5_Ectoderm",
  # "E7.5_Primitive_Streak",  # PLOTS HAVE BEEN DONE WITH THIS UNCOMMENTED
  "E7.5_Endoderm",
  "E7.5_Mesoderm"
)

# Update sample metadata
sample_metadata <- sample_metadata %>% 
  # .[pass_rnaQC==T & stage_lineage%in%opts$stage_lineage]
  .[pass_rnaQC==T & !is.na(id_met) & stage_lineage%in%opts$stage_lineage]
table(sample_metadata$stage)

###############
## Load data ##
###############

# SingleCellExperiment object
sce <- readRDS(io$rna)[,sample_metadata$id_rna]

####################################
## Calculate variability per gene ##
####################################

# Gaussian variance
var <- apply(logcounts(sce),1,var)
dt.var <- data.table(ens_id=names(var), var=var)

# Biological overdispersion using scran
decomp <- modelGeneVar(sce)
dt.biodisp <- data.table(ens_id=rownames(decomp), bio.disp=decomp$bio)

##########
## Plot ##
##########

to.plot <- data.table(
  gene = rowData(sce)$symbol,
  mean = apply(logcounts(sce),1,mean),
  sd = sqrt(var),
  overdispersion = decomp$bio,
  hvg = decomp$p.value < 0.01
) %>% .[is.na(hvg),hvg:=F] %>% setorder(-overdispersion)

to.plot <- to.plot[mean>0.1]

p1 <- ggplot(to.plot, aes(x=mean, y=sd)) +
  geom_point(aes(fill=hvg, alpha=hvg, size=hvg), shape=21, stroke=0.1) +
  stat_smooth() +
  scale_fill_manual(values=c("black","red")) +
  scale_alpha_manual(values=c(0.75,1.5)) +
  scale_size_manual(values=c(0.8,1.6)) +
  guides(fill = guide_legend(override.aes = list(size=2.5))) +
  xlab('Mean expression (RNA)') + ylab('Standard deviation') +
  coord_cartesian(xlim=c(0,10)) +
  ggrepel::geom_text_repel(data=head(to.plot[hvg==TRUE],n=15), aes(x=mean, y=sd, label=gene), size=3.5, color="black") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(color="black", size=rel(0.9)),
    axis.title = element_text(color="black", size=rel(1.1))
  )

p2 <- ggplot(to.plot, aes(x=mean, y=overdispersion)) +
  geom_point(aes(fill=hvg, alpha=hvg, size=hvg), shape=21, stroke=0.1) +
  # stat_smooth() +
  scale_fill_manual(values=c("gray80","red")) +
  scale_alpha_manual(values=c(0.7,0.9)) +
  scale_size_manual(values=c(0.5,2)) +
  guides(fill = guide_legend(override.aes = list(size=2.5))) +
  xlab('Mean expression (RNA)') + ylab('Residual overdispersion (RNA)') +
  coord_cartesian(xlim=c(0,10)) +
  ggrepel::geom_text_repel(data=head(to.plot[hvg==TRUE],n=15), aes(x=mean, y=overdispersion, label=gene), size=5, color="black") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(color="black", size=rel(0.9)),
    axis.title = element_text(color="black", size=rel(1.1))
  )

# p <- cowplot::plot_grid(plotlist=list(p1,p2))

pdf(paste0(io$outdir,"/rna_overdispersion.pdf"), width = 9, height = 5, useDingbats = FALSE)
print(p2)
dev.off()

##########
## Save ##
##########

dt <- merge(dt.var,dt.biodisp, by="ens_id")

# Save results
fwrite(dt, paste0(io$outdir,"/gene_variability.txt.gz"), sep="\t", quote=F)
