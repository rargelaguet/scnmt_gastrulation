
library(data.table)
library(purrr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)


############################
## Define I/O and options ##
############################

## Define I/O 
io <- list()
io$in.metadata <- "/Users/ricard/data/gastrulation/sample_metadata.txt"
io$in.met.stats <- "/Users/ricard/data/gastrulation/met/stats/samples/sample_stats.txt"
io$in.acc.stats <- "/Users/ricard/data/gastrulation/acc/stats/samples/sample_stats.txt"
io$in.rna.stats <- "/Users/ricard/data/gastrulation/rna/stats/sample_stats.txt"
# io$annos_dir  <- "/Users/ricard/data/gastrulation/features/filt"
io$outdir <- "/Users/ricard/gastrulation/metaccrna/qc_stats/out"

## Define options
opts <- list()

# Define which stage and lineages to look at 
opts$stage_lineage10x <- c(
  "E4.5_Epiblast",
  "E4.5_Primitive_endoderm",
  "E5.5_Epiblast",
  "E5.5_Visceral_endoderm",
  "E6.5_Epiblast",
  "E6.5_Primitive_Streak",
  "E6.5_Mesoderm",
  "E6.5_Visceral_endoderm",
  "E7.5_Primitive_Streak",
  "E7.5_Endoderm",
  "E7.5_Mesoderm",
  "E7.5_Epiblast",
  "E7.5_Ectoderm"
)
opts$stage_lineage10x <- NULL

# Define genomic contexts (use NULL for no genomic context filtering)
opts$annos <- NULL

# Define which cells to use
tmp <- fread(io$in.metadata)
if (!is.null(opts$stage_lineage10x) ) {
  tmp <- tmp %>% 
    .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")] %>% 
    .[stage_lineage%in%opts$stage_lineage]
}
opts$rna.cells <- tmp[pass_rnaQC==T,id_rna]
opts$met.cells <- tmp[pass_metQC==T,id_met]
opts$acc.cells <- tmp[pass_accQC==T,id_acc]

###############
## Load data ##
###############

# Load sample metadata
metadata <- fread(io$in.metadata) %>% 
  # .[,c("sample","id_rna","id_met","id_acc","stage","plate")] %>%
  .[,stage_plate:=paste(stage,plate,sep="_")] %>%
  .[id_rna%in%opts$rna.cells | id_met%in%opts$met.cells | id_acc%in%opts$acc.cells  ]

# Load stats
rna.stats <- fread(io$in.rna.stats) %>% 
  merge(metadata[pass_rnaQC==T,c("sample","id_rna")], by="id_rna") %>% .[,id_rna:=NULL] %>%
  .[,total_counts:=log10(total_counts)]

met.stats <- fread(io$in.met.stats) %>% 
  merge(metadata[pass_metQC==T,c("sample","id_met")], by="id_met") %>% .[,id_met:=NULL] %>%
  .[,coverage:=log10(coverage)]

acc.stats <- fread(io$in.acc.stats) %>% 
  merge(metadata[pass_accQC==T,c("sample","id_acc")], by="id_acc") %>% .[,id_acc:=NULL] %>%
  .[,coverage:=log10(coverage)]

# Merge stats
metrna.stats <- merge(rna.stats, met.stats, by="sample") %>% merge(metadata,by="sample")
accrna.stats <- merge(rna.stats, acc.stats, by="sample") %>% merge(metadata,by="sample")
metacc.stats <- merge(met.stats, acc.stats, by="sample") %>% merge(metadata,by="sample")

################
## Histograms ##
################

p.rna <- gghistogram(rna.stats, x="num_genes", y="..count..", fill="#3CB54E") +
  labs(x="Number of expressed genes", y="Count")

pdf(file=sprintf("%s/histogram_coverage_rna.pdf",io$outdir), width=7, height=5)
print(p.rna)
dev.off()

p.met <- gghistogram(met.stats, x="coverage", y="..count..", fill="#F37A71") +
  labs(x="Number of observed CpGs (log)", y="Count")

pdf(file=sprintf("%s/histogram_coverage_met.pdf",io$outdir), width=7, height=5)
print(p.met)
dev.off()

p.acc <- gghistogram(acc.stats, x="coverage", y="..count..", fill="#6691CB") +
  labs(x="Number of observed GpCs (log)", y="Count")

pdf(file=sprintf("%s/histogram_coverage_acc.pdf",io$outdir), width=7, height=5)
print(p.acc)
dev.off()

###############
## Box plots ##
###############

theme_boxplot <- function() {
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
}

# Box plot of RNA complexity per plate
to.plot <- rna.stats %>% merge(metadata, by="sample")
p.rna2 <- ggboxplot(to.plot, x="plate", y="total_counts", fill="#3CB54E", outlier.shape=NA) +
  facet_wrap(~stage, nrow=1, scales="free_x") +
  labs(x="", y="Library size (log)") +
  theme_boxplot()
# pdf(file=sprintf("%s/boxplot_coverage_rna.pdf",io$outdir), width=7, height=5)
# p.rna2
# dev.off()

# Box plot of DNA methylation complexity per plate
to.plot <- met.stats %>% merge(metadata, by="sample") %>% .[,N:=.N, by=c("stage","plate")] %>% .[N>10]
p.met2 <- ggboxplot(to.plot, x="plate", y="coverage", fill="#F37A71", outlier.shape=NA) +
  facet_wrap(~stage, nrow=1, scales="free_x") +
  labs(x="", y="Library size (log)") +
  theme_boxplot()
# pdf(file=sprintf("%s/boxplot_coverage_met.pdf",io$outdir), width=7, height=5)
# p.met2
# dev.off()

# Box plot of chromatin accessibility complexity per plate
to.plot <- acc.stats %>% merge(metadata, by="sample") %>% .[,N:=.N, by=c("stage","plate")] %>% .[N>10]
p.acc2 <- ggboxplot(to.plot, x="plate", y="coverage", fill="#00BFC4", outlier.shape=NA) +
  facet_wrap(~stage, nrow=1, scales="free_x") +
  labs(x="", y="Library size (log)") +
  theme_boxplot()
# pdf(file=sprintf("%s/boxplot_coverage_acc.pdf",io$outdir), width=7, height=5)
# p.acc2
# dev.off()

###################
## Scatter plots ##
###################

# Plot RNA complexity versus DNA methylation coverage
p.metrna <- ggscatter(metrna.stats, x="coverage", y="num_genes", add="reg.line", add.params = list(color="blue", fill="lightgray"), conf.int=TRUE) +
  stat_cor(method = "pearson") +
  labs(y="Number of expressed genes", x="Number of observed CpGs (log)")

pdf(file=sprintf("%s/scatterplot_coverage_metrna.pdf",io$outdir), width=7, height=5, useDingbats=F)
print(p.metrna)
dev.off()

# Plot RNA complexity versus chromatin accessibility coverage
p.accrna <- ggscatter(accrna.stats, x="coverage", y="num_genes", add="reg.line", add.params = list(color="blue", fill="lightgray"), conf.int=TRUE) +
  stat_cor(method = "pearson") +
  labs(y="Number of expressed genes", x="Number of observed GpCs (log)")

pdf(file=sprintf("%s/scatterplot_coverage_accrna.pdf",io$outdir), width=7, height=5, useDingbats=F)
print(p.accrna)
dev.off()

# Plot DNA methylation coverage versus chromatin accessibility coverage
p.metacc <- ggscatter(metacc.stats, x="coverage.x", y="coverage.y", add="reg.line", add.params = list(color="blue", fill="lightgray"), conf.int=TRUE) +
  stat_cor(method = "pearson") +
  labs(x="Number of observed CpGs (log)", y="Number of observed GpCs (log)")

pdf(file=sprintf("%s/scatterplot_coverage_metacc.pdf",io$outdir), width=7, height=5, useDingbats=F)
print(p.metacc)
dev.off()

#######################################################
## RNA library size versus number of expressed genes ##
#######################################################

to.plot <- rna.stats %>% merge(metadata, by="sample")
p.rna3 <- ggscatter(to.plot, x="total_counts", y="num_genes", color="#3CB54E", add="reg.line", add.params = list(color="black", fill="lightgray"), conf.int=TRUE) +
  stat_cor(method = "pearson", size=5) +
  facet_wrap(~stage, nrow=1, scales="free_x") +
  labs(x="Library size (log)", y="Number of expressed genes")
print(p.rna3)

pdf(file=sprintf("%s/scatterplot_coverage_rna_per_stage.pdf",io$outdir), width=12, height=5, useDingbats=F)
print(p.rna3)
dev.off()

model <- lm(num_genes~total_counts, data=to.plot[stage=="E4.5"])
predict.lm(model, data.frame(total_counts=median(to.plot[stage=="E4.5",total_counts])))
predict.lm(model, data.frame(total_counts=median(to.plot[stage=="E7.5",total_counts])))




##################
# Combine plots ##
##################

pdf(file=sprintf("%s/scatterplot_coverage_all.pdf",io$outdir), width=14, height=10, useDingbats=F)
cowplot::plot_grid(plotlist=list(p.metrna, p.accrna, p.metacc, p.rna, p.met, p.acc), nrow=2)
dev.off()

pdf(file=sprintf("%s/boxplot_coverage_all.pdf",io$outdir), width=14, height=5, useDingbats=F)
cowplot::plot_grid(plotlist=list(p.rna2, p.met2, p.acc2), nrow=1)
dev.off()