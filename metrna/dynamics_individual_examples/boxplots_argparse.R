suppressPackageStartupMessages(library(argparse))

source("/Users/ricard/scnmt_gastrulation/settings.R")
source("/Users/ricard/scnmt_gastrulation/metrna/dynamics_individual_examples/load_data.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--gene',    type="character",   help='Feature name for RNA expression (in ENSEMBL ID)')
p$add_argument('--met.id',    type="character",   help='Feature name for DNA methylation')
p$add_argument('--met.anno',    type="character",   help='Genomic context for DNA methylation')
p$add_argument('--stage_lineage',    type="character",   nargs='+',   help='stage and lineages to plot')
p$add_argument('--outdir',  type="character",              help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$gene <- "Rhox4e"
# args$met.id <- "ENSMUSG00000071770"
# args$met.anno <- "prom_2000_2000"
# args$stage_lineage <- c("E4.5_Epiblast", "E5.5_Epiblast", "E6.5_Epiblast", "E6.5_Primitive_Streak", "E6.5_Mesoderm", "E7.5_Epiblast", "E7.5_Ectoderm", "E7.5_Endoderm", "E7.5_Primitive_Streak", "E7.5_Mesoderm")
# args$outdir <- "/Users/ricard/data/scnmt_gastrulation/metrna/plot_individual_examples"
## END TEST ##

####################
## Define options ##
####################

# Define stages to plot

# Define colors for the omics
opts$color <- c(
  # "RNA expression"="#3CB54E",
  # "DNA methylation"="#F37A71"
  "RNA expression"="gray50",
  "DNA methylation"="gray50"
)

# Define minimum coverage
opts$min.cpg <- 3

# Define cells to use
opts$met_cells <- sample_metadata %>% .[pass_metQC==T & stage_lineage%in%args$stage_lineage,id_met]
opts$rna_cells <- sample_metadata %>% .[pass_rnaQC==T & stage_lineage%in%args$stage_lineage,id_rna]

###############
## Load data ##
###############

# Update sample metadata
sample_metadata <- sample_metadata %>%
  .[,c("sample","id_met","id_rna","stage","stage_lineage","lineage10x_2")] %>%
  .[!is.na(id_met) & !is.na(id_rna)] %>%  # optional
  .[id_met%in%opts$met_cells | id_rna %in% opts$rna_cells]

# Load the three omics
dt <- load_data(io, args$gene, args$met.id, args$met.anno, opts$min.cpg)

# Merge data with sample metadata
dt$met <- merge(dt$met, sample_metadata, by="id_met")
dt$rna <- merge(dt$rna, sample_metadata, by="id_rna")

# bind in a single data table
to.plot <- do.call("rbind",list(
  dt$rna[,c("sample","id","stage","stage_lineage","lineage10x_2","value")] %>% .[,assay:="RNA expression"],
  dt$met[,c("sample","id","stage","stage_lineage","lineage10x_2","value")] %>% .[,assay:="DNA methylation"]
)) %>% .[,stage_lineage:=gsub("_"," ",stage_lineage)] %>% 
  .[,assay:=factor(assay,levels=c("RNA expression","DNA methylation"))]
  
##############
## Boxplots ##
##############


p <- ggplot(to.plot, aes(x=stage, y=value)) +
  facet_wrap(~assay, ncol=1, scales="free_y") +
  geom_jitter(aes(fill=assay), size=0.5, stroke=0.2, shape=21, color="black", width=0.1) +
  geom_violin(aes(fill=assay), alpha=0.6, size=0.25) +
  geom_boxplot(aes(fill=assay), alpha=0.6, outlier.shape=NA, width=0.3, size=0.25) +
  scale_fill_manual(values=opts$color) +
  scale_x_discrete(expand = c(0,-1)) +
  labs(x="", y="", title="") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size=rel(1.0), color="black"),
    axis.text.y = element_text(colour="black",size=rel(0.7)),
    axis.line = element_line(colour="black", size=rel(0.7)),
    axis.ticks.y = element_line(size=rel(0.7)),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.position="none"
  )

##########
## Save ##
##########

outfile <- sprintf("%s/boxplot_rna_%s_met_%s.pdf",args$outdir,args$gene,args$met.id)

pdf(outfile, useDingbats=F, width=2, height=3.5)
print(p)
dev.off()
