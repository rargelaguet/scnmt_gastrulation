source("/Users/argelagr/scnmt_gastrulation/settings.R")
source("/Users/argelagr/scnmt_gastrulation/utils.R")
# source("/Users/argelagr/scnmt_gastrulation/metaccrna/dynamics_individual_examples/load_data.R")

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--gene',    type="character",   help='Feature name for RNA expression (in ENSEMBL ID)')
p$add_argument('--met.id',    type="character",   help='Feature name for DNA methylation')
p$add_argument('--met.anno',    type="character",   help='Genomic context for DNA methylation')
p$add_argument('--acc.id',    type="character",   help='Feature name for chromatin accessibility')
p$add_argument('--acc.anno',    type="character",   help='Genomic context for chromatin accessibility')
p$add_argument('--stage_lineage',    type="character",   nargs='+',   help='stage and lineages to plot')
p$add_argument('--outdir',  type="character",              help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
args$gene <- "Rhox5"
args$met.id <- "ENSMUSG00000095180"
args$met.anno <- "prom_2000_2000"
args$acc.id <- "ENSMUSG00000095180"
args$acc.anno <- "prom_2000_2000"
# args$stage_lineage <- c("E4.5_Epiblast", "E5.5_Epiblast", "E6.5_Epiblast", "E6.5_Primitive_Streak", "E6.5_Mesoderm", "E7.5_Epiblast", "E7.5_Ectoderm", "E7.5_Endoderm", "E7.5_Primitive_Streak", "E7.5_Mesoderm")
args$stage_lineage <- c(
  # "E3.5_ICM",
  "E4.5_Epiblast",
  "E4.5_Primitive_endoderm",
  "E5.5_Epiblast",
  "E5.5_Visceral_endoderm",
  "E6.5_Epiblast",
  "E6.5_Primitive_Streak",
  "E6.5_Visceral_endoderm",
  "E7.5_Epiblast",
  "E7.5_Ectoderm",
  "E7.5_Primitive_Streak",
  "E7.5_Endoderm",
  "E7.5_Mesoderm"
)
args$outdir <- "/Users/argelagr/data/scnmt_gastrulation/metaccrna/plot_individual_examples"
## END TEST ##

####################
## Define options ##
####################

# Define stages to plot

# Define colors for the omics
opts$color <- c(
  "RNA expression"="#3CB54E",
  "Chromatin accessibility"="#00BFC4",
  "DNA methylation"="#F37A71"
)

# Define minimum coverage
opts$min.cpg <- 1
opts$min.gpc <- 10

######################
## Load sample data ##
######################

sample_metadata <- fread(io$metadata) %>%
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")] %>%
  .[stage_lineage%in%opts$stage_lineage] %>% .[,stage_lineage:=factor(stage_lineage, levels=opts$stage_lineage)] %>%
  .[,c("sample","id_met","id_rna","id_acc","stage","stage_lineage","lineage10x_2","pass_rnaQC","pass_metQC","pass_accQC")]

# Define cells to use
opts$met_cells <- sample_metadata %>% .[pass_metQC==T,id_met]
opts$rna_cells <- sample_metadata %>% .[pass_rnaQC==T,id_rna]
opts$acc_cells <- sample_metadata %>% .[pass_accQC==T,id_acc]

sample_metadata[,c("pass_rnaQC","pass_metQC","pass_accQC"):=NULL]

###############
## Load data ##
###############

# Load DNA methylation data
met_dt <- fread(sprintf("%s/%s.tsv.gz",io$met_data_parsed,args$met.anno)) %>%
  setnames(c("id_met","id","anno","Nmet","N","value")) %>%
  .[id%in%args$met.id] %>% .[N>=opts$min.cpg]

# Load DNA accessibility data
acc_dt <- fread(sprintf("%s/%s.tsv.gz",io$acc_data_parsed,args$acc.anno)) %>%
  setnames(c("id_acc","id","anno","Nmet","N","value")) %>%
  .[id%in%args$acc.id] %>% .[N>=opts$min.gpc]

# Load RNA data
sce <- load_SingleCellExperiment(io$rna.sce, normalise = TRUE, cells = opts$rna_cells)

# Rename genes
gene_metadata <- fread(io$gene.metadata) %>% .[ens_id%in%rownames(sce) & symbol!=""]
foo <- gene_metadata$symbol
names(foo) <- gene_metadata$ens_id
sce <- sce[rownames(sce) %in% names(foo),]
rownames(sce) <- foo[rownames(sce)]

# Extract data.table
rna_dt <- logcounts(sce[args$gene,]) %>% t %>% as.data.table(keep.rownames = "id_rna") %>% 
  melt(id.vars = "id_rna", value.name = "value", variable.name = "gene")

# Merge data with sample metadata
acc_dt <- merge(acc_dt, sample_metadata[,c("id_acc","sample")], by="id_acc") %>% .[,id_acc:=NULL] %>% .[,assay:="Chromatin accessibility"]
met_dt <- merge(met_dt, sample_metadata[,c("id_met","sample")], by="id_met") %>% .[,id_met:=NULL] %>% .[,assay:="DNA methylation"]
rna_dt <- merge(rna_dt, sample_metadata[,c("id_rna","sample")], by="id_rna") %>% .[,id_rna:=NULL] %>% .[,assay:="RNA expression"]

# bind in a single data table
to.plot <- do.call("rbind",list(rna_dt, met_dt, acc_dt)) %>% 
  merge(sample_metadata[,c("sample","stage","stage_lineage","lineage10x_2")], by="sample") %>% #.[,stage_lineage:=gsub("_"," ",stage_lineage)] %>% 
  .[,assay:=factor(assay,levels=c("RNA expression","DNA methylation","Chromatin accessibility"))]
  
##############
## Boxplots ##
##############

ggplot(to.plot, aes(stage_lineage, y=value)) +
  facet_wrap(~assay, ncol=1, scales="free_y") +
  geom_jitter(aes(color=assay), size=0.5) +
  geom_violin(aes(fill=assay), alpha=0.5, size=0.25) +
  geom_boxplot(aes(fill=assay), alpha=0.5, outlier.shape=NA, width=0.15, size=0.25) +
  scale_fill_manual(values=opts$color) +
  scale_color_manual(values=opts$color) +
  labs(x="", y="", title="") +
  guides(x = guide_axis(angle = 90)) +
  theme_classic() +
  theme(
    axis.title.y = element_text(colour="black", size=rel(1.1), vjust=1.5),
    axis.text.x = element_text(size=rel(1.2), color="black"),
    axis.text.y = element_text(colour="black",size=rel(1.0)),
    axis.line = element_line(colour="black", size=rel(0.7)),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.position="none"
  )

##########
## Save ##
##########

# outfile <- sprintf("%s/boxplot_rna%s_met%s_acc%s.pdf",args$outdir,args$gene,args$met.id,args$acc.id)
# 
# pdf(outfile, useDingbats=F, width=6, height=5)
# print(p)
# dev.off()
