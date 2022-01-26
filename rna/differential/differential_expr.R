here::i_am("rna/differential/differential_expr.R")

suppressMessages(library(edgeR))

source(here::here("settings.R"))
source(here::here("utils.R"))
source(here::here("rna/differential/utils.R"))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--groupA', type="character",  nargs='+',  help='group A')
p$add_argument('--groupB', type="character",  nargs='+',  help='group B')
p$add_argument('--group_label',    type="character",    help='Group label')
p$add_argument('--outfile',        type="character",              help='Output file')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
args$groupA <- ("E3.5_ICM")
args$groupB <- ("E4.5_Epiblast")
args$group_label <- "stage_lineage3"
## END TEST ##

#####################
## Define settings ##
#####################

io$metadata <- file.path(io$basedir,"results/rna/celltype_assignment/sample_metadata_after_celltype_rename.txt.gz")

# Define groups
opts$groups <- c(args$groupA, args$groupB)

# Define FDR threshold
opts$threshold_fdr <- 0.01

# Define minimum logFC for significance
opts$min.logFC <- 1.0

# For a given gene, the minimum fraction of cells that must express it in at least one group
opts$min_detection_rate_per_group <- 0.40

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>% 
  .[pass_rnaQC==T] %>%
  .[,stage_lineage:=paste(stage,celltype,sep="_")] %>%
  .[,stage_lineage2:=paste(stage,celltype2,sep="_")] %>%
  .[,stage_lineage3:=paste(stage,celltype3,sep="_")]

stopifnot(args$group_label%in%colnames(sample_metadata))

sample_metadata <- sample_metadata %>%
  setnames(args$group_label,"group") %>%
  # .[,group:=eval(as.name(args$group_label))] %>%
  .[group%in%c(args$groupA,args$groupB)] %>%
  .[,group:=factor(group,levels=opts$groups)] %>% setorder(group) # Sort cells so that groupA comes before groupB

table(sample_metadata$group)

###############
## Load data ##
###############

# Load SingleCellExperiment object
sce <- load_SingleCellExperiment(
  file = io$rna.sce, 
  normalise = TRUE, 
  cells = sample_metadata$id_rna
)
sce$group <- sample_metadata$group

# Load gene metadata
gene_metadata <- fread(io$gene_metadata) %>%
  .[symbol%in%rownames(sce)] %>%
  .[,c("symbol","ens_id")] %>%
  setnames("symbol","gene")

################
## Parse data ##
################

# calculate detection rate per gene
cdr.dt <- data.table(
  rownames(sce),
  rowMeans(logcounts(sce[,sce$group==opts$groups[1]])>0) %>% round(2),
  rowMeans(logcounts(sce[,sce$group==opts$groups[2]])>0) %>% round(2)
) %>% setnames(c("gene",sprintf("detection_rate_%s",opts$groups[1]),sprintf("detection_rate_%s",opts$groups[2])))
# .[,cdr_diff:=abs(out[,(sprintf("detection_rate_%s",opts$groups[1])),with=F][[1]] - out[,(sprintf("detection_rate_%s",opts$groups[2])),with=F][[1]])] %>%

# Filter genes
sce <- sce[rownames(sce)%in%gene_metadata$gene,]

################################################
## Differential expression testing with edgeR ##
################################################

out <- doDiffExpr(sce, opts$groups, opts$min_detection_rate_per_group) %>%
  # Add sample statistics
  .[,c("groupA_N","groupB_N"):=list(table(sample_metadata$group)[1],table(sample_metadata$group)[2])]%>% 
  # setnames(c("groupA_N","groupB_N"),c(sprintf("N_%s",opts$groups[1]),sprintf("N_%s",opts$groups[2]))) %>%
  # Add gene statistics
  merge(cdr.dt, all.y=T, by="gene") %>%
  merge(gene_metadata, all.y=T, by="gene") %>%
  # Calculate statistical significance
  # .[, sig := (padj_fdr<=opts$threshold_fdr & abs(logFC)>=opts$min.logFC)] %>%
  # .[is.na(sig),sig:=FALSE] %>%
  setorder(padj_fdr, na.last=T)

# Parse columns
out[,c("p.value","padj_fdr","logFC","log_padj_fdr"):=list(signif(p.value,digits=3), signif(padj_fdr,digits=3), round(logFC,3),round(log_padj_fdr,3))]

##################
## Save results ##
##################

fwrite(out, args$outfile, sep="\t", na="NA", quote=F)
