here::i_am("rna/differential/differential_expr.R")

suppressMessages(library(edgeR))

source(here::here("settings.R"))
source(here::here("utils.R"))
source(here::here("rna/differential/utils.R"))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--name_groupA', type="character",  help='name for group A')
p$add_argument('--name_groupB', type="character",  help='name for group B')
p$add_argument('--groupA', type="character",  nargs='+',  help='group A')
p$add_argument('--groupB', type="character",  nargs='+',  help='group B')
p$add_argument('--group_label',    type="character",    help='Group label')
p$add_argument('--outfile',        type="character",              help='Output file')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
args$groupA <- c("E4.5_Epiblast","E5.5_Epiblast")
args$groupB <- c("E6.5_Epiblast","E7.5_Epiblast")
args$name_groupA <- "E4.5_E5.5_Epiblast"
args$name_groupB <- "E6.5_E7.5_Epiblast"
args$group_label <- "stage_lineage3"
## END TEST ##

#####################
## Define settings ##
#####################

io$metadata <- file.path(io$basedir,"results/rna/celltype_assignment/sample_metadata_after_celltype_rename.txt.gz")

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
  setnames(args$group_label,"group")  %>%
  .[group%in%c(args$groupA,args$groupB)]
table(sample_metadata$group)

sample_metadata <- sample_metadata %>%
  .[,group:=ifelse(group%in%args$groupA,args$name_groupA,args$name_groupB)] %>%
  .[,group:=factor(group,levels=c(args$name_groupA,args$name_groupB))]
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
  rowMeans(logcounts(sce[,sce$group==args$name_groupA])>0) %>% round(2),
  rowMeans(logcounts(sce[,sce$group==args$name_groupB])>0) %>% round(2)
# ) %>% setnames(c("gene",sprintf("detection_rate_%s",args$name_groupA),sprintf("detection_rate_%s",args$name_groupB)))
) %>% setnames(c("gene","detection_rate_groupA","detection_rate_groupB"))

# Filter genes
sce <- sce[rownames(sce)%in%gene_metadata$gene,]

################################################
## Differential expression testing with edgeR ##
################################################

out <- doDiffExpr(sce, groups = c(args$name_groupA,args$name_groupB), min_detection_rate_per_group = opts$min_detection_rate_per_group) %>%
  # Add sample statistics
  .[,c("groupA_N","groupB_N"):=list(table(sample_metadata$group)[1],table(sample_metadata$group)[2])]%>% 
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
