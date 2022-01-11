here::i_am("sex/sex_assignment.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--samples',   type="character",   nargs='+',  help='samples')
p$add_argument('--sce',       type="character",               help='SingleCellExperiment file')
p$add_argument('--metadata',  type="character",               help='metadata file')
p$add_argument('--chrY_ratio_threshold',  type="double", default=1e-3,              help='ChrY/chr10 counts ratio threshold')
p$add_argument('--outdir',          type="character",               help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
# args$samples <- opts$samples
# args$sce <- paste0(io$basedir,"/processed_all/SingleCellExperiment.rds")
# args$metadata <- file.path(io$basedir,"results_all/mapping/sample_metadata_after_mapping.txt.gz")
# args$outdir <- paste0(io$basedir,"/results_all/sex_assignment")
# args$chrY_ratio_threshold <- 0.1
## END TEST ##

# I/O
dir.create(args$outdir, showWarnings=F)

###############
## Load data ##
###############

# Load sample metadata
sample_metadata <- fread(args$metadata) %>% .[pass_rnaQC==T & sample%in%args$samples]
print(table(sample_metadata$class))

# Load SingleCellExperiment
sce <- load_SingleCellExperiment(args$sce, cells = sample_metadata$cell, normalise = T)
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

# Load gene metadata
gene_metadata <- fread(io$gene_metadata) %>% 
  .[!grepl("^[Gm|Rik]",symbol)] %>%
  .[symbol%in%rownames(sce)]

################
## Parse data ##
################

# Group genes by chromosome
genes.chrY <- gene_metadata[chr=="chrY",symbol]
genes.chrX <- gene_metadata[chr=="chrX",symbol]
genes.chr10 <- gene_metadata[chr=="chr10",symbol]

# Manual filtering
# For some reason Erdr1 is predicted as Ychr, but the last version of ENSEMBL is in the Xchr
# genes.chrY <- genes.chrY[!genes.chrY=="ENSMUSG00000096768"]

# Create data.table with the expression values of selected genes
dt <- args$samples %>% map(function(i) {
  sce.filt <- sce[,sce$sample==i] %>% .[c(genes.chrY,genes.chrX,genes.chr10),]
  data.table(
    sample = i,
    symbol = rownames(sce.filt),
    expr = logcounts(sce.filt) %>% Matrix::rowMeans()
  )
}) %>% rbindlist %>% merge(gene_metadata[,c("chr","ens_id","symbol")], by="symbol")


###################################################
## Barplots of chrXY/chr1 count ratio per embryo ##
###################################################

# Agregate counts over all genes
sex_assignment.dt <- dt %>% 
  .[,.(expr=mean(expr)),by=c("sample","chr")] %>%
  dcast(sample~chr, value.var="expr") %>%
  .[,ratioY:=round(chrY/chr10,3)] %>% .[,ratioX:=round(chrX/chr10,3)] %>%
  .[,sex:=c("female","male")[as.numeric(ratioY>=args$chrY_ratio_threshold)+1]]

p <- ggbarplot(sex_assignment.dt, x="sample", y="ratioY", fill="sex", sort.val = "asc", palette="Dark2") +
  labs(x="", y="chrY/chr1 expr ratio") +
  guides(x = guide_axis(angle = 90)) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(colour="black",size=rel(0.6)),
    axis.text.y = element_text(colour="black",size=rel(0.8))
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank()
  )

pdf(file.path(args$outdir,"sex_ychr_expr_aggregated.pdf"))
print(p)
dev.off()

#########################################
## Plot expression of individual genes ##
#########################################

to.plot <- dt[chr=="chrY"] %>% 
  merge(sex_assignment.dt[,c("sex","sample")], by="sample") %>%
 .[,foo:=mean(expr),by="symbol"] %>% .[foo>0] %>% .[,foo:=NULL] 

p <- ggbarplot(to.plot, x="symbol", y="expr", facet="sample", fill="sex") +
# p <- ggbarplot(to.plot, x="ens_id", y="counts", facet="sample", fill="gray70") +
  labs(x="", y="Expression levels") +
  guides(x = guide_axis(angle = 90)) +
  theme(
    axis.text.x = element_text(colour="black",size=rel(0.75)),
    axis.ticks.x = element_line(size=rel(0.5)),
    axis.text.y = element_text(colour="black",size=rel(0.5)),
    strip.background = element_blank(),
    strip.text = element_text(size=rel(0.5), color="black")
  )

pdf(file.path(args$outdir,"sex_ychr_expr_per_gene.pdf"), width=14, height=10)
print(p)
dev.off()

#############################
## Plot expression of Xist ##
#############################

to.plot <- dt[symbol=="Xist"] %>% 
  merge(sex_assignment.dt[,c("sex","sample")], by="sample")

p <- ggbarplot(to.plot, x="sample", y="expr", fill="sex") +
  labs(x="", y="Xist expression") +
  guides(x = guide_axis(angle = 90)) +
  theme(
    axis.text.x = element_text(colour="black",size=rel(0.5)),
    axis.text.y = element_text(colour="black",size=rel(0.8))
  )

pdf(file.path(args$outdir,"Xist_expr.pdf"))
print(p)
dev.off()

###########################################
## Scatterplot of chrY/chr1 vs Xist expr ##
###########################################

to.plot <- dt[symbol=="Xist"] %>% 
  setnames("expr","Xist_expr") %>%
  merge(sex_assignment.dt, by="sample")

ggscatter(to.plot, x="ratioY", y="Xist_expr", shape=21, fill="sex", size=2.5) +
  labs(x="chrX/chr1 ratio", y="Xist expression") +
  theme(
    axis.text = element_text(colour="black",size=rel(0.8)),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )

############################
## Update sample metadata ##
############################

sample_metadata_after_sex_assignment <- fread(args$metadata) %>% 
  merge(sex_assignment.dt[,c("sample","sex")], by="sample", all.x = TRUE)

##########
## Save ##
##########

fwrite(sample_metadata_after_sex_assignment, file.path(args$outdir,"sample_metadata_after_sex_assignment.txt.gz"), sep="\t", quote=F, na="NA")
fwrite(sex_assignment.dt, file.path(args$outdir,"sex_assignment.txt.gz"))
