source(here::here("settings.R"))
source(here::here("utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--sce',             type="character",                               help='SingleCellExperiment file')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--outdir',          type="character",                               help='Output file')

args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
# args$sce <- paste0(io$basedir,"/processed/rna/SingleCellExperiment.rds")
# args$metadata <- file.path(io$basedir,"results/rna/qc/sample_metadata_after_qc.txt.gz")# io$metadata
# args$outdir <- paste0(io$basedir,"/results/rna/celltype_assignment/E5.5")
## END TEST ##

# I/O
dir.create(args$outdir, showWarnings=F, recursive=T)

# Options
opts$marker_genes <- list(
  "Visceral_endoderm" = c("Apoa1","Apob","Apoa4","Apom","Muc13"), 
  "Epiblast" = c("Dnmt3b","Tuba1a","Mkrn1","Pim2","Tdgf1")
)

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_rnaQC==TRUE & stage=="E5.5"]

#########################
## Load RNA expression ##
#########################

# Load RNA expression data as SingleCellExperiment object
sce <- load_SingleCellExperiment(args$sce, cells=sample_metadata$cell, normalise = TRUE)

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

#######################
## Feature selection ##
#######################

sce_filt <- sce[unlist(opts$marker_genes),]

#########
## PCA ##
#########

sce_filt <- runPCA(sce_filt, ncomponents = 2) 

####################################
## Assign cell types based on PC1 ##
####################################

pc1 <- reducedDim(sce_filt,"PCA")[,1]

clusters <- c("cluster1","cluster2")[as.numeric(pc1>0)+1]

# this is a bit ad hoc but will do the work
cluster1_epi_expr <- mean(logcounts(sce_filt[opts$marker_genes[["Epiblast"]],clusters=="cluster1"]))
cluster2_epi_expr <- mean(logcounts(sce_filt[opts$marker_genes[["Epiblast"]],clusters=="cluster2"]))

if (mean(cluster1_epi_expr)>mean(cluster2_epi_expr)) {
  cluster2celltype <- c("cluster1"="Epiblast", "cluster2"="Visceral_endoderm")  
} else {
  cluster2celltype <- c("cluster2"="Epiblast", "cluster1"="Visceral_endoderm")  
}

sce_filt$celltype <- stringr::str_replace_all(clusters,cluster2celltype)

###############################
## Plot celltype assignments ##
###############################

to.plot <- data.table(
  cell = colnames(sce_filt),
  celltype = sce_filt$celltype,
  V1 = reducedDim(sce_filt,"PCA")[,1],
  V2 = reducedDim(sce_filt,"PCA")[,2]
)

p <- ggplot(to.plot, aes(x=V1, y=V2, fill=celltype)) +
  geom_point(size=3, shape=21, stroke=0.15) +
  labs(x="PC1", y="PC2") +
  theme_classic() +
  theme(
    legend.position = "top",
    legend.title = element_blank()
  )

pdf(file.path(args$outdir,"pca_celltype.pdf"), width=7, height=5)
print(p)
dev.off()

##########################
## Plot gene expression ##
##########################

# i <- "Epiblast"; j <- "Lrpap1"
for (i in names(opts$marker_genes)) {
  for (j in opts$marker_genes[[i]]) {
    
    to.plot <- data.table(
      gene = j,
      cell = colnames(sce_filt),
      expr = logcounts(sce_filt[j,])[1,],
      V1 = reducedDim(sce_filt,"PCA")[,1],
      V2 = reducedDim(sce_filt,"PCA")[,2]
    )
      
    p <- ggplot(to.plot, aes(x=V1, y=V2, fill=expr)) +
      geom_point(size=3, shape=21, stroke=0.15) +
      scale_fill_gradient2(low = "gray50", mid="gray90", high = "red") +
      labs(x="PC1", y="PC2", title=sprintf("%s (%s marker)",j,i)) +
      theme_classic() +
      theme(
        plot.title = element_text(hjust=0.5),
        legend.position = "none"
      )
    
    outfile <- file.path(args$outdir,sprintf("pca_%s_%s.pdf",i,j))
    pdf(outfile, width=7, height=5)
    print(p)
    dev.off()
  }
}


##########
## Save ##
##########

to_save <- data.table(
  id_rna = colnames(sce_filt),
  celltype = sce_filt$celltype
)

fwrite(to_save, file.path(args$outdir,"celltype_assignment_E55.txt.gz"), sep="\t")
