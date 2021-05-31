library(data.table)
library(purrr)
library(SingleCellExperiment)
library(scater)
library(scran)

#######
## A ##
#######

a <- readRDS("/Users/ricard/data/gastrulation/rna/SingleCellExperiment.rds")
a_counts <- counts(a)

#######
## B ##
#######

b <- readRDS("/Users/ricard/data/scnmt_gastrulation_E3.5/rna/SingleCellExperiment.rds")

# Subset
b <- b[rownames(b)%in%rownames(a)]
b <- b[,!colnames(b)%in%colnames(a)]

b_counts <- counts(b)

#######################################################
## Rename genes and take intersect between data sets ##
#######################################################

# rownames(a_counts) <- stringr::str_split(rownames(a_counts), pattern="_") %>% map_chr(1)

genes_a <- rownames(a_counts)
genes_b <- rownames(b_counts)

genes <- intersect(genes_a,genes_b)

a_counts <- a_counts[genes,]
b_counts <- b_counts[genes,]

#######################
## Create merged SCE ##
#######################

counts <- cbind(a_counts,b_counts)
sce <- SingleCellExperiment(assays = list(counts = counts))

############################
## Calculate size factors ##
############################

clusts = as.numeric(quickCluster(sce, method = "igraph", min.size = 100, BPPARAM = mcparam))
min.clust = min(table(clusts))/2
new_sizes = c(floor(min.clust/3), floor(min.clust/2), floor(min.clust))
sce = computeSumFactors(sce, clusters = clusts, sizes = new_sizes, max.cluster.size = 3000)

to.plot <- data.frame(X = Matrix::colSums(counts(sce)), Y = sizeFactors(sce))
ggplot(to.plot, mapping = aes(x = X, y = Y)) +
  geom_point() +
  labs(x = "Number of counts", y = "Size Factor") +
  theme_classic()

###############
## Normalise ##
###############

# sce <- logNormCounts(sce)

#########################
## Add sample metadata ##
#########################

metadata <- fread("/Users/ricard/data/gastrulation/sample_metadata.txt") %>%
  .[id_rna%in%colnames(sce)]

all(colnames(sce)%in%metadata$id_rna)
all(metadata$id_rna%in%colnames(sce))

metadata <- metadata %>% setkey(id_rna) %>% .[colnames(sce)]

colData(sce) <- metadata %>% as.data.frame %>% tibble::remove_rownames() %>%
  tibble::column_to_rownames("id_rna") %>% .[colnames(counts),] %>% DataFrame

#######################
## Add gene metadata ##
#######################

rowData(sce)$symbol <- rowData(a)[rownames(sce),]$gene_name
# rowData(sce)$ens_id <- rownames(sce)

sum(duplicated(rowData(sce)$symbol))
sum(duplicated(rownames(sce)))

##########
## Save ##
##########

saveRDS(sce, "/Users/ricard/data/gastrulation/rna/SingleCellExperiment.rds")
