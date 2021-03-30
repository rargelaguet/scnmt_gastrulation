suppressMessages(library(data.table))
suppressMessages(library(purrr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(MOFA2))

matrix.please <- function(x) { m<-as.matrix(x[,-1]); rownames(m)<-x[[1]]; m }

############################
## Define I/O and options ##
############################

io <- list()
io$data <- "/Users/ricard/data/mofa2_vignettes/gastrulation_scnmt_mofa_analysis/data.txt.gz"
io$outdir <- "/Users/ricard/data/mofa2_vignettes/gastrulation_scnmt_mofa_analysis/pdf/pca"
io$mofa.model <- "/Users/ricard/data/gastrulation/mofa2/hdf5/test_1.hdf5"

dir.create(io$outdir, showWarnings = F)

###############
## Load data ##
###############

# Load mofa model
mofa <- load_model(io$mofa.model)

# fetch RNA expression data from the model
rna_mtx <- do.call("cbind",mofa@data[["RNA"]])
dim(rna_mtx)

#############
## Run PCA ##
#############

pca <- irlba::prcomp_irlba(rna_mtx, center = FALSE, n = 25)

pca.var.explained <- data.table(
  factor = as.character(1:ncol(pca$x)),
  r2 = (pca$sdev**2)/pca$totalvar,
  model = "PCA"
)

###########################################################
## Correlate with the variance explained by MOFA factors ##
###########################################################

# Extract variance explained estimates 
mofa.var.explained <- data.table(
  factor = as.character(1:mofa@dimensions$K),
  r2 = rowMeans(sapply(mofa@cache$variance_explained$r2_per_factor, function(x) x[,"RNA"])),
  model = "MOFA"
)

# Plot
to.plot <- rbind(pca.var.explained, mofa.var.explained) %>%
  .[factor%in%as.character(1:mofa@dimensions$K)] %>%
  .[,r2:=100*r2]

p <- ggline(to.plot, x="factor", y="r2", color="model") +
  labs(x="Factor index", y="Variance explained (%)") +
  theme(
    legend.title = element_blank()
  )


pdf(paste0(io$outdir,"/r2_mofa_vs_pca.pdf"), width=9, height=5, useDingbats = F)
print(p)
dev.off()