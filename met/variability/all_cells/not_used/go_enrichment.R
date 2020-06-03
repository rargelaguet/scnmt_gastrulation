
#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/scnmt_gastrulation/settings.R")
  io$dir.lib <- "/Users/ricard/bbreg/lib/"
  io$script <- "/Users/ricard/scnmt_gastrulation/metaccrna/dynamics_individual_examples/boxplots_argparse.R"
} else {
  stop("Computer not recognised")
}

R.utils::sourceDirectory(io$dir.lib, modifiedOnly = FALSE)

io$dir.bayes <- paste0(io$scmet,"/betabinomial/bayes/stan")
# io$out.dir <- paste0(io$scmet,"/betabinomial/bayes/pdf")


opts$anno <- c(
  "prom_2000_2000"
)
###############
## Load data ##
###############

# Load HVFs
hvf_res <- list()
for (i in 1:length(opts$anno)) {
  stan <- readRDS(sprintf("%s/allcells_%s_vb.rds", io$dir.bayes, opts$anno[i]))
  hvf_res[[opts$anno[i]]] <- detect_hvf(stan, delta_e = NULL, delta_g = 0.25, efdr = 0.1)
  # hvf_res[[opts$anno[i]]] <- detect_hvf(stan, delta_e = 0.9, delta_g = NULL, efdr = 0.1)
}
rm(stan)

# Load gene metadata
gene_metadata <- fread(io$gene.metadata) %>% 
  .[,c("ens_id","symbol")]

# Load gene set
# io$geneset <- "/Users/ricard/data/MSigDB/v6.0/mus_musculus/C5/bp_binary_matrix_ensembl.rds"
io$geneset <- "/Users/ricard/data/MSigDB/v6.0/mus_musculus/C2/bp_binary_matrix_ensembl.rds"
gene.sets <- readRDS(io$msigFile)

# io$geneset <- "/Users/ricard/data/reactome/v59/mus_musculus/out/mouse_v75_reactome.rds"
# gene.sets <- readRDS(io$geneset)

colnames(gene.sets) <- toupper(colnames(gene.sets))



################
## Parse data ##
################

dt <- hvf_res$prom_2000_2000$summary %>% as.data.table %>%
  setnames("feature_name","ens_id") %>%
  merge(gene_metadata, by="ens_id") %>%
  .[,symbol:=toupper(symbol)]

#########################
## Fisher's exact test ##
#########################

enrichment.dt <- rownames(gene.sets) %>% map(function(i) {
  
  genes <- names(which(gene.sets[i,]==1))
  
  foo = matrix(data = c(
    dt[symbol%in%genes & hvf==TRUE,.N], dt[symbol%in%genes & hvf==FALSE,.N],
    dt[(!symbol%in%genes) & hvf==TRUE,.N], dt[(!symbol%in%genes) & hvf==FALSE,.N]),
    nrow = 2, ncol = 2, dimnames = list(GeneSet = c("True", "False"), HVF = c("True", "False")))

  bar <- fisher.test(foo, alternative = "greater")
  return(data.table(geneset=i, pval=bar$p.value))
}) %>% rbindlist





# Rename pathways
# tmp <- read.table("/Users/ricard/data/reactome/v59/mus_musculus/AllPathways.txt", header=F, quote="", sep="\t", stringsAsFactors=F)[,c(1,2)]
# reactome_meta <- tmp[,2]; names(reactome_meta) <- tmp[,1]
# enrichment.dt$geneset <- stringr::str_replace_all(enrichment.dt$geneset, reactome_meta)
