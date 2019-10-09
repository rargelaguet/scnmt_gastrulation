##########################################################################################
## Functions for mapping scRNAseq experiments to the 10X mouse gastrulation atlas 
##
## Code adapted from: 
## https://github.com/MarioniLab/EmbryoTimecourse2018/blob/master/analysis_scripts/atlas/core_functions.R
##
## Author: Ivan Imaz
##
## Packages to be installed:
## scran, scater, biomaRt, Matrix, irlba, BiocNeighbors, SingleCellExperiment
##########################################################################################

ac <- function(x, alpha=1, ...) { y <- adjustcolor(x, alpha.f=alpha, ...); names(y) <- names(x); return(y)}

celltype_colours <- c("Epiblast" = "#635547",
                     "Primitive Streak" = "#DABE99",
                     "Caudal epiblast" = "#9e6762",
                     
                     "PGC" = "#FACB12",
                     
                     "Anterior Primitive Streak" = "#c19f70",
                     "Notochord" = "#0F4A9C",
                     "Def. endoderm" = "#F397C0",
                     "Gut" = "#EF5A9D",
                     
                     "Nascent mesoderm" = "#C594BF",
                     "Mixed mesoderm" = "#DFCDE4",
                     "Intermediate mesoderm" = "#139992",
                     "Caudal Mesoderm" = "#3F84AA",
                     "Paraxial mesoderm" = "#8DB5CE",
                     "Somitic mesoderm" = "#005579",
                     "Pharyngeal mesoderm" = "#C9EBFB",
                     "Cardiomyocytes" = "#B51D8D",
                     "Allantois" = "#532C8A",
                     "ExE mesoderm" = "#8870ad",
                     "Mesenchyme" = "#cc7818",
                     
                     "Haematoendothelial progenitors" = "#FBBE92",
                     "Endothelium" = "#ff891c",
                     "Blood progenitors 1" = "#f9decf",
                     "Blood progenitors 2" = "#c9a997",
                     "Erythroid1" = "#C72228",
                     "Erythroid2" = "#f79083",
                     "Erythroid3" = "#EF4E22",
                     
                     "NMP" = "#8EC792",
                     
                     "Rostral neurectoderm" = "#65A83E",
                     "Caudal neurectoderm" = "#354E23",
                     "Neural crest" = "#C3C388",
                     "Forebrain/Midbrain/Hindbrain" = "#647a4f",
                     "Spinal cord" = "#CDE088",
                     
                     "Surface ectoderm" = "#f7f79e",
                     
                     "Visceral endoderm" = "#F6BFCB",
                     "ExE endoderm" = "#7F6874",
                     "ExE ectoderm" = "#989898",
                     "Parietal endoderm" = "#1A1A1A"
                     
)

haem_colours <- c(
  "Mes1"= "#c4a6b2",
  "Mes2"= "#ca728c",
  
  "Cardiomyocytes" =  "#B51D8D",  
  
  "BP1" = "#6460c5",
  "BP2" = "#96b8e4",
  "Haem3"= "#02f9ff",
  "BP3" = "#07499f",
  "BP4" = "#036ef9",
  
  "Haem1"= "#bb22a7",
  "Haem2" = "#f695e9",
  "Haem4" = "#4c4a81",
  
  "EC1"= "#006737",
  
  "EC2" = "#5eba46",
  "EC3" = "#818068",
  "EC4"="#d6de22",
  "EC5"="#5f7238",
  "EC6"="#2ab572",
  "EC7"="#000000",
  "EC8"="#a0cb3b",
  
  "Ery1"="#f67a58",
  "Ery2" ="#a26724",
  "Ery3"="#cdaf7f",
  "Ery4"= "#625218",
  
  "My" = "#c62127",
  "Mk"= "#f6931d"
)

stage_colours <- c("E6.5" = "#D53E4F",
                  "E6.75" = "#F46D43",
                  "E7.0" = "#FDAE61",
                  "E7.25" = "#FEE08B",
                  "E7.5" = "#FFFFBF",
                  "E7.75" = "#E6F598",
                  "E8.0" = "#ABDDA4",
                  "E8.25" = "#66C2A5",
                  "E8.5" = "#3288BD",
                  "mixed_gastrulation" = "#A9A9A9")

stage_labels <- names(stage_colours)
names(stage_labels) <- names(stage_colours)
stage_labels[10]    <- "Mixed"

load_data <- function(normalise = TRUE, remove_doublets = FALSE, remove_stripped = FALSE, load_corrected = FALSE, path2atlas=NULL){
  
  if(is.null(path2atlas)){
    path2atlas <- "/nfs/research1/marioni/iimaz/embryo_integ/atlas/"
  }
  
  if(load_corrected & (!remove_doublets | !remove_stripped)){
    message("Using corrected PCs, also removing doublets + stripped now.")
    remove_doublets <- TRUE
    remove_stripped <- TRUE
  }
  

  #counts  <- readRDS(paste0(path2atlas, "raw_counts.rds"))  
  counts  <- Matrix::readMM(paste0(path2atlas, "raw_counts.mtx"))
  genes   <- read.table(paste0(path2atlas, "genes.tsv"), stringsAsFactors = F)
  meta    <- read.table(paste0(path2atlas, "meta.tab"), header = TRUE, 
             sep = "\t", stringsAsFactors = FALSE, comment.char = "$")
  
  rownames(counts) <- genes[,1] #ensembl
  colnames(counts) <- meta$cell

  sce <- SingleCellExperiment::SingleCellExperiment(assays = list("counts" = counts))
  
  if(normalise){
    sfs <- read.table(paste0(path2atlas, "sizefactors.tab"), stringsAsFactors = F)[,1]
    SingleCellExperiment::sizeFactors(sce) <- sfs
    sce <- scater::normalize(sce)
  }
  
  if(remove_doublets){
    sce  <- scater::normalize(sce[,!meta$doublet])
    meta <- meta[!meta$doublet,]
  }
  
  if(remove_stripped){
    sce  <- scater::normalize(sce[,!meta$stripped])
    meta <- meta[!meta$stripped, ]
  }
  
  if(load_corrected){
    corrected <- readRDS(paste0(path2atlas, "corrected_pcas.rds"))
    assign("corrected", corrected, envir = .GlobalEnv)
    
  }
  
  assign("genes", genes, envir = .GlobalEnv)
  assign("meta", meta, envir = .GlobalEnv)
  assign("sce", sce, envir = .GlobalEnv)
  
  
  invisible(0)
}


#ensure counts has columns names for the cells
#match timepoints,samples to the count table
#timepoint_order, sample_order should contain each sample/timepoint ONCE, in correct order for correction
doBatchCorrect <- function(counts, timepoints, samples, timepoint_order, sample_order, npc = 50, pc_override = NULL, BPPARAM = SerialParam()){
  require(BiocParallel)
  
  if(!is.null(pc_override)){
    pca = pc_override
  } else {
    pca = irlba::prcomp_irlba(t(counts), n = npc)$x
    rownames(pca) = colnames(counts)
  }
  
  if(length(unique(samples)) == 1){
    return(pca)
  }
  
  #create nested list
  pc_list    <- lapply(unique(timepoints), function(tp){
    sub_pc   <- pca[timepoints == tp, , drop = FALSE]
    sub_samp <- samples[timepoints == tp]
    list     <- lapply(unique(sub_samp), function(samp){
      sub_pc[sub_samp == samp, , drop = FALSE]
    })
    names(list) <- unique(sub_samp)
    return(list)
  })
  
  names(pc_list) <- unique(timepoints)
  
  #arrange to match timepoint order
  pc_list <- pc_list[order(match(names(pc_list), timepoint_order))]
  pc_list <- lapply(pc_list, function(x){
    x[order(match(names(x), sample_order))]
  })
  
  #perform corrections within list elements (i.e. within stages)
  correct_list <- lapply(pc_list, function(x){
    if(length(x) > 1){
      return(do.call(batchelor::fastMNN, c(x, "pc.input" = TRUE, BPPARAM = BPPARAM))$corrected)
    } else {
      return(x[[1]])
    }
  })
  
  #perform correction over list
  if(length(correct_list)>1){
    correct <- do.call(batchelor::fastMNN, c(correct_list, "pc.input" = TRUE, BPPARAM = BPPARAM))$corrected
  } else {
    correct <- correct_list[[1]]
  }
  
  correct <- correct[match(colnames(counts), rownames(correct)),]
  
  return(correct)
  
}

getHVGs <- function(sce, min.mean = 1e-3, chrY.genes.file = NULL){
  trend  <- scran::trendVar(sce, use.spikes = FALSE, loess.args = list(span = 0.05))
  decomp <- scran::decomposeVar(sce, fit = trend)
  decomp <- decomp[decomp$mean > min.mean,]
  
  #exclude sex genes
  xist <- "ENSMUSG00000086503"
  if(!is.null(chrY.genes.file)){
    #ychr <- read.table("/nfs/research1/marioni/jonny/embryos/data/ygenes.tab", stringsAsFactors = FALSE)[,1]
    ychr <- read.table(chrY.genes.file, stringsAsFactors = FALSE)[,1]
  }else{
    mouse_ensembl <- biomaRt::useMart("ensembl")
    mouse_ensembl <- biomaRt::useDataset("mmusculus_gene_ensembl", mart = mouse_ensembl)
    gene_map <- biomaRt::getBM(attributes=c("ensembl_gene_id", "chromosome_name"),
      filters = "ensembl_gene_id", values = rownames(decomp), mart = mouse_ensembl)
    ychr <- gene_map[gene_map[,2] == "Y", 1]  
  }
  
  other = c("tomato-td") #for the chimera
  decomp = decomp[!rownames(decomp) %in% c(xist, ychr, other),]
  
  decomp$FDR = p.adjust(decomp$p.value, method = "fdr")
  return(rownames(decomp)[decomp$p.value < 0.05])
}

#MAPPING FUNCTIONS

getmode <- function(v, dist) {
  tab <- table(v)
  #if tie, break to shortest distance
  if(sum(tab==max(tab)) > 1){
    tied <- names(tab)[tab == max(tab)]
    sub  <- dist[v %in% tied]
    names(sub) <- v[v %in% tied]
    return(names(sub)[which.min(sub)])
  } else {
    return(names(tab)[which.max(tab)])
  }
}

getcelltypes <- function(v, dist) {
  tab <- table(v)
  #if tie, break to shortest distance
  if(sum(tab == max(tab)) > 1){
    tied <- names(tab)[tab == max(tab)]
    sub  <- dist[v %in% tied]
    names(sub) <- v[v %in% tied]
    return(names(sub)[which.min(sub)])
  } else {
    return(names(tab)[which.max(tab)])
  }
}

getMappingScore <- function(mapping){
  celltypes_accrossK <- matrix(unlist(mapping$celltypes.mapped), 
    nrow=length(mapping$celltypes.mapped[[1]]),
    ncol=length(mapping$celltypes.mapped))
  P <- NULL
  for (i in 1:nrow(celltypes_accrossK)){
    p <- max(table(celltypes_accrossK[i,]))
    index <- which(table(celltypes_accrossK[i,]) == p)
    p <- p/length(mapping$celltypes.mapped)
    P <- c(P,p) 
  }
  return(P)  
}

mnnMap <- function(atlas_pca, atlas_meta, map_pca, map_meta, k_map = 10){
  correct <- scran::fastMNN(atlas_pca, map_pca, pc.input = TRUE)$corrected
  atlas   <- 1:nrow(atlas_pca)
  correct_atlas <- correct[atlas,]
  correct_map   <- correct[-atlas,]
  
  knns <- BiocNeighbors::queryKNN(correct_atlas, correct_map, k = k_map, get.index = TRUE,
    get.distance = FALSE)
  
  #get closest k matching cells
  k.mapped  <- t(apply(knns$index, 1, function(x) atlas_meta$cell[x]))
  celltypes <- t(apply(k.mapped, 1, function(x) atlas_meta$celltype[match(x, atlas_meta$cell)]))
  stages    <- t(apply(k.mapped, 1, function(x) atlas_meta$stage[match(x, atlas_meta$cell)]))
  celltype.mapped <- apply(celltypes, 1, function(x) getmode(x, 1:length(x)))
  stage.mapped    <- apply(stages, 1, function(x) getmode(x, 1:length(x)))
  
  out <- lapply(1:length(celltype.mapped), function(x){
    list(cells.mapped     = k.mapped[x,],
         celltype.mapped  = celltype.mapped[x],
         stage.mapped     = stage.mapped[x],
         celltypes.mapped = celltypes[x,],
         stages.mapped    = stages[x,])
  })
  
  names(out) <- map_meta$cell
  
  return(out)  
  
}

#atlas_meta MUST have a cell column, celltype column and a stage column, spelled exactly like that
#map_meta MUST have a cell column and a stage column, spelled exactly like that
mapWrap <- function(atlas_sce, atlas_meta, map_sce, map_meta, k = 30, mapstage_x=NULL, map2stage_x=NULL, return.list = FALSE){
  
  #prevent duplicate rownames
  # colnames(map_sce) <- paste0("map_", colnames(map_sce))
  # map_meta$cell     <- paste0("map_", map_meta$cell)

  message("Normalizing joint dataset...")
  #easier to avoid directly binding sce objects as it is more likely to have issues
  sce_all <- SingleCellExperiment::SingleCellExperiment(
    list(counts=Matrix::Matrix(cbind(counts(atlas_sce),counts(map_sce)),sparse=TRUE)))
  big_sce <- scater::normalize(sce_all)
  message("Done\n")

  message("Computing highly variable genes...")
  hvgs    <- getHVGs(big_sce)
  message("Done\n")
  
  message("Performing PCA...")
  big_pca <- irlba::prcomp_irlba(t(logcounts(big_sce[hvgs,])), n = 50)$x
  rownames(big_pca) <- colnames(big_sce) 
  atlas_pca <- big_pca[1:ncol(atlas_sce),]
  map_pca   <- big_pca[-(1:ncol(atlas_sce)),]
  message("Done\n")
  
  message("Batch effect correction...")  
  #correct the atlas first
  order_df        <- atlas_meta[!duplicated(atlas_meta$sample), c("stage", "sample")]
  order_df$ncells <- sapply(order_df$sample, function(x) sum(atlas_meta$sample == x))
  order_df$stage  <- factor(order_df$stage, 
                          levels = rev(c("E8.5", 
                                         "E8.25", 
                                         "E8.0", 
                                         "E7.75", 
                                         "E7.5", 
                                         "E7.25", 
                                         "mixed_gastrulation", 
                                         "E7.0", 
                                         "E6.75", 
                                         "E6.5")))
  order_df       <- order_df[order(order_df$stage, order_df$ncells, decreasing = TRUE),]
  order_df$stage <- as.character(order_df$stage)
  
  set.seed(42)
  atlas_corrected <- doBatchCorrect(counts         = logcounts(atlas_sce[hvgs,]), 
                                   timepoints      = atlas_meta$stage, 
                                   samples         = atlas_meta$sample, 
                                   timepoint_order = order_df$stage, 
                                   sample_order    = order_df$sample, 
                                   pc_override     = atlas_pca)
  message("Done\n")
           
  message("MNN mapping...")                        
  if(!is.null(map2stage_x)){
    atlas_stage_index <- NULL
    for (i in seq(from = 1, to = length(map2stage_x))){
      index1 <- which(atlas_meta$stage==map2stage_x[i])
      atlas_stage_index <- c(atlas_stage_index,index1)
    } 
  }
  if(!is.null(mapstage_x)){
    map_stage_index   <- NULL
    for (i in seq(from = 1, to = length(mapstage_x))){    
      index2 <- which(map_meta$stage==mapstage_x[i])
      map_stage_index <- c(map_stage_index,index2)
    } 
    map_meta$stage <- map_meta$stage[map_stage_index]
    map_meta$cell  <- map_meta$cells[map_stage_index]
  }
  
  if(!is.null(map2stage_x) & !is.null(mapstage_x)){
    mapping <- mnnMap(atlas_pca = atlas_corrected[atlas_stage_index,],
                   atlas_meta = atlas_meta[atlas_stage_index,],
                   map_pca = map_pca[map_stage_index,],
                   map_meta = map_meta,k_map = k) 
  }else if(!is.null(map2stage_x) & is.null(mapstage_x)){
    mapping <- mnnMap(atlas_pca = atlas_corrected[atlas_stage_index,],
                   atlas_meta = atlas_meta[atlas_stage_index,],
                   map_pca = map_pca,
                   map_meta = map_meta,k_map = k)   
  }else if(is.null(map2stage_x) & !is.null(mapstage_x)){
    mapping <- mnnMap(atlas_pca = atlas_corrected,
                   atlas_meta = atlas_meta,
                   map_pca = map_pca[map_stage_index,],
                   map_meta = map_meta,k_map = k) 
  }else{
    mapping <- mnnMap(atlas_pca = atlas_corrected,
                   atlas_meta = atlas_meta,
                   map_pca = map_pca,
                   map_meta = map_meta,k_map = k)
  }
  message("Done\n")
             
  if(return.list){
    return(mapping)
  }
  
  message("Computing mapping scores...") 
  out <- list()
  for (i in seq(from = 1, to = k)){
   out$closest.cells[[i]]     <- sapply(mapping, function(x) x$cells.mapped[i])
   out$celltypes.mapped[[i]]  <- sapply(mapping, function(x) x$celltypes.mapped[i])
   out$cellstages.mapped[[i]] <- sapply(mapping, function(x) x$stages.mapped[i])
  }  
  celltype.multinomial.prob <- getMappingScore(out)
  message("Done\n")
  
  message("Writing output...") 
  out$mapping <- data.frame(
  cell            = names(mapping), 
  celltype.mapped = sapply(mapping, function(x) x$celltype.mapped),
  stage.mapped    = sapply(mapping, function(x) x$stage.mapped),
  closest.cell    = sapply(mapping, function(x) x$cells.mapped[1]))

  out$mapping <- cbind(out$mapping,celltype.multinomial.prob)
  message("Done\n")
  
  return(out)
  
}
