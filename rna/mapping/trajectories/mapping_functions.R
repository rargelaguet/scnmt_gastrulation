getHVGs <- function(sce, block, min.mean = 1e-3){
  decomp <- modelGeneVar(sce, block=block)
  decomp <- decomp[decomp$mean > min.mean,]
  decomp$FDR <- p.adjust(decomp$p.value, method = "fdr")
  return(rownames(decomp)[decomp$p.value < 0.10])
}

getmode <- function(v, dist) {
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
  out <- list()
  celltypes_accrossK <- matrix(unlist(mapping$celltypes.mapped),
                               nrow=length(mapping$celltypes.mapped[[1]]),
                               ncol=length(mapping$celltypes.mapped))
  cellstages_accrossK <- matrix(unlist(mapping$cellstages.mapped),
                                nrow=length(mapping$cellstages.mapped[[1]]),
                                ncol=length(mapping$cellstages.mapped))
  out$celltype.score <- NULL
  for (i in 1:nrow(celltypes_accrossK)){
    p <- max(table(celltypes_accrossK[i,]))
    index <- which(table(celltypes_accrossK[i,]) == p)
    p <- p/length(mapping$celltypes.mapped)
    out$celltype.score <- c(out$celltype.score,p)
  }
  out$cellstage.score <- NULL
  for (i in 1:nrow(cellstages_accrossK)){
    p <- max(table(cellstages_accrossK[i,]))
    index <- which(table(cellstages_accrossK[i,]) == p)
    p <- p/length(mapping$cellstages.mapped)
    out$cellstage.score <- c(out$cellstage.score,p)
  }
  return(out)  
}

get_meta <- function(pca_atlas, meta_atlas, pca_query, meta_query, k = 10){
  
  # find knn
  knns <- BiocNeighbors::queryKNN(pca_atlas, pca_query, k = k, 
                                  get.index = TRUE, get.distance = FALSE)
  
  #get closest k matching cells
  k.mapped  <- t(apply(knns$index, 1, function(x) meta_atlas$cell[x]))
  celltypes <- t(apply(k.mapped, 1, function(x) meta_atlas$celltype[match(x, meta_atlas$cell)]))
  stages    <- t(apply(k.mapped, 1, function(x) meta_atlas$stage[match(x, meta_atlas$cell)]))
  celltype.mapped <- apply(celltypes, 1, function(x) getmode(x, 1:length(x)))
  stage.mapped    <- apply(stages, 1, function(x) getmode(x, 1:length(x)))
  out <- lapply(1:length(celltype.mapped), function(x){
    list(cells.mapped     = k.mapped[x,],
         celltype.mapped  = celltype.mapped[x],
         stage.mapped     = stage.mapped[x],
         celltypes.mapped = celltypes[x,],
         stages.mapped    = stages[x,])
  })
  names(out) <- meta_query$cell
  return(out)  
}


joint.normalisation <- function(sce_query, sce_atlas, cosineNorm = FALSE) {
  
  # stopifnot("sample"%in%colnames(colData(sce_atlas)))
  # block <- c(sce_atlas$sample,sce_query$sample) %>% as.factor
  block <- c(rep("atlas",ncol(sce_atlas)),rep("query",ncol(sce_query)))
  
  if (isTRUE(cosineNorm)) {
    assay <- "cosineNorm"
    
    # Log normalisation per data set
    sce_atlas <- multiBatchNorm(sce_atlas, batch=as.factor(sce_atlas$sample))
    if (length(unique(sce_atlas$sample))>1) {
      sce_atlas <- multiBatchNorm(sce_atlas, batch=as.factor(sce_atlas$sample))
    } else {
      sce_atlas <- logNormCounts(sce_atlas)
    }
    
    if (length(unique(sce_query$sample))>1) {
      sce_query <- multiBatchNorm(sce_query, batch=as.factor(sce_query$sample))
    } else {
      sce_query <- logNormCounts(sce_query)
    }
    
    # Cosine normalisation
    assay(sce_atlas, assay) <- cosineNorm(assay(sce_atlas, "logcounts"))
    assay(sce_query, assay) <- cosineNorm(assay(sce_query, "logcounts"))
    
    # Concatenate
    sce_all <- SingleCellExperiment(
      list(cosineNorm=cbind(assay(sce_atlas,assay), assay(sce_query,assay)))
    )
    
  } else {
    assay <- "logcounts"
    
    # Concatenate
    sce_all <- SingleCellExperiment(
      list(counts=Matrix::Matrix(cbind(counts(sce_atlas),counts(sce_query)),sparse=TRUE))
    )
    
    # Log Normalise
    sce_all <- multiBatchNorm(sce_all, batch=block)
  }
  
  # Create block vector
  # sce_all$block_extended <- c(sce_atlas$sample,sce_query$sample) %>% as.factor
  sce_all$block <- c(rep("atlas",ncol(sce_atlas)),rep("query",ncol(sce_query))) %>% as.factor
  
  return(sce_all)
}

mapWrap <- function(sce_atlas, meta_atlas, sce_query, meta_query, k = 30, npcs = 50, genes = NULL, cosineNorm = FALSE){
  
  # Normalisation
  message(sprintf("Normalizing joint dataset using cosineNorm=%s...",cosineNorm))
  sce_all <- joint.normalisation(sce_query, sce_atlas, cosineNorm)
  message("Done\n")
  
  # Feature selection
  if (is.null(genes)) {
    message("Genes not provided. Computing highly variable genes...")
    # hvgs <- getHVGs(sce_all, block=c(meta_atlas$sample, meta_query$sample))
    genes <- getHVGs(sce_all, block=sce_all$block)
    message("Done\n")
  } else {
    message(sprintf("%d Genes provided...",length(genes)))
  }
  
  # Dimensionality reduction
  message("Performing PCA...")
  big_pca <- multiBatchPCA(
    sce_all,
    batch = sce_all$block,
    subset.row = genes,
    d = npcs,
    preserve.single = TRUE,
    assay.type = if (cosineNorm) "cosineNorm" else "logcounts"
  )[[1]]
  rownames(big_pca) <- colnames(sce_all) 
  atlas_pca <- big_pca[1:ncol(sce_atlas),]
  query_pca   <- big_pca[-(1:ncol(sce_atlas)),]
  message("Done\n")
  
  # Mapping query to batch-corrected atlas
  message("MNN mapping...")              
  correct <- reducedMNN(rbind(atlas_pca, query_pca), batch = as.character(sce_all$block))$corrected
  pca_atlas_corrected <- correct[1:nrow(atlas_pca),]
  pca_query_corrected   <- correct[-(1:nrow(atlas_pca)),]
  
  mapping <- get_meta(pca_atlas = pca_atlas_corrected,
                      meta_atlas = meta_atlas,
                      pca_query = pca_query_corrected,
                      meta_query = meta_query,
                      k = k)
  message("Done\n")
  
  # Mapping scores
  message("Computing mapping scores...") 
  out <- list()
  for (i in seq(from = 1, to = k)) {
    out$closest.cells[[i]]     <- sapply(mapping, function(x) x$cells.mapped[i])
    out$celltypes.mapped[[i]]  <- sapply(mapping, function(x) x$celltypes.mapped[i])
    out$cellstages.mapped[[i]] <- sapply(mapping, function(x) x$stages.mapped[i])
  }  
  multinomial.prob <- getMappingScore(out)
  message("Done\n")
  
  # Prepare output
  message("Writing output...") 
  out$pca_atlas_corrected <- pca_atlas_corrected
  out$pca_query_corrected <- pca_query_corrected
  ct <- sapply(mapping, function(x) x$celltype.mapped); is.na(ct) <- lengths(ct) == 0
  st <- sapply(mapping, function(x) x$stage.mapped); is.na(st) <- lengths(st) == 0
  cm <- sapply(mapping, function(x) x$cells.mapped[1]); is.na(cm) <- lengths(cm) == 0
  out$mapping <- data.frame(
    cell            = names(mapping), 
    celltype.mapped = unlist(ct),
    stage.mapped    = unlist(st),
    closest.cell    = unlist(cm))
  
  out$mapping <- cbind(out$mapping,multinomial.prob)
  out$pca <- big_pca
  message("Done\n")
  
  return(out)
  
}
