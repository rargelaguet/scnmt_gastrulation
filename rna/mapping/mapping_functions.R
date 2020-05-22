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
        #return(do.call(scran::fastMNN, c(x, "pc.input" = TRUE, BPPARAM = BPPARAM))$corrected)
        return(do.call(reducedMNN, c(x, BPPARAM = BPPARAM))$corrected) # edited 09.02.2020 because of "Error: 'fastMNN' is not an exported object from 'namespace:scran'", 17.02.2020 changed to reducedMNN because otherwise it thinks PCA space is logcounts which would be utter bullcrap
    } else {
      return(x[[1]])
    }
  })
  
  #perform correction over list
  if(length(correct_list)>1){
      #correct <- do.call(scran::fastMNN, c(correct_list, "pc.input" = TRUE, BPPARAM = BPPARAM))$corrected
      correct <- do.call(reducedMNN, c(correct_list, BPPARAM = BPPARAM))$corrected # edited 09.02.2020 because of "Error: 'fastMNN' is not an exported object from 'namespace:scran'", 17.02.2020 changed to reducedMNN because otherwise it thinks PCA space is logcounts which would be utter bullcrap
  } else {
    correct <- correct_list[[1]]
  }
    
  correct <- correct[match(colnames(counts), rownames(correct)),]
  
  return(correct)
  
}

getHVGs <- function(sce, block, min.mean = 1e-3){
  decomp <- modelGeneVar(sce, block=block)
  decomp <- decomp[decomp$mean > min.mean,]
  decomp$FDR <- p.adjust(decomp$p.value, method = "fdr")
  return(rownames(decomp)[decomp$p.value < 0.01])
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

get_meta <- function(correct_atlas, atlas_meta, correct_map, query_meta, k_map = 10){
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
  names(out) <- query_meta$cell
  return(out)  
}

mapWrap <- function(atlas_sce, atlas_meta, query_sce, query_meta, k = 30, npcs = 50, genes = NULL){
    
  message("Normalizing joint dataset...")
  
  #easier to avoid directly binding sce objects as it is a lot more likely to have issues
  sce_all <- SingleCellExperiment(
    list(counts=Matrix::Matrix(cbind(counts(atlas_sce),counts(query_sce)),sparse=TRUE)))
  
  # Define block structure
  # block <- c(atlas_meta$sample, rep("query",ncol(query_sce)))
  block <- c(rep("atlas",ncol(atlas_sce)), rep("query",ncol(query_sce)))

  # big_sce <- scater::logNormCounts(sce_all) 
  # big_sce <- multiBatchNorm(sce_all, batch=c(atlas_meta$sample, query_meta$batch))
  big_sce <- multiBatchNorm(sce_all, batch=block)
  
  if (is.null(genes)) {
    message("Genes not provided. Computing highly variable genes...")
    hvgs <- getHVGs(big_sce, block=block)
  } else {
    hvgs <- genes
    message(sprintf("%d Genes provided...",length(genes)))
  }
  
  message("Performing PCA...")
  big_pca <- multiBatchPCA(big_sce, batch = block, subset.row = hvgs, d = npcs, preserve.single = TRUE, assay.type = "logcounts")[[1]]
  rownames(big_pca) <- colnames(big_sce) 
  atlas_pca <- big_pca[1:ncol(atlas_sce),]
  query_pca   <- big_pca[-(1:ncol(atlas_sce)),]
  
  message("Batch effect correction for the atlas...")  
  order_df        <- atlas_meta[!duplicated(atlas_meta$sample), c("stage", "sample")]
  order_df$ncells <- sapply(order_df$sample, function(x) sum(atlas_meta$sample == x))
  order_df$stage  <- factor(order_df$stage, levels = rev(c("E8.5","E8.25","E8.0","E7.75","E7.5","E7.25","mixed_gastrulation","E7.0","E6.75","E6.5")))
  order_df       <- order_df[order(order_df$stage, order_df$ncells, decreasing = TRUE),]
  order_df$stage <- as.character(order_df$stage)
  
  set.seed(42)
  atlas_corrected <- doBatchCorrect(counts         = logcounts(atlas_sce[hvgs,]), 
                                    timepoints      = atlas_meta$stage, 
                                    samples         = atlas_meta$sample, 
                                    timepoint_order = order_df$stage, 
                                    sample_order    = order_df$sample, 
                                    pc_override     = atlas_pca,
                                    npc             = npcs)
  
  message("MNN mapping...")                        
  batch <- c(rep("ATLAS",ncol(atlas_sce)), rep("QUERY",ncol(query_sce)))
  correct <- reducedMNN(rbind(atlas_corrected, query_pca), batch=batch)$corrected
  atlas   <- 1:nrow(atlas_pca)
  correct_atlas <- correct[atlas,]
  correct_map   <- correct[-atlas,]

  mapping <- get_meta(correct_atlas = correct_atlas, atlas_meta = atlas_meta,
                      correct_map = correct_map, query_meta = query_meta,
                      k_map = k)
  
  message("Computing mapping scores...") 
  out <- list()
  for (i in seq(from = 1, to = k)) {
    out$closest.cells[[i]]     <- sapply(mapping, function(x) x$cells.mapped[i])
    out$celltypes.mapped[[i]]  <- sapply(mapping, function(x) x$celltypes.mapped[i])
    # out$cellstages.mapped[[i]] <- sapply(mapping, function(x) x$stages.mapped[i])
  }  
  multinomial.prob <- getMappingScore(out)
  
  out$correct_atlas <- correct_atlas
  out$correct_map <- correct_map
  ct <- sapply(mapping, function(x) x$celltype.mapped); is.na(ct) <- lengths(ct) == 0
  # st <- sapply(mapping, function(x) x$stage.mapped); is.na(st) <- lengths(st) == 0
  cm <- sapply(mapping, function(x) x$cells.mapped[1]); is.na(cm) <- lengths(cm) == 0
  
  message("Writing output...") 
  out$mapping.dt <- data.table(
      id_rna          = names(mapping), 
      celltype.mapped = unlist(ct),
      # stage.mapped    = unlist(st),
      closest.cell    = unlist(cm),
      # stage.score = multinomial.prob$cellstage.score,
      celltype.score = multinomial.prob$celltype.score
    )
  
  # out$mapping <- cbind(out$mapping,multinomial.prob)
  out$pca <- big_pca
  
  return(out)
  
}