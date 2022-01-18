
matrix.please<-function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[[1]]
  m
}
# common parameters for fwriting and ensure folder exists
fwrite_tsv <- function(x, file = "", ...){
  dir <- dirname(file)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  fwrite(x, file, sep = "\t", na = "NA", quote = FALSE, ...)
}


#############################
## Commonly-used functions ##
#############################

load_SingleCellExperiment <- function(file, normalise = FALSE, features = NULL, cells = NULL, remove_non_expressed_genes = FALSE) {
  library(SingleCellExperiment); library(scran); library(scater);
  sce <- readRDS(file)
  if (!is.null(cells)) sce <- sce[,cells]
  if (!is.null(features)) sce <- sce[features,]
  if (remove_non_expressed_genes) sce <- sce[which(Matrix::rowSums(counts(sce))>15),]
  if (normalise) sce <- logNormCounts(sce)
  return(sce)
}

load_Seurat <- function(file, assay = "RNA", normalise = FALSE, features = NULL, cells = NULL, remove_non_expressed_genes = FALSE, ...) {
  library(Seurat)
  seurat <- readRDS(file)
  # if (assay%in%Seurat::Assays(seurat)) seurat <- seurat[[assay]]
  if (!is.null(cells)) seurat <- seurat[,cells]
  if (!is.null(features)) seurat <- seurat[features,]
  if (normalise) {
    seurat <- NormalizeData(seurat, normalization.method = "LogNormalize")
    seurat <- ScaleData(seurat, ...)
  }
  if (remove_non_expressed_genes) seurat <- seurat[which(Matrix::rowMeans(seurat@assays[[assay]]@counts)>1e-4),]
  return(seurat)
}

matrix.please<-function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[[1]]
  m
}

pdist <- function(tmat){
  # @param tmat A non-negative matrix with samples by features
  # @reference http://r.789695.n4.nabble.com/dist-function-in-R-is-very-slow-td4738317.html
  mtm <- Matrix::tcrossprod(tmat)
  sq <- rowSums(tmat^2)
  out0 <- outer(sq, sq, "+") - 2 * mtm
  out0[out0 < 0] <- 0
  
  sqrt(out0)
}

smoother_aggregate_nearest_nb <- function(mat, D, k){
  # @param mat A matrix in a shape of #genes x #samples.
  # @param D A predefined distance matrix in a shape of #samples x #samples.
  # @param k An integer to choose \code{k} nearest samples (self-inclusive) to
  #  aggregate based on the distance matrix \code{D}. If \code{k} is greater than
  #  #samples, \code{k} is forced to be #samples to continue aggregation.
  sapply(seq_len(ncol(mat)), function(cid){
    nb_cid <- head(order(D[cid, ]), k)
    closest_mat <- mat[, nb_cid, drop=FALSE]
    # return(Matrix::rowSums(closest_mat))
    return(Matrix::rowMeans(closest_mat,na.rm=T))
  })
}


# TO-FINISH.....
smoother_aggregate_nearest_nb_parallel <- function(mat, D, k, cores=1){
  # @param mat A matrix in a shape of #genes x #samples.
  # @param D A predefined distance matrix in a shape of #samples x #samples.
  # @param k An integer to choose \code{k} nearest samples (self-inclusive) to
  #  aggregate based on the distance matrix \code{D}.

  # library(future)
  library(future.apply)
  plan("multiprocess", workers = ncores)


  # sapply(seq_len(ncol(mat)), function(cid){
  future_sapply(seq_len(ncol(mat)), function(cid){
    nb_cid <- head(order(D[cid, ]), k)
    closest_mat <- mat[, nb_cid, drop=FALSE]
    # return(Matrix::rowSums(closest_mat))
    return(Matrix::rowMeans(closest_mat))
  })
}

# regress_covariates <- function(mtx, vars.to.regress) {
#   data <- scale(t(logcounts(sce_filt)), center = T, scale = F)
#   data_regressed <- apply(data, 2, function(x) {
#     lm.out <- lm(formula=expr~covariate, data=data.frame(expr=x, covariate=factor(sce_filt$stage)));
#     residuals <- lm.out[["residuals"]]+lm.out[["coefficients"]][1]
#   })
# }

# Remove unwanted effects from a matrix
#
# @parm mtx An expression matrix to regress the effects of covariates out
# of should be the complete expression matrix in genes x cells
# @param covariates A matrix or data.frame of latent variables, should be cells
# x covariates, the colnames should be the variables to regress
# @param features_idx An integer vector representing the indices of the
# genes to run regression on
# @param model.use Model to use, one of 'linear', 'poisson', or 'negbinom'; pass
# NULL to simply return mtx
# @param verbose Display a progress bar
#' @importFrom stats as.formula lm
#' @importFrom utils txtProgressBar setTxtProgressBar
#
RegressOutMatrix_parallel <- function(mtx, covariates = NULL, features_idx = NULL, split.by = NULL, block.size = 1000, min.cells.to.block = 3000, ncores = 1, verbose = TRUE) {
  
  library(future)
  library(future.apply)
  plan("multiprocess", workers = ncores)
  
  # Check features_idx
  if (is.null(features_idx)) {
    features_idx <- 1:nrow(mtx)
  }
  if (is.character(features_idx)) {
    features_idx <- intersect(features_idx, rownames(mtx))
    if (length(features_idx) == 0) {
      stop("Cannot use features that are beyond the scope of mtx")
    }
  } else if (max(features_idx) > nrow(mtx)) {
    stop("Cannot use features that are beyond the scope of mtx")
  }
  
  # Check dataset dimensions
  if (nrow(covariates) != ncol(mtx)) {
    stop("Uneven number of cells between latent data and expression data")
  }
  
  # Subset
  mtx <- mtx[features_idx,]
  mtx.dimnames <- dimnames(mtx)
  
  # Define chunck points
  chunk.points <- ChunkPoints(dsize = nrow(mtx), csize = block.size)
  
  # Define cell splitting
  split.cells <- split(colnames(mtx), f = split.by %||% TRUE)

   if (nbrOfWorkers() > 1) {

    # Define chuncks
      chunks <- expand.grid(
        names(split.cells),
        1:ncol(chunk.points),
        stringsAsFactors = FALSE
      )

      # Run RegressOutMatrix in parallel
      mtx.resid <- future_lapply(
        X = 1:nrow(chunks),
        FUN = function(i) {
          row <- chunks[i, ]
          group <- row[[1]]
          index <- as.numeric(row[[2]])
          return(RegressOutMatrix(
            mtx = mtx[chunk.points[1, index]:chunk.points[2, index], split.cells[[group]], drop = FALSE],
            covariates = covariates[split.cells[[group]], , drop = FALSE],
            # features_idx = features_idx[chunk.points[1, index]:chunk.points[2, index]],
            verbose = FALSE
          ))
        }
      )

      # Merge splitted cells
      if (length(split.cells) > 1) {
        merge.indices <- lapply(
          X = 1:length(x = split.cells),
          FUN = seq.int,
          to = length(mtx.resid),
          by = length(split.cells)
        )
        mtx.resid <- lapply(
          X = merge.indices,
          FUN = function(x) {
            return(do.call( 'rbind', mtx.resid[x]))
          }
        )
        mtx.resid <- do.call('cbind', mtx.resid)
      } else {
        mtx.resid <- do.call( 'rbind', mtx.resid)
      }
    } else {
      
      mtx.resid <- lapply(
        X = names(split.cells),
        FUN = function(x) {
          if (verbose && length(split.cells) > 1) {
            message("Regressing out variables from split ", x)
          }
          return(RegressOutMatrix(
            mtx = mtx[, split.cells[[x]], drop = FALSE],
            covariates = covariates[split.cells[[x]], , drop = FALSE],
            features_idx = features_idx,
            verbose = verbose
          ))
        }
      )
      mtx.resid <- do.call('cbind', mtx.resid)
    }
    # dimnames(mtx.resid) <- dimnames(mtx)
    return(mtx.resid)
  }

RegressOutMatrix <- function(mtx, covariates = NULL, features_idx = NULL, verbose = TRUE) {
  
  # Check features_idx
  if (is.null(features_idx)) {
    features_idx <- 1:nrow(mtx)
  }
  if (is.character(features_idx)) {
    features_idx <- intersect(features_idx, rownames(mtx))
    if (length(features_idx) == 0) {
      stop("Cannot use features that are beyond the scope of mtx")
    }
  } else if (max(features_idx) > nrow(mtx)) {
    stop("Cannot use features that are beyond the scope of mtx")
  }
  
  # Check dataset dimensions
  if (nrow(covariates) != ncol(mtx)) {
    stop("Uneven number of cells between latent data and expression data")
  }
  
  # Subset
  mtx <- mtx[features_idx,]
  mtx.dimnames <- dimnames(mtx)
  
  # Create formula for regression
  vars.to.regress <- colnames(covariates)
  fmla <- paste('GENE ~', paste(vars.to.regress, collapse = '+')) %>% as.formula

  # In this code, we'll repeatedly regress different Y against the same X
  # (covariates) in order to calculate residuals.  Rather that repeatedly
  # call lm to do this, we'll avoid recalculating the QR decomposition for the
  # covariates matrix each time by reusing it after calculating it once
  regression.mat <- cbind(covariates, mtx[1,])
  colnames(regression.mat) <- c(colnames(covariates), "GENE")
  qr <- lm(fmla, data = regression.mat, qr = TRUE)$qr
  rm(regression.mat)

  # Make results matrix
  data.resid <- matrix(
    nrow = nrow(mtx),
    ncol = ncol(mtx)
  )

  if (verbose) pb <- txtProgressBar(char = '=', style = 3, file = stderr())

  # Extract residuals from each feature by using the pre-computed QR decomposition
  for (i in 1:length(features_idx)) {
    regression.mat <- cbind(covariates, mtx[features_idx[i], ])
    colnames(regression.mat) <- c(vars.to.regress, 'GENE')
    regression.mat <- qr.resid(qr = qr, y = mtx[features_idx[i],])  # The function qr.resid returns the residuals when fitting y to the matrix with QR decomposition.
    data.resid[i, ] <- regression.mat
    if (verbose) {
      setTxtProgressBar(pb = pb, value = i / length(features_idx))
    }
  }

  if (verbose) close(con = pb)

  dimnames(data.resid) <- mtx.dimnames
  
  return(data.resid)
}



# Generate chunk points
#
# @param dsize How big is the data being chunked
# @param csize How big should each chunk be
#
# @return A matrix where each column is a chunk, row 1 is start points, row 2 is end points
#
ChunkPoints <- function(dsize, csize) {
  return(vapply(
    X = 1L:ceiling(dsize / csize),
    FUN = function(i) {
      return(c(
        start = (csize * (i - 1L)) + 1L,
        end = min(csize * i, dsize)
      ))
    },
    FUN.VALUE = numeric(length = 2L)
  ))
}


"%ni%" <- Negate("%in%")

ggplot_theme_NoAxes <- function() {
  theme(
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
}

minmax.normalisation <- function(x)
{
    return((x-min(x,na.rm=T)) /(max(x,na.rm=T)-min(x,na.rm=T)))
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

GRangesToString <- function(grange, sep = c("-", "-")) {
  regions <- paste0(
    as.character(x = seqnames(x = grange)),
    sep[[1]],
    start(x = grange),
    sep[[2]],
    end(x = grange)
  )
  return(regions)
}


give.n <- function(x){
  return(c(y = mean(x), label = length(x)))
}

sort.abs <- function(dt, sort.field) dt[order(-abs(dt[[sort.field]]))]


# function to pseudobulk a SingleCellExperiment object
pseudobulk_sce_fn <- function(x, assay = NULL, by, fun = c("sum", "mean", "median"), scale = FALSE) {
  
  # check validity of input arguments
  fun <- match.arg(fun)
  if (is.null(assay))  assay <- assayNames(x)[1] 
  
  # store aggregation parameters &
  # nb. of cells that went into aggregation
  md <- metadata(x)
  md$agg_pars <- list(assay = assay, by = by, fun = fun, scale = scale)
  
  # get aggregation function
  # fun <- switch(fun, sum = "rowSums", mean = "rowMeans", median = "rowMedians")
  
  # drop missing factor levels & tabulate number of cells
  cd <- dplyr::mutate_if(as.data.frame(colData(x)), is.factor, droplevels)
  colData(x) <- DataFrame(cd, row.names = colnames(x),check.names = FALSE)
  md$n_cells <- table(as.data.frame(colData(x)[, by]))
  
  # assure 'by' colData columns are factors so that missing combinations aren't dropped
  for (i in by) 
    if (!is.factor(x[[i]])) 
      x[[i]] <- factor(x[[i]])
  
  # split cells & compute pseudo-bulks
  cs <- .split_cells(x, by)
  # pb <- .pb(x, cs, assay, fun)
  pb <- .pb(x=x, by=by, fun=fun)
  if (scale & length(by) == 2) {
    ls <- lapply(.pb(x, cs, "counts", "rowSums"), colSums)
    pb <- lapply(seq_along(pb), function(i) pb[[i]] / 1e6 * ls[[i]])
    names(pb) <- names(ls)
  }
  
  # construct SCE
  pb <- SingleCellExperiment(pb, metadata = md)
  
  # propagate 'colData' columns that are unique across 2nd 'by'
  if (length(by) == 2) {
    cd <- colData(x)
    ids <- colnames(pb)
    counts <- vapply(ids, function(u) {
      m <- as.logical(match(cd[, by[2]], u, nomatch = 0))
      vapply(cd[m, ], function(u) length(unique(u)), numeric(1))
    }, numeric(ncol(colData(x))))
    cd_keep <- apply(counts, 1, function(u) all(u == 1))
    cd_keep <- setdiff(names(which(cd_keep)), by)
    if (length(cd_keep) != 0) {
      m <- match(ids, cd[, by[2]], nomatch = 0)
      cd <- cd[m, cd_keep, drop = FALSE]
      rownames(cd) <- ids
      colData(pb) <- cd
    }
  }
  return(pb)
}


# split cells by cluster-sample
# auxiliary function to pseudobulk a SingleCellExperiment object
#   - by: character vector specifying colData column(s) to split by. If length(by) == 1, a list of length nlevels(colData$by), else,
#          a nested list with 2nd level of length nlevels(colData$by[2])
.split_cells <- function(x, by) {
  if (is(x, "SingleCellExperiment")) x <- colData(x)
  cd <- data.frame(x[by], check.names = FALSE)
  cd <- data.table(cd, cell = rownames(x)) %>% split(by = by, sorted = TRUE, flatten = FALSE)
  purrr::map_depth(cd, length(by), "cell")
}


# auxiliary function to pseudobulk a SingleCellExperiment object
.pb <- function(x, by, fun) {
  
  # compute pseudobulks
  # y <- scuttle::summarizeAssayByGroup(x, assay.type = assay, ids = (ids <- colData(x)[by]), statistics = fun, BPPARAM = BiocParallel::SerialParam())
  y <- scuttle::summarizeAssayByGroup(x, ids = colData(x)[by], statistics = fun)
  colnames(y) <- y[[by[length(by)]]]
  
  if (length(by) == 1)  return(assay(y))
  
  # reformat into one assay per 'by[1]'
  if (is.factor(ids <- y[[by[1]]]))
    ids <- droplevels(ids)
  is <- split(seq_len(ncol(y)), ids)
  ys <- map(is, ~assay(y)[, .])
  
  # fill in missing combinations
  for (i in seq_along(ys)) {
    fill <- setdiff(unique(y[[by[2]]]), colnames(ys[[i]]))
    if (length(fill != 0)) {
      foo <- matrix(0, nrow(x), length(fill))
      colnames(foo) <- fill
      foo <- cbind(ys[[i]], foo)
      o <- paste(sort(unique(y[[by[2]]])))
      ys[[i]] <- foo[, o]
    }
  }
  return(ys)
}


mean_sd <- function(x) data.frame(y=mean(x), ymin=mean(x)-sd(x)/2, ymax=mean(x)+sd(x)/2) 