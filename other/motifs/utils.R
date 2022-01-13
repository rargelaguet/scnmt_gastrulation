
load_motif_annotation <- function(motifSet = "cisbp", species = "Homo sapiens", collection = "CORE", version = 2) {

  if(tolower(motifSet)=="jaspar2020"){
    
    .requirePackage("JASPAR2020",installInfo='BiocManager::install("JASPAR2020")')
    motifs <- TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020, list(species = species, collection = collection))
    obj <- .summarizeJASPARMotifs(motifs)

  }else if(tolower(motifSet)=="jaspar2018"){

    .requirePackage("JASPAR2018",installInfo='BiocManager::install("JASPAR2018")')
    args <- list(species = species, collection = collection, ...)
    motifs <- TFBSTools::getMatrixSet(JASPAR2018::JASPAR2018, args)
    obj <- .summarizeJASPARMotifs(motifs)

  }else if(tolower(motifSet)=="jaspar2016"){

    .requirePackage("JASPAR2016",installInfo='BiocManager::install("JASPAR2018")')
    args <- list(species = species, collection = collection, ...)
    motifs <- TFBSTools::getMatrixSet(JASPAR2016::JASPAR2016, args)
    obj <- .summarizeJASPARMotifs(motifs)

  }else if(tolower(motifSet)=="cisbp"){

    .requirePackage("chromVARmotifs",installInfo='devtools::install_github("GreenleafLab/chromVARmotifs")')
    if(tolower(species) == "mus musculus"){
      if(version == 1){
        message("Using version 1 motifs!")
        data("mouse_pwms_v1")
        motifs <- mouse_pwms_v1        
      }else if(version == 2){
        message("Using version 2 motifs!")
        data("mouse_pwms_v2")
        motifs <- mouse_pwms_v2
      }else{
        stop("Only versions 1 and 2 exist!")
      }
      obj <- .summarizeChromVARMotifs(motifs)
    }else if(tolower(species) == "homo sapiens"){
      if(version == 1){
        message("Using version 1 motifs!")
        data("human_pwms_v1")
        motifs <- human_pwms_v1        
      }else if(version == 2){
        message("Using version 2 motifs!")
        data("human_pwms_v2")
        motifs <- human_pwms_v2
      }else{
        stop("Only versions 1 and 2 exist!")
      }
      obj <- .summarizeChromVARMotifs(motifs)
    }else{
      stop("Species not recognized homo sapiens, mus musculus supported by CisBP!")
    }

  }else if(tolower(motifSet)=="encode"){

    .requirePackage("chromVARmotifs",installInfo='devtools::install_github("GreenleafLab/chromVARmotifs")')
    data("encode_pwms")
    obj <- .summarizeChromVARMotifs(encode_pwms)

  }else if(tolower(motifSet)=="homer"){

    .requirePackage("chromVARmotifs",installInfo='devtools::install_github("GreenleafLab/chromVARmotifs")')
    data("homer_pwms")
    obj <- .summarizeChromVARMotifs(homer_pwms)

  }else{

    stop("Error MotifSet Not Recognized!")

  }

  return(obj)
}


create_motif_overlap_matrix <- function(pwms, subject, genome = "BSgenome.Mmusculus.UCSC.mm10", cutOff = 5e-05, width = 7) {
  
  # sanity checks
  # stopifnot(inherits(pwms,"PWMatrixList"))
  # stopifnot(inherits(subject,"GRanges"))
  
  # Get BSgenome
  .requirePackage(genome)
  BSgenome <- eval(parse(text = genome))
  BSgenome <- validBSgenome(BSgenome)
  seqlevels(subject) <- paste0("chr",levels(seqnames(subject)))
  
  # Calculate Motif Positions
  motifPositions <- motifmatchr::matchMotifs(
    pwms = pwms,
    subject = subject,
    genome = BSgenome, 
    out = "positions", 
    p.cutoff = cutOff, 
    w = width
  )
  
  # Create Motif Overlap Matrix
  allPositions <- unlist(motifPositions)
  overlapMotifs <- findOverlaps(subject, allPositions, ignore.strand=TRUE)
  motifMat <- Matrix::sparseMatrix(
    i = queryHits(overlapMotifs),
    j = match(names(allPositions),names(motifPositions))[subjectHits(overlapMotifs)],
    x = rep(TRUE, length(overlapMotifs)),
    dims = c(length(subject), length(motifPositions))
  )
  colnames(motifMat) <- names(motifPositions)
  motifMat <- SummarizedExperiment::SummarizedExperiment(assays=SimpleList(matches = motifMat), rowRanges = subject)
  
  out <- SimpleList(
    motifMatches = motifMat,
    motifPositions = motifPositions
  )
  return(out)
}


.summarizeJASPARMotifs <- function(motifs = NULL){

  motifNames <- lapply(seq_along(motifs), function(x){
    namex <- make.names(motifs[[x]]@name)
    if(substr(namex,nchar(namex),nchar(namex))=="."){
      namex <- substr(namex,1,nchar(namex)-1)
    }
    namex <- paste0(namex, "_", x)
    namex
  }) %>% unlist(.)

  motifDF <- lapply(seq_along(motifs), function(x){
    data.frame(
      row.names = motifNames[x],
      name = motifs[[x]]@name[[1]],
      ID = motifs[[x]]@ID,
      strand = motifs[[x]]@strand,
      symbol = ifelse(!is.null(motifs[[x]]@tags$symbol[1]), motifs[[x]]@tags$symbol[1], NA) ,
      family = ifelse(!is.null(motifs[[x]]@tags$family[1]), motifs[[x]]@tags$family[1], NA),
      alias = ifelse(!is.null(motifs[[x]]@tags$alias[1]), motifs[[x]]@tags$alias[1], NA),
      stringsAsFactors = FALSE
    )
  }) %>% Reduce("rbind", .) %>% DataFrame
  
  names(motifs) <- motifNames

  out <- list(motifs = motifs, motifSummary = motifDF)

  return(out)
  
}

.summarizeChromVARMotifs <- function(motifs = NULL){

  motifNames <- lapply(seq_along(motifs), function(x){
    namex <- make.names(motifs[[x]]@name)
    if(grepl("LINE", namex)){
      splitNamex <- stringr::str_split(motifs[[x]]@ID, pattern="\\_", simplify = TRUE)
      namex <- splitNamex[1, grep("LINE",splitNamex[1,]) + 1]
    }
    if(substr(namex,nchar(namex),nchar(namex))=="."){
      namex <- substr(namex,1,nchar(namex)-1)
    }
    namex <- paste0(namex, "_", x)
    namex
  }) %>% unlist(.)

  motifNames2 <- lapply(seq_along(motifs), function(x){
    namex <- make.names(motifs[[x]]@name)
    if(grepl("LINE", namex)){
      splitNamex <- stringr::str_split(motifs[[x]]@ID, pattern="\\_", simplify = TRUE)
      namex <- splitNamex[1, grep("LINE",splitNamex[1,]) + 1]
    }
    if(substr(namex,nchar(namex),nchar(namex))=="."){
      namex <- substr(namex,1,nchar(namex)-1)
    }
    namex
  }) %>% unlist(.)

  motifDF <- lapply(seq_along(motifs), function(x){
    df <- data.frame(
      row.names = motifNames[x],
      name = motifNames2[[x]],
      ID = motifs[[x]]@ID,
      strand = motifs[[x]]@strand,
      stringsAsFactors = FALSE
    )
  }) %>% Reduce("rbind", .) %>% DataFrame

  names(motifs) <- motifNames

  out <- list(motifs = motifs, motifSummary = motifDF)

  return(out)

}


.requirePackage <- function(x = NULL, load = TRUE, installInfo = NULL, source = NULL){
  if(x %in% rownames(installed.packages())){
    if(load){
      suppressPackageStartupMessages(require(x, character.only = TRUE))
    }else{
      return(0)
    }
  }else{
    if(!is.null(source) & is.null(installInfo)){
      if(tolower(source) == "cran"){
        installInfo <- paste0('install.packages("',x,'")')
      }else if(tolower(source) == "bioc"){
        installInfo <- paste0('BiocManager::install("',x,'")')
      }else{
        stop("Unrecognized package source, available are cran/bioc!")
      }
    }
    if(!is.null(installInfo)){
      stop(paste0("Required package : ", x, " is not installed/found!\n  Package Can Be Installed : ", installInfo))
    }else{
      stop(paste0("Required package : ", x, " is not installed/found!"))
    }
  }
}

.validGRanges <- function(gr = NULL){
  stopifnot(!is.null(gr))
  if(inherits(gr, "GRanges")){
    return(gr)
  }else{
    stop("Error cannot validate genomic range!")
  }
}



#' @param regions A named `list` of `GRanges` that are to be overlapped with the `peakSet` in the `ArchRProject`.
load_annotation_as_granges_list <- function(regions) {

  if(inherits(regions, "GRanges")) {

    regionPositions <- GRangesList(region = regions)

  } else {

    if(is.null(names(regions))){
      names(regions) <- paste0("Region_", seq_along(regions))
    }

    regionPositions <- lapply(seq_along(regions), function(x) {
      
      if (inherits(regions[[x]], "GRanges")){

          gr <- .validGRanges(regions[[x]])
          
      } else if(is.character(regions[[x]])){

        gr <- tryCatch({
          makeGRangesFromDataFrame(
            df = data.frame(data.table::fread(regions[[x]])), 
            keep.extra.columns = TRUE,
            seqnames.field = "V1",
            start.field = "V2",
            end.field = "V3",
            strand.field="V4"
          )
        }, error = function(y){ print(paste0("Could not successfully get region : ", regions[[x]])) })

      } else{
        print("Unrecognized input in regions please input GRanges, GRangesList, or Paths to bed files!")
      }

      gr

    }) %>% GRangesList

    names(regionPositions) <- names(regions)

    return(regionPositions)
  }
}



#' Get/Validate BSgenome
#' 
#' This function will attempt to get or validate an input as a BSgenome.
#' 
#' @param genome This option must be one of the following: (i) the name of a valid ArchR-supported genome ("hg38", "hg19", or "mm10"),
#' (ii) the name of a `BSgenome` package (for ex. "BSgenome.Hsapiens.UCSC.hg19"), or (iii) a `BSgenome` object.
#' @param masked A boolean describing whether or not to access the masked version of the selected genome. See `BSgenome::getBSgenome()`.
#' @export
validBSgenome <- function(genome = NULL, masked = FALSE){
  
  stopifnot(!is.null(genome))
  if(inherits(genome, "BSgenome")){
    return(genome)
  }else if(is.character(genome)){
    genome <- tryCatch({
      .requirePackage(genome)
      bsg <- eval(parse(text = genome))
      if(inherits(bsg, "BSgenome")){
        return(bsg)
      }else{
        stop("genome is not a BSgenome valid class!")
      }
    }, error = function(x){
      BSgenome::getBSgenome(genome, masked = masked)
    })  
    return(genome)
  }else{
    stop("Cannot validate BSgenome options are a valid BSgenome or character for getBSgenome")
  }  
}
