#' getSignaturesRegions
#'
#' Get peaks found in region signatures given as a bed files and/or dataframes
#' @param object Initialized cisTopic object.
#' @param signature Path to bed file or data frame containing signature regions (or a vector)
#' @param labels A character vector with the names to be given to the signatures (in the same order as the signature paths or named)
#' @param minOverlap Minimum overlap between the regions to consider them as overlapping regions
#' @param ... See \code{findOverlaps} from GenomicRanges
#'
#' @return The regions from the dataset that overlap with the signature are returned as a slot to object@@signatures
#'
#' @examples
#' regions <- c('regions_1.bed', 'regions_2.bed', 'regions_3.bed')
#' cisTopicObject <- getSignaturesRegions(cisTopicObject, regions)
#' cisTopicObject
#' @import GenomicRanges
#' @import S4vectors
#' @export

getSignaturesRegions <- function(
  object,
  signatures,
  labels,
  minOverlap=0.4,
  ...
){
  if (length(object@signatures) > 0){
    if(sum(labels %in% names(object@signatures) > 0)){
      error <- labels[which(labels %in% names(object@signatures))]
      error <- paste(error, collapse=', ')
      stop('There is at least a signature with the same label:', error ,'. Please, rename it.')
    }
  }
  
  regionsSignatures <- llply(signatures, function(signature) .getSignatureRegions(object, signature, ...))
  
  names(regionsSignatures) <- labels
  
  if (length(object@signatures) < 1){
    object@signatures <- regionsSignatures
  } else {
    object@signatures <- c(object@signatures, regionsSignatures)
  }
  
  n.ids <- sapply(regionsSignatures,length)
  regions <-  object@region.names
  idx <- lapply(1:length(regionsSignatures), function (i) match(regionsSignatures[[i]], regions))
  vals <- unlist(idx)
  mat <- as.matrix(sparseMatrix(vals, rep(seq_along(n.ids), n.ids)))
  mat <- apply(mat, 2, function(x) as.factor(x))
  colnames(mat) <- names(regionsSignatures)
  if(nrow(mat) < length(regions)){
    number.rows <- length(regions) - nrow(mat)
    add <- matrix(FALSE, number.rows, length(regionsSignatures))
    mat <- rbind(mat, add)
  }
  rownames(mat) <- regions
  object <- addRegionMetadata(object, as.data.frame(mat))
  return(object)
}


# Helper functions

.getSignatureRegions <- function(
  object,
  signature,
  minOverlap=0.4,
  ...
){
  if (!is.data.frame(signature)){
    regions <- read.table(signature)
    colnames(regions)[1:3] <- c('seqnames', 'start', 'end')
    regions <- makeGRangesFromDataFrame(regions)
  }
  else if (is.data.frame(signature)){
    regions <- makeGRangesFromDataFrame(signature)
  }
  else{
    stop('The signature is not an existing file or a dataframe.')
  }

  coordinates <- object@region.ranges
  regionsSignature <- unique(as.vector(.getOverlapRegionsFromCoordinates(coordinates, regions, minOverlap=minOverlap, ...)[,1]))
  message(paste('The signature contains', length(regions), 'of which', length(regionsSignature), 'overlap the regions in the set.'))
  return(regionsSignature)
}


.getOverlapRegionsFromCoordinates <- function(
  coordinates,
  regions,
  minOverlap=0.4,
  overlapping=TRUE,
  ...)
  {
  dbRegionsOverlap <- findOverlaps(regions, coordinates, type='any', select="all", ignore.strand=TRUE, ...)

  if(minOverlap>0){
      overlaps <- pintersect(regions[queryHits(dbRegionsOverlap)], coordinates[subjectHits(dbRegionsOverlap)])
      percentOverlapCoordinates <- width(overlaps) / width(regions[queryHits(dbRegionsOverlap)])
      percentOverlapRegions <- width(overlaps) / width(coordinates[subjectHits(dbRegionsOverlap)])
      maxOverlap <- apply(cbind(percentOverlapCoordinates, percentOverlapRegions), 1, max)
      dbRegionsOverlap <- dbRegionsOverlap[maxOverlap > minOverlap]
      maxOverlap <- maxOverlap[which(maxOverlap > minOverlap)]
  }

  selectedRegions <- regions[queryHits(dbRegionsOverlap)]
  selectedRegions <- paste(as.vector(seqnames(selectedRegions)), ':', as.vector(start(selectedRegions)), '-', as.vector(end(selectedRegions)), sep='')
  selectedCoordinates <- names(coordinates[subjectHits(dbRegionsOverlap)])

  selectedMapped <- data.frame(selectedCoordinates, selectedRegions, maxOverlap, row.names=NULL)

  if (overlapping != TRUE){
    if(any(duplicated(selectedRegions))){
      selectedMapped <- selectedMapped[order(as.vector(selectedMapped$selectedRegions), -abs(selectedMapped$maxOverlap)), ]
      selectedMapped <- selectedMapped[!duplicated(as.vector(selectedMapped$selectedRegions)), ]
    }
  }
  return(selectedMapped)
}

#' signaturesHeatmap
#'
#' Heatmap containing the row normalised AUC values for the signatures in the topics.
#' @param object Initialized cisTopic object, after regions scores have been calculated with getRegionScores and object@@signatures
#' is filled.
#' @param topics By default all topics will be used, but topics can be selected based on index.
#' @param selected.signatures By default all signatures will be used, but signatures can be selected based on index or name.
#' Alternatively, 'annotation' can be selected to use as signatures the region type labels (e.g. promoter, distal intergenic, ...)
#' For this, the function annotateRegions() must be run first.
#' @param nCores Number of cores to be used for AUCell
#' @param aucMaxRank Threshold to calculate the AUC
#' @param col.low Color to use for lowest signature enrichment
#' @param col.mid Color to use for medium signature enrichment
#' @param col.high Color to use for high signature enrichment
#' @param scale Whether AUC enrichment should be normalized
#' @param ... See \code{Heatmap} from ComplexHeatmap
#'
#' @return Heatmap showing the enrichment per topic per signature
#'
#' @import AUCell
#' @export

signaturesHeatmap <- function(
  object,
  topics = 'all',
  selected.signatures = 'all',
  nCores = 4,
  aucMaxRank = 0.03*nrow(aucellRankings),
  col.low = "dodgerblue",
  col.mid = "floralwhite",
  col.high = "brown1",
  scale=TRUE,
  ...){
  
  # Check info
  if(length(cisTopicObject@signatures) < 1){
    stop('Please, run getSignaturesRegions() first.')
  }
  
  # Check dependencies
  if(! "fastcluster" %in% installed.packages()){
    stop('Please, install fastcluster: \n install.packages("fastcluster")')
  } else {
    require(fastcluster)
  }
  
  if(! "ComplexHeatmap" %in% installed.packages()){
    stop('Please, install ComplexHeatmap: source("https://bioconductor.org/biocLite.R") \nbiocLite("ComplexHeatmap")')
  } else {
    require(ComplexHeatmap)
  }

  # Get scores
  scores <- .getScores(object)
  if (selected.signatures[1] != 'annotation'){
    signatures <- object@signatures
    if (is.null(signatures)){
      stop('Please run getSignaturesRegions() first.')
    }
    if(selected.signatures[1] != 'all'){
      if (sum(selected.signatures %in% names(signatures)) != length(selected.signatures)){
        stop('Check whether the selected signatures have been stored in object@signatures.')
      }
      signatures <- signatures[selected.signatures]
    }
  }
  else{
    signatures <- split(object@region.data, object@region.data$annotation)    
    if (is.null(signatures)){
      stop('Please run annotateRegions() first.')
    }
    signatures <- lapply(signatures, function(x) rownames(x))
  }

  if(topics[1] != 'all'){
    scores <- scores[,topics]
  }

  aucellRankings <- AUCell_buildRankings(as.matrix(scores), nCores=nCores, plotStats=FALSE, verbose = FALSE)
  modulesAUC <- AUCell_calcAUC(signatures, aucellRankings, nCores=nCores, aucMaxRank=aucMaxRank, verbose=FALSE)
  enrichMatrix <- getAUC(modulesAUC)
  
  if (scale){
    enrichMatrix <- t(scale(t(enrichMatrix)))
    name_heatmap <- 'Normalised AUC score'
  } else {
    name_heatmap <- 'AUC score'
  }
  

  cl.topics <- fastcluster::hclust.vector(t(enrichMatrix), method="ward", metric="euclidean")
  dd.col <- as.dendrogram(cl.topics)

  colorPal <- grDevices::colorRampPalette(c(col.low, col.mid, col.high))

  heatmap <- ComplexHeatmap::Heatmap(data.matrix(enrichMatrix), col=colorPal(20), cluster_columns=dd.col,
                                     show_column_names=TRUE, show_row_names = TRUE, 
                                     heatmap_legend_param = list(legend_direction = "horizontal", legend_width = unit(5, "cm"), title_position='topcenter'),
                                     name = name_heatmap, column_title_gp = gpar(fontface = 'bold'), ...)
  ComplexHeatmap::draw(heatmap, heatmap_legend_side = "bottom")
}

#' signatureCellEnrichment
#'
#' Determine signatures enrichment in the cells. If specified, signature enrichment can be plotted.
#' @param object Initialized cisTopic object, after object@@cistromes.regions have been filled (see \code{getCistromes()}).
#' @param selected.signatures By default all signatures will be used, but signatures can be selected based on index or name.
#' Alternatively, 'annotation' can be selected to use as signatures the region type labels (e.g. promoter, distal intergenic, ...)
#' For this, the function annotateRegions() must be run first.
#' @param aucellRankings Precomputed aucellRankings using \code{cisTopic_buildRankings()}. These rankings are not stored in the cisTopicObject due to their size.
#' @param nCores Number of cores to be used for AUCell
#' @param aucMaxRank Threshold to calculate the AUC
#' @param plot Whether enrichment plot should be done. If yes, parameters for plotFeatures will not be ignored.
#' @param ... See \code{plotFeatures()}.
#'
#' @return AUC enrichment values for the signature are stored as a column in object@@cell.data. If specified, cells coloured by
#' their AUC enrichment values will be plotted.
#' @examples
#'
#' cisTopicObject <- getCistromeEnrichment(cisTopicObject, annotation = 'Both', nCores=1)
#' cisTopicObject
#'
#' @import AUCell
#' @export

signatureCellEnrichment <- function(
  object,
  aucellRankings,
  selected.signatures='all',
  nCores = 1,
  aucMaxRank = 0.03*nrow(aucellRankings),
  plot=TRUE,
  ...
){
  # Check info
  if(length(cisTopicObject@signatures) < 1){
    stop('Please, run getSignaturesRegions() first.')
  }
  
  if (selected.signatures != 'annotation'){
    signatures <- object@signatures
    if (is.null(signatures)){
      stop('Please run getSignaturesRegions() first.')
    }
    if(selected.signatures != 'all'){
      if (sum(selected.signatures %in% names(signatures)) != length(selected.signatures)){
        stop('Check whether the selected signatures have been stored in object@signatures.')
      }
      signatures <- signatures[selected.signatures]
    }
  }
  else{
    signatures <- split(object@region.data, object@region.data$annotation)    
    if (is.null(signatures)){
      stop('Please run annotateRegions() first.')
    }
    signatures <- lapply(signatures, function(x) rownames(x))
  }
  
  modulesAUC <- AUCell_calcAUC(signatures, aucellRankings, nCores=nCores, aucMaxRank=aucMaxRank)
  enrichMatrix <- t(getAUC(modulesAUC))
  rownames(enrichMatrix) <- object@cell.names
  object <- addCellMetadata(object, as.data.frame(enrichMatrix))
  
  if (plot){
    plotFeatures(object, target='cell', colorBy=colnames(enrichMatrix), ...)
  }
  
  return(object)
}
