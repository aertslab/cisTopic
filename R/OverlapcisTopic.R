#' getSignaturesRegions
#'
#' Get peaks found in region signatures given as a bed files and/or dataframes
#' @param object Initialized cisTopic object.
#' @param signature Path to bed file or data frame containing signature regions (or a vector)
#' @param labels A character vector with the names to be given to the signatures (in the same order as the signature paths or named)
#' @param minOverlap Minimum overlap between the regions to consider them as overlapping regions
#' @param ... See findOverlaps from GenomicRanges
#'
#' @return The regions from the dataset that overlap with the signature are returned as a slot to object@@signatures
#'
#' @examples
#' regions <- c('regions_1.bed', 'regions_2.bed', 'regions_3.bed')
#' cisTopicObject <- getSignaturesRegions(cisTopicObject, regions)
#' cisTopicObject
#' @importFrom GenomicRanges findOverlaps
#' @importFrom GenomicRanges pintersect
#' @export

getSignaturesRegions <- function(
  object,
  signatures,
  labels,
  minOverlap=0.4,
  ...
){
  regionsSignatures <- llply(signatures, function(signature) .getSignatureRegions(object, signature, ...))
  if (is.null(regionsSignatures)){
    names(regionsSignatures) <- as.character(signatures)
  }
  else{
    if (!is.null(names(labels))){
      labels <- as.vector(labels[as.character(signatures)])
    }
    names(regionsSignatures) <- labels
  }
  
  object@signatures <- regionsSignatures
  return(object)
}


# Helper functions

.getSignatureRegions <- function(
  object,
  signature,
  minOverlap=0.4,
  ...
){
  if (file.exists(signature)){
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
#' For this, the function AnnotateRegions must be run first.
#' @param nCores Number of cores to be used for AUCell
#' @param aucMaxRank Threshold to calculate the AUC
#' @param col.low Color to use for lowest signature enrichment
#' @param col.mid Color to use for medium signature enrichment
#' @param col.high Color to use for high signature enrichment
#' @param ... See \code{AUCell_calcAUC} from AUCell
#'
#' @return Heatmap showing the enrichment per topic per signature
#'
#' @import AUCell
#' @importFrom NMF aheatmap
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
  ...){

  # Get scores
  scores <- .getScores(object)
  if (selected.signatures != 'annotation'){
    signatures <- object@signatures
    if(selected.signatures != 'all'){
      signatures <- signatures[selected.signatures]
    }
  }
  else{
    signatures <- split(object@region.data, object@region.data$annotation)
    signatures <- lapply(signatures, function(x) rownames(x))
  }

  if(topics[1] != 'all'){
    scores <- scores[,topics]
  }

  aucellRankings <- AUCell_buildRankings(as.matrix(scores), nCores=nCores, plotStats=FALSE, verbose = FALSE)
  modulesAUC <- AUCell_calcAUC(signatures, aucellRankings, nCores=nCores, aucMaxRank=aucMaxRank)
  enrichMatrix <- getAUC(modulesAUC)
  enrichMatrix <- t(scale(t(enrichMatrix)))

  cl.topics <- hclust.vector(t(enrichMatrix), method="ward", metric="euclidean")
  dd.col <- as.dendrogram(cl.topics)

  colorPal <- grDevices::colorRampPalette(c("dodgerblue", "floralwhite", "brown1"))

  nmf.options(grid.patch=TRUE)
  NMF::aheatmap(enrichMatrix, scale="none", revC=TRUE, main='Signatures heatmap', sub='Row normalized AUC scores', Colv=dd.col, color = colorPal(20))
}

