#' Get region scores
#'
#' Get region scores per topic
#' @param object Initialized cisTopic object, after the object@@selected.model has been filled.
#' @param method Select the normalization method to use for calculating the region score: 'Zscore' or 'probability'.
#' @param rescale Whether feature scaling should be performed within the topics. By default, it is set to true.
#'
#' @return Returns the region scores per topic in object@@region.data
#'
#' @export

getRegionsScores <- function(
  object,
  method='Zscore',
  scaled=TRUE)
{
  # Restart data frame
  indexes <- grep('Scores_Topic', colnames(object@region.data))
  if (length(indexes) > 0){
    object@region.data <- object@region.data[,-indexes]
  }
  
  # Get scores
  topic.mat <- object@selected.model$topics
  if (method == 'Zscore'){
    normalizedTopics <- topic.mat/(rowSums(topic.mat) + 1e-05)
    scores <- apply(normalizedTopics, 2, function(x) x * (log(x + 1e-05) - sum(log(x + 1e-05))/length(x)))
  }
  else {
    beta <- object@calc.params[['RunModels']]$beta
    scores <- (topic.mat + beta)/rowSums(topic.mat + beta)
  }

  # Scale
  if (scaled == TRUE) {
    scores <- apply(scores, 1, function(x) (x-min(x))/(max(x)-min(x)))
  }
  else {
    scores <- t(scores)
  }

  colnames(scores) <- paste("Scores_Topic", 1:ncol(scores), sep='')
  object@region.data[, colnames(scores)] <- scores

  return(object)
}

#' Binarize cisTopic
#'
#' cisTopic binarization based on cisTopic scores
#' @param object Initialized cisTopic object, after the object@@region.data has been filled with the cisTopic scores.
#' @param method Method for the binarization, 'GammaFit' or 'Predefined'. If GammaFit, a gamma distribution will be fitted
#' to the region scores distributions (a probability threshold must be provided); if threshold, a vector with the top number
#' of regions to be taken per topic must be provided.
#' @param thrP Probability threshold to use as cutoff on the probability distribution when using GammaFit as method.
#' @param plot If TRUE, plots of the gamma fit per topic will be produced and a barplot with the number of regions taken
#' per topic.
#' @param cutoffs Vector (or scalar) with the top number of regions per topic to take per topic if Predefined is selected
#' as method. If scalar, the same number of regions is taken per topic.
#'
#' @return A list with the selected regions per topic stored in object@@binarized.cisTopics
#'
#' @importFrom fitdistrplus fitdist
#' @importFrom plyr llply
#' @export


binarizecisTopics  <- function(
  object,
  method = 'GammaFit',
  thrP = 0.99,
  plot = TRUE,
  cutoffs = NULL
){
  # Get scores
  scores <- .getScores(object)
  object.binarized.cisTopics <- list()

  # Check input
  if (method == 'GammaFit'){
    if (is.null(thrP)){
      stop('If GammaFit is selected as method, a probability threshold for cutoff in the distribution must be provided.')
    }

    for (i in 1:ncol(scores)){
      distr <- fitdist(scores[,i],"gamma", method="mme")
      cutoff <- as.numeric(unlist(quantile(distr, probs = thrP))[1])
      if (plot == TRUE){
        par(mfrow=c(1,1))
        par(las=1)
        hist(scores[,i], breaks=100, prob=TRUE, main=colnames(scores)[i], col=adjustcolor( "dodgerblue", alpha.f = 0.8), xlab='Score')
        curve(dgamma(x, rate = as.numeric(unlist(distr[1]$estimate[2])), shape = as.numeric(unlist(distr[1]$estimate[1]))), add=TRUE, col="magenta", lwd=2)
        abline(v=cutoff, lty=3, lwd=2, col='grey')
        title(sub=paste("Regions with score > ", signif(cutoff, 2), sep=""))
      }
      object.binarized.cisTopics[[colnames(scores)[i]]] <- scores[which(scores[,i] > cutoff), i, drop = FALSE]
    }
  }

  if (method == 'Predefined'){
    if (is.null(cutoffs)){
      stop('If Predefined is selected as method, the threshold cutoffs must be provided.')
    }
    # Correct cutoffs if scalar
    if(length(cutoffs) == 1){
      cutoffs <- rep(cutoffs, ncol(scores))
    }
    else if (length(cutoffs) != ncol(scores)){
      stop('Thresholds for all topics are not provided. Check the length of the cutoff vector.')
    }

    if (plot == TRUE){
      for(i in 1:ncol(scores)){
        par(mfrow=c(1,1))
        par(las=1)
        hist(scores[,i], breaks=100, prob=TRUE, main=colnames(scores)[i], col=adjustcolor( "dodgerblue", alpha.f = 0.8), xlab='Score')
        abline(v=cutoffs[i], lty=3, lwd=2, col='grey')
        title(sub=paste("Regions with score > ", signif(cutoffs[i], 2), sep=""))
      }
    }
    object.binarized.cisTopics <- llply(1:ncol(scores), function (i) scores[which(scores[,i] > cutoff), i, drop = FALSE])
  }

  object.binarized.cisTopics <- llply(1:ncol(scores), function (i) object.binarized.cisTopics[[i]][order(object.binarized.cisTopics[[i]][,1], decreasing = TRUE), , drop=FALSE])
  names(object.binarized.cisTopics) <- colnames(scores)
  object@binarized.cisTopics <- object.binarized.cisTopics
  par(mfrow=c(1,1))
  par(las=2)
  barplot(sapply(object.binarized.cisTopics, nrow), col=adjustcolor( "dodgerblue", alpha.f = 0.8), main='Number of regions selected per topic')

  return(object)
}

# Helper function

.getScores <- function(
  object
){
  scores <- object@region.data[, grep('Scores_Topic', colnames(object@region.data))]
  colnames(scores) <- paste('Topic', 1:ncol(scores), sep='')
  return(scores)
}

#' getBedFiles
#'
#' Saves topics as bed files
#' @param object Initialized cisTopic object, after object@@binarized.cisTopics has been filled.
#' @param path Path to save the coordinates of the cisTopic specific regions as bed files.
#' @param seqlengths Annotation from where to retrieve the seqlengths
#'
#' @export

getBedFiles <- function(
  object,
  path
){
  dir.create(path, showWarnings = FALSE)
  object.binarized.cisTopics <- object@binarized.cisTopics
  coordinates <- object@region.data[ , c('seqnames', 'start', 'end')]
  for (i in 1:length(object.binarized.cisTopics)){
    write.table(coordinates[rownames(object.binarized.cisTopics[[i]]),], file=paste(path, '/Topic_', i, '.bed', sep=''), row.names=FALSE, col.names = FALSE, quote=FALSE,  sep = "\t", eol = "\n")
  }
}

#' getBigwigFiles
#'
#' Saves topics as bigwigs files
#' @param object Initialized cisTopic object, after running getRegionScores
#' @param path Path to save bigwig files
#'
#' @importFrom rtracklayer export.bw
#' @export

getBigwigFiles <- function(
  object,
  path,
  seqlengths
){
  dir.create(path, showWarnings = FALSE)
  scores <- .getScores(object)
  granges <- object@region.ranges
  seqlengths(granges) <- seqlengths[names(seqlengths(granges))]
  lapply(1:ncol(scores), function(i) .getOneBigWig(granges, scores, path, i))
}

.getOneBigWig <- function(
  granges,
  scores,
  path,
  topic
){
  coord <- granges
  column <- paste('Topic', topic, sep='')
  elementMetadata(coord)[['score']] <- scores[,column]
  con <- paste(path, 'Topic_', topic, '.bw', sep='')
  export.bw(coord, con)
}

