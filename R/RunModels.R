#' Run Latent Dirichlet Allocation in a cisTopic object
#'
#' Run Latent Dirichlet Allocation in a given cisTopic object.
#' @param object Initialized cisTopic object.
#' @param topic Integer or vector of integers indicating the number of topics in the model/s (by default it is 25). We recommend to try several values
#' if possible, and select the best model based on the highest likelihood.
#' @param nCores Number of cores to use. By default it is 1, but if several models with distinct number of topics are being tested; it is
#' recommended to increase it to the number of models tested (or capacity of the machine). Parellelization is done with snow.
#' @param seed Seed for the assignment initialization for making results reproducible.
#' @param iterations Number of iterations over the data set. By default, 500 iterations are taken. However, we advise to use logLikelihoodByIter to check whether the
#' log likelihood of the model is stabilized with this parameters.
#' @param burnin Number of iterations to discard from the assingment counting. By default, 250 iterations are discarded. This number has to be
#' lower than the number of iterations.
#' @param alpha Scalar value indicating the (symmetric) Dirichlet hyperparameter for topic proportions. By default, it is set to 50.
#' @param alphaByTopic Logical indicating whether the scalar given in alpha has to be divided by the number of topics. By default, it is set to true.
#' @param beta Scalar value indicating the (symmetric) Dirichlet hyperparameter for topic multinomilas. By default, it is set to 0.1.
#' @param returnType Defines what has to be returned to the cisTopic object: either 'allModels' or 'selectedModel'. 'allModels' will return
#' a list with all the fitted models (as lists) to object@@models, while 'selectedModel' will return the model with the best log likelihood to
#' object@@selected.model, and a dataframe with the log likelihood of the other models to object@@log.lik. By default, this function will return all models for allowing posterior selection; however, note that if the number
#' of models and the size of the data is considerably big, returning all models may be memory expensive.
#' @param ... See \code{lda.collapsed.gibbs.sampler} from the package lda.
#'
#' @return Returns a cisTopic object with the models stored in object@@models.
#' If specified, only the best model based on log likelihood is returned in object@@selected.model, and the rest of log likelihood values are stored in object@@log.lik.
#'
#' @details The selected parameters are adapted from Griffiths & Steyvers (2004).
#' @importFrom lda lda.collapsed.gibbs.sampler
#' @importFrom parallel makeCluster
#' @import doSNOW
#' @importFrom plyr llply
#' @export
#'
#' @examples
#' bamfiles <- c('example_1.bam', 'example_2.bam', 'example_3.bam')
#' regions <- 'example.bed'
#' cisTopicObject <- CreatecisTopicObjectfromBAM(bamfiles, regions)
#' cisTopicObject <- RunModel(cisTopicObject)

runModels <- function(
  object,
  topic=25,
  nCores=1,
  seed=123,
  iterations = 500,
  burnin = 250,
  alpha = 50,
  alphaByTopic = TRUE,
  beta=0.1,
  returnType='allModels',
  ...
) {
  # Take binary count matrix
  object.binary.count.matrix <- object@binary.count.matrix

  # Prepare data
  cellList <- lapply(seq_len(ncol(object.binary.count.matrix)), function(i) rbind(as.integer(seq(0,nrow(object.binary.count.matrix)-1,1)), as.integer(object.binary.count.matrix[,i])))
  names(cellList) <- colnames(object.binary.count.matrix)
  cellList_nozeros <- lapply(cellList, function(x) {x[,x[2,]!=0]})
  cellList <- lapply(cellList_nozeros, function(x) {colnames(x) <- rownames(object.binary.count.matrix)[x[1,]+1];x})
  regionList <- rownames(object.binary.count.matrix)

  if (length(topic) > 1){
    if (length(topic) < nCores){
      print(paste('The number of cores (', nCores, ') is higher than the number of topics (', topic,').', sep=''))
    }

    if (nCores > 1){
      # Run models with SNOW
      cl <- makeCluster(nCores, type = "SOCK")
      registerDoSNOW(cl)
      clusterEvalQ(cl, library(lda))
      clusterExport(cl, c("cellList", "regionList", "topic", "iterations", "burnin", "alpha", "beta"), envir=environment())
      opts <- list(preschedule=TRUE)
      clusterSetRNGStream(cl, seed)
      if (alphaByTopic==TRUE){
        models <- suppressWarnings(llply(.data=topic, .fun=function(t) lda.collapsed.gibbs.sampler(cellList, t, regionList, num.iterations=iterations, alpha=alpha/t, eta=beta, compute.log.likelihood = TRUE, burnin=burnin, ...) , .parallel = TRUE, .paropts = list(.options.snow=opts), .inform=FALSE))
      }
      else{
        models <- suppressWarnings(llply(.data=topic, .fun=function(t) lda.collapsed.gibbs.sampler(cellList, t, regionList, num.iterations=iterations, alpha=alpha, eta=beta, compute.log.likelihood = TRUE, burnin=burnin, ...) , .parallel = TRUE, .paropts = list(.options.snow=opts), .inform=FALSE))
      }
    }
    else{
      if (alphaByTopic==TRUE){
        models <- suppressWarnings(llply(.data=topic, .fun=function(t) lda.collapsed.gibbs.sampler(cellList, t, regionList, num.iterations=iterations, alpha=alpha/t, eta=beta, compute.log.likelihood = TRUE, burnin=burnin, ...), .progress = progress_text(char = ".")))
      }
      else{
        models <- suppressWarnings(llply(.data=topic, .fun=function(t) lda.collapsed.gibbs.sampler(cellList, t, regionList, num.iterations=iterations, alpha=alpha, eta=beta, compute.log.likelihood = TRUE, burnin=burnin, ...), .progress = progress_text(char = ".")))
      }
    }
    stopCluster(cl)
  }

  else{
    if (alphaByTopic==TRUE){
      models <- llply(lda.collapsed.gibbs.sampler(cellList, topic, regionList, num.iterations=iterations, alpha=alpha/topic, eta=beta, compute.log.likelihood = TRUE, burnin=burnin, ...), .progress = progress_text(char = "."))
    }
    else{
      models <- llply(lda.collapsed.gibbs.sampler(cellList, topic, regionList, num.iterations=iterations, alpha=alpha, eta=beta, compute.log.likelihood = TRUE, burnin=burnin, ...), .progress = progress_text(char = "."))
    }
  }

  if(returnType=='allModels'){
    object@models <- models
  }
  if(returnType=='selectedModel') {
    object@selected.model <- selectModel(models)
  }

  object@calc.params[['runModels']] <- c(as.list(environment(), all = TRUE)[names(formals("runModels"))[-1]], list(...))

  return(object)
}

#' Model selection based on log likelihood
#'
#' Plots log likelihood of the different models and selects the best one based on the maximum likelihood (or specified by the user).
#' @param object Initialized cisTopic object or list of LDA models.
#' @param select Number of topics of the selected model. If NULL, the model with the best log likelihood is picked.
#' @param keepBinaryMatrix Whether to keep the binary accessibility matrix within the cisTopic object.
#' @param keepModels Whether to keep all the models within the cisTopic object.
#' @param ... Ignored.
#'
#' @return Returns a cisTopic object (when the input is a cisTopic object) with the selected model stored in object@@selected.model, and the log likelihoods of the models in object@@log.lik.
#' Otherwise, it returns a list with the selected models.
#'
#' @export
#'
#' @examples
#' bamfiles <- c('example_1.bam', 'example_2.bam', 'example_3.bam')
#' regions <- 'example.bed'
#' cisTopicObject <- CreatecisTopicObjectfromBAM(bamfiles, regions)
#' cisTopicObject <- RunModel(cisTopicObject)
#' cisTopicObject <- selectModel(cisTopicObject)

selectModel <- function(
  object,
  select=NULL,
  keepBinaryMatrix=TRUE,
  keepModels=TRUE,
  ...
){
  # Checking input
  if (as.vector(class(object)) == 'cisTopic'){
    models <- object@models
  }
  else if (is.list(object)){
    models <- object
  }
  else{
    stop('Incorrect input. Is it a cisTopic object or a list of models?')
  }

  # Make data frame
  loglikelihood <- sapply(seq_along(models), FUN=function(i) models[[i]]$log.likelihood[2,ncol(models[[i]]$log.likelihood)])
  topics <-  sapply(seq_along(models), FUN=function(i) nrow(models[[i]]$topics))
  object.log.lik <- data.frame(topics=topics, LL=loglikelihood)
  object.log.lik <- object.log.lik[order(object.log.lik$topics),]

  par(bty = 'n')
  plot(object.log.lik$topics, object.log.lik$LL, xlab="Number of topics", ylab="log P(D|M,T)", type='o', pch=16, col='black', main='Model selection')

  if (is.null(select)){
    points(object.log.lik$topics[which(object.log.lik$LL == max(object.log.lik$LL))], max(object.log.lik$LL), pch=4, col='red', lwd = 7)
    title(sub=paste("Best model:", object.log.lik$topics[which(object.log.lik$LL == max(object.log.lik$LL))], 'topics'))
    selected.model <- models[[which(object.log.lik$LL == max(object.log.lik$LL))]]
  }
  else{
    points(object.log.lik$topics[which(object.log.lik$topics == select)], object.log.lik$LL[which(object.log.lik$topics == select)], pch=4, col='red', lwd = 7)
    title(sub=paste("Selected model:", object.log.lik$topics[which(object.log.lik$topics == select)], 'topics'))
    selected.model <- models[[which(object.log.lik$topics == select)]]
  }

  if (as.vector(class(object)) == 'cisTopic'){
    object@log.lik <- object.log.lik
    object@selected.model <- selected.model
    if (keepBinaryMatrix != TRUE){
      object@binary.count.matrix <- NULL
    }
    if (keepModels != TRUE){
      object@models <- NULL
    }

    return(object)
  }
  else if (is.list(object)){
    return(selected.model)
  }
}

#' Plots log likelihood in the models per iteration.
#'
#' Plots the log likelihood of the different models through iterations. Log likelihood after burnin should be stable, otherwise more iterations
#' may be required (together with an increase in the burnin).
#' @param object Initialized cisTopic object.
#' @param ... Ignored.
#'
#' @return Plots the log likelihood of the models thorugh the different iterations after burnin.
#'
#' @export
#'
#' @examples
#' bamfiles <- c('example_1.bam', 'example_2.bam', 'example_3.bam')
#' regions <- 'example.bed'
#' cisTopicObject <- CreatecisTopicObjectfromBAM(bamfiles, regions)
#' cisTopicObject <- RunModel(cisTopicObject)
#' cisTopicObject <- logLikelihoodByIter(cisTopicObject)

logLikelihoodByIter <- function(
  object,
  ...
){
  models <- object@models
  burnin <- object@calc.params[['runModels']][['burnin']]
  iterations <- object@calc.params[['runModels']][['iterations']]

  loglikelihood_iterations <- sapply(seq_along(models), FUN=function(i) models[[i]]$log.likelihood[2,])
  topics <-  sapply(seq_along(models), FUN=function(i) nrow(models[[i]]$topics))
  colnames(loglikelihood_iterations) <- paste(topics, 'topics')

  col <- rainbow(length(topics),s = 0.5)
  par(bty = 'n')
  matplot(1:iterations,loglikelihood_iterations,  type = 'l', lty=1, lwd=4, col = col, xlab="Iteration number", ylab="log P(D|M,T)", main='Likelihood stabilization')
  legend("bottomright", legend = colnames(loglikelihood_iterations), fill=col)
}


#' Probability of each region in each cell in the data set.
#'
#' Calculates the probability of each region in each cell.
#'
#' @param object Initialized cisTopic object.
#' @param ... Ignored.
#'
#' @return Returns a matrix where the rows are the regions, the columns the cells, and the values are the probabilities of seeing a
#' region in a matrix.
#'
#' @importFrom lda predictive.distribution
#' @export
#'
#' @examples
#' bamfiles <- c('example_1.bam', 'example_2.bam', 'example_3.bam')
#' regions <- 'example.bed'
#' cisTopicObject <- CreatecisTopicObjectfromBAM(bamfiles, regions)
#' cisTopicObject <- RunModel(cisTopicObject)
#' cisTopicObject <- selectModel(cisTopicObject)
#' cisTopicObject <- predictiveDistributioncisTopicObject)

predictiveDistribution <- function(
  object,
  ...
){
  document_expects <- object@selected.model$document_expects
  topics <- object@selected.model$topics
  alpha <- object@calc.params[['runModels']]$alpha
  beta <- object@calc.params[['runModels']]$beta
  object@selected.model[['predictive.distribution']] <- predictive.distribution(document_expects, topics, alpha, beta)
  return(object)
}
