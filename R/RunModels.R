#' Run Latent Dirichlet Allocation in a cisTopic object
#'
#' Run Latent Dirichlet Allocation in a given cisTopic object.
#' @param object Initialized cisTopic object.
#' @param topic Integer or vector of integers indicating the number of topics in the model/s (by default it is a vector with 2, 10, 20, 30, 40 and 50 topics). We recommend to try several values
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
#' @param addModels Whether models should be added if there is a pre-existing list of models or should be overwritten by new models. If TRUE, parameters are setted to match the existing models.
#' @param ... See \code{lda.collapsed.gibbs.sampler} from the package lda.
#'
#' @return Returns a cisTopic object with the models stored in object@@models.
#' If specified, only the best model based on log likelihood is returned in object@@selected.model, and the rest of log likelihood values are stored in object@@log.lik.
#'
#' @details The selected parameters are adapted from Griffiths & Steyvers (2004).
#' @importFrom lda lda.collapsed.gibbs.sampler
#' @import parallel
#' @import doSNOW
#' @import Matrix
#' @importFrom plyr llply
#' @export
#'
#' @examples
#' bamfiles <- c('example_1.bam', 'example_2.bam', 'example_3.bam')
#' regions <- 'example.bed'
#' cisTopicObject <- createcisTopicObjectfromBAM(bamfiles, regions)
#' cisTopicObject <- runModels(cisTopicObject)

runModels <- function(
  object,
  topic=c(2, 10, 20, 30, 40, 50),
  nCores=1,
  seed=123,
  iterations = 500,
  burnin = 250,
  alpha = 50,
  alphaByTopic = TRUE,
  beta=0.1,
  returnType='allModels',
  addModels = TRUE,
  ...
) {
  if (burnin >= iterations){
    stop('The number of iterations must be higher than the burnin!')
  } 
  
  if (!is.null(object@calc.params[['runModels']])){
    if (!addModels){
      print('Resetting models...')
      object@calc.params[['runModels']] <- NULL
      object@models <- list()
      object@calc.params[['runModels']] <- c(as.list(environment(), all = TRUE)[names(formals("runModels"))[-1]], list(...))
    }
    else {
      seed <- object@calc.params[['runModels']]$seed
      iterations <- object@calc.params[['runModels']]$iterations
      burnin <- object@calc.params[['runModels']]$burnin
      alpha <- object@calc.params[['runModels']]$alpha
      alphaByTopic <- object@calc.params[['runModels']]$alphaByTopic
      beta <- object@calc.params[['runModels']]$beta
      cat(paste('In order to compare these models with previous models, these parameters will take the values from the previous run: \n seed:', seed, '\n iterations:', iterations, '\n burnin:', burnin, '\n alpha:', alpha, '\n alphaByTopic:', alphaByTopic, '\n beta:', beta))
    }
  }
  else{
    object@calc.params[['runModels']] <- c(as.list(environment(), all = TRUE)[names(formals("runModels"))[-1]], list(...))
  }
  
  # Take binary count matrix
  object.binary.count.matrix <- object@binary.count.matrix
  cellnames <- colnames(object.binary.count.matrix)
  regionnames <- rownames(object.binary.count.matrix)
  
  # Prepare data
  print('Formatting data...')
  cellList <- split(as.integer(object.binary.count.matrix@i), rep(seq_along(object.binary.count.matrix@p+1), times=diff(c(object.binary.count.matrix@p+1, length(object.binary.count.matrix@i) + 1))))
  rm(object.binary.count.matrix)
  cellList <- lapply(cellList, function(x) {x <- rbind(x, rep(as.integer(1), length(x)))}) 
  names(cellList) <- cellnames
  cellList <- lapply(cellList, function(x) {colnames(x) <- regionnames[x[1,]+1];x})
  regionList <- regionnames
  
  if (length(topic) > 1){
    if (length(topic) < nCores){
      print(paste('The number of cores (', nCores, ') is higher than the number of models (', length(topic),').', sep=''))
    }
    
    if (nCores > 1){
      # Run models with SNOW
      print('Exporting data...')
      cl <- makeCluster(nCores, type = "SOCK")
      registerDoSNOW(cl)
      clusterEvalQ(cl, library(lda))
      clusterExport(cl, c("cellList", "regionList", "topic", "iterations", "burnin", "alpha", "beta"), envir=environment())
      opts <- list(preschedule=TRUE)
      clusterSetRNGStream(cl, seed)
      if (alphaByTopic==TRUE){
        print('Running models...')
        models <- suppressWarnings(llply(.data=topic, .fun=function(t) lda.collapsed.gibbs.sampler(cellList, t, regionList, num.iterations=iterations, alpha=alpha/t, eta=beta, compute.log.likelihood = TRUE, burnin=burnin, ...)[-1] , .parallel = TRUE, .paropts = list(.options.snow=opts), .inform=FALSE))
      }
      else{
        print('Running models...')
        models <- suppressWarnings(llply(.data=topic, .fun=function(t) lda.collapsed.gibbs.sampler(cellList, t, regionList, num.iterations=iterations, alpha=alpha, eta=beta, compute.log.likelihood = TRUE, burnin=burnin, ...)[-1] , .parallel = TRUE, .paropts = list(.options.snow=opts), .inform=FALSE))
      }
      stopCluster(cl)
    }
    else{
      if (alphaByTopic==TRUE){
        print('Running models...')
        models <- suppressWarnings(llply(.data=topic, .fun=function(t) lda.collapsed.gibbs.sampler(cellList, t, regionList, num.iterations=iterations, alpha=alpha/t, eta=beta, compute.log.likelihood = TRUE, burnin=burnin, ...)[-1], .progress = progress_text(char = ".")))
      }
      else{
        print('Running models...')
        models <- suppressWarnings(llply(.data=topic, .fun=function(t) lda.collapsed.gibbs.sampler(cellList, t, regionList, num.iterations=iterations, alpha=alpha, eta=beta, compute.log.likelihood = TRUE, burnin=burnin, ...)[-1], .progress = progress_text(char = ".")))
      }
    }

  }
  
  else{
    set.seed(seed)
    if (alphaByTopic==TRUE){
      print('Running models...')
      models <- llply(lda.collapsed.gibbs.sampler(cellList, topic, regionList, num.iterations=iterations, alpha=alpha/topic, eta=beta, compute.log.likelihood = TRUE, burnin=burnin, ...)[-1], .progress = progress_text(char = "."))
    }
    else{
      print('Running models...')
      models <- llply(lda.collapsed.gibbs.sampler(cellList, topic, regionList, num.iterations=iterations, alpha=alpha, eta=beta, compute.log.likelihood = TRUE, burnin=burnin, ...)[-1], .progress = progress_text(char = "."))
    }
  }
  
  if (!is.null(object@models)){
    if(length(topic) == 1){
      models <- .addModels(c(object@models, list(models))) 
    } else {
      models <- .addModels(c(object@models, models))
    }
  } else {
    names(models) <- laply(1:length(models), function(x) sapply(models[x], function(y) nrow(y$topic_sums)))
  }
  
  if(returnType=='allModels'){
    object@models <- models
  }
  if(returnType=='selectedModel') {
    object@selected.model <- selectModel(models)
  }
  
  return(object)
}

# Helper function

.addModels <- function(
  modelList
){
  names(modelList) <- laply(1:length(modelList), function(x) sapply(modelList[x], function(y) nrow(y$topic_sums)))
  modelList <- modelList[as.character(sort(as.numeric(names(modelList))))]
  modelList <- modelList[unique(names(modelList))]
  return(modelList)
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
#' The unnormalized cell assignments throughtout the sampling iterations are stored in \code{cisTopicObject@@selected.model$document_expects}; 
#' while the corresposding unnormalized region assignments are stored in \code{cisTopicObject@@selected.model$topics}.
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
    if (length(models) < 1){
      stop('Please, run runModels() first.')
    }
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
      object@models <- list()
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
#' @param select Plot log-likelihood per iteration for selected topics (as a vector with the number of topics per model).
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
  select=NULL,
  ...
){
  models <- object@models
  
  if (length(models) < 1){
    stop('Please, run runModels() first.')
  }
  
  burnin <- object@calc.params[['runModels']][['burnin']]
  iterations <- object@calc.params[['runModels']][['iterations']]

  loglikelihood_iterations <- sapply(seq_along(models), FUN=function(i) models[[i]]$log.likelihood[2,])
  topics <-  sapply(seq_along(models), FUN=function(i) nrow(models[[i]]$topics))
  colnames(loglikelihood_iterations) <- paste(topics, 'topics')
  
  if (!is.null(select)){
    loglikelihood_iterations <- loglikelihood_iterations[,paste0(select, ' topics')]
  }

  col <- .distinctColorPalette(ncol(loglikelihood_iterations))
  par(bty = 'n')
  matplot(1:iterations,loglikelihood_iterations,  type = 'l', lty=1, lwd=4, col = col, xlab="Iteration number", ylab="log P(D|M,T)", main='Likelihood stabilization')
  abline(v = burnin, lty=2, col='grey')
  legend("bottomright", legend = colnames(loglikelihood_iterations), fill=col)
}


#' Probability of each region in each cell in the data set.
#'
#' Calculates the probability of each region in each cell.
#'
#' @param object Initialized cisTopic object.
#' @param big.matrix If having big data, we recommend to use the bigmemory package for the calculations.
#' @param ... Ignored.
#'
#' @return Returns a matrix where the rows are the regions, the columns the cells, and the values are the probabilities of seeing a
#' region in a matrix.
#'
#' @import Matrix
#' @export
#'
#' @examples
#' bamfiles <- c('example_1.bam', 'example_2.bam', 'example_3.bam')
#' regions <- 'example.bed'
#' cisTopicObject <- createcisTopicObjectfromBAM(bamfiles, regions)
#' cisTopicObject <- runModels(cisTopicObject)
#' cisTopicObject <- selectModel(cisTopicObject)
#' cisTopicObject <- predictiveDistribution(cisTopicObject)

predictiveDistribution <- function(
  object,
  big.matrix=FALSE,
  ...
){
  
  if (length(object@selected.model) < 1){
    stop('Please, run selectModel() first.')
  }
  
  document_expects <- object@selected.model$document_expects
  topics <- object@selected.model$topics
  alpha <- object@calc.params[['runModels']]$alpha/length(object@selected.model$topic_sums)
  beta <- object@calc.params[['runModels']]$beta
  pred.matrix <- as.matrix(.predictive.distribution(document_expects, topics, alpha, beta, big.matrix))
  return(pred.matrix)
}


# Helper function
.predictive.distribution <- function(
  document_expects,
  topics,
  alpha,
  beta,
  big.matrix
){
  if (big.matrix){
    if(! "bigmemory" %in% installed.packages()){
      stop('Please, install bigmemory: \n install.packages("bigmemory")')
    } else {
      require(bigmemory)
    }
    smoothed.topics <- as.big.matrix(t((topics + beta)/Matrix::rowSums(topics + beta)))
    props_mat <- as.big.matrix(apply(document_expects, 2, function(x) {(x + alpha)/sum(x + alpha)}))
  } else {
    smoothed.topics <- t((topics + beta)/Matrix::rowSums(topics + beta))
    props_mat <- apply(document_expects, 2, function(x) {(x + alpha)/sum(x + alpha)})
  }
  mat <- smoothed.topics %*% props_mat 
  return(mat)
}