#' Run tSNE on the cell-cisTopic/region-cisTopic distributions
#'
#' Run tSNE on the cell-cisTopic/region-cisTopic distributions.
#' @param object Initialized cisTopic object, after the object@@selected.model has been filled.
#' @param target Whether dimensionality reduction should be applied on cells ('cell') or regions (region). Note that for speed and clarity
#' reasons, dimesionality reduction on regions will only be done using the regions assigned to topics with high confidence 
#' (see binarizecisTopics()).
#' @param method Select the method for processing the cell assignments: 'Z-score' and 'Probability'. In the case of regions, 
#' an additional method, 'NormTop' is available (see getRegionScores()).
#' @param seed Integer for making the results reproducible.
#' @param ... See \code{Rtsne} from the package Rtsne.
#'
#' @return Returns a cisTopic object with the tSNE coordinates stored in object@@dr in object@@dr$cell$tSNE or object@@dr$region$tSNE 
#' depending on the target.
#'
#' @details 'Z-score' computes the Z-score for each topic assingment per cell/region. 'Probability' divides the topic assignments by the total number
#' of assignments in the cell/region in the last iteration plus alpha. If using 'NormTop', regions are given an score defined by: \eqn{\beta_{w, k} (\log
#' \beta_{w,k} - 1 / K \sum_{k'} \log \beta_{w,k'})}.
#'
#' 
#' @export
#' 
#' @examples
#' cisTopicObject <- runtSNE(cisTopicobject, target='cell', method='Z-score')
#' cisTopicObject

runtSNE <- function(
  object,
  target,
  method='Z-score',
  seed=123,
  ...
){
  
  # Check info
  if (length(object@selected.model) < 1){
    stop('Please, run selectModel() first.')
  }
  
  # Check dependencies
  if(! "Rtsne" %in% installed.packages()){
    stop('Please, install Rtsne: \n install.packages("Rtsne")')
  } else {
    require(Rtsne)
  }
  
  modelMat <- modelMatSelection(object, target, method)

  set.seed(seed)
  tSNE <- Rtsne::Rtsne(t(modelMat), ...)
  rownames(tSNE$Y) <- colnames(modelMat)
  colnames(tSNE$Y) <- paste0('tSNE', 1:ncol(tSNE$Y))

  object@dr[[target]][['tSNE']] <- tSNE$Y
  object@calc.params[[target]][['tSNE']] <- c(as.list(environment(), all = TRUE)[names(formals("runtSNE"))[c(2,3)]], list(...))
  return(object)
}


#' Run umap on the cell-cisTopic/region-cisTopic distributions
#'
#' Run umap on the cell-cisTopic/region-cisTopic distributions.
#' @param object Initialized cisTopic object, after the object@@selected.model has been filled.
#' @param target Whether dimensionality reduction should be applied on cells ('cell') or regions (region). Note that for speed and clarity
#' reasons, dimesionality reduction on regions will only be done using the regions assigned to topics with high confidence 
#' (see binarizecisTopics()).
#' @param method Select the method for processing the cell assignments: 'Z-score' and 'Probability'. In the case of regions, 
#' an additional method, 'NormTop' is available (see getRegionScores()).
#' @param seed Integer for making the results reproducible.
#' @param ... See \code{umap} from the package umap.
#'
#' @return Returns a cisTopic object with the umap coordinates stored in object@@dr in object@@dr$cell$umap or object@@dr$region$umap 
#' depending on the target.
#'
#' @details 'Z-score' computes the Z-score for each topic assingment per cell/region. 'Probability' divides the topic assignments by the total number
#' of assignments in the cell/region in the last iteration plus alpha. If using 'NormTop', regions are given an score defined by: \eqn{\beta_{w, k} (\log
#' \beta_{w,k} - 1 / K \sum_{k'} \log \beta_{w,k'})}.
#'
#' 
#' @export
#' 
#' @examples
#' cisTopicObject <- runtSNE(cisTopicobject, target='cell', method='Z-score')
#' cisTopicObject

runUmap <- function(
  object,
  target,
  method='Z-score',
  seed=123,
  ...
){
  
  # Check info
  if (length(object@selected.model) < 1){
    stop('Please, run selectModel() first.')
  }
  
  # Check dependencies
  if(! "umap" %in% installed.packages()){
    stop('Please, install umap: \n install.packages("umap")')
  } else {
    require(umap)
  }
  
  modelMat <- modelMatSelection(object, target, method)
  
  set.seed(seed)
  Umap <- umap::umap(t(modelMat), ...)
  rownames(Umap$layout) <- colnames(modelMat)
  colnames(Umap$layout) <- paste0('UMAP', 1:ncol(Umap$layout))
  
  object@dr[[target]][['Umap']] <- Umap$layout
  object@calc.params[[target]][['Umap']] <- c(as.list(environment(), all = TRUE)[names(formals("runUmap"))[c(2,3)]], list(...))
  return(object)
}

#' Calculate Diffusion Components on the cell-cisTopic distributions
#'
#' Calculate Diffusion Components on the cell-cisTopic distributions
#' @param object Initialized cisTopic object, after the object@@selected.model has been filled.
#' @param target Whether dimensionality reduction should be applied on cells ('cell') or regions (region). Note that for speed and clarity
#' reasons, dimesionality reduction on regions will only be done using the regions assigned to topics with high confidence 
#' (see binarizecisTopics()).
#' @param method Select the method for processing the cell assignments: 'Z-score' and 'Probability'. In the case of regions, 
#' an additional method, 'NormTop' is available (see getRegionScores()).
#' @param seed Integer for making the results reproducible.
#' @param ... See \code{DiffusionMap} from the package destiny.
#'
#' @return Returns a cisTopic object with the tSNE coordinates stored in object@@dr$DM.
#'
#' @details 'Z-score' computes the Z-score for each topic assingment per cell/region. 'Probability' divides the topic assignments by the total number
#' of assignments in the cell/region in the last iteration plus alpha. If using 'NormTop', regions are given an score defined by: \eqn{\beta_{w, k} (\log
#' \beta_{w,k} - 1 / K \sum_{k'} \log \beta_{w,k'})}.
#'
#' @export

runDM <- function(
  object,
  target,
  method='Z-score',
  seed=123,
  ...
){
  
  # Check info
  if (length(object@selected.model) < 1){
    stop('Please, run selectModel() first.')
  }
  
  # Check dependencies
  if(! "destiny" %in% installed.packages()){
    stop('Please, install destiny: \n source("https://bioconductor.org/biocLite.R") \n biocLite("destiny")')
  } else {
    require(destiny)
  }
  
  modelMat <- modelMatSelection(object, target, method)

  set.seed(seed)
  dif <- destiny::DiffusionMap(t(modelMat), ...)
  coord <- dif@eigenvectors
  rownames(coord) <- colnames(modelMat)
  colnames(coord) <- paste0('DC', seq(1:ncol(coord)))
  object@dr[[target]][['DiffusionMap']] <- coord
  object@calc.params[[target]][['runDM']] <- c(as.list(environment(), all = TRUE)[names(formals("runDM"))[c(2,3)]], list(...))
  return(object)
}

#' Calculate Principal Components on the cell-cisTopic distributions
#'
#' Calculate Principal Components (PCs) on the cell-cisTopic distributions
#' @param object Initialized cisTopic object, after the object@@selected.model has been filled.
#' @param target Whether dimensionality reduction should be applied on cells ('cell') or regions (region). Note that for speed and clarity
#' reasons, dimesionality reduction on regions will only be done using the regions assigned to topics with high confidence 
#' (see binarizecisTopics()).
#' @param method Select the method for processing the cell assignments: 'Z-score' and 'Probability'. In the case of regions, 
#' an additional method, 'NormTop' is available (see getRegionScores()).
#' @param ... See \code{prcomp} from the package stats.
#'
#' @return Returns a cisTopic object with a list of PCA information stored in object@@dr$cell$PCA or object@@dr$region$PCA.
#'
#' @slot loadings Matrix whose columns contain eigenvectors
#' @slot sdev Standard deviations of the PCs
#' @slot var.coord Coordinates of the variables (correlation between the variables and the PCs)
#' @slot var.cos2 Cos2 of the variables. Measures their representation quality.
#' @slot var.contrib Contributions of the variables to the PCs
#' @slot ind.coord Coordinates of individuals
#' @slot ind.cos2 Cos2 of the individuals
#' @slot ind.contrib Contributions of the individuals to the PCs
#' @slot eigs Eigenvalues, which measure the variability retained per PC
#' @slot variance.explained Percentage of variance explained by each component
#'
#' @details 'Z-score' computes the Z-score for each topic assingment per cell/region. 'Probability' divides the topic assignments by the total number
#' of assignments in the cell/region in the last iteration plus alpha. If using 'NormTop', regions are given an score defined by: \eqn{\beta_{w, k} (\log
#' \beta_{w,k} - 1 / K \sum_{k'} \log \beta_{w,k'})}.
#'
#' @importFrom stats prcomp
#' @export

runPCA <- function(
  object,
  target,
  method='Z-score',
  seed=123,
  ...
){
  
  # Check info
  if (length(object@selected.model) < 1){
    stop('Please, run selectModel() first.')
  }
  
  inimodelMat <- modelMatSelection(object, target, method)

  set.seed(seed)
  modelMat <- scale(t(inimodelMat), center=TRUE, scale=FALSE, ...)
  rownames(modelMat) <- colnames(inimodelMat)
  rm(inimodelMat)
  colnames(modelMat) <- paste('Topic', 1:ncol(modelMat), sep='')
  pca <- prcomp(modelMat)

  # Compute Coordinates
  loadings <- pca$rotation
  sdev <- pca$sdev
  var.coord <- t(t(loadings) * sdev)

  # Compute Cos2
  var.cos2 <- var.coord^2

  # Compute variable contributions
  var.contrib <- t(100*(t(var.coord^2))/colSums(var.coord^2))

  # Coordinates of individuals
  ind.coord <- pca$x

  # Cos2 of individuals
  pca.scale <- rep(1, ncol(modelMat))
  pca.center <- rep(1, ncol(modelMat))
  distance <- apply(modelMat, 1, function(row) sum(((row-pca.center)/pca.scale)^2))
  ind.cos2 <- apply(ind.coord, 2, function(ind) ind^2/distance)

  # Eigenvalues and variance explained
  eigs <- pca$sdev^2
  variance.explained <- round(eigs / sum(eigs) * 100, digits = 2)

  # Individual contributions
  ind.contrib <- t(apply(ind.coord, 1,  function(ind) 100*(1/nrow(modelMat))*(ind^2/eigs)))

  object@dr[[target]][['PCA']] <- list(loadings=loadings, sdev=sdev, var.coord=var.coord, var.cos2=var.cos2, var.contrib=var.contrib,
                             ind.coord=ind.coord, ind.cos2=ind.cos2, ind.contrib=ind.contrib, eigs=eigs, variance.explained=variance.explained)
  object@calc.params[[target]][['runPCA']] <- c(as.list(environment(), all = TRUE)[names(formals("runPCA"))[c(2,3)]], list(...))
  return(object)
}

#' Retrieve normalised topic-cell and region-topic assignments
#'
#' Retrieve topic-cell and region-topic assignments
#' @param object Initialized cisTopic object, after the object@@selected.model has been filled.
#' @param target Whether dimensionality reduction should be applied on cells ('cell') or regions ('region'). Note that for speed and clarity
#' reasons, dimesionality reduction on regions will only be done using the regions assigned to topics with high confidence 
#' (see binarizecisTopics()).
#' @param method Select the method for processing the cell assignments: 'Z-score' and 'Probability'. In the case of regions, 
#' an additional method, 'NormTop' is available (see getRegionScores()).
#' @param all.regions If target is region, whether to return a matrix with all regions or only regions belonging to binarized 
#' topics (see binarizecisTopics()).
#'
#' @details 'Z-score' computes the Z-score for each topic assingment per cell/region. 'Probability' divides the topic assignments by the total number
#' of assignments in the cell/region in the last iteration plus alpha. If using 'NormTop', regions are given an score defined by: \eqn{\beta_{w, k} (\log
#' \beta_{w,k} - 1 / K \sum_{k'} \log \beta_{w,k'})}.
#'
#' @importFrom stats prcomp
#' @export

modelMatSelection <- function(
  object,
  target,
  method,
  all.regions=FALSE
){
  # Check info
  if (length(object@selected.model) < 1){
    stop('Please, run selectModel() first.')
  }
  
  if (target == 'cell'){
    if (method == 'Z-score'){
      modelMat <- scale(object@selected.model$document_expects, center=TRUE, scale=TRUE)
    }
    
    else if (method == 'Probability'){
      alpha <- object@calc.params[['runModels']]$alpha/length(object@selected.model$topic_sums)
      modelMat <- apply(object@selected.model$document_sums, 2, function(x) {(x + alpha)/sum(x + alpha)})
    }
    else{
      stop('Incorrect method selected. Chose method between "Z-score" and "Probability".')
    }
    colnames(modelMat) <- object@cell.names
    rownames(modelMat) <- paste0('Topic', 1:nrow(modelMat))
  }
  
  else if (target == 'region'){
    if (!all.regions){
      if (length(object@binarized.cisTopics) < 1){
        stop('Please, use binarizecisTopics() first for defining the high confidence regions for dimensionality reduction!')
      }
      else {
        regions <- unique(unlist(lapply(object@binarized.cisTopics, rownames)))
      }
    }
    
    topic.mat <- object@selected.model$topics
    
    if (method == 'NormTop'){
      normalizedTopics <- topic.mat/(rowSums(topic.mat) + 1e-05)
      modelMat <- apply(normalizedTopics, 2, function(x) x * (log(x + 1e-05) - sum(log(x + 1e-05))/length(x)))
    }
    
    else if (method == 'Z-score'){
      modelMat <- scale(object@selected.model$topics, center=TRUE, scale=TRUE)
    }
    
    else if (method == 'Probability'){
      beta <- object@calc.params[['runModels']]$beta
      topic.mat <- object@selected.model$topics
      modelMat <-  (topic.mat + beta)/rowSums(topic.mat + beta)
    }
    
    else{
      stop('Incorrect method selected. Chose "NormTop", "Z-score" and "Probability".')
    }
    
    colnames(modelMat) <- object@region.names
    rownames(modelMat) <- paste0('Topic', 1:nrow(modelMat))
    
    if (!all.regions){
      modelMat <- modelMat[,regions]
    }
  }
  
  else{
    stop('Please, provide target="cell" or "region".')
  }
  
  return(modelMat)
}


#' Plot features (cells or regions) based on dimensionality reduction over cell-topic or topic-region distributions
#'
#' Plot features (cells or regions) based on dimensionality reduction over cell-topic or topic-region distributions
#'
#' @param object Initialized cisTopic object, after the object@@dr has been filled.
#' @param target Whether dimensionality reduction should be applied on cells ('cell') or regions (region). Note that for speed and clarity
#' reasons, dimesionality reduction on regions will only be done using the regions assigned to topics with high confidence 
#' (see binarizecisTopics()).
#' @param method Select the dimensionality reduction method to use for plotting: 'tSNE', 'Umap',  'PCA',  'Biplot', 'DM' (for Diffusion Map).
#' @param colorBy Select the cell metadata used to colour the plots with. By default, all features are used.
#' @param dim Dimensions to use in the plot (2 or 3). Biplot is only doable in 2D.
#' @param intervals Intervals to apply on the color palette for coloring continuous variables. By default, it is 10.
#' @param topic_contr Color by topic distribution ('Z-score' or 'Probability'). 
#' @param topics Vector containing the numbers of the topics for which the plot has to be made if topic_contr is not null.
#' @param col.low Color to use for lowest topic enrichment
#' @param col.mid Color to use for medium topic enrichment
#' @param col.high Color to use for high topic enrichment
#' @param labels Labels for the Z-score in the continuous variables plots
#' @param cex.legend Size of the legend
#' @param cex.dot Size of the dot
#' @param colsVar List specifying the colors to use for each label in each colouring level for cell metadata
#' @param plot_ly Whether plot_ly should be used for the plots
#' @param legend Whether plots should be given with a legend. If FALSE, the legend is given in an independent following plot.
#' @param ... Ignored.
#'
#' @return Plots cell states based on the dimensionality reduction method selected, coloured by the given metadata (one plot per feature).
#'
#' @export

plotFeatures <- function(
  object,
  target,
  method='tSNE',
  colorBy=NULL,
  dim=2,
  intervals=10,
  topic_contr=NULL,
  topics=NULL,
  col.low = "pink",
  col.mid = "red",
  col.high = "darkred",
  labels=3,
  colVars=list(),
  plot_ly=FALSE,
  legend=TRUE,
  cex.legend=0.7,
  factor.min=0.05,
  factor.max=0.75,
  cex.dot=1,
  ...
){
  if (is.null(colorBy) && is.null(topic_contr)){
    stop('Please select whether to color by some feature or topic contributions.')
  }
  # Check dependencies
  if(! "plotly" %in% installed.packages() && plot_ly ){
    stop('Please, install plotly: \n install.packages("plotly")')
  }  else {
    require(plotly)
  }
  
  if(! "scatterplot3d" %in% installed.packages() && !plot_ly && dim == 3 ){
    stop('Please, install scatterplot3d: \n install.packages("scatterplot3d")')
  } else {
    require(scatterplot3d)
  }
  
  # Select coordinates to plot
  coordinates <- .selectCoordinates(object, target, method, dim)  
  
  # Select cell/region data
  if (target == 'cell'){
    feature.data <- object@cell.data[rownames(coordinates),]
  }
  else if (target == 'region'){
    feature.data <- object@region.data[rownames(coordinates),]
  }
  feature.names <- rownames(feature.data)
  
  # Get par values
  par.opts <- par()

  # Variables to color by
  if (!is.null(colorBy)){
    
    if (sum(colorBy %in% colnames(feature.data)) != length(colorBy)){
      stop(paste('The variable', colorBy[-which(colorBy %in% colnames(feature.data))], 'is not included in the', target, 'data. Please check and re-run.'))
    }
    
    for (columnName in colorBy){
      variable <- setNames(feature.data[,columnName], feature.names)
      
      # Continuous plotting
      if(is.numeric(variable)){
        colorPal <- grDevices::colorRampPalette(c(col.low, col.mid, col.high))
        if (method != 'Biplot'){
          .plotContinuous(coordinates, variable, feature.names, colorPal, main=columnName, intervals=intervals, dim=dim, plot_ly=plot_ly, legend=legend, factor.min = factor.min, factor.max=factor.max, cex.dot=cex.dot)
        }
        else{
          coordinates <- .rescaleVector(coordinates[,c(1:2)])
          var.coord <- .rescaleVector(object@dr[[target]][['PCA']]$var.coord[,c(1:2)])
          .plotContinuous(coordinates, variable, feature.names, colorPal, main=columnName, intervals=intervals, dim=2, var.coord=var.coord, plot_ly=plot_ly, legend=legend, factor.min = factor.min, factor.max=factor.max, cex.dot=cex.dot)
        }
        
      }
      
      else{
        if (method != 'Biplot'){
          .plotFactor(coordinates, variable, feature.names, main=columnName, dim=dim, colVars=colVars, plot_ly=plot_ly, legend=legend, cex.legend = cex.legend, factor.min = factor.min, factor.max=factor.max, cex.dot=cex.dot)
        }
        else{
          coordinates <- .rescaleVector(coordinates[,c(1:2)])
          var.coord <- .rescaleVector(object@dr[[target]][['PCA']]$var.coord[,c(1:2)])
          .plotFactor(coordinates, variable, feature.names, main=columnName, dim=2, var.coord=var.coord, colVars=colVars, plot_ly=plot_ly, legend=legend, cex.legend = cex.legend, factor.min = factor.min, factor.max = factor.max, cex.dot=cex.dot)
        }
      }
    }
  }

    
  if (!is.null(topic_contr)){
    topic.mat <- modelMatSelection(object, target, topic_contr)
    
    if (!is.null(topics)){
      topic.mat <- topic.mat[topics,,drop=FALSE]
    }
    
    colorPal <- grDevices::colorRampPalette(c(col.low, col.mid, col.high))
    for (i in 1:nrow(topic.mat)){
      if(method != 'Biplot'){
        .plotContinuous(coordinates, topic.mat[i,], feature.names, colorPal, main=rownames(topic.mat)[i], intervals=intervals, dim=dim, labels=labels, plot_ly=plot_ly, legend=legend, factor.min = factor.min, factor.max=factor.max, cex.dot=cex.dot)
      }
      else{
        coordinates <- .rescaleVector(coordinates[,c(1:2)])
        var.coord <- .rescaleVector(object@dr[[target]][['PCA']]$var.coord[,c(1:2)])
        .plotContinuous(coordinates, topic.mat[i,], feature.names, colorPal, main=rownames(topic.mat)[i], intervals=intervals, dim=2, var.coord = var.coord, labels=labels, plot_ly=plot_ly, legend=legend, factor.min = factor.min, factor.max=factor.max, cex.dot=cex.dot)
      }
    }
  }
  # Restore par
  suppressWarnings(par(par.opts))
}



# Helper Functions

.selectCoordinates <- function(
  object,
  target,
  method,
  dim=dim
){
  if (!target %in% c('cell', 'region')){
    stop('Please, provide target="cell" or "region".')
  } 
  
  if (method == 'tSNE'){
    if (is.null(object@dr[[target]][['tSNE']])){
      stop(paste0('Please, run first: cisTopicObject <- runtSNE(cisTopicObject, target="', target,'", ...).'))
    }
    
    coordinates <- object@dr[[target]][['tSNE']]
    
    if (ncol(coordinates) == 2 && dim > 2){
      stop(paste0('Please, run first: cisTopicObject <- runtSNE(cisTopicObject, target="', target,'", dim=3 ...)for a 3D tSNE.'))
    }
  }
  
  else if (method == 'Umap'){
    if (is.null(object@dr[[target]][['Umap']])){
      stop(paste0('Please, run first: cisTopicObject <- runUmap(cisTopicObject, target="', target,'", ...).'))
    }
    
    coordinates <- object@dr[[target]][['Umap']]
    
    if (ncol(coordinates) == 2 && dim > 2){
      stop(paste0('Please, run first: cisTopicObject <- runUmap(cisTopicObject, target="', target,'", n_components=3 ...)for a 3D tSNE.'))
    }
  }
  
  else if (method == 'PCA' || method == 'Biplot'){
    if (is.null(object@dr[[target]][['PCA']])){
      stop(paste0('Please, run first: cisTopicObject <- runPCA(cisTopicObject, target="', target,'", ...).'))
    }
    coordinates <- object@dr[[target]][['PCA']]$ind.coord
  }
  
  else if (method == 'DM'){
    if (is.null(object@dr[[target]][['DiffusionMap']])){
      stop(paste0('Please, run first: cisTopicObject <- runDM(cisTopicObject, target="', target,'", ...).'))
    }
    coordinates <- object@dr[[target]][['DiffusionMap']]
  }
  return(coordinates)
}
  

.plotContinuous <- function(
  coordinates,
  variable,
  names,
  colorPal,
  main,
  intervals,
  dim=2,
  var.coord=NULL,
  labels=3,
  plot_ly=FALSE,
  legend=TRUE,
  factor.min=0.05,
  factor.max=0.2,
  cex.dot=1,
  ...
){
  cellColor <- setNames(adjustcolor(colorPal(intervals), alpha=.8)[as.numeric(cut(variable,breaks=10, right=F,include.lowest=T))], names)
  par(bty='n')
  
  if (!is.null(var.coord)){
    # Plot the correlation circle
    a <- seq(0, 2*pi, length = 100)
    plot(cos(a), sin(a), type = 'l', col="gray", xlab = "PC1",  ylab = "PC2", xlim = c(-1-factor.min, 1+factor.max), ylim= c(-1-factor.min, 1+factor.max))
    abline(h = 0, v = 0, lty = 2)

    # Plot dots
    points(coordinates, col=cellColor[rownames(coordinates)], pch=16, cex=cex.dot)

    # Add active variables
    arrows(0, 0, var.coord[, 1], var.coord[, 2], length = 0.1, angle = 15, code = 2)

    # Add labels
    text(var.coord, labels=rownames(var.coord), cex = 0.75, adj=1, font=2)
  }
  else{
    if (dim == 2){
      if (!plot_ly){
        plot(coordinates, col=cellColor[rownames(coordinates)], pch=16, xlab=colnames(coordinates)[1], ylab=colnames(coordinates)[2], xlim = c(min(coordinates[,1])-factor.min*abs(min(coordinates[,1])), max(coordinates[,1])+factor.max*abs(max(coordinates[,1]))), 
             ylim=c(min(coordinates[,2])-factor.min*abs(min(coordinates[,2])), max(coordinates[,2])+factor.max*abs(max(coordinates[,2]))), main=main, cex=cex.dot)
      } else {
        p <- plot_ly(x = coordinates[,1], y = coordinates[,2], color = variable, colors=adjustcolor(colorPal(intervals), alpha=.8)) %>%
          add_markers() %>%
          layout(title = main,
                 scene = list(xaxis = list(title = colnames(coordinates)[1]),
                 yaxis = list(title = colnames(coordinates)[2])))
        print(p)
      }
    }
    if (dim == 3) {
      if (!plot_ly){
        scatterplot3d(coordinates[,1], coordinates[,2], coordinates[,3], color=cellColor[rownames(coordinates)], pch=16, xlab=colnames(coordinates)[1], ylab=colnames(coordinates)[2], zlab = colnames(coordinates)[3], main=main,
                      xlim = c(min(coordinates[,1])-factor.min*abs(min(coordinates[,1])), max(coordinates[,1])+factor.max*abs(max(coordinates[,1]))), ylim = c(min(coordinates[,2])-factor.min*abs(min(coordinates[,2])), max(coordinates[,2])+factor.max*abs(max(coordinates[,2]))),
                      zlim = c(min(coordinates[,3])-factor.min*abs(min(coordinates[,3])), max(coordinates[,3])+factor.max*abs(max(coordinates[,3]))), cex.symbols=cex.dot)
      } else {
        p <- plot_ly(x = coordinates[,1], y = coordinates[,2], z = coordinates[,3], color = variable, colors=adjustcolor(colorPal(intervals), alpha=.8)) %>%
          add_markers() %>%
          layout(title = main,
                 scene = list(xaxis = list(title = colnames(coordinates)[1]),
                 yaxis = list(title = colnames(coordinates)[2]),
                 zaxis = list(title = colnames(coordinates)[3])))
        print(p)
      }
    }
  }
  
  if (!plot_ly){
    if (legend){
      .vertical.image.legend(c(min(variable), max(variable)), colorPal(intervals))
    } else {
      plot.new()
      .vertical.image.legend(c(min(variable), max(variable)), colorPal(intervals))
    }
  }
}

.distinctColorPalette <-function(k) {
  set.seed(123)
  if(packageVersion("scales") >= '1.1.0'){
    ColorSpace <- t(unique(col2rgb(scales::hue_pal(l=85)(2e3))))
  } else {
    ColorSpace <- t(unique(col2rgb(scales::hue_pal(l=60:100)(2e3))))
  }
  km <- kmeans(ColorSpace, k, iter.max=20)
  colors <- rgb(round(km$centers), maxColorValue=255)
  return(colors)
}

# Adapted from aqfig


.vertical.image.legend <- function(zlim, col){
  starting.par.settings <- par(no.readonly=TRUE)
  mai <- par("mai")
  fin <- par("fin")
  x.legend.fig <- c( 1.0-(mai[4]/fin[1]), 1.0 )
  y.legend.fig <- c( mai[1]/fin[2], 1.0-(mai[3]/fin[2]) )
  x.legend.plt <- c( x.legend.fig[1]+(0.08*(x.legend.fig[2]-x.legend.fig[1])),
                     x.legend.fig[2]-(0.6*(x.legend.fig[2]-x.legend.fig[1])) )
  y.legend.plt <- y.legend.fig
  cut.pts <- seq(zlim[1], zlim[2], length=length(col)+1)
  z <- ( cut.pts[1:length(col)] + cut.pts[2:(length(col)+1)] ) / 2
  par(new=TRUE, pty="m", plt=c(x.legend.plt, y.legend.plt), bty='o')
  image(x=1, y=z, z=matrix(z, nrow=1, ncol=length(col)),
        col=col, xlab="", ylab="", xaxt="n", yaxt="n")
  axis(4, mgp = c(3, 0.2, 0), las = 2, cex.axis=0.6, tcl=-0.1)
  box()
  mfg.settings <- par()$mfg
  par(starting.par.settings)
  par(mfg=mfg.settings, new=FALSE)
}

.plotFactor <- function(
  coordinates,
  variable,
  names,
  main,
  dim=2,
  var.coord=NULL,
  colVars=list(),
  plot_ly = FALSE,
  legend = TRUE, 
  cex.legend=0.7,
  factor.min=0.05,
  factor.max=0.5,
  cex.dot=1,
  ...
){
  levels <- as.vector(sort(unique(variable)))
  
  if(is.null(colVars[[main]])) {
    colVars[[main]] <- setNames(.distinctColorPalette(k=length(levels)), levels)
  }
  
  cellColor <- setNames(colVars[[main]][variable], names)
  
  par(bty = 'n')
  
  if (!is.null(var.coord)){
    # Plot the correlation circle
    a <- seq(0, 2*pi, length = 100)
    plot(cos(a), sin(a), type = 'l', col="gray", xlab = "PC1",  ylab = "PC2", xlim = c(-1-factor.min, 1+factor.max), ylim = c(-1-factor.min, 1+factor.max))
    abline(h = 0, v = 0, lty = 2)

    # Plot dots
    points(coordinates, col=cellColor[rownames(coordinates)], pch=16, cex=cex.dot)

    # Add active variables
    arrows(0, 0, var.coord[, 1], var.coord[, 2], length = 0.1, angle = 15, code = 2)

    # Add labels
    text(var.coord, labels=rownames(var.coord), cex = 0.75, adj=1, font=2)
  }
  else{
    if (dim == 2){
      if (!plot_ly){
        plot(coordinates, col=cellColor[rownames(coordinates)], pch=16, xlab=colnames(coordinates)[1], ylab=colnames(coordinates)[2], xlim = c(min(coordinates[,1])-factor.min*abs(min(coordinates[,1])), max(coordinates[,1])+factor.max*abs(max(coordinates[,1]))), 
             ylim=c(min(coordinates[,2])-factor.min*abs(min(coordinates[,2])), max(coordinates[,2])+factor.max*abs(max(coordinates[,2]))), main=main, cex=cex.dot)
      } else {
        p <- plot_ly(x = coordinates[,1], y = coordinates[,2], color = variable, colors= unique(colVars[[main]][variable])) %>%
          add_markers() %>%
          layout(title = main,
                 scene = list(xaxis = list(title = colnames(coordinates)[1]),
                 yaxis = list(title = colnames(coordinates)[2])))
        print(p)
      }
    }
    if (dim == 3) {
      if (!plot_ly){
        s3d <-  scatterplot3d(coordinates[,1], coordinates[,2], coordinates[,3], color=cellColor[rownames(coordinates)], pch=16, xlab=colnames(coordinates)[1], ylab=colnames(coordinates)[2], zlab = colnames(coordinates)[3], main=main, 
                              xlim = c(min(coordinates[,1])-factor.min*abs(min(coordinates[,1])), max(coordinates[,1])+factor.max*abs(max(coordinates[,1]))), ylim = c(min(coordinates[,2])-factor.min*abs(min(coordinates[,2])), max(coordinates[,2])+factor.max*abs(max(coordinates[,2]))),
                              zlim = c(min(coordinates[,3])-factor.min*abs(min(coordinates[,3])), max(coordinates[,3])+factor.max*abs(max(coordinates[,3]))), cex.symbols=cex.dot)
      } else {
        p <- plot_ly(x = coordinates[,1], y = coordinates[,2], z = coordinates[,3], color = variable, colors = unique(colVars[[main]][variable])) %>%
          add_markers() %>%
          layout(title = main,
                 scene = list(xaxis = list(title = colnames(coordinates)[1]),
                 yaxis = list(title = colnames(coordinates)[2]),
                 zaxis = list(title = colnames(coordinates)[3])))
        print(p)
      }
    }
  }
  
  if (!plot_ly){
    if (legend){
      if (dim == 2){
        legend("topright", legend =  names(colVars[[main]]), fill=colVars[[main]], cex=cex.legend)
      } else {
        legend(s3d$xyz.convert(max(coordinates[,1]), max(coordinates[,2]), max(coordinates[,3])), fill = colVars[[main]], legend = names(colVars[[main]]), cex = cex.legend)
      }
    }
    else {
      plot.new()
      legend(x="top", legend = names(colVars[[main]]), fill=colVars[[main]], cex=cex.legend, title=as.expression(bquote(bold(.(main)))))
    }
  }
}

.rescaleVector<- function(coordinates){
  vectors.length <- apply(coordinates, 1, function(x) sqrt(x[1]^2+x[2]^2))
  target.vector.length <- vectors.length/max(vectors.length)
  new.rescaled <- t(sapply(1:nrow(coordinates), function(i) (coordinates[i,]/sqrt(sum(coordinates[i,]^2))*target.vector.length[i])))
  rownames(new.rescaled) <- rownames(coordinates)
  return(new.rescaled)
}

#' Cell-cisTopic distribution heatmap
#'
#' Plot cell states based on dimensionality reduction over cell-cisTopic distributions
#' @param object Initialized cisTopic object, after the object@@selected.model has been filled.
#' @param method Select the normalization method to use for plotting: 'Z-score' or 'Probability'.
#' @param colorBy Select the cell metadata used to colour the plots with. By default, all categorical features are used.
#' @param colVars List specifying the colors to use for each label in each colouring level
#' @param col.low Color to use for lowest topic enrichment
#' @param col.mid Color to use for medium topic enrichment
#' @param col.high Color to use for high topic enrichment
#' @param select.cells If a subset of cells want to be used for making the heatmap, selected cell names can be provided (as a vector).
#' @param ... See \code{Heatmap} from ComplexHeatmap
#'
#' @return Heatmap clustering cells based on their cell-cisTopic distributions.
#' 
#' @details 'Z-score' computes the Z-score for each topic assingment per cell/region and 'Probability' divides the topic assignments by the total number
#' of assignments in the cell/region in the last iteration plus alpha.
#'
#' @export

cellTopicHeatmap <- function(
  object,
  method ='Z-score',
  colorBy=NULL,
  colVars=NULL,
  col.low = "floralwhite",
  col.mid = "pink",
  col.high = "red",
  select.cells=NULL,
  ...)
  {
  
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
  
  # Check info
  if (length(object@selected.model) < 1){
    stop('Please, run selectModel() first.')
  }
  
  topic.mat <- modelMatSelection(object, 'cell', method)

  rownames(topic.mat) <- paste('Topic', seq(1,nrow(topic.mat)))
  colnames(topic.mat) <- object@cell.names
  
  if (!is.null(select.cells)){
    if (length(select.cells) > 1){
      topic.mat <- topic.mat[ ,select.cells]
      object.cell.data <- object@cell.data[select.cells,]
    }else{
      stop('Only a cell name has been provided. Please select a group of cells and provide them as a vector.')
    }
  } else {
    object.cell.data <- object@cell.data
  }

  cl.cells <- fastcluster::hclust.vector(t(topic.mat), method="ward", metric="euclidean")
  dd.cells <- as.dendrogram(cl.cells)
  colorPal <- grDevices::colorRampPalette(c(col.low, col.mid, col.high))
  
  if (is.null(colorBy)){
    heatmap <- ComplexHeatmap::Heatmap(data.matrix(topic.mat), col=colorPal(20), cluster_columns=dd.cells, name=method,
                                       show_column_names=FALSE, show_row_names = TRUE, 
                                       heatmap_legend_param = list(legend_direction = "horizontal", legend_width = unit(5, "cm"), title_position='topcenter'),
                                       column_title = "Topic contribution per cell", column_title_gp = gpar(fontface = 'bold'), ...)
    ComplexHeatmap::draw(heatmap, heatmap_legend_side = "bottom")
  } else {
    for (variable in colorBy){
      if(is.null(colVars[[variable]])) {
        colVars[[variable]] <- setNames(.distinctColorPalette(length(unique(object@cell.data[,variable]))), as.vector(sort(unique(object@cell.data[,variable]))))
        cellColor <- setNames(colVars[[variable]][object.cell.data[,variable]], rownames(object.cell.data))
      }
    }
    annotation <- ComplexHeatmap::HeatmapAnnotation(df = object.cell.data[,colorBy,drop=FALSE], col = colVars, which='column', width = unit(5, "mm"))
    heatmap <- ComplexHeatmap::Heatmap(data.matrix(topic.mat), col=colorPal(20), cluster_columns=dd.cells, name=method,
                                       show_column_names=FALSE, show_row_names = TRUE, top_annotation = annotation, 
                                       heatmap_legend_param = list(legend_direction = "horizontal", legend_width = unit(5, "cm"), title_position='topcenter'),
                                       column_title = "Topic contribution per cell", column_title_gp = gpar(fontface = 'bold'), ...)
    ComplexHeatmap::draw(heatmap, heatmap_legend_side = "bottom", annotation_legend_side = "right")
  }
}
