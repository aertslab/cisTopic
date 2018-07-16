#' Run tSNE on the cell-cisTopic distributions
#'
#' Run tSNE on the cell-cisTopic distributions.
#' @param object Initialized cisTopic object, after the object@@selected.model has been filled.
#' @param method Select the method for processing the cell assignments: 'Zscore',  'probability', 'predictive.distribution'.
#' @param seed Integer for making the results reproducible.
#' @param ... See Rtsne from the package Rtsne.
#'
#' @return Returns a cisTopic object with the tSNE coordinates stored in object@@dr$tSNE.
#'
#' @details Zscore computes the Z-score for each topic assingment per cell, probability divides the topic assignments by the total number
#' of assignments in the cell in the last iteration plus alpha and predictive distribution calculates the probability of the seeing
#' the regions per cell.
#'
#' @importFrom Rtsne Rtsne
#' @importFrom lda predictive.distribution
#' @export

runtSNE <- function(
  object,
  method='Zscore',
  seed=123,
  ...
){
  if (method == 'Zscore'){
    modelMat <- scale(object@selected.model$document_expects, center=TRUE, scale=TRUE)
  }

  if (method == 'probability'){
    alpha <- 50/length(cisTopicObject@selected.model$topic_sums)
    modelMat <- apply(object@selected.model$document_sums, 2, function(x) {(x + alpha)/sum(x + alpha)})
  }

  if (method == 'predictive.distribution'){
    if (is.null(object@selected.model$predictive.distribution)){
      object <- predictiveDistribution(object)
    }
    modelMat <- object@selected.model$predictive.distribution
  }

  set.seed(seed)
  tSNE <- Rtsne(t(modelMat), ...)
  rownames(tSNE$Y) <- object@cell.names
  if (ncol(tSNE$Y) == 2){
    colnames(tSNE$Y) <- c("tsne1", "tsne2")
  }
  else if (ncol(tSNE$Y) == 3){
    colnames(tSNE$Y) <- c("tsne1", "tsne2", "tsne3")
  }

  object@dr[['tSNE']] <- tSNE$Y
  object@calc.params[['tSNE']] <- c(as.list(environment(), all = TRUE)[names(formals("runtSNE"))[c(2,3)]], list(...))
  return(object)
}


#' Calculate Diffusion Components on the cell-cisTopic distributions
#'
#' Calculate Diffusion Components on the cell-cisTopic distributions
#' @param object Initialized cisTopic object, after the object@@selected.model has been filled.
#' @param method Select the method for processing the cell assignments: 'Zscore',  'probability', 'predictive.distribution'.
#' @param seed Integer for making the results reproducible.
#' @param ... See DiffusionMap from the package destiny.
#'
#' @return Returns a cisTopic object with the tSNE coordinates stored in object@@dr$DM.
#'
#' @details Zscore computes the Z-score for each topic assingment per cell, probability divides the topic assignments by the total number
#' of assignments in the cell in the last iteration plus alpha and predictive distribution calculates the probability of the seeing
#' the regions per cell.
#'
#' @importFrom destiny DiffusionMap
#' @importFrom lda predictive.distribution
#' @export

runDM <- function(
  object,
  method='Zscore',
  seed=123,
  ...
){
  if (method == 'Zscore'){
    modelMat <- scale(object@selected.model$document_expects, center=TRUE, scale=TRUE)
  }

  if (method == 'probability'){
    alpha <- 50/length(cisTopicObject@selected.model$topic_sums)
    modelMat <- apply(object@selected.model$document_sums, 2, function(x) {(x + alpha)/sum(x + alpha)})
  }

  if (method == 'predictive.distribution'){
    if (is.null(object@selected.model$predictive.distribution)){
      object <- predictiveDistribution(object)
    }
    modelMat <- object@selected.model$predictive.distribution
  }

  set.seed(seed)
  dif <- DiffusionMap(t(modelMat), ...)
  coord <- dif@eigenvectors
  rownames(coord) <- object@cell.names
  colnames(coord) <- paste('DC', seq(1:ncol(coord)))
  object@dr[['DiffusionMap']] <- coord
  object@calc.params[['runDM']] <- c(as.list(environment(), all = TRUE)[names(formals("runDM"))[c(2,3)]], list(...))
  return(object)

}

#' Calculate Principal Components on the cell-cisTopic distributions
#'
#' Calculate Principal Components (PCs) on the cell-cisTopic distributions
#' @param object Initialized cisTopic object, after the object@@selected.model has been filled.
#' @param method Select the method for processing the cell assignments: 'Zscore',  'probability', 'predictive.distribution'.
#' @param ... See prcomp from the package stats.
#'
#' @return Returns a cisTopic object with a list of PCA information stored in object@@dr$PCA.
#'
#' @slot loadings Matrix whose columns contain eigenvectors
#' @slot sdev Standard deviations of the PCs
#' @slot var.coord Coordinates of the variables (correlation between the variables and the PCs)
#' @slot var.cos2 Cos2 of the variables. Measures their representation quality.
#' @slot var.contrib Contributions of the variables to the PCs
#' @slot ind.coord Coordinates of individuals
#' @slot ind.cos2 Cos2 of the individuals.
#' @slot ind.contrib Contributions of the individuals to the PCs
#' @slot eigs Eigenvalues, which measure the variability retained per PC
#' @slot variance.explained Percentage of variance explained by each component
#'
#' @details Zscore computes the Z-score for each topic assingment per cell, probability divides the topic assignments by the total number
#' of assignments in the cell in the last iteration plus alpha and predictive distribution calculates the probability of the seeing
#' the regions per cell.
#'
#' @importFrom stats prcomp
#' @importFrom lda predictive.distribution
#' @export

runPCA <- function(
  object,
  method='Zscore',
  seed=123,
  ...
){
  if (method == 'Zscore'){
    modelMat <- scale(object@selected.model$document_expects, center=TRUE, scale=TRUE)
  }

  if (method == 'probability'){
    alpha <- 50/length(cisTopicObject@selected.model$topic_sums)
    modelMat <- apply(object@selected.model$document_sums, 2, function(x) {(x + alpha)/sum(x + alpha)})
  }

  if (method == 'predictive.distribution'){
    if (is.null(object@selected.model$predictive.distribution)){
      object <- predictiveDistribution(object)
    }
    modelMat <- object@selected.model$predictive.distribution
  }

  set.seed(seed)
  modelMat <- scale(t(modelMat), center=TRUE, scale=FALSE, ...)
  rownames(modelMat) <- object@cell.names
  colnames(modelMat) <-  paste('Topic', 1:ncol(modelMat), sep='')
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

  object@dr[['PCA']] <- list(loadings=loadings, sdev=sdev, var.coord=var.coord, var.cos2=var.cos2, var.contrib=var.contrib,
                             ind.coord=ind.coord, ind.cos2=ind.cos2, ind.contrib=ind.contrib, eigs=eigs, variance.explained=variance.explained)
  object@calc.params[['runPCA']] <- c(as.list(environment(), all = TRUE)[names(formals("runPCA"))[c(2,3)]], list(...))
  return(object)
}

#' Plot cell states based on dimensionality reduction over cell-cisTopic distributions
#'
#' Plot cell states based on dimensionality reduction over cell-cisTopic distributions
#' @param object Initialized cisTopic object, after the object@@selected.model has been filled.
#' @param method Select the dimensionality reduction method to use for plotting: 'tSNE',  'PCA',  'Biplot', 'DiffusionMap'.
#' @param colorBy Select the cell metadata used to colour the plots with. By default, all features are used.
#' @param dim Dimensions to use in the plot (2 or 3). Biplot is only doable in 2D.
#' @param intervals Intervals to apply on the color palette for coloring continuous variables. By default, it is 10.
#' @param topic_contr Color by topic distribution ('Zscore' or 'probability'). 
#' @param topics Vector containing the numbers of the topics for which the plot has to be made if topic_contr is not null.
#' @param col.low Color to use for lowest topic enrichment
#' @param col.mid Color to use for medium topic enrichment
#' @param col.high Color to use for high topic enrichment
#' @param labels Labels for the Zscore in the continuous variables plots.
#' @param colsVar List specifying the colors to use for each label in each colouring level for cell metadata
#'
#' @return Plots cell states based on the dimensionality reduction method selected, coloured by the given metadata (one plot per feature).
#'
#' @import scatterplot3d
#' @export

plotCellStates <- function(
  object,
  method='tSNE',
  colorBy=NULL,
  dim=2,
  intervals=10,
  topic_contr='Zscore',
  topics='all',
  col.low = "dodgerblue",
  col.mid = "floralwhite",
  col.high = "brown1",
  labels=3,
  colVars=list(),
  ...
){
  if (method == 'tSNE'){
    if (is.null(object@dr[['tSNE']])){
      stop('Run the function runtSNE first!')
    }
    coordinates <- object@dr[['tSNE']]
    if (ncol(coordinates) == 2 && dim > 2){
      stop('If you want to plot tSNE in 3D, please use runtSNE with dims = 3 as argument.')
    }
  }
  else if (method == 'PCA' || method == 'Biplot'){
    if (is.null(object@dr[['PCA']])){
      stop('Run the function runPCA first!')
    }
    coordinates <- object@dr[['PCA']]$ind.coord
  }
  else if (method == 'DM'){
    if (is.null(object@dr[['DiffusionMap']])){
      stop('Run the function runDM first!')
    }
    coordinates <- object@dr[['DiffusionMap']]
  }

  if (is.null(colorBy)){
    colorBy <- colnames(object@cell.data)
  }

  for (var in colorBy){
    variable <- object@cell.data[,var]
    names(variable) <- object@cell.names
    if(is.numeric(variable)){
      colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))
      if (method != 'Biplot'){
        .plotContinuous(coordinates, variable, object@cell.names, colorPal, main=var, intervals=intervals, dim=dim)
      }
      else{
        coordinates <- .rescaleVector(coordinates[,c(1:2)])
        var.coord <- .rescaleVector(object@dr[['PCA']]$var.coord[,c(1:2)])
        .plotContinuous(coordinates, variable, object@cell.names, colorPal, main=var, intervals=intervals, dim=2, var.coord=var.coord)
      }

    }
    else{
      if (method != 'Biplot'){
        .plotFactor(coordinates, variable, object@cell.names, main=var, dim=dim, colVars=colVars)
      }
      else{
        coordinates <- .rescaleVector(coordinates[,c(1:2)])
        var.coord <- .rescaleVector(object@dr[['PCA']]$var.coord[,c(1:2)])
        .plotFactor(coordinates, variable, object@cell.names, main=var, dim=2, var.coord=var.coord, colVars=colVars)
      }
    }
  }

  if (!is.null(topic_contr)){
    if (topic_contr == 'Zscore'){
      topic.mat <- scale(object@selected.model$document_expects, center=TRUE, scale=TRUE)
    }
    if (topic_contr == 'probability'){
      alpha <- 50/length(cisTopicObject@selected.model$topic_sums)
      topic.mat <- apply(object@selected.model$document_sums, 2, function(x) {(x + alpha)/sum(x + alpha)})
    }
    rownames(topic.mat) <- paste('Topic', 1:nrow(topic.mat))
   
    if(topics != 'all'){
      topic.mat <- topic.mat[topics,,drop=FALSE]
    }
    
    colorPal <- grDevices::colorRampPalette(c(col.low, col.mid, col.high))
    for (i in 1:nrow(topic.mat)){
      if(method != 'Biplot'){
        .plotContinuous(coordinates, topic.mat[i,], object@cell.names, colorPal, main=rownames(topic.mat)[i], intervals=intervals, dim=dim, labels=labels)
      }
      else{
        coordinates <- .rescaleVector(coordinates[,c(1:2)])
        var.coord <- .rescaleVector(object@dr[['PCA']]$var.coord[,c(1:2)])
        .plotContinuous(coordinates, topic.mat[i,], object@cell.names, colorPal, main=rownames(topic.mat)[i], intervals=intervals, dim=2, var.coord = var.coord, labels=labels)
      }
    }
  }
}


#' Helper Functions

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
  ...
){
  par.opts <- par()
  cellColor <- setNames(adjustcolor(colorPal(intervals), alpha=.8)[as.numeric(cut(variable,breaks=10, right=F,include.lowest=T))], names)
  layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
  par(bty = 'n')
  if (!is.null(var.coord)){
    # Plot the correlation circle
    a <- seq(0, 2*pi, length = 100)
    plot( cos(a), sin(a), type = 'l', col="gray", xlab = "PC1",  ylab = "PC2")
    abline(h = 0, v = 0, lty = 2)

    # Plot dots
    points(coordinates, col=cellColor[rownames(coordinates)], pch=16)

    # Add active variables
    arrows(0, 0, var.coord[, 1], var.coord[, 2], length = 0.1, angle = 15, code = 2)

    # Add labels
    text(var.coord, labels=rownames(var.coord), cex = 0.75, adj=1, font=2)
  }
  else{
    if (dim == 2){plot(coordinates, col=cellColor[rownames(coordinates)], pch=16, xlab=colnames(coordinates)[1], ylab=colnames(coordinates)[2])}
    if (dim == 3) {scatterplot3d(coordinates[,1], coordinates[,2], coordinates[,3], color=cellColor[rownames(coordinates)], pch=16, xlab=colnames(coordinates)[1], ylab=colnames(coordinates)[2], zlab = colnames(coordinates)[3], main=main)}
  }

  legend_image <- as.raster(matrix(rev(colorPal(intervals)), ncol=1))
  plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = main)
  labels.name <- seq(from = min(variable), to = max(variable), length.out = labels)
  labels.name <- round(labels.name, 2)
  text(x=1.5, y = seq(0,1,l=labels), cex=1, labels = labels.name)
  rasterImage(legend_image, 0, 0, 1, 1)
  suppressWarnings(par(par.opts))
}

.plotFactor <- function(
  coordinates,
  variable,
  names,
  main,
  dim=2,
  var.coord=NULL,
  colVars=list(),
  ...
){
  par.opts <- par()
  levels <- as.vector(sort(unique(variable)))
  if(is.null(colVars[[main]])) {
    colVars[[main]] <- setNames(rainbow(length(unique(variable)), s=0.5), levels)
  }
  cellColor <- setNames(colVars[[main]][variable], names)

  layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
  par(bty = 'n')
  if (!is.null(var.coord)){
    # Plot the correlation circle
    a <- seq(0, 2*pi, length = 100)
    plot( cos(a), sin(a), type = 'l', col="gray", xlab = "PC1",  ylab = "PC2")
    abline(h = 0, v = 0, lty = 2)

    # Plot dots
    points(coordinates, col=cellColor[rownames(coordinates)], pch=16)

    # Add active variables
    arrows(0, 0, var.coord[, 1], var.coord[, 2], length = 0.1, angle = 15, code = 2)

    # Add labels
    text(var.coord, labels=rownames(var.coord), cex = 0.75, adj=1, font=2)
  }
  else {
    if (dim == 2){plot(coordinates, col=cellColor[rownames(coordinates)], pch=16, xlab=colnames(coordinates)[1], ylab=colnames(coordinates)[2])}
    if (dim == 3) {scatterplot3d(coordinates[,1], coordinates[,2], coordinates[,3], color=cellColor[rownames(coordinates)], pch=16, xlab=colnames(coordinates)[1], ylab=colnames(coordinates)[2], zlab = colnames(coordinates)[3])}
  }

  plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = main)
  legend(x=0, y=1, legend=names(colVars[[main]]), pch=16, col=colVars[[main]], bty='n', x.intersp = 0.2)
  suppressWarnings(par(par.opts))
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
#' @param method Select the normalization method to use for plotting: 'Zscore' or 'probability'.
#' @param colorBy Select the cell metadata used to colour the plots with. By default, all categorical features are used.
#' @param colVars List specifying the colors to use for each label in each colouring level
#' @param col.low Color to use for lowest topic enrichment
#' @param col.mid Color to use for medium topic enrichment
#' @param col.high Color to use for high topic enrichment
#' @param ... See \code{aheatmap} from NMF
#'
#' @return Heatmap clustering cells based on their cell-cisTopic distributions.
#'
#' @importFrom fastcluster hclust.vector
#' @importFrom NMF aheatmap
#' @export

cellTopicHeatmap <- function(
  object,
  method ='Zscore',
  colorBy=NULL,
  colVars=list(),
  col.low = "dodgerblue",
  col.mid = "floralwhite",
  col.high = "brown1",
  ...)
  {
  if (method == 'Zscore'){
    topic.mat <- scale(object@selected.model$document_expects, center=TRUE, scale=TRUE)
  }
  else if (method == 'probability'){
    alpha <- 50/length(cisTopicObject@selected.model$topic_sums)
    topic.mat <- apply(object@selected.model$document_sums, 2, function(x) {(x + alpha)/sum(x + alpha)})
  }

  rownames(topic.mat) <- paste('Topic', seq(1,nrow(topic.mat)))
  colnames(topic.mat) <- object@cell.names

  cl.cells <- hclust.vector(t(topic.mat), method="ward", metric="euclidean")
  dd.cells <- as.dendrogram(cl.cells)
  colorPal <- grDevices::colorRampPalette(c(col.low, col.mid, col.high))

  object.cell.data <- object@cell.data
  if (is.null(colorBy)){
    colorBy <- colnames(object.cell.data)[-which(as.vector(sapply(object.cell.data, is.numeric)))]
  }

  for (variable in colorBy){
    if(is.null(colVars[[variable]])) {
      colVars[[variable]] <- setNames(rainbow(length(unique(object.cell.data[,variable])), s=0.5), unique(object.cell.data[,variable]))
      cellColor <- setNames(colVars[[variable]][object.cell.data[,variable]], object@cell.names)
    }
  }

  nmf.options(grid.patch=TRUE)
  NMF::aheatmap(topic.mat, scale="none", revC=TRUE, main='cisTopic contributions per cell', sub='Column normalized topic contribution',
                Colv=dd.cells, annCol=object.cell.data[object@cell.names, colorBy, drop=FALSE],
                annColor=colVars, labCol=NA,
                color = colorPal(20))
}
