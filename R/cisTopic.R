#' The cisTopic Class
#'
#' The cisTopic object contains the initial data and the rest of the objects (e.g. models, distributions),
#' created during the cisTopic workflow. cisTopic is ideally suited for single cell epigenomics experiments.
#'
#' Each cisTopic object contains several objects, which are listed below:
#'
#' @slot count.matrix Matrix with the number of counts per region (rows) in each cell (columns). When dealing with single cell methylation data,
#' we use beta values instead of counts.
#' @slot binary.count.matrix Matrix with regions as rows and cells as columns, where 1 represents an accessible region and 0 a
#' non-accessible region. By default, a region is considered accessible if there is at least one count within the region.
#' @slot is.acc Counts threshold to determine if a region is accessible (0 by default)
#' @slot models List with all the LDA models trained. Each model is represented as an object of the class cisTopicModel. This slot will only be
#' functional when the number of topics and the size of the data set allows it.
#' @slot log.lik Data frame containing the log likelihood of the models built.
#' @slot selected.model List with the information regards the model selected. The unnormalized cell assignments throughtout the sampling iterations are
#' stored in \code{cisTopicObject@@selected.model$document_expects}; while the corresposding unnormalized region assignments are stored in \code{cisTopicObject@@selected.model$topics}.
#' @slot dr Dimensionality reduction outputs, which can be applied on 'cell' or 'region'. Within each slot ('cell' or 'regions'), there are
#' sub-slots named by technique: 'Umap', 'tSNE', 'DiffussionMap', 'PCA', 'Biplot'.
#' @slot calc.params Named list to store all calculation-related parameter choices.
#' @slot cell.names Vector with the names of all the single cells (column names of the counts/accessibility matrix).
#' @slot cell.data Data frame that contains the meta-information about each cell, starting with number of counts detected in the defined regions (nCounts), and accessible
#' regions (nAcc); more information is added using \code{AddCellMetadata}.
#' @slot region.names Vector with the position coordinates of all the  regions (row names of the counts/accessibility matrix).
#' @slot region.ranges GRanges object with the regions coordinates.
#' @slot region.data Data frame that contains the meta-information about each region, starting with the region coordinates ('seqnames', 'start', 'end') and width (width), number of
#' counts mapped to that region accross the data set (nCounts), number of cells in which the region is accessible (nCells); more information is added using \code{AddRegionMetadata}
#' @slot binarized.cisTopics List containing the regions that are considered as part of a cisTopic.
#' @slot signatures List containing regions belonging to given region signatures.
#' @slot binarized.regions.to.Rct List containing the RcisTarget regions that map to the top regions in each of the cisTopics.
#' @slot binarized.RcisTarget List containing objects of the class RcisTarget produced by running RcisTarget in the top binarized regions
#' for each cisTopic.
#' @slot binarized.rGREAT List with rGREAT results when using the top binarized regions.
#' @slot cistromes.ctx List with cistromes based on ctx regions.
#' @slot cistromes.regions List with cistromes based on data regions.
#' @slot cistromes.genes List with cistromes based on data regions converted to genes.
#' @slot project.name Vector with the name of the project (for record keeping).
#' @slot other List to store any kind of related data (e.g. region lists).
#' @slot version Version of package used in object creation.
#'
#' @name cisTopic
#' @rdname cisTopic
#' @aliases cisTopic-class
#' @exportClass cisTopic
#' @importFrom Rcpp evalCpp
#' @import methods
#' @importFrom utils globalVariables

cisTopic <- methods::setClass(
  "cisTopic",
  slots = c(
    count.matrix = "ANY",
    binary.count.matrix = "ANY",
    is.acc = 'numeric',
    models = "list",
    selected.model = "ANY",
    log.lik = "data.frame",
    dr = "list",
    calc.params = "list",
    cell.names = "character",
    cell.data = "data.frame",
    region.names = "character",
    region.ranges = "ANY",
    region.data = "data.frame",
    binarized.cisTopics = "list",
    signatures = "list",
    binarized.regions.to.Rct = "list",
    binarized.RcisTarget="list",
    binarized.rGREAT="list",
    cistromes.ctx="list",
    cistromes.regions="list",
    cistromes.genes="list",
    project.name="character",
    other="ANY",
    version="ANY"
  )
)

#' Show method for cisTopic
#'
#' @param object A cisTopic object
#' @name show
#'

setMethod(
  f = "show",
  signature = "cisTopic",
  definition = function(object) {
    cat(
      "An object of class",
      class(object),
      "in project",
      object@project.name,
      "\n",
      length(x = object@region.names),
      "regions across",
      length(x = object@cell.names),
      "samples.\n"
    )
    invisible(x = NULL)
  }
)
