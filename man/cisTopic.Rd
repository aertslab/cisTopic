% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cisTopic.R
\docType{class}
\name{cisTopic}
\alias{cisTopic}
\alias{cisTopic-class}
\title{The cisTopic Class}
\description{
The cisTopic object contains the initial data and the rest of the objects (e.g. models, distributions),
created during the cisTopic workflow. cisTopic is ideally suited for single cell epigenomics experiments.
}
\details{
Each cisTopic object contains several objects, which are listed below:
}
\section{Slots}{

\describe{
\item{\code{count.matrix}}{Matrix with the number of counts per region (rows) in each cell (columns). When dealing with single cell methylation data,
we use beta values instead of counts.}

\item{\code{binary.count.matrix}}{Matrix with regions as rows and cells as columns, where 1 represents an accessible region and 0 a
non-accessible region. By default, a region is considered accessible if there is at least one count within the region.}

\item{\code{is.acc}}{Counts threshold to determine if a region is accessible (0 by default)}

\item{\code{models}}{List with all the LDA models trained. Each model is represented as an object of the class cisTopicModel. This slot will only be
functional when the number of topics and the size of the data set allows it.}

\item{\code{log.lik}}{Data frame containing the log likelihood of the models built.}

\item{\code{selected.model}}{List with the information regards the model selected. The unnormalized cell assignments throughtout the sampling iterations are
stored in \code{cisTopicObject@selected.model$document_expects}; while the corresponding unnormalized region assignments are stored in \code{cisTopicObject@selected.model$topics}.}

\item{\code{dr}}{Dimensionality reduction outputs, which can be applied on 'cell' or 'region'. Within each slot ('cell' or 'regions'), there are
sub-slots named by technique: 'Umap', 'tSNE', 'DiffussionMap', 'PCA', 'Biplot'.}

\item{\code{calc.params}}{Named list to store all calculation-related parameter choices.}

\item{\code{cell.names}}{Vector with the names of all the single cells (column names of the counts/accessibility matrix).}

\item{\code{cell.data}}{Data frame that contains the meta-information about each cell, starting with number of counts detected in the defined regions (nCounts), and accessible
regions (nAcc); more information is added using \code{AddCellMetadata}.}

\item{\code{region.names}}{Vector with the position coordinates of all the  regions (row names of the counts/accessibility matrix).}

\item{\code{region.ranges}}{GRanges object with the regions coordinates.}

\item{\code{region.data}}{Data frame that contains the meta-information about each region, starting with the region coordinates ('seqnames', 'start', 'end') and width (width), number of
counts mapped to that region accross the data set (nCounts), number of cells in which the region is accessible (nCells); more information is added using \code{AddRegionMetadata}}

\item{\code{binarized.cisTopics}}{List containing the regions that are considered as part of a cisTopic.}

\item{\code{signatures}}{List containing regions belonging to given region signatures.}

\item{\code{binarized.regions.to.Rct}}{List containing the RcisTarget regions that map to the top regions in each of the cisTopics.}

\item{\code{binarized.RcisTarget}}{List containing objects of the class RcisTarget produced by running RcisTarget in the top binarized regions
for each cisTopic.}

\item{\code{binarized.rGREAT}}{List with rGREAT results when using the top binarized regions.}

\item{\code{cistromes.ctx}}{List with cistromes based on ctx regions.}

\item{\code{cistromes.regions}}{List with cistromes based on data regions.}

\item{\code{cistromes.genes}}{List with cistromes based on data regions converted to genes.}

\item{\code{project.name}}{Vector with the name of the project (for record keeping).}

\item{\code{other}}{List to store any kind of related data (e.g. region lists).}

\item{\code{version}}{Version of package used in object creation.}
}}

