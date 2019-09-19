#' Initialize and setup the cisTopic object starting from CellRanger ATAC output matrix.
#'
#' Initializes the cisTopic object from CellRanger ATAC output files and defined regions. 
#' @param data_folder Matrix folder from CellRanger ATAC (/outs/filtered_peak_bc_matrix/, after untar file).''
#' @param metrics List of CellRanger barcodes and their metrics (/outs/singlecell.csv)
#' @param project.name Project name (string).
#' @param min.cells Minimal number of cells in which the region has to be accessible. By default, all regions accessible in at least one cell are kept.
#' @param min.regions Minimal number of regions that have to be accessible within a cell to be kept. By default, all cells with at least one region accessible are kept.
#' @param is.acc Number of counts necessary to consider a region as accessible.
#' @param keepCountsMatrix Whether to keep the counts matrix or not inside the object. For large matrices, we recommend to set this
#' to FALSE.
#' @return Returns a cisTopic object with the counts data stored in object@@count.matrix.
#' object@@binary.count.matrix, object@@cell.names, object@@cell.data (including counting statistics), object@@regions.ranges, object@@regions.data are also initialized.
#'
#' @import Matrix
#' 
#' @export
#'
#' @examples
#' data_folder <- '/outs/filtered_peak_bc_matrix'
#' cisTopicObject <- createcisTopicObjectFrom10Xmatrix(data_folder, metrics)
#' cisTopicObject

createcisTopicObjectFrom10Xmatrix <- function(
  data_folder,
  metrics,
  project.name = "cisTopicProject",
  min.cells = 1,
  min.regions = 1,
  is.acc = 1,
  keepCountsMatrix = TRUE,
  ...
) {
  # Check dependencies
  
  if(! "Matrix" %in% installed.packages()){
    stop('Please, install Matrix: \n install.packages("Matrix")')
  } else {
    suppressMessages(require(Matrix))
  }
  
  # Prepare barcodes and metrics
  CR_metrics <- fread(metrics, header = T, sep = ',')
  CR_metrics <- CR_metrics[which(unlist(as.vector(CR_metrics[,'is__cell_barcode'] == 1))),]
  CR_metrics <- as.data.frame(CR_metrics)
  rownames(CR_metrics) <- unlist(as.vector(CR_metrics[,1]))
  CR_metrics <- CR_metrics[,-1]
  
  barcodes <- rownames(CR_metrics) # default (e.g. V2)
  barcodesFile <- paste0(data_folder, '/barcodes.tsv') # overwrite if exists
  if(file.exists(barcodesFile)) barcodes <- read.table(paste0(data_folder, '/barcodes.tsv'))[,1]
  
  # Read peaks
  regions <- paste0(data_folder, '/peaks.bed')
  regions_frame <- read.table(regions)
  if (ncol(regions_frame) >= 5){
    regions_frame <- regions_frame[,c(1:3, 5)]
  } else {
    regions_frame <- regions_frame[,c(1:3)]
    regions_frame <- cbind(regions_frame, rep('*', nrow(regions_frame)))
  }
  
  colnames(regions_frame) <- c('Chr', 'Start', 'End', 'Strand')
  
  
  m <- readMM(paste0(data_folder, '/matrix.mtx'))
  colnames(m) <- barcodes
  rownames(m) <- paste0(regions_frame[,1], ':', regions_frame[,2], '-', regions_frame[,3])
  
  # Prepare cell data
  print('Creating cisTopic object...')
  count.matrix <- m
  
  if(any(!barcodes %in% rownames(CR_metrics))) warning("Some barcodes are missing from the metrics, the count matrix might include empty droplets (e.g. non-cells).")
  CR_metrics <- CR_metrics[barcodes,,drop=FALSE]
  rownames(CR_metrics) <- barcodes
  cell.data <- CR_metrics
  colnames(cell.data) <- paste0('10X_', colnames(cell.data))
  
  object <- createcisTopicObject(count.matrix = count.matrix, 
                                 project.name = project.name,
                                 min.cells = min.cells,
                                 min.regions = min.regions,
                                 is.acc = is.acc,
                                 keepCountsMatrix=keepCountsMatrix)
  
  object <- addCellMetadata(object, cell.data = as.data.frame(cell.data))
  return(object)
}


#' Initialize and setup the cisTopic object starting from CellRanger ATAC output files and defined regions
#'
#' Initializes the cisTopic object from CellRanger ATAC output files and defined regions. 
#' @param fragments Fragment file from CellRanger ATAC (/outs/fragments.tsv.gz)
#' @param regions Path to the bed file with the defined regions.
#' @param metrics List of CellRanger barcodes and their metrics (/outs/singlecell.csv)
#' @param project.name Project name (string).
#' @param min.cells Minimal number of cells in which the region has to be accessible. By default, all regions accessible in at least one cell are kept.
#' @param min.regions Minimal number of regions that have to be accessible within a cell to be kept. By default, all cells with at least one region accessible are kept.
#' @param is.acc Number of counts necessary to consider a region as accessible.
#' @param keepCountsMatrix Whether to keep the counts matrix or not inside the object. For large matrices, we recommend to set this
#' to FALSE.
#'
#' @return Returns a cisTopic object with the counts data stored in object@@count.matrix.
#' object@@binary.count.matrix, object@@cell.names, object@@cell.data (including counting statistics), object@@regions.ranges, object@@regions.data are also initialized.
#'
#' @import Matrix
#' @import data.table
#' @import GenomicRanges
#' @import dplyr
#' 
#' @export
#'
#' @examples
#' fragments <- '/outs/fragments.tsv.gz'
#' regions <- '/outs/peaks.bed'
#' metrics <- '/outs/singlecell.csv'
#' cisTopicObject <- createcisTopicObjectFrom10X(fragments, regions, metrics)
#' cisTopicObject

createcisTopicObjectFrom10X <- function(
  fragments,
  regions,
  metrics,
  project.name = "cisTopicProject",
  min.cells = 1,
  min.regions = 1,
  is.acc = 1,
  keepCountsMatrix = TRUE,
  ...
) {
  # Check dependencies
  
  if(! "Matrix" %in% installed.packages()){
    stop('Please, install Matrix: \n install.packages("Matrix")')
  } else {
    suppressMessages(require(Matrix))
  }
  
  if(! "GenomicRanges" %in% installed.packages()){
    stop('Please, install Matrix: \n install.packages("GenomicRanges")')
  } else {
    suppressMessages(require(GenomicRanges))
  }
  
  if(! "data.table" %in% installed.packages()){
    stop('Please, install Matrix: \n install.packages("data.table")')
  } else {
    suppressMessages(require(data.table))
  }
  if(! "dplyr" %in% installed.packages()){
    stop('Please, install Matrix: \n install.packages("dplyr")')
  } else {
    suppressMessages(require(dplyr))
  }
  
  print('Counting...')
  # Prepare annotation
  regions_frame <- read.table(regions)
  if (ncol(regions_frame) >= 5){
    regions_frame <- regions_frame[,c(1:3, 5)]
  } else {
    regions_frame <- regions_frame[,c(1:3)]
    regions_frame <- cbind(regions_frame, rep('*', nrow(regions_frame)))
  }
  
  colnames(regions_frame) <- c('Chr', 'Start', 'End', 'Strand')
  peak_ranges <- makeGRangesFromDataFrame(as.data.frame(regions_frame))
  
  # Prepare barcodes and metrics
  CR_metrics <- fread(metrics, header = T, sep = ',')
  CR_metrics <- CR_metrics[which(unlist(as.vector(CR_metrics[,'is__cell_barcode'] == 1))),]
  CR_metrics <- as.data.frame(CR_metrics)
  rownames(CR_metrics) <- unlist(as.vector(CR_metrics[,1]))
  CR_metrics <- CR_metrics[,-1]
  barcodes <- rownames(CR_metrics)
  
  # Select fragments from good cells
  good_fragments <- data.table::fread(cmd=paste0("zcat < ", fragments)) %>% 
    data.frame() %>% filter(V4 %in% barcodes) %>%  # filter for barcodes in our search set
    GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "V1", start.field = "V2", end.field = "V3", keep.extra.columns = TRUE)
  
  # Get the overlaps with peaks
  Frag2peak <- GenomicRanges::findOverlaps(peak_ranges, good_fragments)
  
  # Establish a numeric index for the barcodes for sparse matrix purposes
  id <- factor(as.character(GenomicRanges::mcols(good_fragments)$V4), levels = barcodes)
  
  # Make sparse matrix with counts with peaks by  unique barcode
  countdf <- data.frame(peaks = S4Vectors::queryHits(Frag2peak),
                        sample = as.numeric(id)[S4Vectors::subjectHits(Frag2peak)]) %>%
    dplyr::group_by(peaks, sample) %>% dplyr::summarise(count = n()) %>% data.matrix()
  
  m <- Matrix::sparseMatrix(i = c(countdf[,1], length(peak_ranges)),
                            j = c(countdf[,2], length(barcodes)),
                            x = c(countdf[,3],0))
  colnames(m) <- barcodes
  rownames(m) <- paste0(regions_frame[,1], ':', regions_frame[,2], '-', regions_frame[,3])
  
  # Prepare cell data
  print('Creating cisTopic object...')
  count.matrix <- m
  
  cell.data <- CR_metrics
  colnames(cell.data) <- paste0('10X_', colnames(cell.data))

  object <- createcisTopicObject(count.matrix = count.matrix, 
                                 project.name = project.name,
                                 min.cells = min.cells,
                                 min.regions = min.regions,
                                 is.acc = is.acc,
                                 keepCountsMatrix=keepCountsMatrix)
  
  object <- addCellMetadata(object, cell.data = as.data.frame(cell.data))
  return(object)
}


#' Initialize and setup the cisTopic object starting from bam files and defined regions
#'
#' Initializes the cisTopic object from the raw bam files and defined regions
#' @param bamfiles List with the paths to the bam files.
#' @param regions Path to the bed file with the defined regions.
#' @param project.name Project name (string).
#' @param min.cells Minimal number of cells in which the region has to be accessible. By default, all regions accessible in at least one cell are kept.
#' @param min.regions Minimal number of regions that have to be accessible within a cell to be kept. By default, all cells with at least one region accessible are kept.
#' @param is.acc Number of counts necessary to consider a region as accessible.
#' @param keepCountsMatrix Whether to keep the counts matrix or not inside the object. For large matrices, we recommend to set this
#' to FALSE.
#' @param paired Whether data should be treated as paired end or not. If it is FALSE, we count a read if its 5' end falls within
#' the region, if false, fragments will be counted instead of individual reads.
#' @param ... See \code{featureCounts} function from Rsubread.
#'
#' @return Returns a cisTopic object with the counts data stored in object@@count.matrix.
#' object@@binary.count.matrix, object@@cell.names, object@@cell.data (including counting statistics), object@@regions.ranges, object@@regions.data are also initialized.
#'
#' @import Matrix
#' 
#' @export
#'
#' @examples
#' bamfiles <- c('example_1.bam', 'example_2.bam', 'example_3.bam')
#' regions <- 'example.bed'
#' cisTopicObject <- createcisTopicObjectFromBAM(bamfiles, regions)
#' cisTopicObject

createcisTopicObjectFromBAM <- function(
  bamfiles,
  regions,
  project.name = "cisTopicProject",
  min.cells = 1,
  min.regions = 1,
  is.acc = 1,
  keepCountsMatrix = TRUE,
  paired = FALSE,
  ...
) {
  # Check dependencies
  if(! "Rsubread" %in% installed.packages()){
    if (Sys.info()['sysname'] != 'Windows'){
      stop('Please, install Rsubread: \n source("https://bioconductor.org/biocLite.R") \n biocLite("Rsubread")')
    }
    else{
      stop('The Rsubread package is not available for Windows. See our FAQ webpages for possible alternatives at: ')
    }
  } else {
    suppressMessages(require(Rsubread))
  }
  
  if(! "Matrix" %in% installed.packages()){
    stop('Please, install Matrix: \n install.packages("Matrix")')
  } else {
    suppressMessages(require(Matrix))
  }
  
  # Prepare annotation
  regions_frame <- read.table(regions)
  if (ncol(regions_frame) >= 5){
    regions_frame <- regions_frame[,c(1:3, 5)]
  }
  else {
    regions_frame <- regions_frame[,c(1:3)]
    regions_frame <- cbind(regions_frame, rep('*', nrow(regions_frame)))
  }
  
  GeneID <- paste(regions_frame[,1], ':', regions_frame[,2], '-', regions_frame[,3], sep='')
  regions_frame <- cbind(GeneID, regions_frame)
  colnames(regions_frame) <- c('GeneID', 'Chr', 'Start', 'End', 'Strand')

  # Count reads
  count.data <- Rsubread::featureCounts(bamfiles, annot.ext=regions_frame, read2pos=5, isPairedEnd=paired, ...)

  # Prepare cell data
  print('Creating cisTopic object...')
  count.matrix <- Matrix(count.data$counts, sparse=TRUE)
  cell.data <- t(count.data$stat)
  colnames(cell.data) <- cell.data[1,]
  cell.data <- cell.data[-1,]
  cell.data <- cell.data[,c('Assigned', 'Unassigned_NoFeatures')]
  row_names <- rownames(cell.data)
  cell.data <- apply(cell.data, 2, function(x) as.numeric(as.character(x)))
  rownames(cell.data) <- row_names

  Total_reads <- cell.data[,'Assigned'] + cell.data[,'Unassigned_NoFeatures']
  pct_ReadsInPeaks <- cell.data[,'Assigned']/Total_reads
  cell.data <- cbind(Total_reads, pct_ReadsInPeaks, cell.data)
  

  
  object <- createcisTopicObject(count.matrix = count.matrix, 
                                 project.name = project.name,
                                 min.cells = min.cells,
                                 min.regions = min.regions,
                                 is.acc = is.acc,
                                 keepCountsMatrix=keepCountsMatrix)
  
  object <- addCellMetadata(object, cell.data = as.data.frame(cell.data))
  return(object)
}

#' Initialize and setup the cisTopic object starting from methylation files
#'
#' Initializes the cisTopic object from DNA methylation calls. These files must be tab separated files containing chromosome, position, number of methylated
#' reads and total number of reads.
#' @param methfiles List with the paths to the tab separated files containing the methylation information.
#' @param regions Path to the bed file with the defined regions.
#' @param project.name Project name (string).
#' @param min.cells Minimal number of cells in which the region has to be accessible. By default, all regions accessible in at least one cell are kept.
#' @param min.regions Minimal number of regions that have to be accessible within a cell to be kept. By default, all cells with at least one region accessible are kept.
#' @param is.acc Minimal beta value to be considered methylated. As default, it is 0.5.
#' @param ... Ignored
#'
#' @return Returns a cisTopic object with the counts data (in this case, beta values data) stored in object@@count.matrix. NA's are filled as 0's.
#' object@@binary.count.matrix, object@@cell.names, object@@cell.data (including counting statistics), object@@regions.ranges, object@@regions.data are also initialized.
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomicRanges findOverlaps
#' @importFrom stats aggregate
#' @export
#'
#' @examples
#' methfiles <- c('meth_1.txt', 'meth_2.txt', 'meth_3.txt')
#' regions <- 'example.bed'
#' cisTopicObject <- createcisTopicObjectFromBAM(methfiles, regions)
#' cisTopicObject

createcisTopicObjectFromMeth <- function(
  methfiles,
  regions,
  project.name = "cisTopicProject",
  min.cells = 1,
  min.regions = 1,
  is.acc = 0.5,
  ...
) {
  # Prepare annotation
  regions_frame <- read.table(regions)[,c(1:3)]
  colnames(regions_frame) <- c('seqnames', 'start', 'end')
  rownames(regions_frame) <- paste(regions_frame$seqnames, ':', regions_frame$start, '-', regions_frame$end, sep='')
  regions_granges <- makeGRangesFromDataFrame(as.data.frame(regions_frame))

  # Prepare beta matrix, cell data and region data
  beta.matrix <- matrix(0, nrow(regions_frame), length(methfiles), sparse=TRUE)
  rownames(beta.matrix) <- rownames(regions_frame)
  colnames(beta.matrix) <- methfiles

  cell.data <- matrix(0, length(methfiles), 2)
  rownames(cell.data) <- methfiles
  colnames(cell.data) <- c('Methylated reads', 'Total reads')

  region.data <- matrix(0, nrow(regions_frame), 2)
  rownames(region.data) <- rownames(regions_frame)
  colnames(region.data) <- c('Methylated reads', 'Total reads')

  for (file in methfiles){
    # Read data
    print(paste('Reading file ', file, '...', sep=''))
    sample <- read.table(file)
    sample <- sample[, c(1,2,2,3,4)]
    colnames(sample) <- c('seqnames', 'start', 'end', 'Methylated', 'Total_reads')
    sample_granges <- makeGRangesFromDataFrame(as.data.frame(sample), keep.extra.columns=TRUE)

    # Fill cell data
    cell.data[file, 1] <- sum(sample_granges$Methylated)
    cell.data[file, 2] <- sum(sample_granges$Total_reads)

    # Calculate aggregate reads per region (methylated and total)
    overlaps <- findOverlaps(sample_granges, regions_granges)
    sites <- sample_granges[queryHits(overlaps)]
    sumMethylated <- aggregate(sites$Methylated, list(subjectHits(overlaps)), sum)
    sumReads <- aggregate(sites$Total_reads, list(subjectHits(overlaps)), sum)

    # Fill region data
    rownames(sumMethylated) <- rownames(regions_frame)[sumMethylated[,1]]
    region.data[rownames(sumMethylated),1] <- region.data[rownames(sumMethylated),1] + sumMethylated[,2]

    rownames(sumReads) <- rownames(regions_frame)[sumReads[,1]]
    region.data[rownames(sumReads),2] <- region.data[rownames(sumReads),2] + sumReads[,2]

    # Calculate beta
    aggrBeta <- sumMethylated[,2]/sumReads[,2]
    names(aggrBeta) <- rownames(regions_frame)[sumMethylated[,1]]
    beta.matrix[names(aggrBeta), file] <- aggrBeta
  }
  object <- createcisTopicObject(count.matrix = beta.matrix, project.name = project.name, min.cells = min.cells, min.regions = min.regions, is.acc = is.acc, ...)
  object <- addCellMetadata(object, cell.data = as.data.frame(cell.data))
  object <- addRegionMetadata(object, region.data = as.data.frame(region.data))
  return(object)
}


#' Initialize and setup the cisTopic object
#'
#' Initializes the cisTopic object from a counts matrix
#' @param count.matrix Count matrix containing cells as rows and regions as columns. The row names must be the coordinates in position format (e.g. chr1:110-610).
#' We recommend to use as input sparse matrices.
#' @param project.name Project name (string).
#' @param min.cells Minimal number of cells in which the region has to be accessible. By default, all regions accessible in at least one cell are kept.
#' @param min.regions Minimal number of regions that have to be accessible within a cell to be kept. By default, all cells with at least one region accessible are kept.
#' @param is.acc Number of counts necessary to consider a region as accessible. When using single cell methylation data, this threshold
#' represents the beta value from which above regions will be considered methylated.
#' @param keepCountsMatrix Whether to keep the counts matrix or not inside the object. For large matrices, we recommend to set this
#' to FALSE.
#' @param ... Ignored
#'
#' @return Returns a cisTopic object with the counts data stored in object@@count.matrix.
#' object@@binary.count.matrix, object@@cell.names, object@@cell.data, object@@regions.ranges, object@@regions.data are also initialized.
#'
#' @details For single cell methylation data, the input matrix should contain columns as cells, regions as rows and beta values as values instead of counts.
#' @importFrom utils packageVersion
#' @importFrom Matrix colSums rowSums
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @export
#'
#' @examples
#' cisTopic_mel <- createcisTopicObject(count.matrix = count.matrix)
#' cisTopic_mel
#'

createcisTopicObject <- function(
  count.matrix,
  project.name = "cisTopicProject",
  min.cells = 1,
  min.regions = 1,
  is.acc = 1,
  keepCountsMatrix=TRUE,
  ...
) {
  cisTopic.version <- packageVersion("cisTopic")
  object <- new(
    Class = "cisTopic",
    is.acc = is.acc,
    project.name = project.name,
    version = cisTopic.version
  )
  
  # Binarize and filter matrix
  object.binary.count.matrix <- Matrix(1*(count.matrix >= is.acc), sparse=TRUE)
  num.acc.cells <- Matrix::rowSums(object.binary.count.matrix)
  num.acc.regions <- Matrix::colSums(object.binary.count.matrix)
  cells.use <- which(num.acc.regions >= min.regions)
  regions.use <- which(num.acc.cells >= min.cells)
  
  count.matrix <- count.matrix[regions.use, cells.use]
  
  # Take values from the count matrix
  
  # Set cell names
  object@cell.names <- colnames(count.matrix)
  
  # Set region names
  object@region.names <- rownames(count.matrix)
  
  # Set region ranges
  seqnames <- sapply(strsplit(rownames(count.matrix), split = ":"), "[", 1)
  coord <- sapply(strsplit(rownames(count.matrix), split = ":"), "[", 2)
  start <- sapply(strsplit(coord, split = "-"), "[", 1)
  end <- sapply(strsplit(coord, split = "-"), "[", 2)
  bed_coord <- cbind(seqnames, start, end)
  rownames(bed_coord) <- rownames(count.matrix)
  object@region.ranges <- makeGRangesFromDataFrame(as.data.frame(bed_coord))
  
  # Cell data
  nCounts_celldata <- Matrix::colSums(count.matrix)
  
  # Region data
  nCounts_regiondata <- Matrix::rowSums(count.matrix)
  
  
  if (keepCountsMatrix == TRUE){
    object@count.matrix <- count.matrix
  } 
  
  rm(count.matrix)
  
  object.binary.count.matrix <- object.binary.count.matrix[regions.use, cells.use]
  
  # Set cells data
  nCounts <- nCounts_celldata
  rm(nCounts_celldata)
  nAcc <- Matrix::colSums(object.binary.count.matrix)
  object.cell.data <- cbind(nCounts, nAcc)
  object.cell.data <- apply(object.cell.data, 2, function(x) as.numeric(as.character(x)))
  rownames(object.cell.data) <- object@cell.names
  object@cell.data <- as.data.frame(object.cell.data)
  
  # Set regions data
  nCounts <- nCounts_regiondata
  rm(nCounts_regiondata)
  nCells <- as.numeric(Matrix::rowSums(object.binary.count.matrix))
  width <- abs(as.numeric(end)-as.numeric(start))
  object.region.data <- cbind(seqnames, start, end, width, nCounts, nCells)
  rownames(object.region.data) <- object@region.names
  object@region.data <- as.data.frame(object.region.data)
  object@region.data[,2:ncol(object@region.data)] <- apply(object@region.data[,2:ncol(object@region.data)], 2, function(x) as.numeric(as.character(x)))
  
  object@binary.count.matrix <- object.binary.count.matrix
  rm(object.binary.count.matrix)
  
  object@calc.params[['createcisTopicObject']] <- list(min.cells = min.cells, min.regions = min.regions, is.acc = is.acc)
  
  return(object)
}

#' Add cells metadata to a cisTopic object
#'
#' Add cells metadata to a given cisTopic object
#' @param object cisTopic object
#' @param cell.data Additional metadata to add to the cisTopic object. It should be a data frame were the rows are cells and the columns data fields.
#' All cells included in the cisTopic object must be included in this dataframe (additional cells will be ignored).
#' @param ... Ignored
#'
#' @return Returns a cisTopic object with the new data added in object@@cell.data.
#'
#' @export
#'
#' @examples
#' cisTopic_mel <- createcisTopicObject(count.matrix = count.matrix)
#' cisTopic_mel <- addCellMetadata(cisTopic_mel, cell.data = cell.data)
#' cisTopic_mel
#'

addCellMetadata <- function(
  object,
  cell.data,
  ...
) {

  object.cell.data <- object@cell.data
  if(is.null(cell.data)){
    stop('Please, provide data.')
  }
  else{
    if (is.null(rownames(cell.data))){
      stop('Please, set row names for the cell data.')
    }
    if (is.null(colnames(cell.data))){
      stop('Please, provide column names for the cell data.')
    }
    if (sum(rownames(cell.data) %in% rownames(object.cell.data)) < nrow(object.cell.data)){
      stop('Are all the cells included in the new metadata?')
    }
    cell.data <- as.data.frame(cell.data)
    cell.data <- cell.data[rownames(object.cell.data),,drop=FALSE]
    cell.data <- droplevels(cell.data)
    if (sum(colnames(object.cell.data) %in% colnames(cell.data)) > 0){
      object.cell.data <- object.cell.data[,-which(colnames(object.cell.data) %in% colnames(cell.data))]
    }
    column_names<- c(colnames(object.cell.data), colnames(cell.data))
    object.cell.data <- cbind(object.cell.data, cell.data)
    colnames(object.cell.data) <- column_names
    object@cell.data <- object.cell.data

    return(object)
  }
}

#' Add regions metadata to a cisTopic object
#'
#' Add regions metadata to a given cisTopic object
#' @param object cisTopic object
#' @param region.data Additional metadata to add to the cisTopic object. It should be a data frame were the rows are regions and the columns data fields.
#' All regions included in the cisTopic object must be included in this dataframe (additional cells will be ignored), named by their position coordinates.
#' @param ... Ignored
#'
#' @return Returns a cisTopic object with the new data added in object@@region.data.
#'
#' @export
#'
#' @examples
#' cisTopic_mel <- createcisTopicObject(count.matrix = count.matrix)
#' cisTopic_mel <- addRegionMetadata(cisTopic_mel, region.data = region.data)
#' cisTopic_mel
#'

addRegionMetadata <- function(
  object,
  region.data,
  ...
) {

  object.region.data <- object@region.data
  if(is.null(region.data)){
    stop('Please, provide data.')
  }
  else{
    if (is.null(rownames(region.data))){
      stop('Please, set row names for the region data.')
    }
    if (is.null(colnames(region.data))){
      stop('Please, set column name for the region data.')
    }
    if (sum(rownames(region.data) %in% rownames(object.region.data)) < nrow(object.region.data)){
      stop('Are all the regions included in the new metadata?')
    }
    region.data <- as.data.frame(region.data)
    region.data <- region.data[rownames(object.region.data),,drop=FALSE]
    region.data <- droplevels(region.data)
    if (sum(colnames(object.region.data) %in% colnames(region.data)) > 0){
      object.region.data <- object.region.data[,-which(colnames(object.region.data) %in% colnames(region.data))]
    }
    column_names <- c(colnames(object.region.data), colnames(region.data))
    object.region.data <- cbind(object.region.data, region.data)
    colnames(object.region.data) <- column_names
    object@region.data <- object.region.data

    return(object)
  }
}

#' Rename cells in a cisTopic object
#'
#' Rename cells in a cisTopic object
#' @param object cisTopic object
#' @param names Vector with the new names to be given to the cells.
#' @param ... Ignored
#'
#' @return Returns a cisTopic object with updated cell names in the relevant slots.
#'
#' @export
#'
#' @examples
#' cisTopic_mel <- createcisTopicObject(count.matrix = count.matrix)
#' cisTopic_mel <- renameCells(cisTopic_mel, names)
#'
#' cisTopic_mel
#'

renameCells <- function(
  object,
  names,
  ...
) {
  object@cell.names <- names
  rownames(object@cell.data) <- names
  
  if (!is.null(object@count.matrix)){
    colnames(object@count.matrix) <- names
  }
  
  if (!is.null(object@binary.count.matrix)){
    colnames(object@binary.count.matrix) <- names
  }
  
  if (!is.null(object@dr[['cell']][['tSNE']])){
    rownames(object@dr[['cell']][['tSNE']]) <- names
  }
  
  if (!is.null(object@dr[['cell']][['Umap']])){
    rownames(object@dr[['cell']][['Umap']]) <- names
  }
  
  if (!is.null(object@dr[['cell']][['DiffusionMap']])){
    rownames(object@dr[['cell']][['DiffusionMap']]) <- names
  }
  
  if (!is.null(object@dr[['cell']][['PCA']])){
    rownames(object@dr[['cell']][['PCA']]$ind.coord) <- names
  }

  return(object)
}
 
