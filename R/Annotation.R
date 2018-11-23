#' annotateRegions
#'
#' Add peak annotation to the regions (e.g. type of region, closest gene, ...)
#' @param object Initialized cisTopic object.
#' @param txdb TxDb object matching the genome used for mapping
#' @param annoDb Annotation package matchng the organism of origin
#' @param ... See annotatePeak from ChIPseeker
#'
#' @return Peak annotation is added as new columns to object@@region.data
#'
#' @details See \code{SignaturesHeatmap} for estimating the region types per topic.
#' @export

annotateRegions <- function(
  object,
  txdb,
  annoDb,
  ...
){
  # Check depedencies
  if(! "ChIPseeker" %in% installed.packages()){
    stop('Please, install ChIPseeker: \n source("https://bioconductor.org/biocLite.R") \nbiocLite("ChIPseeker")')
  } else {
    require(ChIPseeker)
  }
  
  regions <- object@region.ranges
  elementMetadata(regions)[["regionNames"]] <- names(regions)
  object.region.data <- object@region.data
  regionAnno <- as.data.frame(annotatePeak(regions, TxDb=txdb, annoDb=annoDb, ...))
  regionAnno[grep("Exon", regionAnno$annotation), "annotation"] <- "Exon"
  regionAnno[grep("Intron", regionAnno$annotation), "annotation"] <- "Intron"
  rownames(regionAnno) <- regionAnno$regionNames
  regionAnno <- regionAnno[rownames(object.region.data), !(colnames(regionAnno) %in% c("regionNames", "strand"))]
  object.region.data[, colnames(regionAnno)] <- regionAnno
  object@region.data <- object.region.data
  return(object)
}

#' GREAT
#'
#' Run rGREAT in the binarized cisTopics
#' @param object Initialized cisTopic object, after object@@binarized.cisTopics has been filled.
#' @param genome Genome to which the data was aligned (see rGREAT for options).
#' @param liftOver GRangesList object containing the original coordinates (in the same format as
#' object@@region.names) in the data set as slot names and the corresponding mapping regions as a GRanges object in the slot.
#' @param fold_enrichment Minimum binomial fold enrichment to keep term.
#' @param geneHits Minimum number of genes associated to keep term.
#' @param sign Maximum adjusted p-value to keep term.
#' @param request_interval Time interval for two requests.
#' @param ... See \code{submitGreatJob} from rGREAT
#'
#' @return Non-empty GREAT results per topic are return as a list to object@@binarized.rGREAT
#'
#' @details This function works with regions annotated to hg19, hg18, mm10, mm9 and danRer7. For other genomes, a liftOver step to one of the available genomes is required.
#'
#'
#' @examples
#' object <- GREAT(object, request_interval = 10)
#'
#' @importFrom plyr llply
#' @export

GREAT <- function(
  object,
  genome='hg19',
  liftOver=NULL,
  fold_enrichment=2,
  geneHits=1,
  sign=0.05,
  request_interval=20,
  ...
  ){
  # Check info
  if (length(object@binarized.cisTopics) < 1){
    stop('Please, run binarizecisTopics() first.')
  }
  
  # Check dependencies
  if(! "rGREAT" %in% installed.packages()){
    stop('Please, install rGREAT: \n source("https://bioconductor.org/biocLite.R") \nbiocLite("rGREAT")')
  } else {
    require(rGREAT)
  }
  
  if (is.null(liftOver)){
    object.binarized.rGREAT <- llply(1:length(object@binarized.cisTopics), function(i) .doGREAT(object@region.ranges[rownames(object@binarized.cisTopics[[i]])], genome, fold_enrichment, geneHits, sign, request_interval, ...))
  } else {
    object.binarized.rGREAT <- llply(1:length(object@binarized.cisTopics), function(i) .doGREAT(unlist(liftOver[rownames(object@binarized.cisTopics[[i]])], recursive = TRUE, use.names = TRUE), genome, fold_enrichment, geneHits, sign, request_interval, ...))
  }
  names(object.binarized.rGREAT) <- names(object@binarized.cisTopics)
  object@binarized.rGREAT <- object.binarized.rGREAT
  return(object)
}

# Helper functions

.doGREAT <- function(
  coord,
  genome,
  fold_enrichment,
  geneHits,
  sign,
  request_interval,
  ...
  )
  {
  coord <- sortSeqlevels(coord)
  coord <- sort(coord)
  job <- submitGreatJob(coord, species=genome, request_interval = request_interval, ...)
  tb <- getEnrichmentTables(job, ontology=availableOntologies(job), request_interval = request_interval)
  for (i in 1:length(tb)){
    tb[[i]] <- tb[[i]][-which(tb[[i]][,'Binom_Fold_Enrichment'] < fold_enrichment),]
    if (length(which(tb[[i]][,'Hyper_Observed_Gene_Hits'] < geneHits)) > 1){
      tb[[i]] <- tb[[i]][-which(tb[[i]][,'Hyper_Observed_Gene_Hits'] < geneHits),]
    }
    tb[[i]] <- tb[[i]][-which(tb[[i]][,'Binom_Adjp_BH'] > sign),]
    tb[[i]] <- tb[[i]][-which(tb[[i]][,'Hyper_Adjp_BH'] > sign),]
  }
  tb <- tb[sapply(tb, function(x) dim(x)[1]) > 0]
  return(tb)
}

#' ontologyDotPlot
#'
#' Dot plots for visualizing GREAT results
#' @param object Initialized cisTopic object, after object@@binarized.rGREAT has been filled.
#' @param top Number of top terms per topic to visualize
#' @param topics Topics to plot. By default, all of them; if not a vector with the topic numbers must be provided.
#' @param ontology Ontology from which terms are selected to visualize.
#' @param var.y Varible (factor) to visualize in the y axis. By default it is ID, but it can be also name.
#' @param var.col Variable whose values will be visualized with the color scale
#' @param var.size Variable whose values will be visualized with the size scale
#' @param order.by Order terms in the y axis (descendent order) based on given variable
#' @param col.low Lower color in the color scale
#' @param col.mid Middle color in the color scale
#' @param col.high Higher color in the color scale
#' @param min.size Minimum dot size
#' @param max.size Maximum dot size
#'
#' @return Dot plot with the topics, terms and the selected enrichments to visualize
#'
#'
#' @import ggplot2
#' @import plyr
#' @importFrom data.table rbindlist
#'
#'
#' @export

ontologyDotPlot <- function(
  object,
  top=5,
  topics='all',
  ontology='GO Biological Process',
  var.y='ID',
  var.col='Binom_Adjp_BH',
  var.size='Binom_Fold_Enrichment',
  order.by='Binom_Fold_Enrichment',
  col.low = "brown1",
  col.mid = "floralwhite",
  col.high = "dodgerblue",
  min.size=2,
  max.size=10
)
{
  # Check dependencies
  if(! "data.table" %in% installed.packages()){
    stop('Please, install data.table: \n install.packages("data.table")')
  } else {
    require(data.table)
  }
  
  if(! "ggplot2" %in% installed.packages()){
    stop('Please, install data.table: \n install.packages("ggplot2")')
  } else {
    require(ggplot2)
  }
  
  if(order.by %in% c('Binom_Raw_PValue', 'Binom_Adjp_BH', 'Hyper_Raw_PValue', 'Hyper_Adjp_BH')){
    decreasing <- FALSE
  }
  else{
    decreasing <- TRUE
  }

  if (topics[1] == 'all'){
    topics <- which(laply(1:length(object@binarized.rGREAT), function(i) !is.null(object@binarized.rGREAT[[i]][[ontology]])) == TRUE)
    GOdata <- llply(topics, function(i) object@binarized.rGREAT[[i]][[ontology]])
    GOdata <- llply(1:length(GOdata), function(i) GOdata[[i]][order(GOdata[[i]][order.by], decreasing = decreasing),])
    GOdata <- llply(1:length(GOdata), function(i) GOdata[[i]][1:top,])
  }
  else{
    GOdata <- llply(topics, function(i) object@binarized.rGREAT[[i]][[ontology]])
    GOdata <- llply(1:length(GOdata), function(i) GOdata[[i]][order(GOdata[[i]][order.by], decreasing = decreasing),])
    GOdata <- llply(1:length(GOdata), function(i) GOdata[[i]][1:top,])
    names(GOdata) <- names(object@binarized.rGREAT)[topics]
  }

  names(GOdata) <- names(object@binarized.rGREAT)[topics]
  TopicID <- llply(1:length(GOdata), function(i) names(GOdata)[i])
  GOdata <- Map(cbind, GOdata, TopicID = TopicID)
  GOdata <- llply(1:length(GOdata), function(i) GOdata[[i]][!is.na(GOdata[[i]][,1]),])
  GOdata <- rbindlist(GOdata)
  GOdata <- GOdata[order(GOdata[[order.by]], decreasing = !decreasing),]
  GOdata[[var.y]] <- factor(GOdata[[var.y]], levels = unique(GOdata[[var.y]]))

  colorPal <- grDevices::colorRampPalette(c(col.low, col.mid, col.high))
  p <- ggplot(data = GOdata, mapping = aes_string(x = 'TopicID', y = var.y)) +
    geom_point(mapping = aes_string(size = var.size, color = var.col)) +
    scale_radius(range = c(min.size, max.size)) +
    scale_colour_gradientn(colors=colorPal(10)) +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))

  print(p)
}
