#' binarizedcisTopicsToCtx
#'
#' Convert binarized cisTopics to the corresponding cisTarget coordinates
#' @param object Initialized cisTopic object, after the object@@binarized.cisTopics has been filled.
#' @param genome Genome to which the data has been mapped. If you are using liftOver, provide the genome to which
#' you data has been lift-overed. The available genomes are hg19, dm3, dm6 and mm9.
#' @param liftOver GRangesList object containing the original coordinates (in the same format as
#' object@@region.names) in the data set as slot names and the corresponding mapping regions as a GRanges object in the slot.
#' @param ... See \code{findOverlap} from GenomicRanges
#'
#' @return Returns the cisTarget regions that correspond to the binarized topics in object@@binarized.regions.to.Rct
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @export

binarizedcisTopicsToCtx <- function(
  object,
  genome,
  liftOver = NULL,
  minOverlap = 0.4,
  ...
){
  # Check input
  if(length(object@binarized.cisTopics) < 1){
    stop('Please, run binarizecisTopics() first.')
  }
  
  if (genome == 'hg19'){
    data(hg19_CtxRegions)
    CtxRegions <- makeGRangesFromDataFrame(hg19_CtxRegions, keep.extra.columns = TRUE)
  }
  else if (genome == 'dm6'){
    data(dm6_CtxRegions)
    CtxRegions <- makeGRangesFromDataFrame(dm6_CtxRegions, keep.extra.columns = TRUE)
  }
  else if (genome == 'dm3'){
    data(dm3_CtxRegions)
    CtxRegions <- makeGRangesFromDataFrame(dm3_CtxRegions, keep.extra.columns = TRUE)
  }
  else if (genome == 'mm9'){
    data(mm9_CtxRegions)
    CtxLabel <- paste0('mm9_r70__', mm9_CtxRegions$seqnames, ':', mm9_CtxRegions$start+1, '-', mm9_CtxRegions$end)
    mm9_CtxRegions$CtxLabel <- CtxLabel
    CtxRegions <- makeGRangesFromDataFrame(mm9_CtxRegions, keep.extra.columns = TRUE)
  }
  else{
    stop('The genome required is not available! Try using the liftover option.')
  }

  if (is.null(liftOver)){
    object.binarized.regions.to.Rct <- llply(1:length(object@binarized.cisTopics), function(i) .getCtxRegions(object@region.ranges[rownames(object@binarized.cisTopics[[i]]),], CtxRegions, minOverlap = minOverlap, ...))
  } else {
    object.binarized.regions.to.Rct <- llply(1:length(object@binarized.cisTopics), function(i) .getCtxRegions(unlist(liftOver[rownames(object@binarized.cisTopics[[i]])], recursive = TRUE, use.names = TRUE), CtxRegions, minOverlap = minOverlap, ...))
  }
  names(object.binarized.regions.to.Rct) <- names(object@binarized.cisTopics)
  object@binarized.regions.to.Rct <- object.binarized.regions.to.Rct
  return(object)
}

# Helper functions

.getCtxRegions <- function(
  coordinates,
  CtxRegions,
  minOverlap=0.4,
  ...
){
  selectedCtxRegions <- unique(as.vector(.getOverlapRegionsFromCoordinates(coordinates, CtxRegions, minOverlap=minOverlap, ...)[,2]))
  selectedCtxLabels <- as.vector(CtxRegions[selectedCtxRegions]@elementMetadata[,"CtxLabel"])
  message(paste('Number of regions selected:', length(selectedCtxLabels)))
  return(selectedCtxLabels)
}

#' scoredRegionsToCtx
#'
#' Map regions to the most overlapping cisTarget region
#' @param object Initialized cisTopic object, after the object@@binarized.cisTopics has been filled.
#' @param genome Genome to which the data has been mapped. The available genomes are hg19, dm3, dm6 and mm9.
#' @param liftOver GRangesList object containing the original coordinates (in the same format as
#' object@@region.names) in the data set as slot names and the corresponding mapping regions as a GRanges object in the slot.
#' @param ... See \code{findOverlap} from GenomicRanges
#'
#' @return Returns the region scores per topic in object@@region.data
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @export

scoredRegionsToCtx <- function(
  object,
  genome,
  liftOver = NULL, 
  minOverlap = 0.4,
  ...
){
  if (genome == 'hg19'){
    data(hg19_CtxRegions)
    CtxRegions <- makeGRangesFromDataFrame(hg19_CtxRegions, keep.extra.columns = TRUE)
  }
  else if (genome == 'dm6'){
    data(dm6_CtxRegions)
    CtxRegions <- makeGRangesFromDataFrame(dm6_CtxRegions, keep.extra.columns = TRUE)
  }
  else if (genome == 'dm3'){
    data(dm3_CtxRegions)
    CtxRegions <- makeGRangesFromDataFrame(dm3_CtxRegions, keep.extra.columns = TRUE)
  }
  else if (genome == 'mm9'){
    data(mm9_CtxRegions)
    CtxLabel <- paste0('mm9_r70__', mm9_CtxRegions$seqnames, ':', mm9_CtxRegions$start+1, '-', mm9_CtxRegions$end)
    mm9_CtxRegions$CtxLabel <- CtxLabel
    CtxRegions <- makeGRangesFromDataFrame(mm9_CtxRegions, keep.extra.columns = TRUE)
  }

  if (is.null(liftOver)){
    coordinates <- object@region.ranges
  } else {
    coordinates <- unlist(liftOver, recursive = TRUE, use.names = TRUE)
  }

  selectedRegions <- .getOverlapRegionsFromCoordinates(coordinates, CtxRegions, minOverlap=minOverlap, overlapping = FALSE, ...)
  selectedCtxLabels <- as.vector(CtxRegions[as.vector(selectedRegions[,2])]@elementMetadata[,"CtxLabel"])
  names(selectedCtxLabels) <- as.vector(selectedRegions[,1])
  object@region.data$CtxLabels <- selectedCtxLabels[object@region.names]
  message(paste('Number of regions selected:', length(selectedCtxLabels)))
  return(object)
}

#' topicsRcisTarget
#'
#' Run RcisTarget in the binarized topics
#' @param object Initialized cisTopic object, after the object@@binarized.regions.to.Rct has been filled.
#' @param pathToFeather Path to the feather database to use. Note that this database has to match the genome used for mapping.
#' @param genome Genome to which the data was aligned or liftovered (hg19, mm9, dm3 or dm6).
#' @param motifRankings Feather database corresponding to the genome
#' @param reduced_database Whether a reduced version of the database (e.g. background database) is being used or not (by default, it is set to false).
#' @param nesThreshold Minimum enrichment score that a motif should have to be included in the table
#' @param rocthr ROC threshold for AUC calculation. For mouse and human, we recommend 0.005; for fly, 0.01.
#' @param maxRank Number of regions consider for visualization. For mouse and human, we recommend 20000; for fly, 5000.
#' @param ... See RcisTarget
#'
#' @return Motif enrichment table is stored in object@@binarized.RcisTarget
#' @examples
#' pathToFeather <- "hg19-regions-1M-9species.all_regions.mc9nr.feather"
#' cisTopicObject <- topicRcisTarget(cisTopicObject, genome='hg19', pathToFeather)
#' cisTopicObject
#' @import RcisTarget
#' @importFrom parallel makeCluster
#' @import doSNOW
#' @import feather
#' @importFrom plyr llply
#' @export

topicsRcisTarget <- function(
  object,
  genome,
  pathToFeather,
  reduced_database = FALSE,
  nesThreshold = 3,
  rocthr = 0.005,
  maxRank = 20000,
  nCores = 1,
  ...
){
  # Check input
  if(length(object@binarized.regions.to.Rct) < 1){
    stop('Please, run binarizedcisTopicsToCtx() first.')
  }
  
  # Check dependencies
  if(! "RcisTarget" %in% installed.packages()){
    stop('Please, install RcisTaregt: \n install_github("aertslab/RcisTarget")')
  } else {
    require(RcisTarget)
  }
  
  if (genome == 'hg19'){
    data(motifAnnotations_hgnc)
    motifAnnot <- motifAnnotations_hgnc
    if(rocthr!=0.005 | maxRank!=20000){
      warning("For Homo sapiens the recommended settings are: rocthr=0.005, maxRank=20000")
    } 
  }
  else if (genome == 'mm9'){
    data(motifAnnotations_mgi)
    motifAnnot <- motifAnnotations_mgi
    if(rocthr!=0.005 | maxRank!=20000){
      warning("For Mus musculus the recommended settings are: rocthr=0.005, maxRank=20000")
    } 
  }
  else if (genome == 'dm3'){
    data(motifAnnotations_dmel)
    motifAnnot <- motifAnnotations_dmel
    if(rocthr!=0.01 | maxRank!=5000){
      warning("For Drosophila melanogaster the recommended settings are: rocthr=0.01, maxRank=5000")
    } 
  }
  else if (genome == 'dm6'){
    data(motifAnnotations_dmel)
    motifAnnot <- motifAnnotations_dmel
    if(rocthr!=0.01 | maxRank!=5000){
      warning("For Drosophila melanogaster the recommended settings are: rocthr=0.01, maxRank=5000")
    } 
  } else {
    stop('The genome required is not available! Try using the liftover option.')
  }

  topicsList <- object@binarized.regions.to.Rct

  if (reduced_database == FALSE){
    ctxreg <- unique(as.vector(unlist(object@binarized.regions.to.Rct)))
    motifRankings <- importRankings(pathToFeather, columns = c('features', ctxreg))
  }
  else{
    motifRankings <- importRankings(pathToFeather)
    ctxregions <- colnames(getRanking(motifRankings))[-1]
    topicsList <- llply(1:length(topicsList), function(i) topicsList[[i]][which(topicsList[[i]] %in% ctxregions)])
    names(topicsList) <- names(object@binarized.regions.to.Rct)
  }
  
    columnsinRanking <- feather_metadata(pathToFeather)[[2]][2]
  
  if (length(topicsList) < nCores){
    print(paste('The number of cores (', nCores, ') is higher than the number of topics (', topic,').', sep=''))
  }

  if(nCores > 1){
    cl <- makeCluster(nCores, type = "SOCK")
    registerDoSNOW(cl)
    print(paste('Exporting data to clusters...'))
    clusterEvalQ(cl, library(RcisTarget))
    clusterExport(cl, c("topicsList", "motifRankings", "motifAnnot", "nesThreshold", "rocthr", "columnsinRanking", "maxRank"), envir=environment())
    print(paste('Running RcisTarget...'))
    cisTopic.cisTarget <- suppressWarnings(llply(1:length(topicsList), function (i) cisTarget(topicsList[[i]],
                                    motifRankings,
                                    motifAnnot = motifAnnot,
                                    nesThreshold = nesThreshold,
                                    aucMaxRank = rocthr * columnsinRanking,
                                    geneErnMmaxRank = maxRank,
                                    nCores=1
    ), .parallel = TRUE))
    stopCluster(cl)
  }
  else{
    cisTopic.cisTarget <- suppressWarnings(llply(1:length(topicsList), function (i) cisTarget(topicsList[[i]],
                                                                                             motifRankings,
                                                                                             motifAnnot = motifAnnot,
                                                                                             nesThreshold = nesThreshold,
                                                                                             aucMaxRank = rocthr * columnsinRanking,
                                                                                             geneErnMmaxRank = maxRank,
                                                                                             nCores=1
    )))
  }

  object.binarized.RcisTarget <- list()

  for (i in 1:length(cisTopic.cisTarget)){
    if(nrow(cisTopic.cisTarget[[i]]) > 0){
      colnames(cisTopic.cisTarget[[i]])[c(1, 7, 9)] <- c('cisTopic', 'nEnrRegions', 'enrichedRegions')
      cisTopic.cisTarget[[i]]$cisTopic <- rep(paste('Topic', i, sep='_'), nrow(cisTopic.cisTarget[[i]]))
      object.binarized.RcisTarget[[i]] <- addLogo(cisTopic.cisTarget[[i]])
    } else {
      cisTopic.cisTarget[[i]] <- NULL
    }
    
  }

  object@binarized.RcisTarget <- object.binarized.RcisTarget
  return(object)
}

#' makeBackgroundFeather
#'
#' Create a background cisTarget ranking out of the genome-wide cisTarget regions based on given regions.
#' @param pathToFeather Path to the feather database to use. Note that this database has to match the genome used for mapping.
#' @param subsetRegions Subset of regions to be re-ranked (optional).
#' @param subsetMotifs Subset of motifs to be used from the database (optional).
#' @param pathToSave Path to save background database.
#' @param ... See RcisTarget
#'
#' @return A file containing the feather database.
#' 
#' @details The created database can be used as input for \code{topicsRcisTarget()}; using reduced_database = TRUE. Subsetting by 
#' regions and motifs is optional, but at least one of them as to be filled.
#' @examples
#' pathToFeather <- "hg19-regions-1M-9species.all_regions.mc9nr.feather"
#' subsetRegion <- cisTopicObject@@region.data$CtxLabels
#' pathToSave <- "hg19-regions-subsetMelanoma.feather"
#' makeBackgroundFeather(pathToFeather, subsetRegion, pathToSave)
#' @import feather
#' @export

makeBackgroundFeather <- function(
  pathToFeather,
  subsetRegions = NULL,
  subsetMotifs = NULL,
  pathToSave,
  ...
){
  # Check dependencies
  if(! "tibble" %in% installed.packages()){
    stop('Please, install tibble: \n install.packages("tibble")')
  }
  
  
  if (is.null(subsetRegions) & is.null(subsetMotifs)){
    stop('Please, introduce regions and/or motifs to subset for.')
  }
  
  motifRankings <- feather(pathToFeather)
  motifRankings <- motifRankings[which(motifRankings$features %in% subsetMotifs), c(1, which(colnames(motifRankings) %in% subsetRegions))]
  reRankedMotifRankings <- tibble::as_tibble(t(apply(motifRankings[,2:ncol(motifRankings)], 1, function(x) as.integer(rank(x)))))
  colnames(reRankedMotifRankings) <- colnames(motifRankings)[-1]
  reRankedMotifRankings <- tibble::add_column(reRankedMotifRankings, features=motifRankings$features, .before=1)
  write_feather(reRankedMotifRankings, pathToSave)
}



#' getCistromes
#'
#' Get cistromes formed based on motif enrichment with RcisTarget
#' @param object Initialized cisTopic object, after object@@binarized.RcisTarget and object@@region.data$CtxLabels have been filled.
#' @param annotation Annotations to be used ('TF_highConf', 'Both'). By default, only the both high confidence and indirect annotation is used.
#' @param gene.cistrome Whether the cistromes with gene linked gene symbols have to be formed (based on the closest gene). It requires to 
#' run first annotateRegions().
#' @param region.cistrome Whether the cistromes with regions in the data set linked to the ctx regions have to be formed (based on maximum overlap). It requires to 
#' run first scoredRegionstoCtx().
#' @param nCores Number of cores to be used (by default 1).
#' @param ... Ignored
#'
#' @return Cistromes are stored as ctx regions (object@@cistromes.ctx), regions (object@@cistromes.regions) and if specified, as gene symbols
#' (object@@cistromes.genes). Cistromes containing extended in their name are formed by both the high confidence and indirect annotation; otherwise,
#' only by the high confidence features.
#' 
#' @examples
#'
#' cisTopicObject <- getCistromes(cisTopicObject, annotation = 'Both', gene.cistrome=FALSE, nCores=1)
#' cisTopicObject
#'
#' @importFrom parallel makeCluster
#' @import doSNOW
#' @importFrom data.table rbindlist
#' @importFrom plyr llply
#' @export

getCistromes <- function(
  object,
  annotation = 'Both',
  gene.cistrome = FALSE,
  region.cistrome=TRUE,
  nCores = 1,
  ...
){
  # Check input
  if(length(object@binarized.RcisTarget) < 1){
    stop('Please, run topicsRcisTarget() first.')
  }
  
  object.binarized.RcisTarget <- object@binarized.RcisTarget
  if(nCores > 1){
    cl <- makeCluster(nCores, type = "SOCK")
    registerDoSNOW(cl)
    print(paste('Exporting data to clusters...'))
    clusterEvalQ(cl, library("data.table"))
    clusterExport(cl, c("object.binarized.RcisTarget", "annotation", ".onecisTopicGetCtxCistromes"), envir=environment())
    print(paste('Annotating...'))
    object.cistromes.ctx <- suppressWarnings(llply(1:length(object.binarized.RcisTarget), function (i) .onecisTopicGetCtxCistromes(object.binarized.RcisTarget[[i]], annotation = annotation
    ), .parallel = TRUE))
    stopCluster(cl)
  }
  else{
    object.cistromes.ctx <- suppressWarnings(llply(1:length(object.binarized.RcisTarget), function (i) .onecisTopicGetCtxCistromes(object.binarized.RcisTarget[[i]], annotation = annotation
    )))
  }

  object@cistromes.ctx <- object.cistromes.ctx
  
  if(region.cistrome){
    object.cistromes.regions <- object.cistromes.ctx
    
    for (i in 1:length(object.binarized.RcisTarget)){
      for (j in 1:length(object.cistromes.ctx[[i]])){
        object.cistromes.regions[[i]][[j]] <- as.vector(unlist(object.cistromes.ctx[[i]][[j]]))
        object.cistromes.regions[[i]][[j]] <- unique(as.vector(unlist(rownames(object@region.data[which(object@region.data$CtxLabels %in% object.cistromes.ctx[[i]][[j]]),]))))
        object.cistromes.regions[[i]][[j]] <- object.cistromes.regions[[i]][[j]][!is.na(object.cistromes.regions[[i]][[j]])]
      }
      names(object.cistromes.regions[[i]]) <- sapply(strsplit(names(object.cistromes.regions[[i]]), split = " (",  fixed = TRUE), "[", 1)
      names(object.cistromes.regions[[i]]) <- paste(names(object.cistromes.regions[[i]]), " (",lengths(object.cistromes.regions[[i]]), "p)", sep="")
    }
    
    object@cistromes.regions <- object.cistromes.regions
  } 
  else{
    print('Region cistromes have not been computed. Please, run scoredRegionstoCtx() first.')
  }
  
  if (gene.cistrome){
    if ('SYMBOL' %in% colnames(object@region.data)){
      object.cistromes.genes <- object.cistromes.regions
      
      for (i in 1:length(object.binarized.RcisTarget)){
        for (j in 1:length(object.cistromes.regions[[i]])){
          object.cistromes.genes[[i]][[j]] <- as.vector(unlist(object.cistromes.regions[[i]][[j]]))
          object.cistromes.genes[[i]][[j]] <- unique(as.vector(unlist(object@region.data[object.cistromes.regions[[i]][[j]], 'SYMBOL'])))
          object.cistromes.genes[[i]][[j]] <- object.cistromes.genes[[i]][[j]][!is.na(object.cistromes.genes[[i]][[j]])]
        }
        names(object.cistromes.genes[[i]]) <- sapply(strsplit(names(object.cistromes.genes[[i]]), split = " (",  fixed = TRUE), "[", 1)
        names(object.cistromes.genes[[i]]) <- paste(names(object.cistromes.genes[[i]]), " (",lengths(object.cistromes.genes[[i]]), "g)", sep="")
      }
      
      object@cistromes.genes <- object.cistromes.genes
    }
    else{
      print('Gene cistromes have not been computed. Please, run annotateRegions() first.')
    }
  }
  return(object)
}


# Helper functions.

# Return cistromes based on i-cisTarget regions.
.onecisTopicGetCtxCistromes <- function(
  motifEnrichmentTable,
  annotation='TF_highConf'
  ){
  if (annotation == 'TF_highConf'){
    motifEnrichment.asIncidList <- apply(motifEnrichmentTable, 1, function(oneMotifRow) {
      peaks <- strsplit(oneMotifRow["enrichedRegions"], ";")[[1]]
      TFs <- strsplit(oneMotifRow["TF_highConf"], "; ")[[1]]
      oneMotifRow <- data.frame(rbind(oneMotifRow), stringsAsFactors=FALSE)
      oneMotifRow <-  data.frame(oneMotifRow[rep(1, length(TFs)),c("NES", "motif")],  TFs, stringsAsFactors = FALSE)
      frame <- data.frame(oneMotifRow[rep(1, length(peaks)),c("NES", "motif", "TFs")], rep("TF_highConf", length(peaks)), peaks, stringsAsFactors = FALSE)
      if(nrow(oneMotifRow) > 1){
        for (i in 2:nrow(oneMotifRow)){
          frame <- rbind(frame, data.frame(oneMotifRow[rep(i, length(peaks)),c("NES", "motif", "TFs")], rep("TF_highConf", length(peaks)), peaks, stringsAsFactors = FALSE))
        }
      }
      frame
    })
  }
  else {
    motifEnrichment.asIncidList_direct <- apply(motifEnrichmentTable, 1, function(oneMotifRow) {
      peaks <- strsplit(oneMotifRow["enrichedRegions"], ";")[[1]]
      TFs <- strsplit(oneMotifRow["TF_highConf"], "; ")[[1]]
      oneMotifRow <- data.frame(rbind(oneMotifRow), stringsAsFactors=FALSE)
      oneMotifRow <-  data.frame(oneMotifRow[rep(1, length(TFs)),c("NES", "motif")],  TFs, stringsAsFactors = FALSE)
      frame <- data.frame(oneMotifRow[rep(1, length(peaks)),c("NES", "motif", "TFs")], rep("TF_highConf", length(peaks)), peaks, stringsAsFactors = FALSE)
      if(nrow(oneMotifRow) > 1){
        for (i in 2:nrow(oneMotifRow)){
          frame <- rbind(frame, data.frame(oneMotifRow[rep(i, length(peaks)),c("NES", "motif", "TFs")], rep("TF_highConf", length(peaks)), peaks, stringsAsFactors = FALSE))
        }
      }
      frame
    })
    motifEnrichment.asIncidList_inferred <- apply(motifEnrichmentTable, 1, function(oneMotifRow) {
      peaks <- strsplit(oneMotifRow["enrichedRegions"], ";")[[1]]
      TFs <- strsplit(oneMotifRow["TF_lowConf"], "; ")[[1]]
      oneMotifRow <- data.frame(rbind(oneMotifRow), stringsAsFactors=FALSE)
      oneMotifRow <-  data.frame(oneMotifRow[rep(1, length(TFs)),c("NES", "motif")],  TFs, stringsAsFactors = FALSE)
      frame <- data.frame(oneMotifRow[rep(1, length(peaks)),c("NES", "motif", "TFs")], rep("TF_lowConf", length(peaks)), peaks, stringsAsFactors = FALSE)
      if(nrow(oneMotifRow) > 1){
        for (i in 2:nrow(oneMotifRow)){
          frame <- rbind(frame, data.frame(oneMotifRow[rep(i, length(peaks)),c("NES", "motif", "TFs")], rep("TF_lowConf", length(peaks)), peaks, stringsAsFactors = FALSE))
        }
      }
      frame
    })
    motifEnrichment.asIncidList <- rbind(motifEnrichment.asIncidList_direct, motifEnrichment.asIncidList_inferred)
  }


  motifEnrichment.asIncidList <- rbindlist(motifEnrichment.asIncidList)
  colnames(motifEnrichment.asIncidList) <- c("NES", "motif", "TF", "annot","region")
  motifEnrichment.asIncidList$TF <- gsub(' (.*).', '', motifEnrichment.asIncidList$TF)
  motifEnrichment.asIncidList <- data.frame(motifEnrichment.asIncidList, stringsAsFactors = FALSE)

  # Get targets for each TF, but keep info about best motif/enrichment
  # (directly annotated motifs are considered better)
  cistromeTargetsInfo <- lapply(split(motifEnrichment.asIncidList, motifEnrichment.asIncidList$TF), function(tfTargets){
    tfTable <- as.data.frame(do.call(rbind, lapply(split(tfTargets, tfTargets$region), function(enrOneGene){
      TF_highConf <- "TF_highConf" %in% enrOneGene$annot
      enrOneGeneByAnnot <- enrOneGene
      if(annotation == 'TF_highConf') enrOneGeneByAnnot <- enrOneGeneByAnnot[which(enrOneGene$annot == "TF_highConf"),]
      bestMotif <- which.max(enrOneGeneByAnnot$NES)

      cbind(TF=unique(enrOneGene$TF), region=unique(enrOneGene$region), nMotifs=nrow(enrOneGene),
            bestMotif=as.character(enrOneGeneByAnnot[bestMotif,"motif"]), NES=as.numeric(enrOneGeneByAnnot[bestMotif,"NES"]), TF_highConf=TF_highConf)
    })), stringsAsFactors=FALSE)
    tfTable[order(tfTable$NES, decreasing = TRUE),]
  })

  rm(motifEnrichment.asIncidList)
  cistromeTargetsInfo <- rbindlist(cistromeTargetsInfo)
  colnames(cistromeTargetsInfo) <- c("TF", "region", "nMotifs", "bestMotif", "NES", "TF_highConf")

  # Split into cistromes
  cistromeTargetsInfo_splitByAnnot <- split(cistromeTargetsInfo, cistromeTargetsInfo$TF_highConf)
  cistromes <- sapply(split(cistromeTargetsInfo_splitByAnnot[["TRUE"]], cistromeTargetsInfo_splitByAnnot[["TRUE"]][,"TF"]), function(x) sort(as.character(unlist(x[,"region"]))))
  if (annotation == 'Both'){
    cistromes_extended <- sapply(split(cistromeTargetsInfo_splitByAnnot[["FALSE"]],cistromeTargetsInfo_splitByAnnot[["FALSE"]][,"TF"]), function(x) unname(x[,"region"]))
    cistromes_extended[which(names(cistromes_extended) %in% names(cistromes))] <- lapply(names(cistromes_extended)[which(names(cistromes_extended) %in% names(cistromes))], function(tf) sort(unique(c(cistromes[[tf]], cistromes_extended[[tf]]))))
    names(cistromes_extended) <- paste(names(cistromes_extended), "_extended", sep="")
    
    if (class(cistromes) == 'matrix'){
      cistromes_extended[[colnames(cistromes)]] <- as.vector(cistromes)
      cistromes <- cistromes_extended
    } else {
      cistromes <- c(cistromes, cistromes_extended)
    }
  }
  names(cistromes) <- paste(names(cistromes), " (",lengths(cistromes), "r)", sep="")
  return(cistromes)
}

#' getCistromeEnrichment
#'
#' Determine topic-specific cistrome enrichment in the cells. If specified, cistrome enrichment can be plotted.
#' @param object Initialized cisTopic object, after object@@cistromes.regions have been filled (see \code{getCistromes()}).
#' @param topic Topic number of the cistrome.
#' @param TFname Name of the transcription factor linked to the cistrome.
#' @param annotation Annotations to be used ('TF_highConf', 'Both'). By default, only the high confidence annotation is used.
#' @param aucellRankings Precomputed aucellRankings using \code{cisTopic_buildRankings()}. These rankings are not stored in the cisTopicObject due to their size.
#' @param nCores Number of cores to be used for AUCell
#' @param aucMaxRank Threshold to calculate the AUC
#' @param plot Whether enrichment plot should be done. If yes, parameters for plotFeatures will not be ignored.
#' @param ... See \code{plotFeatures()}.
#'
#' @return AUC enrichment values for the signature are stored as a column in object@@cell.data. If specified, cells coloured by
#' their AUC enrichment values will be plotted.
#' @examples
#'
#' cisTopicObject <- getCistromeEnrichment(cisTopicObject, annotation = 'Both', nCores=1)
#' cisTopicObject
#'
#' @import AUCell
#' @export

getCistromeEnrichment <- function(
  object,
  topic,
  TFname,
  annotation = 'TF_highConf',
  aucellRankings,
  nCores = 1,
  aucMaxRank = 0.03*nrow(aucellRankings),
  plot=TRUE,
  ...
){
  topicCistromes <- cisTopicObject@cistromes.regions[[topic]]
  
  if (is.null(grep(TFname, names(topicCistromes)))){
    stop(paste0('The specified cistrome cannot be found. Please, check whether there is a cistrome for ', TFname, ' in topic ', topic, ' with ', annotation, '.'))
  }
  else{
    cistromeNames <- names(topicCistromes)[grep(TFname, names(topicCistromes))]
    if (annotation == 'TF_highConf'){
      cistromeNames <- cistromeNames[grep(paste0(TFname,' '), cistromeNames)]
    }
    else if (annotation == 'Both'){
      cistromeNames <- cistromeNames[grep(paste0(TFname,'_'), cistromeNames)]
    }
  }
  
  cistrome <- topicCistromes[cistromeNames]
  
  modulesAUC <- AUCell_calcAUC(cistrome, aucellRankings, nCores=nCores, aucMaxRank=aucMaxRank)
  enrichMatrix <- t(getAUC(modulesAUC))
  rownames(enrichMatrix) <- object@cell.names
  colnames(enrichMatrix) <- paste0('Topic', topic, '_', colnames(enrichMatrix))
  object <- addCellMetadata(object, as.data.frame(enrichMatrix))
  
  if (plot){
    plotFeatures(object, target='cell', colorBy=colnames(enrichMatrix))
  }
  
  return(object)
}
  
