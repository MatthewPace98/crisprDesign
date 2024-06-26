#' @title Add CFD and MIT scores to a \linkS4class{GuideSet} object. 
#' @description Add CFD and MIT off-target scores to a 
#'     \linkS4class{GuideSet} object. 
#'     Both the CFD and MIT methods are available for the SpCas9 nuclease.
#'     The CFD method is also available for the CasRx nuclease.
#'     Other nucleases are currently not supported. 
#' 
#' @param object A \linkS4class{GuideSet} object or a 
#'     \linkS4class{PairedGuideSet} object.
#'     \code{crisprNuclease(object)} must be either using SpCas9 or CasRx.
#' @param max_mm The maximimum number of mismatches between the spacer sequence
#'     and the protospacer off-target sequence to be considered in the 
#'     off-target score calculations. Off-targets with a number of 
#'     mismatches greater than \code{max_mm} will be excluded; this is useful
#'     if one wants to avoid the aggregated off-target scores to be driven by a
#'     large number of off-targets that have low probability of cutting.
#' @param includeDistance Should a distance penalty for the MIT score be
#'     included? TRUE by default. 
#' @param offset Numeric value specifying an offset to add to the denominator
#'     when calcuting the aggregated score (inverse summation formula).
#'     0 by default.
#' @param ... Additional arguments, currently ignored.
#' 
#' @return A \code{GuideSet} or a \code{PairedGuideSet} object with added 
#'     scores. The alignments annotation returned by \code{alignments(object)}
#'     will have additional column storing off-target scores. Those scores
#'     representing the off-target score for each gRNA and off-target pair.
#'     For SpCas9, a column containing an aggregated specificity off-target 
#'     score for each scoring method is added to the metadata columns 
#'     obtained by \code{mcols(object)}.
#' 
#' @details See the \pkg{crisprScore} package for a description of the 
#'     different off-target scoring methods.
#' 
#' @examples
#' 
#' data(guideSetExampleWithAlignments, package="crisprDesign")
#' gs <- guideSetExampleWithAlignments
#' gs <- addOffTargetScores(gs)
#' 
#' @author Jean-Philippe Fortin, Luke Hoberecht
#' 
#' @seealso \code{link{addOnTargetScores}} to add on-target scores.
#' 
#' @export
#' @rdname addOffTargetScores
#' @importFrom S4Vectors split mcols<-
setMethod("addOffTargetScores",
          "GuideSet",
          function(object,
                   max_mm=2,
                   includeDistance=TRUE,
                   offset=0
){
    
    object <- .validateGuideSet(object)
    .checkOffTargetScoresParameters(guideSet=object,
                                    max_mm=max_mm,
                                    includeDistance=includeDistance,
                                    offset=offset)
    if (!.hasUniqueSpacers(object)){
        warning("There are duplicated gRNAs in the GuideSet object, which",
                " can lead to incorrect off-target scores. Consider removing ",
                "duplicated gRNAs first before running addOffTargetScores. ")
    }
    
    object <- .addOffTargetScoresToAlignments(guideSet=object,
                                              includeDistance=includeDistance)
    object <- .addOffTargetScoresToGuideSet(guideSet=object,
                                              max_mm=max_mm,
                                              offset=offset)
    return(object)
})

#' @rdname addOffTargetScores
#' @export
setMethod("addOffTargetScores",
          "PairedGuideSet", 
          function(object,
                   max_mm=2,
                   includeDistance=TRUE,
                   offset=0
){
    object <- .validatePairedGuideSet(object)
    unifiedGuideSet <- .pairedGuideSet2GuideSet(object)
    unifiedGuideSet <- addOffTargetScores(unifiedGuideSet,
                                          max_mm=2,
                                          includeDistance=TRUE,
                                          offset=0)
    out <- .addColumnsFromUnifiedGuideSet(object,
                                          unifiedGuideSet)
    
    return(out)
})


#' @rdname addOffTargetScores
#' @export
setMethod("addOffTargetScores", "NULL", function(object){
    return(NULL)
})


#' @importFrom S4Vectors mcols isTRUEorFALSE
.checkOffTargetScoresParameters <- function(guideSet,
                                            max_mm,
                                            includeDistance,
                                            offset
){
    if (!"alignments" %in% colnames(S4Vectors::mcols(guideSet))){
        stop("Alignments must be added to guideSet prior to calculating ",
             "off-target scores; see ?addSpacerAlignments")
    }
    stopifnot("'max_mm' must be a single non-negative integer value" = {
        is.vector(max_mm, mode="numeric") &&
            length(max_mm) == 1 &&
            max_mm >= 0 &&
            max_mm %% 1 == 0
    })
    stopifnot("'includeDistance' must be TRUE or FALSE" = {
        S4Vectors::isTRUEorFALSE(includeDistance)
    })
    stopifnot("'offset' must be a single non-negative numeric value" = {
        is.vector(offset, mode="numeric") &&
            length(offset) == 1 &&
            offset >= 0
    })
    invisible(NULL)
}


# Add off-target scores for each off-target stored in the alignments data.frame
#' @importFrom utils data
#' @importFrom S4Vectors split mcols<-
#' @importFrom crisprScore getCFDScores
#' @importFrom crisprScore getMITScores
.addOffTargetScoresToAlignments <- function(guideSet,
                                            includeDistance=TRUE
){
    crisprNuclease <- crisprNuclease(guideSet)
    utils::data(SpCas9,
                package="crisprBase",
                envir=environment())
    utils::data(CasRx,
                package="crisprBase",
                envir=environment())
    isCas9  <- .identicalNucleases(crisprNuclease,
                                   SpCas9,
                                   checkSpacerLength=FALSE)
    isCasRx <- .identicalNucleases(crisprNuclease,
                                   CasRx,
                                   checkSpacerLength=FALSE)
    if (!isCas9 & !isCasRx){
        stop("Nuclease must be either SpCas9 or CasRx ",
             "for off-target scoring")
    }
    spacerLen <- spacerLength(crisprNuclease)
    if (isCasRx & spacerLen>27){
        stop("For CasRx, spacer length must be at most 27nt ",
             "for off-target scoring.")
    }
    if (isCas9 & spacerLen>20){
        stop("For SpCas9, spacer length must be at most 20nt ",
             "for off-target scoring.")
    }
   
    aln <- alignments(guideSet)
    spacers      <- as.character(aln$spacer)
    protospacers <- as.character(aln$protospacer)
    pams         <- as.character(aln$pam)
    if (isCasRx){
        protospacers <- DNAStringSet(protospacers)
        protospacers <- reverseComplement(protospacers)
        protospacers <- as.character(protospacers)
    }

    if (isCas9 & spacerLen == 19){
        spacers      <- paste0("G", spacers, recycle0=TRUE)
        protospacers <- paste0("G", protospacers, recycle0=TRUE)
    }
    if (isCas9){
        nuclease <- "SpCas9"
    } else if (isCasRx){
        nuclease <- "CasRx"
    }

    if (isCas9 | isCasRx){
        score_cfd <- crisprScore::getCFDScores(spacers=spacers,
                                               protospacers=protospacers,
                                               pams=pams,
                                               nuclease=nuclease)
        aln$score_cfd <- score_cfd$score
    }
    if (isCas9){
        score_mit <- crisprScore::getMITScores(spacers=spacers,
                                               protospacers=protospacers,
                                               pams=pams,
                                               includeDistance=includeDistance)
        aln$score_mit <- score_mit$score
        
  # CRISTA
  run_crista <- FALSE  # CRISTA is currently too inefficient, and deprecated by its authors 
  if (run_crista == TRUE){
     
  extendedSequences <- .getExtendedSequences(guideSet,
                                             start=-22,
                                             end=6)
  good <- !is.na(extendedSequences)
  scores <- rep(NA, length(extendedSequences))
  seqs <- extendedSequences[good]

  # Generates spacer list to match protospacer list length
  i <- 1  
  crista_spacers <- vector()
  previous_spacer <- spacers[1]  
  for (j in 1:length(spacers)) {
    if (spacers[j] != previous_spacer) {
      i <- i + 1  
    }    
    previous_spacer <- spacers[j]
    crista_spacers[j] <- seqs[i]
  }
  
  crista_protospacers <- as.character(aln$protospacer)
  #results <- crisprScore::getCRISTAScores(protospacer=crista_protospacers, 
  #                                        spacer=crista_spacers)
  
  score_crista <- results$score
  aln$score_crista <- score_crista
  }
}
  guideSetSpacers <- spacers(guideSet, as.character=TRUE)
  aln <- S4Vectors::split(aln,
                          f=factor(aln$spacer,
                          levels=unique(guideSetSpacers)))
  aln <- aln[guideSetSpacers]
  names(aln) <- names(guideSet)
  S4Vectors::mcols(guideSet)[["alignments"]] <- aln
          
    
  return(guideSet)
}
          

.getExtendedSequences <- function(guideSet,
                                  start,
                                  end
){
  guideSet <- .validateGuideSet(guideSet)
  
  gr <- guideSet
  wh_neg <- which(as.character(strand(gr))=="-")
  # The order of resizing IRanges matters
  # to presever the validity of a positive width.
  if (start>0 & end>0){
    end(gr)   <- end(guideSet)+end
    start(gr) <- start(guideSet)+start
    start(gr)[wh_neg] <- start(guideSet)[wh_neg]-end
    end(gr)[wh_neg]   <- end(guideSet)[wh_neg]-start
  } else {
    start(gr) <- start(guideSet)+start
    end(gr)   <- end(guideSet)+end
    end(gr)[wh_neg]   <- end(guideSet)[wh_neg]-start
    start(gr)[wh_neg] <- start(guideSet)[wh_neg]-end
  }
  
  gr <- GenomicRanges::trim(gr) #Taking care of invalid values
  
  good <- which(as.character(strand(gr)) %in% c("+", "-"))
  out <- rep(NA_character_, length(gr))
  names(out) <- names(gr)
  if (length(good)==0){
    return(out)
  } 
  if (targetOrigin(guideSet)=="customSequences"){
    seqs <- getSeq(customSequences(guideSet),gr[good])
  } else {
    seqs <- getSeq(bsgenome(guideSet), gr[good])
  }
  seqs <- as.character(seqs)
  
  #Making sure the sequences are not out of bound:
  len <- end-start+1 # Expected length
  seqs[seqs==""] <- NA
  seqs[nchar(seqs)<len] <- NA
  out[good] <- seqs
  return(out)
}

# Add aggregated off-target scores to the GuideSet
#' @importFrom S4Vectors mcols mcols<-
.addOffTargetScoresToGuideSet <- function(guideSet,
                                          max_mm,
                                          offset
){
    aln <- alignments(guideSet)
    aln <- as.data.frame(S4Vectors::mcols(aln),
                         stringsAsFactors=FALSE)
    validMismatchCount <- aln$n_mismatches <= max_mm
    aln <- aln[validMismatchCount, , drop=FALSE]
    aln <- split(aln, f=aln$spacer)

    .getAggregateScore <- function(score){
        vapply(aln, function(x){
            nmm <- x[["n_mismatches"]]
            x <- x[[score]]
            if (sum(nmm==0)==0){
                x <- sum(x, na.rm=TRUE) + 1
            } else {
                x <- sum(x, na.rm=TRUE) + offset
            }
            return(1/x)
        }, FUN.VALUE=numeric(1))
    }
    
    spacers <- spacers(guideSet, as.character=TRUE)
    crisprNuclease <- crisprNuclease(guideSet)
    utils::data(SpCas9,
                package="crisprBase",
                envir=environment())
    utils::data(CasRx,
                package="crisprBase",
                envir=environment())
    isCas9  <- .identicalNucleases(crisprNuclease,
                                   SpCas9,
                                   checkSpacerLength=FALSE)
    isCasRx <- .identicalNucleases(crisprNuclease,
                                   CasRx,
                                   checkSpacerLength=FALSE)


    if (isCasRx){
        stop("Aggregation is not implemented yet for CasRx")
    }
    if (isCas9){
        cfd <- .getAggregateScore("score_cfd")
        cfdIndices <- match(spacers, names(cfd))
        S4Vectors::mcols(guideSet)[["score_cfd"]] <- cfd[cfdIndices]
    } 
    if (isCas9){
        mit <- .getAggregateScore("score_mit")
        mitIndices <- match(spacers, names(mit))
        S4Vectors::mcols(guideSet)[["score_mit"]] <- mit[mitIndices]
    }
    return(guideSet)
}
