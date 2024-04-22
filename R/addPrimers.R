#' @title Design primers for a \linkS4class{GuideSet} object. 
#' @description Design primers for a 
#'     \linkS4class{GuideSet} object. 
#' 
#' @author Matthew Pace, Jean-Philippe Fortin, Luke Hoberecht
#' 
#' @export
#' @rdname addPrimers
#' @importFrom S4Vectors split mcols<-
setMethod("addPrimers",
          "GuideSet",
          function(object,
                   flank_len=75,
                   name="primer_design"
){
    object <- .validateGuideSet(object)
    object <- .addPrimersToGuideSet(guideSet=object,
                          flank_len=flank_len,
                           s4_colname=name)
    return(object)
})


#' @rdname addPrimers
#' @export
setMethod("addPrimers", 
          "NULL", 
          function(object){
    return(NULL)
})

.addPrimersToGuideSet <- function(guideSet,
                                  flank_len,
                                  s4_colname
){
  start <- -flank_len-20
  end <- flank_len-1
  
  extendedSequences <- .getExtendedSequences(guideSet,
                                             start=start,
                                             end=end)
  good <- !is.na(extendedSequences)
  seqs <- extendedSequences[good]
  
  ideal_len <- flank_len*2+20+(20*2) # flanking regions, spacer, and primer pair
  range_low <- ideal_len-50
  range_high <- ideal_len+50
  df_list <- list()
  n_primer_pairs <- 3
  for (i in seq_along(extendedSequences)) {
    input_str <- paste0("SEQUENCE_ID=", names(seqs)[i], "\n",
                        "SEQUENCE_TEMPLATE=", seqs[[i]], "\n",
                        "PRIMER_OPT_SIZE=22\n", # Optimal Primer Length
                        "PRIMER_MIN_TM=55\n", # Minimum Tm
                        "PRIMER_OPT_TM=65\n", # Optimal Tm
                        "PRIMER_MAX_TM=70\n", # Maximum Tm
                        "PRIMER_OPT_GC_PERCENT=50\n", # Optimal GC%
                        "PRIMER_SALT_MONOVALENT=40\n", # Concentration of Monovalent Cations
                        "PRIMER_SALT_DIVALENT=2.5\n", # Concentration of Divalent Cations
                        "PRIMER_DNA_CONC=50\n", # Annealing Oligo Concentration
                        "PRIMER_DNTP_CONC=0.3\n", # dNTP Concentration
                        "PRIMER_NUM_RETURN=", n_primer_pairs, "\n",
                        "SIZE_RANGE=", range_low, "-", range_high, "\n", 
                        "PRIMER_PRODUCT_SIZE_RANGE=", range_low, "-", range_high, "\n", 
                        "=") 
            
    # write the input string to a temporary file
    input_file <- tempfile()
    writeLines(input_str, input_file)
    
    cmd <- paste("cat", input_file, "| primer3/src/primer3_core")
    
    # Capture the output of the command
    cmd_output <- system(cmd, intern = TRUE)
    
    # Remove empty strings from the output
    cmd_output <- cmd_output[nzchar(cmd_output)]
    
    # Initialize vectors to hold column names and values
    col_names <- c()
    values <- c()
    
    # Loop over each line in the output to separate column names and values
    for(line in cmd_output) {
      split_line <- strsplit(line, "=")[[1]]
      if(length(split_line) == 2) { # If line can be split into name and value
        col_names <- c(col_names, split_line[1])
        values <- c(values, split_line[2])
      }
    }
    
    # Create a dataframe from column names and values and add it to the list
    df_list[[i]] <- as.data.frame(t(values), stringsAsFactors = FALSE)
    colnames(df_list[[i]]) <- col_names
  }

  # Get all unique column names from all data frames in the list
all_colnames <- unique(unlist(lapply(df_list, names)))

# add missing columns
harmonize_df <- function(df, all_colnames) {
  missing_colnames <- setdiff(all_colnames, names(df))
  df[missing_colnames] <- NA  # add missing columns with NA
  df <- df[, all_colnames]  # optional, to order columns similarly
  return(df)
}

harmonized_df_list <- lapply(df_list, harmonize_df, all_colnames = all_colnames)

# rbind all harmonized data frames
primers <- do.call(rbind, harmonized_df_list)

  # Subset the dataframe to include only columns whose names end with "SEQUENCE" or "PRODUCT_SIZE"
  #primers <- primers[, grepl("SEQUENCE$|PRODUCT_SIZE$", colnames(primers))]

  s4_colname <- paste0("primer3_", s4_colname)
  S4Vectors::mcols(guideSet)[[s4_colname]] <- primers
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


.validateGuideSet <- function(obj,
                              eMessage=NULL
){
  if (is.null(eMessage)){
    eMessage <- "guideSet argument must be a GuideSet object."
  }
  isGuideSet <- methods::is(obj, "GuideSet")
  if (!isGuideSet){
    stop(eMessage)
  }
  return(obj)
}
