#' Create transcript models
#'
#' Function for generating binned transcript models given a set of bins
#'
#' @param transcripts A \code{\link[GenomicRanges]{GRanges-class}}
#' object
#' @param bins A \code{\link[GenomicRanges]{GRangesList-class}} object,
#' each element in the list contains the bins for each corresponding group of
#' transcripts
#' @param bin_size width of bins
#' @param transcript_name_column A string that indicates which column in the
#' GRanges object contain the transcript names
#' @return A list of matrices where each row in the matrix corresponds to a
#' bin and each column is a transcript
#' @export
create_transcript_models <- function(transcripts, bins, bin_size,
                                     transcript_name_column) {
  # check input class
  if (!methods::is(transcripts, "GRanges")) {
    stop("transcripts is not a GRanges object")
  }
  # Compute percent overlap of each transcript per bin and cast to matrix
  # First find all pairwise overlaps
  ovr <- IRanges::findOverlapPairs(transcripts, bins)
  # Compute bin percent overlap for all intersections
  ovr_intersect <- IRanges::pintersect(ovr, drop.nohit.ranges = FALSE)
  ovr_val <- GenomicRanges::width(ovr_intersect) / bin_size

  # Compute the dimensions of model matrix for each loci
  d <- S4Vectors::elementNROWS(ovr_val)
  # Pre-extract transcript names for ease of access
  transcript_names <-
    GenomicRanges::values(transcripts)[, transcript_name_column]
  # Reshape intersection output into matricies of the correct dimension
  tx_matrix_models <- lapply(
    S4Vectors::split(d, names(d))[names(bins)], function(x) {
      # Retrive rows in transcripts that each intersection pair came from
      lookup_ind <- which(names(ovr_val) == names(x[1]))
      # Place into matrix of correct dimension
      m <- matrix(data = BiocGenerics::unlist(ovr_val[lookup_ind]),
                  nrow = x[1],
                  ncol = length(x))
      # Set column names to correct transcript names
      colnames(m) <- transcript_names[lookup_ind]
      return(m)
    }
  )
  return(tx_matrix_models)
}
