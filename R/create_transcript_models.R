#' Create transcript models
#'
#' Function for generating binned transcript models given a set of bins
#'
#' @param transcript_groups A \code{\link[GenomicRanges]{GRangesList-class}}
#' object, each element in the list contains transcripts from one gene or
#' several genes which are overlapped.
#' @param bins A \code{\link[GenomicRanges]{GRangesList-class}} object,
#' each element in the list contains the bins for each corresponding group of
#' transcripts
#' @param transcript_name_column A string that indicates which column in the
#' GRanges object contain the transcript names
#' @return A list of matrices where each row in the matrix corresponds to a
#' bin and each column is a transcript
#' @export
create_transcript_models <- function(transcript_groups, bins,
                                     transcript_name_column) {
  # check input class
  if (!methods::is(transcript_groups, "GRangesList")) {
    stop("transcript_groups is not a GRangesList object")
  }
  if (!methods::is(bins, "GRangesList")) {
    stop("bins is not a GRangesList object")
  }
  if (length(transcript_groups) != length(bins)) {
    stop("transcript_groups and bins should be lists of the same length")
  }
  # Compute percent overlap of each transcript per bin and cast to matrix
  tx_matrix_models <- mapply(FUN = function(tx_grp, group_bins) {
    # Get intersection of all bins against all transcripts
    overlaps <- GenomicRanges::pintersect(
      rep(group_bins, length(tx_grp)),
      rep(tx_grp, each = length(group_bins))
    )

    # Compute percent overlap
    percent_overlap <- matrix(IRanges::width(overlaps) /
                        IRanges::width(group_bins)[1],
                        nrow = length(group_bins),
                        byrow = FALSE)
    colnames(percent_overlap) <-
      S4Vectors::elementMetadata(tx_grp)[, transcript_name_column]
    return(percent_overlap)
  }, tx_grp = transcript_groups, group_bins = bins, SIMPLIFY = F)
  return(tx_matrix_models)
}
