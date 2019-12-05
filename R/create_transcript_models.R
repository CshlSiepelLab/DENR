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
#' @return A list of matrices where each row in the matrix corresponds to a
#' bin and each column is a transcript
#' @export
create_transcript_models <- function(transcript_groups, bins) {
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
    # Get which transcripts intersect which bins
    hits <- GenomicRanges::findOverlaps(group_bins, tx_grp)
    # Compute number of bases the intersect the trancript per bin
    overlaps <- GenomicRanges::pintersect(
      group_bins[S4Vectors::queryHits(hits)],
      tx_grp[S4Vectors::subjectHits(hits)]
    )
    # Compute percent overlap
    percent_overlap <- IRanges::width(overlaps) /
                        IRanges::width(group_bins[S4Vectors::queryHits(hits)])
    # Reshape into matrix
    po_matrix <- reshape2::acast(
        cbind(S4Vectors::as.data.frame(hits), percent_overlap),
        formula = queryHits ~ subjectHits,
        value.var = "percent_overlap",
        fill = 0
    )
    colnames(po_matrix) <- tx_grp$TXNAME
    return(po_matrix)
  }, tx_grp = transcript_groups, group_bins = bins, SIMPLIFY = F)
  return(tx_matrix_models)
}
