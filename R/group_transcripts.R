#' @title Group transcripts
#'
#' @description
#' Creates groups of transcripts on the same strand based on
#' proximity. The groups are constructed of connected transcripts
#' such that a pair of transcripts are considered connected if
#' they are within a given distance \emph{d} of each other. The
#' group then consists of all transcripts that are connected to
#' at least one other member.
#'
#' @param transcript_granges \code{\link[GenomicRanges]{GRanges-class}}
#' object
#' @param distance the distance within which two transcripts are
#' considered connected
#' @param threads The number of threads to use
#'
#' @return A \code{\link[GenomicRanges]{GenomicRangesList-class}}
#' object
#'
#' @rdname group_transcripts
#' @export
group_transcripts <- function(transcript_granges, distance = 0, threads = 1) {
  # Check that txdb is an actual transcript database object
  if (!methods::is(transcript_granges, "GRanges")) {
    stop("transcript is not a TxDb object")
  }
  # Ensure everything is sorted, otherwise later assignment assumptions fail
  transcript_granges <- GenomeInfoDb::sortSeqlevels(transcript_granges)
  transcript_granges <- GenomicRanges::sort(transcript_granges)
  # Assign a unique id to each row
  transcript_granges$unique_id <- seq_len(length(transcript_granges))
  # Create expanded granges
  tx_granges_expand <- GenomicRanges::resize(transcript_granges,
                             width = GenomicRanges::width(transcript_granges) +
                               floor(distance / 2),
                             fix = "end")
  tx_granges_expand <- GenomicRanges::resize(tx_granges_expand,
                             width = GenomicRanges::width(tx_granges_expand) +
                               ceiling(distance / 2),
                             fix = "start")
  tx_granges_expand <- GenomicRanges::trim(tx_granges_expand)
  # Reduce granges to get unions of overlapping ranges which will become the
  # transcript groups. The -/+ 1 bit is to assure that perfectly adjoining
  # groups are not merged
  GenomicRanges::end(tx_granges_expand) <-
    GenomicRanges::end(tx_granges_expand) - 1
  tx_reduce <- GenomicRanges::reduce(tx_granges_expand, ignore.strand = FALSE)
  GenomicRanges::end(tx_reduce) <- GenomicRanges::end(tx_reduce) + 1
  # Get the group assignments
  group_assignment <- GenomicRanges::findOverlaps(query = transcript_granges,
                                                  subject = tx_reduce)
  # Create group vector string
  uid <- paste0(S4Vectors::subjectHits(group_assignment), "_",
                GenomicRanges::strand(transcript_granges))

  # Split the transcripting into their groups
  gr_groups <- GenomicRanges::split(transcript_granges,
                                    f = uid)

  return(gr_groups)
}
