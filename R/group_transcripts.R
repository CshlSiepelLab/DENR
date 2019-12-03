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
  # Ensure everything is sorted
  transcript_granges <- GenomeInfoDb::sortSeqlevels(transcript_granges)
  transcript_granges <- GenomicRanges::sort(transcript_granges)
  # Assign a unique id to each row
  transcript_granges$unique_id <- seq_len(length(transcript_granges))
  # Get chromosome and strand strings
  strand <- S4Vectors::decode(GenomicRanges::seqnames(transcript_granges))
  chrom <- S4Vectors::decode(GenomicRanges::strand(transcript_granges))
  # Create iterator over transcript_granges with chromosome_strand information
  tx_gr_iter <- iterators::isplit(transcript_granges,
                                  paste0(strand, "_", chrom))

  if (threads > 1) {
    # Create parallel cluster and register it with the foreach backend
    cluster <- snow::makeCluster(threads)
    doSNOW::registerDoSNOW(cluster)
    `%doloop%` <- foreach::`%dopar%`
  } else {
    `%doloop%` <- foreach::`%do%`
  }
  # Iterate over
  tx_groups <- foreach::foreach(tx = tx_gr_iter, .combine = "rbind") %doloop% {
    # Get distance between end of transcript and start of following trans.
    # Negative values mean they overlap
    dist_to_next <- with(tx, GenomicRanges::start(value)[-1] -
                           GenomicRanges::end(value[1:(length(value) - 1)]))
    # Compute where breaks should be
    break_after <- which(dist_to_next > distance)
    # Compute groupings
    if (length(break_after) > 0) {
      num_groups <- length(break_after) + 1
      members_per_group <- c(diff(c(0, break_after)),
                             length(tx$value) -
                               break_after[length(break_after)])
      group_labels <- paste0(tx$key, ":",
                             rep(x = 1:num_groups, times = members_per_group))
    } else {
      group_labels <- rep(paste0(tx$key, ":", 1), length(tx$value))
    }
    uid_grouping <- data.frame(uid = tx$value$unique_id,
                               group_labels = group_labels)
    return(uid_grouping)
  }
  # Stop cluster
  if (threads > 1) {
    snow::stopCluster(cluster)
  }
  # Split the transcripting into their groups
  tx_groups <- tx_groups[order(tx_groups$uid), ]
  gr_groups <- GenomicRanges::split(transcript_granges,
                                    f = tx_groups$group_labels)
  return(gr_groups)
}

## Appease R CMD check
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("tx"))
}
