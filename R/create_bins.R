#' Create bins for multiple groups of overlapped transcripts
#'
#' Function for generating bins for mutiple groups of overlapped transcripts.
#'
#' @param transcript_groups A \code{\link[GenomicRanges]{GRangesList-class}}
#' object,each element in the list contains transcripts from one gene or several
#' genes which are overlapped.
#' @param bin_size An integer, used to tile the gene region. Default is 250bp.
#' @return A \code{\link[GenomicRanges]{GRangesList-class}} object, containing
#' the binned region.
#' @export
create_bins <- function(transcript_groups, bin_size = 250) {
  # check input class
  if (!methods::is(transcript_groups, "GRangesList")) {
    stop("transcript_groups is not a GRangesList object")
  }
  if (!methods::is(bin_size, "numeric") || bin_size < 0) {
    stop("bin_size is not a positive number")
  }

  # First reduce transcript groups to single range using large gapwidth to
  # ensure single range
  red_tx_grps <- unlist(
    GenomicRanges::reduce(transcript_groups, min.gapwidth = 1e9))
  # Pre-compute number of bins in each range
  bin_count <- ceiling(GenomicRanges::width(red_tx_grps) / bin_size)
  # Pre-build rle for each element
  chrom <-
    S4Vectors::Rle(as.vector(GenomicRanges::seqnames(red_tx_grps)), bin_count)
  strnd <-
    S4Vectors::Rle(as.vector(GenomicRanges::strand(red_tx_grps)), bin_count)
  # Create Granges as single GRanges
  tiles <- GenomicRanges::GRanges(
    chrom,
    IRanges::IRanges(start = unlist(
      seqv(from = GenomicRanges::start(red_tx_grps),
           to = GenomicRanges::end(red_tx_grps),
           by = bin_size)),
      width = bin_size),
    strand = strnd)
  # Return split GRangesList
  return(GenomicRanges::split(
    x = tiles, rep(names(transcript_groups), bin_count))
  )
}

#' Vectorized seq.default
#'
#' Version of the \code{seq} function that takes vectorized arguments
#'
#' @inheritParams base::seq
seqv <- Vectorize(seq.default, vectorize.args = c("from", "to"))
