#' density ratio calculation
#'
#' Calculates the log2 ratio of the polymerase upstream and downstream of the TSS on the
#' sense strand. This ratio is used to determine if a transcript may be active. Used
#' during data addition to determine if a transcript specified as inactive in the
#' \code{fit} call should be overriden.NOTE: ratio is calculated with an additional
#' 1 pseudocount/kb added to the numerator and denominator
#'
#' @inheritParams add_data
#' @param up_width The length of the region upstream of the TSS to calculate the
#' polymerase denisty in
#' @param up_shift how far upstream of the TSS \code{up_width} is shifted
#' @param body_width The length of the region downstream of the TSS to calculate the
#' polymerase denisty in
#' @param body_shift how far downstream of the TSS \code{body_width} is shifted
#'
#' @return a named vector of the polymerase ratios
upstream_polymerase_ratio <- function(tq, bigwig_plus, bigwig_minus,
                                      up_width = 5e3, up_shift = 5e2,
                                      body_width = 5e3, body_shift = 5e2) {
  # Get upstream regions
  up_gr <- GenomicRanges::promoters(tq@transcripts, upstream = up_width - 1,
                                    downstream = 1)
  # Match bigwig and GRanges seqlevels
  up_gr[GenomicRanges::strand(up_gr) == "+"] <-
    GenomicRanges::shift(up_gr[GenomicRanges::strand(up_gr) == "+"],
                         shift = -up_shift - 1)
  up_gr[GenomicRanges::strand(up_gr) == "-"] <-
    GenomicRanges::shift(up_gr[GenomicRanges::strand(up_gr) == "-"],
                         shift = up_shift + 1)

  # Get body regions
  body_gr <- GenomicRanges::promoters(tq@transcripts, upstream = 0,
                                    downstream = body_width)
  body_gr[GenomicRanges::strand(body_gr) == "+"] <-
    GenomicRanges::shift(body_gr[GenomicRanges::strand(body_gr) == "+"],
                         shift = body_shift)
  body_gr[GenomicRanges::strand(body_gr) == "-"] <-
    GenomicRanges::shift(body_gr[GenomicRanges::strand(body_gr) == "-"],
                         shift = -body_shift)

  # Match seqinfo
  up_gr <- apply_bigwig_seqinfo(up_gr, bigwig_file = bigwig_plus)
  body_gr <- apply_bigwig_seqinfo(body_gr, bigwig_file = bigwig_plus)

  # Trim regions
  up_gr <- GenomicRanges::trim(up_gr)
  body_gr <- GenomicRanges::trim(body_gr)

  up_count <- numeric(length(up_gr))
  body_count <- numeric(length(body_gr))
  # Sumarize bigwig
  up_count[S4Vectors::decode(GenomicRanges::strand(up_gr) == "+")] <-
    summarize_bigwig(bigwig_file = bigwig_plus,
                     up_gr[GenomicRanges::strand(up_gr) == "+"],
                     summary_operation = "mean")
  up_count[S4Vectors::decode(GenomicRanges::strand(up_gr) == "-")] <-
    summarize_bigwig(bigwig_file = bigwig_minus,
                     up_gr[GenomicRanges::strand(up_gr) == "-"],
                     summary_operation = "mean")
  body_count[S4Vectors::decode(GenomicRanges::strand(body_gr) == "+")] <-
    summarize_bigwig(bigwig_file = bigwig_plus,
                     body_gr[GenomicRanges::strand(body_gr) == "+"],
                     summary_operation = "mean")
  body_count[S4Vectors::decode(GenomicRanges::strand(body_gr) == "-")] <-
    summarize_bigwig(bigwig_file = bigwig_minus,
                     body_gr[GenomicRanges::strand(body_gr) == "-"],
                     summary_operation = "mean")
  # Calculate log2(ratio) with 1 additional read/kb count in each range
  upr <- log2(abs(body_count) + 1e-3) - log2(abs(up_count) + 1e-3)
  names(upr) <- get_tx_id(tq)
  return(upr)
}
