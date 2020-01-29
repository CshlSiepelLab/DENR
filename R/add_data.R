#' @title Add data
#'
#' @description adds data from bigwig to a \link{transcript_quantifier-class}
#' object
#'
#' @inheritParams summarize_bigwig
#' @param transcript_quantifier A \link{transcript_quantifier-class} object
#' @param bigwig_plus the path to a bigwig for reads on the plus strand
#' @param bigwig_minus the path to a bigwig for reads on the minus strand
#'
#' @name add_data
#' @rdname add_data
#' @include transcript_quantifier-class.R
#' @return A \link{transcript_quantifier-class} object with count data added.
#' @export
methods::setGeneric("add_data",
                    function(transcript_quantifier,
                             bigwig_plus = NULL, bigwig_minus = NULL,
                             summary_operation = "sum") {
                      standardGeneric("add_data")
                    })

#' @rdname add_data
methods::setMethod("add_data",
  signature(transcript_quantifier = "transcript_quantifier"),
  function(transcript_quantifier,
           bigwig_plus = NULL, bigwig_minus = NULL,
           summary_operation = "sum") {
      bins <- transcript_quantifier@bins
      strands <- runValue(GenomicRanges::strand(bins))
      # summarize bigwig files by strands
      bw_counts <-
          c(
              summarize_bigwig(bigwig_plus, bins[unlist(strands == "+")],
                               summary_operation),
              summarize_bigwig(bigwig_minus, bins[unlist(strands == "-")],
                               summary_operation)
          )
      # reorder the counts as the order in bins
      transcript_quantifier@counts <-
          bw_counts[names(transcript_quantifier@models)]
      return(transcript_quantifier)
})
