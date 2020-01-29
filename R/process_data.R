if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(".Generic", "chrom_bins"))
}

#' @title Assert chromosome exists
#'
#' @description  Checks whether a chromosome exists in the target bigwig. Throws
#' error if it does not.
#' @param bigwig_file a string pointing to a bigwig file
#' @param chromosome chromosome name
#' @rdname assert_chromosome_exists
#' @return silent if true, else error
#' @export
assert_chromosome_exists <- function(chromosome, bigwig_file) {
  bw_seqnames <- GenomicRanges::seqnames(
    rtracklayer::seqinfo(rtracklayer::BigWigFile(bigwig_file))
  )
  pass <- all(chromosome %in% bw_seqnames)
  if (!pass) {
    not_incl <- setdiff(chromosome, bw_seqnames)
    stop(paste("chromosome",
               paste(not_incl, collapse = ", "),
               "does not exist in", bigwig_file))
  }
}

#' @title Summarize bigwig
#'
#' @description  Summarizes a bigwig over a set of ranges in a
#' \code{\link[GenomicRanges]{GRangesList-class}}
#' or a \code{\link[GenomicRanges]{GRanges-class}} object
#'
#' @inheritParams assert_chromosome_exists
#' @param bins A \code{\link[GenomicRanges]{GRangesList-class}} or
#' \code{\link[GenomicRanges]{GRanges-class}} object where all elements in a
#' given GRanges (or GRangesList element) are from the same chromosome
#' @param summary_operation the summary opperation to apply per bin (e.g.
#' sum, mean, median, etc.) Defaults to "sum"
#'
#' @return A list of vectors with each one corresponding to one set of bins and
#' each element of a vector corresponding to a bin
#'
#' @name summarize_bigwig
#' @rdname summarize_bigwig
#'
#' @export
methods::setGeneric("summarize_bigwig",
                    function(bigwig_file, bins, summary_operation = "sum") {
                      standardGeneric("summarize_bigwig")
})

#' @rdname summarize_bigwig
methods::setMethod("summarize_bigwig", signature(bins = "GRangesList"),
                   function(bigwig_file, bins, summary_operation = "sum") {
  # Flatten GRangesList
  flat_bins <- BiocGenerics::unlist(bins, use.names = TRUE)

  # Get the values back as flat
  sum_flat <- unlist(summarize_bigwig(
    bigwig_file,
    flat_bins,
    summary_operation
  ))

  # Re-list the values
  sum_list <- split(sum_flat, attr(sum_flat, "names"))
  return(sum_list)
})

#' @rdname summarize_bigwig
methods::setMethod("summarize_bigwig", signature(bins = "GRanges"),
                   function(bigwig_file, bins, summary_operation = "sum") {
  # ** Checks **
  if (!file.exists(bigwig_file)) {
    stop("bigwig_file path does not exist")
  }

  # Get unique chromosome names in bins
  bin_chrom <- unique(S4Vectors::runValue(GenomicRanges::seqnames(bins)))

  # Check that query chromosome is present in bigwig
  assert_chromosome_exists(chromosome = bin_chrom, bigwig_file)

  # ** End checks **
  # Read in the relevant bigwig regions
  import_range <- GenomicRanges::reduce(bins)
  imported_bw <- rtracklayer::import.bw(
    con = bigwig_file, which = import_range,
    as = "RleList")
  # Pre-allocate vector to store scores in
  scores <- numeric(length = length(bins))
  for (chrom in bin_chrom) {
    # Get bins that are relevant for that chromosome
    incl_bins <- as.logical(GenomicRanges::seqnames(bins) == chrom)
    # Lookup correct method
    operation <- methods::selectMethod(summary_operation,
                                       signature = "RleViews")
    # Summarize reads across bins
    tryCatch({
      .Generic <<- summary_operation # nolint
      scores[incl_bins] <- unlist(BiocGenerics::do.call(
        what = operation,
        args = list(
          x = IRanges::Views(
            subject = imported_bw[chrom],
            IRanges::ranges(bins[incl_bins]))
          )
        ))
      .Generic <<- NULL # nolint
    }, error = function(e) {
      err_string <- stringr::str_match(as.character(e), ":\\s*(.*)")[, 2]
      stop(paste0("Bigwig summarization with '", summary_operation,
                  "' failed, ", err_string))
    })
  }
  names(scores) <- names(bins)
  return(scores)
})
