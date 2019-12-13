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
    stop(paste("chromosome", chromosome, "does not exist in", bigwig_file))
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
#' @param threads number of threads to use
#'
#' @name summarize_bigwig
#' @rdname summarize_bigwig
#'
#' @export
methods::setGeneric("summarize_bigwig",
                    function(bigwig_file, bins, summary_operation = "sum",
                             threads = 1) {
                      standardGeneric("summarize_bigwig")
})

#' @rdname summarize_bigwig
methods::setMethod("summarize_bigwig", signature(bins = "GRangesList"),
                   function(bigwig_file, bins, summary_operation = "sum",
                            threads = 1) {
  # Create iterator over bins
  bins_iter <- iterators::iter(bins)

  if (threads > 1) {
    # Create parallel cluster and register it with the foreach backend
    cluster <- snow::makeCluster(threads)
    doSNOW::registerDoSNOW(cluster)
    `%doloop%` <- foreach::`%dopar%`
  } else {
    `%doloop%` <- foreach::`%do%`
  }

  # Iterate over bin groups and get read sums
  sum_list <- foreach::foreach(chrom_bins = bins_iter,
                               .noexport = c("bins"),
                               .errorhandling = "stop") %doloop% {
          summarize_bigwig(bigwig_file, chrom_bins, summary_operation)
  }

  # Set names
  names(sum_list) <- names(bins)
  return(sum_list)
})

#' @rdname summarize_bigwig
methods::setMethod("summarize_bigwig", signature(bins = "GRanges"),
                   function(bigwig_file, bins, summary_operation = "sum",
                            threads = 1) {
  # ** Checks **
  if (!file.exists(bigwig_file)) {
    stop("bigwig_file path does not exist")
  }

  # Check if more than one unique seqname in bins
  bin_chrom <- S4Vectors::runValue(GenomicRanges::seqnames(bins))
  if (length(bin_chrom) > 1) {
    stop("One or more groups of contained multiple chromosomes")
  }

  # Check that query chromosome is present in bigwig
  assert_chromosome_exists(chromosome = bin_chrom, bigwig_file)

  # ** End checks **

  # Read in bigwig region
  import_range <- GenomicRanges::reduce(bins)
  imported_bw <- rtracklayer::import.bw(con = bigwig_file, which = import_range,
                                          as = "RleList")
  # Lookup correct method
  operation <- methods::selectMethod(summary_operation, signature = "RleViews")
  # Summarize reads across bins
  tryCatch({
    .Generic <<- summary_operation # nolint
    read_summary <- BiocGenerics::do.call(what = operation,
      list(x = IRanges::Views(imported_bw, IRanges::ranges(bins))[[1]]))
    .Generic <<- NULL # nolint
  }, error = function(e) {
    err_string <- stringr::str_match(as.character(e), ":\\s*(.*)")[, 2]
    stop(paste0("Bigwig summarization with '", summary_operation, "' failed, ",
               err_string))
  })
  return(list(read_summary))
})

utils::globalVariables(c(".Generic", "chrom_bins"))
