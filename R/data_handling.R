#' @title Add data
#'
#' @description adds data from bigwig to a \link{transcript_quantifier-class}
#' object
#'
#' @param tq A \link{transcript_quantifier-class} object
#' @param bigwig_plus the path to a bigwig for reads on the plus strand
#' @param bigwig_minus the path to a bigwig for reads on the minus strand
#'
#' @name add_data
#' @rdname add_data
#' @include transcript_quantifier-class.R
#' @return A \link{transcript_quantifier-class} object with count data added.
#' @export
methods::setGeneric("add_data",
                    function(tq, bigwig_plus = NULL, bigwig_minus = NULL) {
                      standardGeneric("add_data")
                    })

#' @rdname add_data
methods::setMethod("add_data",
  signature(tq = "transcript_quantifier"),
  function(tq, bigwig_plus = NULL, bigwig_minus = NULL) {
      bins <- tq@bins
      summary_operation <- "mean"
      strands <- S4Vectors::runValue(GenomicRanges::strand(bins))
      # summarize bigwig files by strands
      bw_counts <-
          c(
              lapply(summarize_bigwig(bigwig_plus, bins[unlist(strands == "+")],
                               summary_operation), abs),
              lapply(summarize_bigwig(bigwig_minus, bins[unlist(strands == "-")],
                               summary_operation), abs)
          )
      # reorder the counts as the order in bins
      tq@counts <- bw_counts[names(tq@models)]
      # Add count metadata
      tq@count_metadata$bigwig_plus <- bigwig_plus
      tq@count_metadata$bigwig_plus <- bigwig_minus
      tq@count_metadata$library_size <- abs(total_coverage(bigwig_plus)) +
        abs(total_coverage(bigwig_minus))
      return(tq)
})

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
  bw_seqnames <- GenomicRanges::seqnames(rtracklayer::seqinfo(
    rtracklayer::BigWigFile(bigwig_file)))
  pass <- all(chromosome %in% bw_seqnames)
  if (!pass) {
    not_incl <- setdiff(chromosome, bw_seqnames)
    stop(paste("chromosome",
               paste(not_incl, collapse = ", "),
               "does not exist in", bigwig_file))
  }
}

#' @title Total coverage
#'
#' @description  computes total coverage in bigwig
#' @param x a string pointing to a bigwig file
#'
#' @name total_coverage
#'
#' @export
total_coverage <- function(x) {
  # Check file existance
  if (!file.exists(x)) {
    stop("file does not exist")
  }
  # Check correct suffix
  ext <- tools::file_ext(x)
  if (ext %in% c("bw", "bigWig")) {
    tot_reads <- sum(sapply(rtracklayer::import.bw(con = x, as = "RleList"),
                            function(x) sum(abs(x))))
  } else {
    stop(ext, " is an unsupported file type")
  }
  return(tot_reads)
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
#' @param autostyle_seqlevels Logical. If \code{TRUE} matches the seqlevel style of the
#' \code{bigwig_file} amd \code{bins} object. (Default: TRUE)
#'
#' @return A list of vectors with each one corresponding to one set of bins and
#' each element of a vector corresponding to a bin
#'
#' @name summarize_bigwig
#' @rdname summarize_bigwig
#'
#' @importClassesFrom GenomicRanges GRanges GRangesList
#'
#' @export
methods::setGeneric("summarize_bigwig",
                    function(bigwig_file, bins, summary_operation = "sum",
                             autostyle_seqlevels = TRUE) {
                      standardGeneric("summarize_bigwig")
                    })

#' @rdname summarize_bigwig
methods::setMethod("summarize_bigwig", signature(bins = "GRangesList"),
                   function(bigwig_file, bins, summary_operation = "sum",
                            autostyle_seqlevels = TRUE) {
                     # Flatten GRangesList
                     flat_bins <- BiocGenerics::unlist(bins, use.names = TRUE)

                     # Get the values back as flat
                     sum_flat <- unlist(summarize_bigwig(
                       bigwig_file,
                       flat_bins,
                       summary_operation
                     ))

                     # Re-list the values
                     sum_list <- split(sum_flat, factor(attr(sum_flat, "names"),
                                                        levels = names(bins)))
                     return(sum_list)
                   })

#' @rdname summarize_bigwig
methods::setMethod("summarize_bigwig", signature(bins = "GRanges"),
                   function(bigwig_file, bins, summary_operation = "sum",
                            autostyle_seqlevels = TRUE) {
                     # ** Checks **
                     if (!file.exists(bigwig_file)) {
                       stop("bigwig_file path does not exist")
                     }

                     # Autostyle seqlevels if specified, this avoids issues
                     # with queries failing due to chr1 vs. 1, etc. type
                     # seqnames
                     if (autostyle_seqlevels) {
                       # supress errors complaining about equivalent maps
                       withCallingHandlers({
                         GenomeInfoDb::seqlevelsStyle(bins) <-
                           GenomeInfoDb::seqlevelsStyle(
                             rtracklayer::BigWigFile(bigwig_file))[1]
                       }, warning = function(w) {
                         if (startsWith(conditionMessage(w),
                                        "found more than one best sequence"))
                           invokeRestart("muffleWarning")
                       })

                     }

                     # Restrice seqlevels in query to those in use
                     GenomeInfoDb::seqlevels(bins) <-
                       GenomeInfoDb::seqlevelsInUse(bins)

                     # Get unique chromosome names in bins
                     bin_chrom <-
                       unique(S4Vectors::runValue(GenomicRanges::seqnames(bins)))

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
                         err_string <-
                           stringr::str_match(as.character(e), ":\\s*(.*)")[, 2]
                         stop(paste0("Bigwig summarization with '", summary_operation,
                                     "' failed, ", err_string))
                       })
                     }
                     names(scores) <- names(bins)
                     return(scores)
                   })
