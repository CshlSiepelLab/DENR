#' @title View transcript abundances
#'
#' @description Produces a table with estimated trancript abundances for a given
#' \code{\link{transcript_quantifier-class}} object.
#' @param norm_method type of normalization to use: "tpm" or "tmm". Defaults
#' to "tpm". Described further in details.
#'
#'
#' @details We provide two methods for normalizing transcript abundances, with each
#' having different strengths. Both of the methods start with the polymerase desity (
#' polymerase per bp) estimated by the DENR method but handle sequencing depth
#' differently:
#' \itemize{
#'  \item{"tpm"}{Abundances are returned in the form analogous to TPM
#'  (transcripts per million) where the output units are normalized for sequencing depth
#'  by dividing the polymerase density by the total read count and multiplying by 10^6}
#'  \item{"tmm"}{Normalizes per transcript abundance using the inverse of the trimmed
#'  mean abundance for transcripts with abundance > 0. Uses transcripts in the 20%-80%
#'  quantile range.}
#' }
#'
#' @inheritParams add_data
#'
#' @include transcript_quantifier-class.R
#' @name transcript_abundance
#' @rdname transcript_abundance
#'
#' @export
methods::setGeneric("transcript_abundance",
                    function(tq, norm_method = c("tpm", "tmm")) {
                      standardGeneric("transcript_abundance")
                    })

#' @rdname transcript_abundance
methods::setMethod("transcript_abundance",
                   signature(tq = "transcript_quantifier"),
  function(tq, norm_method = c("tpm", "tmm")) {
    # Unlist abundance table for ease of access
    abundance_vector <- unlist(tq@model_abundance)
    # Compute 1/total_abundance then mutiply by 1e6 so that it becomes the correct
    # scaling factor for a metric analogous to TPM
    norm_method <- norm_method[1]
    if (norm_method == "tpm") {
      # Note: no need to correct for transcript length because abundances are already
      # based on per bp density
      scale_factor <-  1e6 / sum(abundance_vector)
    } else if (norm_method == "tmm") {
      scale_factor <- 1 / mean(abundance_vector[abundance_vector > 0], trim = 0.2)
    } else {
      stop("Invalid normalization method")
    }
    # Compute final index to match transcripts to indices
    abundance_lookup_index <-
      cumsum(!duplicated(with(tq@transcript_model_key,
                              paste0(group, "_", model))))
    # Add abundances to table
    abundance_dt <- data.table::data.table(
      transcript_name = tq@transcript_model_key$tx_name,
      abundance = abundance_vector[abundance_lookup_index] * scale_factor,
      model = with(tq@transcript_model_key, paste0(group, "_", model))
    )
    # Sort to match transcripts
    ord <- base::match(
      S4Vectors::elementMetadata(tq@transcripts)[[tq@column_identifiers[1]]],
      abundance_dt$transcript_name
    )
    abundance_dt <- abundance_dt[ord, ]
    # Add gene names if they exist
    if (!is.na(tq@column_identifiers[2])) {
      abundance_dt[, gene_name := GenomicRanges::values(tq@transcripts)[, tq@column_identifiers[2]]] #nolint
    }
    return(abundance_dt)
  }
)

#' @title View gene abundances
#'
#' @description Produces a table with estimated gene abundances for a given
#' \code{\link{transcript_quantifier-class}} object if gene identifiers are
#' present
#'
#' @inheritParams transcript_abundance
#'
#' @include transcript_quantifier-class.R
#' @name gene_abundance
#' @rdname gene_abundance
#'
#' @export
methods::setGeneric("gene_abundance",
                    function(tq, norm_method = c("tpm", "tmm")) {
                      standardGeneric("gene_abundance")
                    })

#' @rdname gene_abundance
methods::setMethod("gene_abundance",
                   signature(tq = "transcript_quantifier"),
                   function(tq, norm_method = c("tpm", "tmm")) {
                     # If no gene catagory present return error
                     if (is.na(tq@column_identifiers[2])) {
                       stop("No gene level annotation present")
                     }
                     # Get abundance at the per transcript level
                     tx_abund_tab <- transcript_abundance(tq, norm_method)
                     # Remove all duplicated models and sum over transcripts per
                     # gene
                     gene_abund_tab <-
                       tx_abund_tab[!duplicated(model),
                                    .(abundance = sum(abundance)),
                              by = "gene_name"]
                     return(gene_abund_tab)
                   }
)

## Appease R CMD check
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("gene_name", "abundance", "model"))
}
