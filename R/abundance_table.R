#' @title View abundances
#'
#' @description Produces a table with estimated trancript abundances for a given
#' \code{\link{transcript_quantifier-class}} object
#'
#' @inheritParams add_data
#'
#' @include transcript_quantifier-class.R
#' @name abundance_table
#' @rdname abundance_table
#'
#' @export
methods::setGeneric("abundance_table",
                    function(transcript_quantifier) {
                      standardGeneric("abundance_table")
                    })

#' @rdname abundance_table
methods::setMethod("abundance_table",
                   signature(transcript_quantifier = "transcript_quantifier"),
  function(transcript_quantifier) {
    # Alias for ease of use
    tq <- transcript_quantifier
    # Unlist abundance table for ease of access
    abundance_vector <- unlist(tq@model_abundance)
    # Compute final index to match transcripts to indices
    abundance_lookup_index <-
      cumsum(!duplicated(with(tq@transcript_model_key,
                              paste0(group, "_", model))))
    # Add abundances to table
    abundance_dt <- data.table::data.table(
      transcript_name = tq@transcript_model_key$tx_name,
      abundance = abundance_vector[abundance_lookup_index],
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

## Appease R CMD check
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("gene_name"))
}