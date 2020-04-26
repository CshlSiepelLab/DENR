#' @title View transcript abundances
#'
#' @description Produces a table with estimated trancript abundances for a given
#' \code{\link{transcript_quantifier-class}} object
#'
#' @inheritParams add_data
#'
#' @include transcript_quantifier-class.R
#' @name transcript_abundance
#' @rdname transcript_abundance
#'
#' @export
methods::setGeneric("transcript_abundance",
                    function(tq) {
                      standardGeneric("transcript_abundance")
                    })

#' @rdname transcript_abundance
methods::setMethod("transcript_abundance",
                   signature(tq = "transcript_quantifier"),
  function(tq) {
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

#' @title View gene abundances
#'
#' @description Produces a table with estimated gene abundances for a given
#' \code{\link{transcript_quantifier-class}} object if gene identifiers are
#' present
#'
#' @inheritParams add_data
#'
#' @include transcript_quantifier-class.R
#' @name gene_abundance
#' @rdname gene_abundance
#'
#' @export
methods::setGeneric("gene_abundance",
                    function(tq) {
                      standardGeneric("gene_abundance")
                    })

#' @rdname gene_abundance
methods::setMethod("gene_abundance",
                   signature(tq = "transcript_quantifier"),
                   function(tq) {
                     # If no gene catagory present return error
                     if (is.na(tq@column_identifiers[2])) {
                       stop("No gene level annotation present")
                     }
                     # Get abundance at the per transcript level
                     tx_abund_tab <- transcript_abundance(tq)
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
