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
    # Reformat abundance estimates
    abundance_tab <- apply(tq@transcript_model_key, 1, function(x) {
      a <- tq@model_abundance[[as.integer(x[2])]][as.integer(x[3])]
      data.frame(transcript_name = as.character(x[1]),
               abundance = as.numeric(a),
               model = paste0(as.integer(x[2]), "_", as.integer(x[3])))
    })
    abundance_tab <- do.call("rbind",
                            args = c(abundance_tab, stringsAsFactors = FALSE))
    # Clean up the output
    i <- sapply(abundance_tab, is.factor)
    abundance_tab[i] <- lapply(abundance_tab[i], as.character)
    rownames(abundance_tab) <- NULL
    # Sort to match transcripts
    ord <- base::match(
      S4Vectors::elementMetadata(tq@transcripts)[[tq@column_identifiers[1]]],
      abundance_tab$transcript_name
    )
    abundance_tab <- abundance_tab[ord, ]
    # Add gene names if they exist
    if (!is.na(tq@column_identifiers[2])) {
      abundance_tab$gene_name <- # nolint
        GenomicRanges::values(tq@transcripts)[, tq@column_identifiers[2]]
    }
    return(abundance_tab)
  }
)
