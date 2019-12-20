#' transcript_quantifier_valid
#'
#' Checks that a \code{transcript_quantifier} object is valid
#' @param object a \code{transcript_quantifier} object
#'
#' @return TRUE if valid, else errors
transcript_quantifier_valid <- function(object) {
  errors <- c()
  # Check that lists contain only matricies
  if (!is_matrix_list(object@models)) {
    errors <- c(errors, "All elements of models must be matricies")
  }
  if (!is_vector_list(object@masks)) {
    errors <- c(errors, "All elements of masks must be vectors")
  }
  # If counts are specified
  if (length(object@counts) != 0) {
    # Check that it is a list of matricies
    if (!is_matrix_list(object@counts)) {
      errors <- c(errors, "All elements of counts must be matricies")
    }
    # Check that counts has the same number of rows per element as models has
    # columns if counts exist
    if (!all(mapply(nrow, object@counts, SIMPLIFY = T) ==
            mapply(ncol, object@models, SIMPLIFY = T))) {
      errors <- c(errors, "Number of model and count bins differ")
    }
  }
  if (length(errors) == 0) TRUE else errors
}

#' Class transcript_quantifier
#'
#' Class \code{transcript_quantifier} holds encodings for each transcript and
#' can be augmented to include data and transcript abundance estimates using the
#' \code{fit()} function
#'
#' @slot transcripts a  \code{\link[GenomicRanges]{GRanges-class}} that holds
#' all the transcript coordinates
#' @slot column_identifiers a two element character vector that holds the column
#' names in \code{transcripts} for the trancript and gene identifiers
#' respectively. Only a column for the transcript identifier is required,
#' missing gene identifiers are set to NA.
#' @slot bins a \code{\link[GenomicRanges]{GRangesList-class}} object that
#' records the bins for each group of transcripts
#' @slot models a list of matrices with numeric values between 0 and 1
#' representing the fractional overlap of per transcript per bin where the rows
#' are the transcripts and the columns are the bins. Rownames are the
#' corresponding transcripts
#' @slot masks a list of matrices with numeric values between 0 and 1
#' that are used to modify the models
#' @slot transcript_model_key A three column \code{data.frame} that maps
#' transcripts to their group and model
#' @slot counts a list of vectors containing the read counts per bin.
#' Initialized empty.
#' @slot model_abundance A list of vectors corresponding to \code{models}
#' of transcript abundances. Initialized at 0.
#'
#' @name transcript_quantifier-class
#' @rdname transcript_quantifier-class
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom GenomicRanges CompressedGRangesList
#' @exportClass transcript_quantifier
methods::setClass("transcript_quantifier",
                  slots = c(transcripts = "GRanges",
                            column_identifiers = "character",
                            bins = "CompressedGRangesList",
                            models = "list",
                            masks = "list",
                            transcript_model_key = "data.frame",
                            counts = "list",
                            model_abundance = "list"),
                  validity = transcript_quantifier_valid
)

#' transcript_quantifier
#'
#' Contructs an object that holds the transcript models
#' @param transcripts a \link[GenomicRanges]{GRanges-class} object that must
#' contain a metadata column with a transcript id and may contain an additional
#' column with a gene id
#' @param gene_name_column a string that indicates which column in the
#' GRanges object contains the gene names (not required)
#' @inheritParams create_bins
#' @inheritParams group_transcripts
#' @inheritParams create_model_masks
#' @inheritParams reduce_transcript_models
#' @inheritParams create_transcript_models
#' @param threads number of threads that can be used
#'
#' @return an \code{\link{transcript_quantifier-class}} object
#'
#' @export
transcript_quantifier <- function(transcripts, transcript_name_column,
                             gene_name_column = NULL,
                             bin_size = 100, distance = 0,
                             mask_start_bins = NULL, mask_end_bins = NULL,
                             bin_operation = c("round", "floor", "ceiling"),
                             threads = 1) {
  # **Some checks prior to beginning construction**
  # Check for correct GRanges object metadata
  if (!transcript_name_column %in%
      colnames(S4Vectors::elementMetadata(transcripts))) {
    stop(paste("transcripts does not have a column matching",
               transcript_name_column))
  }
  if (!is.character(
    S4Vectors::elementMetadata(transcripts)[, transcript_name_column])) {
    stop(paste("transcripts column", transcript_name_column,
               "must be of class 'character'"))
  }
  if (!is.null(gene_name_column)) {
    # Alias for easier use
    gnc <- gene_name_column
    if (!gnc %in%
       colnames(S4Vectors::elementMetadata(transcripts))) {
      stop(paste("transcripts does not have a column matching",
                 gnc))
    }
    if (!is.character(
      S4Vectors::elementMetadata(transcripts)[, gnc])) {
      stop(paste("transcripts column", gnc, "is class", paste0("'",
                 class(S4Vectors::elementMetadata(transcripts)[, gnc]), "'"),
                 "must be of class 'character'"))
    }
  }

  # Check binsize
  if (bin_size < 1) {
    stop("Binsize must be >= 1")
  }
  # Check distance
  if (distance < 0) {
    stop("Distance cannot be less than 0")
  }

  # **End checks**

  # Group transcripts
  tx_grps <- group_transcripts(transcripts,
                              distance = distance,
                              threads = threads)
  # Get group strands
  grp_strand <- as.character(unlist(lapply(tx_grps, function(x) {
    return(S4Vectors::runValue(GenomicRanges::strand(x))[1])
  })))

  # Bin transcripts regions
  grp_bins <- create_bins(transcript_groups = tx_grps,
                bin_size = bin_size)

  # Create transcript models
  tx_models <- create_transcript_models(tx_grps, grp_bins,
                                        transcript_name_column)

  # Create masks
  model_masks <- create_model_masks(transcript_models = tx_models,
                                    strand = grp_strand,
                                    mask_start_bins = mask_start_bins,
                                    mask_end_bins = mask_end_bins)

  # Reduce transcript models and generate transcript_model_key
  reduced_models <- reduce_transcript_models(
    transcript_models_ls = mask_transcripts(tx_models, model_masks),
    bin_operation = bin_operation)

  if (is.null(gene_name_column)) {
    gene_name_column <- NA
  }

  # Initialize model_abundances
  abundance <- lapply(reduced_models[[1]], function(x) return(numeric(ncol(x))))

  # Return transcript model object
  return(methods::new(Class = "transcript_quantifier",
                transcripts = transcripts,
                column_identifiers = c(transcript_id = transcript_name_column,
                                gene_id = gene_name_column),
                bins = grp_bins,
                models = reduced_models[[1]],
                masks = model_masks,
                transcript_model_key = reduced_models[[2]],
                counts = list(),
                model_abundance = abundance
               ))
}
