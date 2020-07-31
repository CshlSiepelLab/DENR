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
    if (!is_vector_list(object@counts)) {
      errors <- c(errors, "All elements of counts must be matricies")
    }
    # Check that counts has the same number of rows per element as models has
    # columns if counts exist
    data_mod_dim_equal <- mapply(function(counts, models) {
      length(counts) == nrow(models)
    }, counts = object@counts, models = object@models)
    if (!all(data_mod_dim_equal)) {
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
#' @slot transcripts a \code{\link[GenomicRanges]{GRanges-class}} that holds
#' all the transcript coordinates
#' @slot column_identifiers a two element character vector that holds the column
#' names in \code{transcripts} for the trancript and gene identifiers
#' respectively. Only a column for the transcript identifier is required,
#' missing gene identifiers are set to NA.
#' @slot bins a \code{\link[GenomicRanges]{GRangesList-class}} object that
#' records the bins for each group of transcripts
#' @slot bin_size width of bins used to create transcript models
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
#' @slot upstream_polymerase_ratios the log2 ratio of mean counts for each transcript
#' from the region immediately upstream and downstream of each TSS;
#' log2([500bp, 5.5kb] / [-5.5kb, -500bp])
#' @slot count_metadata holds information about the files the count data came from
#' including names and total counts
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
                            bin_size = "integer",
                            models = "list",
                            masks = "list",
                            transcript_model_key = "data.frame",
                            counts = "list",
                            upstream_polymerase_ratios = "numeric",
                            count_metadata = "list",
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
#' @param distance the distance within which two transcripts are
#' considered connected (must be at least bin size, defaults to bin size). The
#' smaller this value is the more efficiently the model can be fit. Only
#' increase this if you are masking large regions at the starts and ends of
#' genes.
#' @inheritParams create_bins
#' @inheritParams group_transcripts
#' @inheritParams create_model_masks
#' @inheritParams reduce_transcript_models
#' @inheritParams create_transcript_models
#'
#' @return an \code{\link{transcript_quantifier-class}} object
#'
#' @export
transcript_quantifier <- function(transcripts, transcript_name_column,
                             gene_name_column = NULL,
                             bin_size = 250, distance = NULL,
                             mask_start_bins = NULL, mask_end_bins = NULL,
                             bin_operation = c("round", "floor", "ceiling")) {
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
  duplicated_transcript_names <-
    any(duplicated(S4Vectors::elementMetadata(transcripts)[, transcript_name_column]))
  if (duplicated_transcript_names) {
    stop("One or more transcript names is duplicated, transcript names must be unique")
  }
  # Checks on gene names
  if (!is.null(gene_name_column)) {
    # Alias for easier use
    gnc <- gene_name_column
    if (!gnc %in%
       colnames(S4Vectors::elementMetadata(transcripts))) {
      stop("transcripts does not have a column matching ", gnc)
    }
    if (!is.character(
      S4Vectors::elementMetadata(transcripts)[, gnc])) {
      stop(paste("transcripts column", gnc, "is class", paste0("'",
                 class(S4Vectors::elementMetadata(transcripts)[, gnc]), "'"),
                 "must be of class 'character'"))
    }
  }
  # Set default distance
  if (is.null(distance)) {
    message("Setting distance to default of ", bin_size)
    distance <- bin_size
  }
  # This is important to ensure no spurious overlaps in all vs. all overlap
  # during transcript model construction
  if (distance < bin_size) {
    stop("distance must be equal to or greater than bin size")
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

  # Force copy of object underlying GRanges to prevent any weird side effects if
  # GRanges is using a data.table or something else that can modify in place
  transcripts <- GenomicRanges::makeGRangesFromDataFrame(
    data.table::copy(data.table::as.data.table(transcripts)),
    keep.extra.columns = T)

  # Group transcripts
  message("Grouping transcripts...")
  tx_grps <- group_transcripts(transcripts, distance = distance)

  # Ensure is sorted
  tx_grps <- GenomeInfoDb::sortSeqlevels(tx_grps)
  tx_grps <- GenomicRanges::sort(tx_grps)

  # Get group strands
  grp_strand <- stringr::str_extract(names(tx_grps), "[+-]")

  # Bin transcripts regions
  message("Binning loci...")
  grp_bins <- create_bins(transcript_groups = tx_grps,
                bin_size = bin_size)

  # Create transcript models
  message("Creating transcript models ...")
  tx_models <- create_transcript_models(transcripts = transcripts,
                                        bins = grp_bins,
                                        bin_size = bin_size,
                                        transcript_name_column)

  # Create masks
  message("Creating masks ...")
  model_masks <- create_model_masks(transcript_models = tx_models,
                                    strand = grp_strand,
                                    mask_start_bins = mask_start_bins,
                                    mask_end_bins = mask_end_bins)

  # Reduce transcript models and generate transcript_model_key
  message("Merging redundant models ...")
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
                bin_size = as.integer(bin_size),
                models = reduced_models[[1]],
                masks = model_masks,
                transcript_model_key = reduced_models[[2]],
                counts = list(),
                upstream_polymerase_ratios = numeric(0),
                count_metadata = list(bigwig_plus = NA_character_,
                                      bigwig_minus = NA_character_,
                                      library_size = NA_real_),
                model_abundance = abundance
               ))
}

#' get_tx_id
#'
#' Convenience function for getting transcript names from a
#' \code{\link{transcript_quantifier-class}} object
#' @param tq \code{\link{transcript_quantifier-class}} object
#'
#' @return a character vector of transcript ids in the same order as they are in the
#' \code{tq@transcripts} slot
get_tx_id <- function(tq) {
  return(GenomicRanges::mcols(tq@transcripts)[, tq@column_identifiers[1]])
}

#' @inherit methods::show
methods::setMethod("show", signature = "transcript_quantifier", function(object) {
  num_transcripts <- length(object@transcripts)
  num_models <- sum(unlist(lapply(object@models, ncol)))
  num_loci <- length(object@models)
  bin_size <- object@bin_size
  bwp <- object@count_metadata$bigwig_plus
  bwm <- object@count_metadata$bigwig_minus

  if (!is.na(object@column_identifiers["gene_id"])) {
    num_genes <- length(unique(
      S4Vectors::mcols(object@transcripts)[[object@column_identifiers["gene_id"]]]))
    gene_string <- paste("Number of Genes:", num_genes)
  } else {
    gene_string <- "No gene id present"
  }

  write("A transcript_quantifier object with:", file = stdout())
  write(paste(num_transcripts, "transcripts converted to", num_models, "models",
              "grouped into", num_loci, "loci"), file = stdout())
  write(gene_string, file = stdout())
  write(paste("bin size:", bin_size), file = stdout())
  write(paste("Bigwig data (plus):", bwp), file = stdout())
  write(paste("Bigwig data (minus):", bwm), file = stdout())
})
