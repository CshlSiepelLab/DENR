#' Create model masks
#'
#' Function for generating masks
#'
#' @param transcript_models A list of matricies containing transcript models
#' @param mask_start_bins the number of bins with non-zero entries at the start
#' of a transcript that should be masked (Used for masking both start and end if
#' strand is \'*\')
#' @param mask_end_bins the number of bins with non-zero entries at the end
#' of a transcript that should be masked (Ignored if strand is \'*\')
#' @param strand a vector containing one entry per group of transcript models
#' indicating the strand
#'
#' @rdname create_model_masks
#'
#' @include type_checkers.R
#' @return A list of vectors, where each vector indicates the bins that should
#' be masked
#' @export
create_model_masks <- function(transcript_models, strand, mask_start_bins = 0,
                               mask_end_bins = 0) {
  # ***Checks***
  if (!is_strand_vector(strand)) {
    stop("strand must be a vector containing only '+','-', or '*'")
  }
  if (!is_matrix_list(transcript_models)) {
    stop("transcript models must be a list of matricies")
  }
  if (mask_start_bins < 0 || mask_end_bins < 0) {
    stop("Number of bins to mask must be > 0")
  }
  if (length(strand) != length(transcript_models)) {
    stop("strand and transcript_models must be the same length")
  }

  # ***End checks***

  # Get indicies of masked bins per transcript group
  all_masks <- mapply(function(tx_models, strand) {
    group_masks <- apply(tx_models, 2, function(tx) {
      first_non_zero <- min(which(tx != 0))
      last_non_zero <- max(which(tx != 0))
      # Orient the first/last axis along the 5' to 3' axis
      if (strand == "+") {
        mask <- unique(c(seq(from = first_non_zero,
                             length.out = mask_start_bins,
                             by = 1),
                         seq(from = last_non_zero,
                             length.out = mask_end_bins,
                             by = -1)))
      } else if (strand == "-") {
        mask <- unique(c(seq(from = first_non_zero,
                             length.out = mask_end_bins,
                             by = 1),
                         seq(from = last_non_zero,
                             length.out = mask_start_bins,
                             by = -1)))
      } else if (strand == "*") {
        mask <- unique(c(seq(from = first_non_zero,
                            length.out = mask_start_bins,
                            by = 1),
                        seq(from = last_non_zero,
                            length.out = mask_start_bins,
                            by = -1)))
      }
      # Ensure all masked bins are between the first and last non-zero entries
      mask <- mask[mask >= first_non_zero & mask <= last_non_zero]
      return(mask)
    })
    # Ensure masks are returned as one 1D vector per transcript group
    return(sort(unique(as.vector(unlist(group_masks)))))
  }, tx_models = transcript_models, strand = as.list(strand), SIMPLIFY = F)
  return(all_masks)
}

#' Mask transcripts
#'
#' Function for generating masks
#'
#' @param transcript_models A list of matricies containing transcript models
#' @param masks a list of
#'
#' @rdname mask_transcripts
#'
#' @include type_checkers.R
#' @return A list of matrices where each row in the matrix corresponds to a
#' transcript and each column is a bin
#' @export
mask_transcripts <- function(transcript_models, masks) {
  # ***Checks***
  if (!is_matrix_list(transcript_models)) {
    stop("transcript models must be a list of matricies")
  }
  if (length(masks) != length(transcript_models)) {
    stop("masks and transcript_models must be the same length")
  }

  # ***End checks***

  # Mask transcript models
  masked_transcripts <- mapply(function(tx_models, masks) {
    tx_models[masks, ] <- 0
    return(tx_models)
  }, tx_models = transcript_models, masks = masks, SIMPLIFY = F)
  return(masked_transcripts)
}
