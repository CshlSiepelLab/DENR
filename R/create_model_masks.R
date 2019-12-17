#' Create model masks
#'
#' Function for generating masks
#'
#' @param transcript_models A list of matricies containing transcript models.
#' @param mask_start_bins A numeric vertor with length two which giving the
#' number of bins should be masked near the start of a transcript. The first
#' integer is the number of bins will be masked inside the transcript, while
#' the second interger is the number outside the transcript. Default c(0, 0).
#' @param mask_end_bins A numeric vertor with length two which giving the
#' number of bins should be masked near the end of a transcript. The first
#' integer is the number of bins will be masked inside the transcript, while
#' the second interger is the number outside the transcript. Default c(0, 0).
#' @param strand a vector containing one entry per group of transcript models
#' indicating the strand.
#'
#' @rdname create_model_masks
#'
#' @include type_checkers.R
#' @return A list of vectors, where each vector indicates the bins that should
#' be masked
#' @export
#'
#'
create_model_masks <- function(transcript_models, strand,
                               mask_start_bins = NULL,
                               mask_end_bins = NULL) {
  if (is.null(mask_start_bins)) mask_start_bins <- c(0, 0)
  if (is.null(mask_end_bins)) mask_end_bins <- c(0, 0)
  # ***Checks***
  if (!is_strand_vector(strand, allow_star = FALSE)) {
    stop("strand must be a vector containing only '+','-'")
  }
  if (!is_matrix_list(transcript_models)) {
    stop("transcript models must be a list of matricies")
  }
  if (!(is_length_two_vector(mask_start_bins) &&
          is_length_two_vector(mask_end_bins))) {
        stop("the length of input vector should be two,
             and the numbers in the vertor should be non-negative integers")
    }
  if (length(strand) != length(transcript_models)) {
    stop("strand and transcript_models must be the same length")
  }
  # ***End checks***
 mask_start_in <- mask_start_bins[1]
 mask_start_out <- mask_start_bins[2]
 mask_end_in <- mask_end_bins[1]
 mask_end_out <- mask_end_bins[2]

  # Get indicies of masked bins per transcript group
  all_masks <- mapply(function(tx_models, strand) {
    group_masks <- apply(tx_models, 2, function(tx) {
      first_non_zero <- min(which(tx != 0))
      last_non_zero <- max(which(tx != 0))
      # Orient the first/last axis along the 5' to 3' axis
      if (strand == "+") {
        mask <- unique(c(seq(from = first_non_zero,
                             length.out = mask_start_in,
                             by = 1),
                         seq(from = first_non_zero - 1,
                             length.out = mask_start_out,
                             by = -1),
                         seq(from = last_non_zero,
                             length.out = mask_end_in,
                             by = -1),
                         seq(from = last_non_zero + 1,
                             length.out = mask_end_out,
                             by = 1)))
      } else if (strand == "-") {
        mask <- unique(c(seq(from = first_non_zero,
                             length.out = mask_end_in,
                             by = 1),
                         seq(from = first_non_zero - 1,
                             length.out = mask_end_out,
                             by = -1),
                         seq(from = last_non_zero,
                             length.out = mask_start_in,
                             by = -1),
                         seq(from = last_non_zero + 1,
                             length.out = mask_start_out,
                             by = 1)))
      }
      # Ensure all masked bins are within the tx range being considered
      mask <- mask[mask >= 1 & mask <= length(tx)]
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
