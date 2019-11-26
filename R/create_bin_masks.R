#' Create data masks
#'
#' Creates a set of bin masks for transcript starts and ends
#'
#' @param txdb a \code{\link[GenomicFeatures]{TxDb-class}} object
#' @param bins a \code{\link[GenomicRanges]{GRangesList-class}} object
#' @param mask_start_bins The number of bins at the start of each transcript to
#' mask
#' @param mask_end_bins The number of bins at the end of each transcript to
#' mask
#' @rdname create_bin_masks
#' @return A list of bit vectors
#' @export
create_bin_masks <- function(txdb, bins, mask_start_bins=0, mask_end_bins=0) {
    print(txdb)
}
