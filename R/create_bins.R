#' Create bins for multiple groups of overlapped transcripts
#'
#' Function for generating bins for mutiple groups of overlapped transcripts.
#'
#' @param gr_ls A \code{\link[GenomicRanges]{GRangesList-class}} object,
#' each element in the list contains transcripts from one gene or several genes
#' which are overlapped.
#' @param bin_size An integer, used to tail the gene region. Default is 50bp.
#' @return A \code{\link[GenomicRanges]{GRangesList-class}} object, containing
#' the binned region.
#' @export
 
create_bins <- function(gr_ls, bin_size = 50) {
    # check input class
    if (!methods::is(gr_ls, "GRangesList")) {
        stop("input is not a GRangesList object")
    } else if (!(methods::is(bin_size, "numeric") & bin_size > 0)) {
        stop("bin_size is not a positive number")
    }
    gr_bin_ls <- methods::as(lapply(gr_ls, get_gr_bin, bin_size = bin_size),
                             "GRangesList")
    return(gr_bin_ls)
}

#' Create bins for a group of overlapped transcripts
#'
#' Function for generating bins for a group of overlapped transcripts.
#'
#' @param gr A \code{\link[GenomicRanges]{GRanges-class}} object,
#' it contains transcripts from one gene or several genes
#' which are overlapped.
#' @param bin_size An integer, used to tail the gene region. Default is 50bp.
#' @importFrom S4Vectors runValue
#' @return A \code{\link[GenomicRanges]{GRanges-class}} object, containing
#' the binned region.

get_gr_bin <- function(gr, bin_size) {
    g_start <- min(GenomicRanges::start(gr))
    g_end <- max(GenomicRanges::end(gr))
    g_starts <- seq(g_start, g_end, by = bin_size)
    gr_binned <-
        data.frame(seqnames = runValue(GenomicRanges::seqnames(gr)),
                                start = g_starts,
                                end = g_starts + bin_size - 1,
                                strand = runValue(GenomicRanges::strand(gr)))
    gr_binned <- GenomicRanges::makeGRangesFromDataFrame(gr_binned)
    return(gr_binned)
}
