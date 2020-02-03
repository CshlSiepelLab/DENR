#' Create bins for multiple groups of overlapped transcripts
#'
#' Function for generating bins for mutiple groups of overlapped transcripts.
#'
#' @param transcript_groups A \code{\link[GenomicRanges]{GRangesList-class}}
#' object,each element in the list contains transcripts from one gene or several
#' genes which are overlapped.
#' @param bin_size An integer, used to tail the gene region. Default is 50bp.
#' @return A \code{\link[GenomicRanges]{GRangesList-class}} object, containing
#' the binned region.
#' @importFrom data.table :=
#' @export
 
create_bins <- function(transcript_groups, bin_size = 50) {
    # check input class
    if (!methods::is(transcript_groups, "GRangesList")) {
        stop("transcript_groups is not a GRangesList object")
    }
    if (!methods::is(bin_size, "numeric") || bin_size < 0) {
        stop("bin_size is not a positive number")
    }
    gr <- BiocGenerics::unlist(transcript_groups, use.names = TRUE)
    gr_dt <- data.table::as.data.table(gr)
    gr_dt[, select_groups := names(gr)]
    gr_bin_ls <- with(
        gr_dt[, bin_seq(seqnames[1],
                       start[1],
                       max(end),
                       strand[1],
                       bin_size), by = "select_groups"],
        GenomicRanges::split(
            GenomicRanges::GRanges(seqname, IRanges::IRanges(start, end),
                                   strand = strand),
            select_groups)
    )
    return(gr_bin_ls)
}

#' Create bins for a group of overlapped transcripts
#'
#' Function for generating bins for a group of overlapped transcripts.
#'
#' @param seqname chromsome name
#' @param start start (single value)
#' @param end end (single value)
#' @param strand strand (single value)
#' @inheritParams create_bins
#' @return A \code{\link[data.table]{data.table}} object, containing
#' the binned region.

bin_seq <- function(seqname, start, end, strand, bin_size) {
    g_starts <- seq(start,
                    length.out = ceiling((end - start) / bin_size),
                    by = bin_size)
    gr_binned <-
        data.table::data.table(seqname,
                               start = g_starts,
                                end = g_starts + bin_size - 1,
                                strand = strand)
    return(gr_binned)
}

## Appease R CMD check
if (getRversion() >= "2.15.1") {
    utils::globalVariables(c("select_groups", "seqnames", "start",
                             "end", "strand"))
}
