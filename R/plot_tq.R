#' @title Plot transcript_abundance object
#'
#' @description plot transcripts and masks
#'
#' @inheritParams add_data
#' @param gene_name Name of your query gene. Only works if gene ids were given
#' to the \link{transcript_quantifier-class} object
#' @param chrom the chromosome of your query region (single value)
#' @param start the start position of your query region (single value)
#' @param end the end position of your query region (single value)
#' @param strand the strand specificity of your query region (default: ANY)
#'
#' @return plotted tracks of transcripts and masks of your query gene
#'
#' @include transcript_quantifier-class.R
#' @name plot_model
#' @export
plot_model <- function(transcript_quantifier, gene_name = NULL,
                       chrom = NULL, start = NULL, end = NULL,
                       strand = NULL) {
  # Alias for easier use
  tq <- transcript_quantifier
  # Check that only gene name or position information is specified
  if (!is.null(gene_name) & (!is.null(chrom) | ! is.null(start) |
                            !is.null(end))) {
    stop("Only gene name OR positional information can be specified")
  }

  # Check for gene name specification
  if (!is.null(gene_name)) {
    gid_col <- tq@column_identifiers[2]
    if (is.na(gid_col)) {
      stop("transcript_quantifier object not built with gene ids")
    }
    if (!gene_name %in%
       S4Vectors::elementMetadata(tq@transcripts)[, gid_col]) {
      stop(paste("Gene", gene_name, "does not exist"))
    }
  } else if (!is.null(chrom) & is.numeric(start) & is.numeric(end)) {
    if (length(chrom) > 1 | length(start) > 1 | length(end) > 1) {
      stop("Only one value may be specified for each positional option")
    }
    # Check that position overlaps with one or more transcripts
    query_range <- GenomicRanges::GRanges(chrom, IRanges::IRanges(start, end),
                                          strand = strand)
    query_inter <- GenomicRanges::intersect(query_range, tq@transcripts,
                                            ignore.strand = is.null(strand))
    if (length(query_inter) == 0) {
      stop("Positional query does not intersect any transcripts")
    }
  } else {
    stop("Incorrect positional specification")
  }

  ## ** End checking **

  # Get target transcripts
  target_tx <- get_transcripts(tq, gene_name, chrom, start, end, strand)

  # Define bounds of plot range
  if (is.null(start)) {
    chrom <- GenomicRanges::seqnames(target_tx)[1]
    start <- min(GenomicRanges::start(target_tx))
    end <- max(GenomicRanges::end(target_tx))
  }

  # Get masks
  tx_col <- tq@column_identifiers[1]
  gid_col <- tq@column_identifiers[2]
  transcripts <- S4Vectors::elementMetadata(target_tx)[, tx_col]
  target_masks <- get_masks(tq, transcripts)
  GenomicRanges::strand(target_masks) <- "*"

  # Get count data
  data <- get_data(tq, chrom, start, end)
  data <- data[data$value > 0]

  # Transcripts track
  if (!is.na(gid_col)) {
    gene_names <- S4Vectors::elementMetadata(target_tx)[, gid_col]
  } else {
    gene_names <- rep("placeholder", length(target_tx))
  }
  unique_gene_names <- unique(unlist(gene_names))

  tx_track <- Gviz::AnnotationTrack(target_tx, name = "transcripts",
                                    id = transcripts,
                                    feature = unlist(tq@transcripts$gene_id))
  Gviz::feature(tx_track)

  # Masks track
  masks_track <- Gviz::AnnotationTrack(target_masks, name = "masks")

  # Data plotting track
  d_track_plus <- NULL
  d_track_minus <- NULL
  if (length(data[GenomicRanges::strand(data) == "+"]) > 0) {
    d_track_plus <- Gviz::DataTrack(data[GenomicRanges::strand(data) == "+"],
                                    type = "h",
                                    name = "Summarized read counts (+)",
                                    col = "red"
                                    )
  }
  if (length(data[GenomicRanges::strand(data) == "-"]) > 0) {
    data[GenomicRanges::strand(data) == "-"]$value <-
      -abs(data[GenomicRanges::strand(data) == "-"]$value)
    d_track_minus <- Gviz::DataTrack(data[GenomicRanges::strand(data) == "-"],
                                     type = "h",
                                     name = "Summarized read counts (-)",
                                     col = "blue"
                                     )
  }

  # Plot tracks
  args <- list(trackList = list(tx_track,
                                masks_track,
                                d_track_plus,
                                d_track_minus),
               from = start, to = end,
               featureAnnotation = "feature")
  args$trackList <- #nolint
    args$trackList[!unlist(lapply(args$trackList, is.null))] #nolint
  gene_color_list <- viridisLite::viridis(length(unique_gene_names))
  names(gene_color_list) <- unique_gene_names
  args <- c(args, gene_color_list)

  return(do.call(Gviz::plotTracks, args))
}

#' @title Get target transcripts
#'
#' @description Retrives transcripts corresponding to a gene or genomic region
#' as a \link[GenomicRanges]{GRanges-class} object
#' @inheritParams plot_model
#'
#' @return a \link[GenomicRanges]{GRanges-class} object
#' @name get_transcripts
get_transcripts <- function(transcript_quantifier, gene_name = NULL,
                            chrom = NULL, start = NULL,
                            end = NULL, strand = NULL) {
  # Alias for ease of use
  tq <- transcript_quantifier
  # Look up transcripts
  if (!is.null(gene_name)) {
    gid_col <- tq@column_identifiers[2]
    genes <- S4Vectors::elementMetadata(tq@transcripts)[, gid_col]
    out <- tq@transcripts[genes == gene_name]
  } else {
    query_range <- GenomicRanges::GRanges(chrom,
                                          IRanges::IRanges(start, end),
                                          strand)
    overlaps <- GenomicRanges::findOverlaps(query_range, tq@transcripts,
                                            ignore.strand = is.null(strand))
    out <- tq@transcripts[S4Vectors::subjectHits(overlaps)]
  }
  return(out)
}

#' @title Get target masks
#'
#' @description Retrives masks corresponding to a gene or genomic region
#' as a \link[GenomicRanges]{GRanges-class} object
#'
#' @inheritParams plot_model
#' @param transcripts a character vector of transcript names
#'
#' @return a \link[GenomicRanges]{GRanges-class} object
#' @name get_masks
get_masks <- function(transcript_quantifier, transcripts) {
  tq <- transcript_quantifier
  key <- tq@transcript_model_key
  target_group <- unique(key[key$tx_name %in% transcripts, ]$group)
  target_masks <- GenomicRanges::GRanges()
  for (g in target_group) {
    target_masks <- c(target_masks, tq@bins[[g]][tq@masks[[g]]])
  }
  return(target_masks)
}

#' @title Get data
#'
#' @description Retrives data for plotting
#'
#' @param data_source a \link{transcript_quantifier-class} object or string
#' giving the path to a bigwig file
#' @inheritParams plot_model
#'
#' @return a \link[GenomicRanges]{GRanges-class} objects with the read counts in
#' the \code{value} column
#'
#' @include transcript_quantifier-class.R
#' @name get_data
#' @rdname get_data
methods::setGeneric("get_data",
                    function(data_source, chrom, start, end) {
                      standardGeneric("get_data")
                    })

#' @rdname get_data
methods::setMethod("get_data",
                   signature(data_source = "transcript_quantifier"),
                   function(data_source, chrom, start, end) {
                     #Get the target transcripts
                     target_tx <- get_transcripts(data_source, chrom = chrom,
                                                  start = start,
                                                  end = end)
                     # Lookup the relevant groups
                     tx_col <- data_source@column_identifiers[1]
                     transcripts <-
                       S4Vectors::elementMetadata(target_tx)[, tx_col]
                     key <- data_source@transcript_model_key
                     target_group <-
                       unique(key[key$tx_name %in% transcripts, ]$group)

                     # Combine the data and the GRanges object
                     value_granges <- BiocGenerics::Reduce("c", mapply(
                       function(bins, counts) {
                         bins$value <- counts
                         return(bins)
                       },
                       data_source@bins[target_group],
                       data_source@counts[target_group],
                       SIMPLIFY = FALSE)
                     )
                     return(value_granges)
                   }
)

#' @rdname get_data
methods::setMethod("get_data",
                   signature(data_source = "character"),
                   function(data_source, chrom, start, end) {
                     stop("Not implemented")
                   }
)
