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
#' @param ymax_bw maximum value of y-axis for bw files when plotting
#' @param ymax_abundance maximum value of y-axis for abundance when plotting
#'
#' @return plotted tracks of transcripts and masks of your query gene
#'
#' @include transcript_quantifier-class.R
#' @name plot_model
#' @export
plot_model <- function(tq,
                       gene_name = NULL,
                       chrom = NULL,
                       start = NULL,
                       end = NULL,
                       strand = NULL,
                       bigwig_plus = NULL,
                       bigwig_minus = NULL,
                       ymax_bw = NULL,
                       ymax_abundance = NULL) {
    # Check that only gene name or position information is specified
    if (!is.null(gene_name) & (!is.null(chrom) | !is.null(start) |
                               !is.null(end))) {
        stop("Only gene name OR positional information can be specified")
    }

    # Check for gene name specification
    if (!is.null(gene_name)) {
        gid_col <- tq@column_identifiers[2]
        if (is.na(gid_col)) {
            stop("transcript_quantifier object not built with gene ids")
        }
        if (!all(gene_name %in%
            S4Vectors::elementMetadata(tq@transcripts)[, gid_col])) {
            stop(paste("Gene", gene_name, "does not exist"))
        }
    } else if (!is.null(chrom) &
               is.numeric(start) & is.numeric(end)) {
        if (length(chrom) > 1 | length(start) > 1 | length(end) > 1) {
            stop("Only one value may be specified for each positional option")
        }
        # Check that position overlaps with one or more transcripts
        query_range <-
            GenomicRanges::GRanges(chrom, IRanges::IRanges(start, end), strand = strand)
        query_inter <-
            GenomicRanges::intersect(query_range, tq@transcripts,
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
        chrom <-
          S4Vectors::runValue(droplevels(GenomicRanges::seqnames(target_tx)[1]))
        start <- min(GenomicRanges::start(target_tx))
        end <- max(GenomicRanges::end(target_tx))
    }

    # genome coordination track
    axis_track <- Gviz::GenomeAxisTrack(target_tx)

    # Get masks
    tx_col <- tq@column_identifiers[1]
    gid_col <- tq@column_identifiers[2]
    transcripts <- S4Vectors::elementMetadata(target_tx)[, tx_col]
    target_masks <- get_masks(tq, transcripts)
    target_masks <- GenomicRanges::reduce(target_masks)

    # Transcripts track
    if (!is.na(gid_col)) {
        gene_names <- S4Vectors::elementMetadata(target_tx)[, gid_col]
    } else {
        gene_names <- rep("placeholder", length(target_tx))
    }
    unique_gene_names <- unique(unlist(gene_names))

    # Add columns for gene region track
    target_tx$gene <- unlist(GenomicRanges::values(target_tx)[[gid_col]])
    target_tx$feature <- target_tx$gene
    target_tx$transcript <- unlist(GenomicRanges::values(target_tx)[[tx_col]])

    tx_track <- Gviz::GeneRegionTrack(
        target_tx,
        name = "transcripts"
    )

    tx_track@dp@pars$shape <- "arrow"

    # Masks track
    masks_track <- Gviz::AnnotationTrack(target_masks,
                                         name = "masks",
                                         feature = Gviz::strand(target_masks),
                                         shape = "box",
                                         chromosome = chrom)
    # Some objects for data retrieval
    data_tracks <- list()
    bw_files <- c(`+` = bigwig_plus, `-` = bigwig_minus)
    strand_col <- c(`+` = "blue", `-` = "red")
    bw_max <- 0
    # Get count data it it exists
    if (any(!is.null(bw_files))) {
        for (s in names(bw_files)) {
            file <- bw_files[s]
            if (!file.exists(file)) {
                stop(paste("File", file, "does not exist"))
            }
            bw <- rtracklayer::import(file,
                       which = GenomicRanges::GRanges(
                           seqnames = chrom,
                           ranges = IRanges::IRanges(start = start, end = end)
                       ))
            bw$score <- abs(bw$score)

            # Get potential ymax
            bw_max <-
                max(bw_max, max(stats::quantile(abs(bw$score), 0.99))) * 1.05
            data_tracks[[s]] <- Gviz::DataTrack(
                                 range = bw,
                                 type = "h",
                                 window = -1,
                                 windowSize = tq@bin_size,
                                 name = paste0("PRO-seq (", s, ")"),
                                 col = strand_col[s],
                                 strand = s,
                                 chromosome = chrom)

        }
    } else {
        if (length(tq@counts) > 0) {
            data <- get_data(tq, chrom, start, end)
            bw_max <- stats::quantile(data$value, 0.99) * 1.05
            # plot datatrack for read count
            for (s in unique(GenomicRanges::strand(data))) {
                data$value[data$value == 0] <- NA
                data_tracks[[as.character(s)]] <-
                    Gviz::DataTrack(
                        data[GenomicRanges::strand(data) == s],
                        type = "h",
                        window = -1,
                        windowSize = tq@bin_size,
                        name = paste0("Summarized read counts (", s, ")"),
                        col = strand_col[as.character(s)],
                        strand = as.character(s),
                        chromosome = chrom
                    )
            }
        }
    }

    abundance_max <- 0
    abundance_tracks <- list()

    if (any(unlist(tq@model_abundance) != 0)) {
      # Get abundance data
      abundance <- get_abundance(tq, chrom, start, end)
      # Only keep transcripts with non-zero expression for visualization
      abundance <- abundance[, apply(mcols(abundance), 2, max) != 0]
      # plot datatrack for abundance
      for (s in unique(GenomicRanges::strand(abundance))) {
          abundance_tracks[[as.character(s)]] <- Gviz::DataTrack(
              abundance[GenomicRanges::strand(abundance) == s],
              type = "histogram",
              name = paste0("Predicted abundance (", s, ")"),
              groups = colnames(
                  S4Vectors::elementMetadata(
                      abundance[GenomicRanges::strand(abundance) == s])),
              legend = TRUE,
              strand = s,
              chromosome = chrom,
              stackedBars = TRUE
        )
      }
      abundance_max <-
          ceiling(abs(max(as.matrix(S4Vectors::mcols(abundance)))))
    }

    # Plot tracks
    args <- list(
        trackList = list(
            axis_track,
            tx_track,
            masks_track,
            data_tracks[["+"]],
            data_tracks[["-"]],
            abundance_tracks[["+"]],
            abundance_tracks[["-"]]
        ),
        from = start,
        to = end,
        chromosome = chrom,
        transcriptAnnotation = "transcript"
    )

    # Override default ymax with user specification if given
    if (!is.null(ymax_bw)) {
        bw_max <- ymax_bw
    }

    if (!is.null(ymax_abundance)) {
        abundance_max <- ymax_abundance
    }

    datatrack_names <- c("Summarized read counts (+)",
                         "Summarized read counts (-)",
                         "PRO-seq (+)", "PRO-seq (-)")

    # Set max for all valid data tracks
    for (track in seq_along(args$trackList)) {
        if (class(args$trackList[[track]]) == "DataTrack") {
            if (args$trackList[[track]]@name %in% datatrack_names) {
                args$trackList[[track]] <-
                    set_datatrack_ylim(args$trackList[[track]], c(0, bw_max))
            }
            if (args$trackList[[track]]@name %in%
                c("Predicted abundance (+)", "Predicted abundance (-)")) {
                args$trackList[[track]] <-
                    set_datatrack_ylim(args$trackList[[track]], c(0, abundance_max))
            }
        }
    }

    args$trackList <- #nolint
        args$trackList[!unlist(lapply(args$trackList, is.null))] #nolint
    # Add colors for genes and masks
    gene_colors <- viridisLite::viridis(length(unique_gene_names))
    names(gene_colors) <- unique_gene_names

    mask_colors <- c("blue", "red")
    names(mask_colors) <- c("+", "-")

    args <- c(args, gene_colors, mask_colors)

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
get_transcripts <- function(tq,
                            gene_name = NULL,
                            chrom = NULL,
                            start = NULL,
                            end = NULL,
                            strand = NULL) {
    # Look up transcripts for given gene name and assign chrom, start and end
    if (!is.null(gene_name)) {
        gid_col <- tq@column_identifiers[2]
        genes <- S4Vectors::elementMetadata(tq@transcripts)[, gid_col]
        tx <- tq@transcripts[genes %in% gene_name]
        chrom <- S4Vectors::runValue(GenomeInfoDb::seqnames(tx))[[1]]
        start <- min(BiocGenerics::start(tx))
        end <- max(BiocGenerics::end(tx))
    }

    query_range <- GenomicRanges::GRanges(chrom,
                                              IRanges::IRanges(start, end),
                                              strand)
    overlaps <- GenomicRanges::findOverlaps(query_range, tq@transcripts,
                                        ignore.strand = is.null(strand))
    out <- tq@transcripts[S4Vectors::subjectHits(overlaps)]

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
get_masks <- function(tq, transcripts) {
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
#' @description Retrives data for plotting. Returned ranges may go a bit outside
#' the queried ones if the overlapping bins do not perfectly align with the
#' user provided boundries
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
                       target_tx <-
                           get_transcripts(data_source,
                                           chrom = chrom,
                                           start = start,
                                           end = end)
                       # Lookup the relevant groups
                       tx_col <- data_source@column_identifiers[1]
                       transcripts <- S4Vectors::elementMetadata(target_tx)[, tx_col]
                       key <- data_source@transcript_model_key
                       target_group <- unique(key[key$tx_name %in% transcripts, ]$group)

                       # Combine the data and the GRanges object
                       value_granges <-
                         BiocGenerics::Reduce("c",
                                              mapply(
                                                function(bins, counts) {
                                                  bins$value <- counts
                                                  return(bins)
                                                },
                                                data_source@bins[target_group],
                                                data_source@counts[target_group],
                                                SIMPLIFY = FALSE
                                              ), init = GenomicRanges::GRanges()
                         )
                       return(value_granges)
                   })

#' @rdname get_data
methods::setMethod("get_data",
                   signature(data_source = "character"),
                   function(data_source, chrom, start, end) {
                       stop("Not implemented")
                   })

#' @title retrieve abundance for plotting
#'
#' @inheritParams plot_model
#'
#' @return a \link[GenomicRanges]{GRanges-class} objects with the adbundance in
#' the metadate column(s).
#'
#' @include transcript_quantifier-class.R
#' @name get_abundance
get_abundance <-
    function(tq, chrom, start, end) {
        target_tx <- get_transcripts(tq, chrom = chrom, start = start, end = end)
        # Lookup the relevant groups
        tx_col <- tq@column_identifiers[1]
        transcripts <- S4Vectors::elementMetadata(target_tx)[, tx_col]
        key <- tq@transcript_model_key
        target_group <- unique(key[key$tx_name %in% transcripts, ]$group)
        tx_models <- tq@models[target_group]
        tx_bins <- tq@bins[target_group]
        tx_abundances <- tq@model_abundance[target_group]

        tx_meta <-
            mapply(function(tx_model, tx_abundance) {
                t(t(tx_model) * tx_abundance)
            },
            tx_models,
            tx_abundances, SIMPLIFY = FALSE)

        value_granges <-
            BiocGenerics::Reduce("c", mapply(function(tx_bin, tx_meta) {
                GenomicRanges::values(tx_bin) <- tx_meta
                return(tx_bin)
            },
            tx_bins,
            tx_meta, SIMPLIFY = FALSE), init = GenomicRanges::GRanges())
        # Set NA to zero
        GenomicRanges::values(value_granges) <-
          apply(GenomicRanges::values(value_granges), 2, function(x) {
            x[is.na(x)] <- 0
            return(x)
          })
        return(value_granges)
    }

#' @title Set datatrack ylim
#'
#' @description Set ylim of the \link[Gviz]{DataTrack-class} object
#' @param data_track a \link[Gviz]{DataTrack-class} object
#' @param ylim a length two numeric vector holding the ylim boundaries.
#'
#' @return a \link[Gviz]{DataTrack-class} object
#' @name set_datatrack_ylim
set_datatrack_ylim <- function(data_track, ylim) {
    s <- Gviz::strand(data_track)
    s <- factor(as.character(s), levels = c("+", "-", "*"))
    scale <- do.call(c(I, rev, I)[as.numeric(s)][[1]], list(ylim))
    Gviz::displayPars(data_track)$ylim <- scale
    return(data_track)
}
