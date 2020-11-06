setOldClass("loess")

## Appease R CMD check, so many becasue of data.table
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("group", "percent_match", ".SD", "bin_start",
                           "bin_end", "percent_transcribed",
                           "expected_bin_start", "tx_start", "tx_end",
                           "expected_bin_end", "y", "region", "."))
}

#' Class transcript_shape_profile
#'
#' Class \code{transcript_shape_profile} holds data for a set of transcripts
#' that are used to generate an empirical shape profile. Can be created by
#' calling \code{transcript_shape_profile}. A shape profile can be used with
#' \code{modify_model} to alter the models in a \code{transcript_quantifier}
#' object.
#'
#' @slot transcripts a \link[GenomicRanges]{GRanges-class} object that contains
#' the ranges of transcripts used to fit the empirical shape profile.
#' @slot bin_size bin size used to create transcript profile functions
#' @slot linear_head_length distance into gene from TSS that is treated as the
#' head region for the purposes of position scaling
#' @slot linear_tail_length distance into gene from the transcript end that is
#' treated as the tail region for the purposes of position scaling
#' @slot head_percent The head region of the gene is rescaled so the position
#' vector that is paired with the count vector is between
#' \code{[0, head_percent]}
#' @slot tail_percent The tail region of the gene is rescaled so the position
#' vector that is paired with the count vector is between
#' \code{[1 - tail_percent, 1]}
#' @slot shape_profile a loess function that is used to predict the shape
#' profile
#' @slot shape_scale_factor a scalar value to multiply the output of the loess
#' function by so that the profile factors predicted by the loess fit have a
#' median value of 1 in the central body of the gene
#'
#' @name transcript_shape_profile-class
#' @rdname transcript_shape_profile-class
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom GenomicRanges CompressedGRangesList
#' @exportClass transcript_quantifier
methods::setClass("transcript_shape_profile",
                  slots = c(transcripts = "GRanges",
                            bin_size = "integer",
                            linear_head_length = "integer",
                            linear_tail_length = "integer",
                            head_percent = "numeric",
                            tail_percent = "numeric",
                            shape_profile = "loess",
                            shape_scale_factor = "numeric")
)


#' Creates transcript_shape_profile
#'
#' Creates object of \code{transcript_shape_profile}, which holds data for a set
#' of transcripts that are used to generate an empirical shape profile. That
#' shape profile can then be used with \code{modify_model} to alter the models
#' in a \code{transcript_quantifier} object. To get the list of transcripts to
#' create the shape profile a series of filters are applied and then the longest
#' of the remaining transcripts for each overlapping set of transcripts is
#' chosen. Then each transcript is binned and ...
#'
#' @param transcripts a \link[GenomicRanges]{GRanges-class} object that contains
#' the ranges of potential transcripts to be used to fit an empirical shape
#' profile. This list will be filtered using the criteria specified by the
#' \code{min_transcript_length}, \code{max_group_size}, and
#' \code{percent_model_match} parameters
#' @param min_percent_transcribed all transcripts with pro-seq signal for less
#' than this percent of the transcript are removed
#' @param min_transcript_length all transcripts shorter than this are removed
#' from consideration (default: 5e3)
#' @param percent_model_match Based on the transcript ranges for each transcript
#' certain regions are expected to contain nascent RNAP signal or not. For each
#' transcript the percent agreement with that expectation is computed and all
#' transcripts below that are removed (default: 0.8).
#'
#' @inheritParams rescale_fixed_width_position
#' @inheritParams create_bins
#' @inheritParams add_data
#' @inheritParams stats::loess
#' @importFrom data.table :=
#'
#' @name transcript_shape_profile
#' @rdname transcript_shape_profile
#' @export
transcript_shape_profile <- function(transcripts,
                                     bigwig_plus = NULL, bigwig_minus = NULL,
                                     bin_size = 250,
                                     min_percent_transcribed = 0.8,
                                     min_transcript_length = 8e3,
                                     percent_model_match = 0.8,
                                     linear_head_length = 10e3,
                                     linear_tail_length = 5e3,
                                     span = 0.4) {
  ## Input checking ##

  if (class(transcripts) != "GRanges") {
    stop("transcripts is not a GRanges object")
  }

  if (bin_size >= linear_head_length || bin_size >= linear_tail_length) {
    stop("bin_size must be less than linear head/tail length")
  }

  ## End input checking ##

  # 0) Remove unused seqlevels and seqlevels that are not in bigwig
  GenomeInfoDb::seqlevels(transcripts) <-
    GenomeInfoDb::sortSeqlevels(GenomeInfoDb::seqlevelsInUse(transcripts))
  bw_seqlevels <- intersect(
    GenomeInfoDb::seqlevels(rtracklayer::BigWigFile(bigwig_minus)),
    GenomeInfoDb::seqlevels(rtracklayer::BigWigFile(bigwig_plus))
  )

  tx_count_0 <- length(transcripts)
  message(paste("Initial transcripts:", tx_count_0))

  drop_lvls <- setdiff(GenomeInfoDb::seqlevels(transcripts), bw_seqlevels)
  if (length(drop_lvls) > 0) {
    transcripts <- GenomeInfoDb::dropSeqlevels(transcripts, drop_lvls,
                                               pruning.mode = "coarse")
    tx_count_1 <- length(transcripts)
    message(
      paste("There are", tx_count_0 - tx_count_1,
             "transcripts from a seqlevel not present in both bigwig files.",
             "Dropping those transcripts")
    )
    message("Transcripts after seqlevels drop: ", tx_count_1)
  }

  # 1) Filter out short transcripts
  transcripts <- transcripts[GenomicRanges::width(transcripts) >= min_transcript_length]
  tx_count_2 <- length(transcripts)
  message(paste("Transcripts after size filter:", tx_count_2))
  if (tx_count_2 == 0) {
    stop("No transcripts remaining")
  }

  # 2) Group transcripts and filter by group size
  tx_grps <- group_transcripts(transcripts, distance = 0)
  group_count_1 <- length(tx_grps)
  message(paste("Groups formed:", group_count_1))

  # 3) Get model agreement for each transcript
  message("Computing GOF per transcript given data ...")
  mod_agree <- model_agreement(tx_grps, bigwig_plus, bigwig_minus)

  # 4) Select final candidate transcript, one per group based on best matched
  #    one
  message("Filtering transcripts by minimum GOF and % transcribed ...")
  ma_dt <- data.table::as.data.table(mod_agree)
  ma_dt[, group := names(mod_agree)]
  ma_dt <- ma_dt[percent_transcribed >= min_percent_transcribed &
                   percent_match >= percent_model_match]
  message("Groups remaining: ", length(unique(ma_dt$group)))
  if (nrow(ma_dt) == 0) {
    stop("No transcripts passed filters to be used for model fitting")
  }
  message("Selecting best fitting transcript per group ...")
  final_tx <- ma_dt[, .SD[which.max(percent_match)], by = "group"]

  # 5) Convert candidates to GRanges
  final_tx_gr <- with(final_tx,
                      GenomicRanges::GRanges(seqnames,
                                             IRanges::IRanges(start, end),
                                             strand = strand,
                                             percent_transcribed,
                                             percent_match
                      )
  )

  # 5) Group and bin transcripts
  message("Collecting data for profile model ...")
  final_tx_grp <- group_transcripts(final_tx_gr)
  final_bins <- create_bins(final_tx_grp, bin_size = bin_size)
  plus_grps <- grep("+", names(final_bins))
  minus_grps <- grep("-", names(final_bins))

  # 6) Extract the relevant scores
  plus_counts <- summarize_bigwig(bigwig_file = bigwig_plus, final_bins[plus_grps])
  minus_counts <- lapply(summarize_bigwig(bigwig_file = bigwig_plus,
                                  final_bins[minus_grps]), rev) # reverse minus
  raw_counts <- c(plus_counts, minus_counts)

  # 7) Rescale positions so that data can be plotted on common x-axis
  head_percent <- 0.2 # Fraction of range dedicated to head region
  tail_percent <- 0.2 # Fraction of range dedicated to tail region
  rescaled_position_data <- data.table::rbindlist(
    lapply(raw_counts, function(x) {
      pos <-
        rescale_fixed_width_position(bin_size = bin_size,
                                     length_out = length(x),
                                     linear_head_length = linear_head_length,
                                     linear_tail_length = linear_tail_length,
                                     head_percent = head_percent,
                                     tail_percent = tail_percent
        )
      # Subtract small value for min/max scaling so that it doesn't fail
      # when min(x) == max(x)
      return(data.table::data.table(
        pos,
        scaled_count = (x - (min(x) - 1e-6)) / (max(x) - (min(x) - 1e-6)))
      )
    })
  )

  # 8) Fit profile model
  message(paste("Fitting profile model using", length(raw_counts), "genes..."))
  ctrl <- stats::loess.control()
  ctrl[["trace.hat"]] <- "approximate"
  profile_function <- stats::loess(scaled_count ~ pos,
                            data = rescaled_position_data, span = span,
                            control = ctrl)
  # Strip out unnecessary bits of loess object to reduce memory footprint
  profile_function <- strip_loess(profile_function)

  # 9) Compute scaling factor for profile model
  tmp <- data.table::data.table(
    pos = seq(head_percent, 1 - tail_percent, length.out = 1e4)
  )
  out <- stats::predict(profile_function, tmp)
  shape_scale_factor <- 1 / stats::median(out)

  return(methods::new("transcript_shape_profile",
                      transcripts = final_tx_gr,
                      bin_size = as.integer(bin_size),
                      linear_head_length = as.integer(linear_head_length),
                      linear_tail_length = as.integer(linear_tail_length),
                      head_percent = head_percent,
                      tail_percent = tail_percent,
                      shape_profile = profile_function,
                      shape_scale_factor = shape_scale_factor))
}

#' Computes GOF statistic for groups of transcripts
#'
#'
#' Annotates transcriptes from the input \link[GenomicRanges]{GRangesList-class}
#' object with two additional columns: \code{percent_match} and
#' \code{percent_transcribed}. Within the genomic interval covered by each group
#' of transcripts, certain bases are predicted to be transcribed or not for each
#' transcript annotation. Based on the supplied bigwigs each base in the
#' relevant region is then annotated as containig polymerase (being transcribed)
#' or not. \code{percent_match} is the fraction of bases whose state, as
#' predicted by the transcript annotations matches the annotation from the data.
#' \code{percent_transcribed} is then the percent of the region predicted by the
#' transcript annotation to be transcribed that actually is.
#'
#' @param tx_grps a \link[GenomicRanges]{GRangesList-class} object that contains
#' groups of overlapping transcripts
#'
#' @inheritParams add_data
#'
#' @name model_agreement
model_agreement <- function(tx_grps, bigwig_plus, bigwig_minus) {
  # Pre-compute regions that will need to be imported from bigwig
  tx <- unlist(tx_grps)
  import_range <- GenomicRanges::reduce(tx)
  import_range <- GenomicRanges::split(import_range,
                                       GenomicRanges::strand(import_range))

  # Get GRanges objects of what is considered transcribed on each strand
  # according to the bigwig file
  transcribed_thresh <- 1 # >= this is considered transcribed, control of this
                          # is not acessible to user for now
  transcribed_plus <- get_transcribed_regions(import_range[["+"]],
                                              bigwig_plus, transcribed_thresh)
  GenomicRanges::strand(transcribed_plus) <- "+"
  transcribed_minus <- get_transcribed_regions(import_range[["-"]],
                                              bigwig_minus, transcribed_thresh)
  GenomicRanges::strand(transcribed_minus) <- "-"
  transcribed <- c(transcribed_plus, transcribed_minus)
  # remove unused seqlevels
  GenomeInfoDb::seqlevels(transcribed) <-
      GenomeInfoDb::seqlevelsInUse(transcribed)
  # use complement to get untranscribed regions
  untranscribed <- GenomicRanges::gaps(transcribed)
  untranscribed <- untranscribed[GenomicRanges::strand(untranscribed) != "*"]

  # Get start-stop of joint ranges for each group
  group_ranges <- GenomicRanges::GRanges(
    seqnames = unlist(S4Vectors::runValue(GenomicRanges::seqnames(tx_grps))),
    IRanges::IRanges(min(GenomicRanges::start(tx_grps)),
                     max(GenomicRanges::end(tx_grps))),
    strand = unlist(S4Vectors::runValue(GenomicRanges::strand(tx_grps)))
  )
  group_sizes <- S4Vectors::elementNROWS(tx_grps)
  group_ranges_repl <- rep(group_ranges, group_sizes)
  names(group_ranges) <-  names(tx_grps)

  # Create ranges of the left and right untranscribed regions for each
  # transcript
  chr <- GenomicRanges::seqnames(tx)
  strands <- GenomicRanges::strand(tx)
  tx_untrans_left <- GenomicRanges::GRanges(chr, IRanges::IRanges(
    GenomicRanges::start(group_ranges_repl), GenomicRanges::start(tx) - 1,
    names = names(tx)
  ), strand = strands)
  tx_untrans_right <- GenomicRanges::GRanges(chr, IRanges::IRanges(
    GenomicRanges::end(tx) + 1, GenomicRanges::end(group_ranges_repl),
    names = names(tx)
  ), strand = strands)

  # Get the number of bins expected to be transcribed that are untranscribed
  unexpected_untrans <- covered_bases(query = tx, subject = untranscribed)

  # Get the number of bins expected to be untranscribed that are transcribed
  unexpected_trans_left <- covered_bases(tx_untrans_left, transcribed)
  unexpected_trans_right <- covered_bases(tx_untrans_right, transcribed)

  # Compute the expected transcribed and untranscribed bits too
  percent_trans <- covered_bases(tx, transcribed) / GenomicRanges::width(tx)

  # Compute the total number of errors for each transcript
  total_err <-
    unexpected_trans_left + unexpected_trans_right + unexpected_untrans
  percent_err <- total_err / GenomicRanges::width(group_ranges[names(tx)])
  percent_match <- 1 - percent_err
  tx$percent_transcribed <- percent_trans
  tx$percent_match <- percent_match
  return(tx)
}

#' Annotates transcribed regions
#'
#'
#' Returns \link[GenomicRanges]{GRanges-class} that contains all transcribed
#' regions within \code{ranges} with respect to \code{bw}
#'
#' @param ranges a \link[GenomicRanges]{GRanges-class}
#' @param bw the path to a bigwig file
#' @param transcribed_thresh the amount of signal required for a region to be
#' considered transcribed
#'
#' @name get_transcribed_regions
get_transcribed_regions <- function(ranges, bw, transcribed_thresh = 1) {
  # Import data for each strand and convert to reduced GRanges that indicates
  # regions that are considered "transcribed"
  imported_bw <- rtracklayer::import.bw(
    con = bw, which = ranges,
    as = "GRanges")
  transcribed <- GenomicRanges::reduce(
    imported_bw[abs(GenomicRanges::score(imported_bw)) >= transcribed_thresh]
  )
  return(transcribed)
}

#' Gets per-query coverage
#'
#'
#' Returns vector with number of bases overlapped by the subject intervals per
#' query interval.
#'
#' @param query a \link[GenomicRanges]{GRanges-class}
#' @param subject a \link[GenomicRanges]{GRanges-class}
#'
#' @name covered_bases
covered_bases <- function(query, subject) {
  # Pre-allocate output vector with intial value set to zero
  covered <- integer(length(query))
  # Reduce subject to guarentee unique overlaps
  subject <- GenomicRanges::reduce(subject)
  # Save initial query rows and split by strand and chromosome
  query$uniq_row_id <- seq_along(query)
  query_l <- GenomicRanges::split(query, GenomicRanges::seqnames(query))
  # Get the coverages to the subject per strand level
  sub_covr_plus <-
    GenomicRanges::coverage(subject[GenomicRanges::strand(subject) == "+"])
  sub_covr_minus <-
    GenomicRanges::coverage(subject[GenomicRanges::strand(subject) == "-"])
  sub_covr_star <-
    GenomicRanges::coverage(subject[GenomicRanges::strand(subject) == "*"])
  # Iterate over the seqlevels where there is coverage
  for (s in GenomeInfoDb::seqlevels(subject)) {
    q_split <- GenomicRanges::split(query_l[[s]], GenomicRanges::strand(query_l[[s]]))
    # Subjects on the +/* strand cover + query items
    covered[q_split[["+"]]$uniq_row_id] <- unlist(
      IRanges::viewSums(IRanges::Views(sub_covr_plus[[s]],
                                       IRanges::ranges(q_split[["+"]]))) +
        IRanges::viewSums(IRanges::Views(sub_covr_star[[s]],
                                         IRanges::ranges(q_split[["+"]])))
    )
    # Subjects on the -/* strand cover - query items
    covered[q_split[["-"]]$uniq_row_id] <- unlist(
      IRanges::viewSums(IRanges::Views(sub_covr_minus[[s]],
                                       IRanges::ranges(q_split[["-"]]))) +
        IRanges::viewSums(IRanges::Views(sub_covr_star[[s]],
                                         IRanges::ranges(q_split[["-"]])))
    )
    # Subjects on the +/-/* strand cover * query items
    covered[q_split[["*"]]$uniq_row_id] <- unlist(
      IRanges::viewSums(IRanges::Views(sub_covr_plus[[s]],
                                       IRanges::ranges(q_split[["*"]]))) +
        IRanges::viewSums(IRanges::Views(sub_covr_minus[[s]],
                                         IRanges::ranges(q_split[["*"]]))) +
        IRanges::viewSums(IRanges::Views(sub_covr_star[[s]],
                                         IRanges::ranges(q_split[["*"]])))
    )
  }
  return(covered)
}

#' Tripartite linear rescaling
#'
#' Given a length, a fixed bin width, and head/tail lengths this function
#' generates a sequence [0, 1] that is composed of three linear sub-scales for
#' the head, mid, and tail region. These scales are created such that positions
#' <= \code{linear_head_length} are converted to a uniform sequence between
#' \code{[0, head_percent]} and values
#' \code{[bin_size * length_out - linear_tail_length, bin_size]} are converted
#' to a uniform sequence \code{[tail_percent, 1]} with the mid values falling
#' in the remaining range.
#'
#' @param bin_size width of each position range that is being converted
#' @param length_out length of position vector to be created
#' @param linear_head_length total width of head region (# of bins * width)
#' @param linear_tail_length total width of tail region (# of bins * width)
#' @param head_percent fraction of [0, 1] to be occupied by head region
#' @param tail_percent fraction of [0, 1] to be occupied by tail region
#' @param short_heuristic whether sequences shorter than the sum of the head and
#' tail length should be accommodated by proportionally rescaling the head and
#' tail lengths. Otherwise an error will be thrown.
#'
#' @name rescale_fixed_width_position
rescale_fixed_width_position <- function(bin_size, length_out,
                                         linear_head_length, linear_tail_length,
                                         head_percent = 0.2,
                                         tail_percent = 0.2,
                                         short_heuristic = TRUE) {
  # total width of all bins
  total_width <- bin_size * length_out

  # Handling of cases where total sequence width is too short
  if (linear_head_length + linear_tail_length > total_width) {
    # If the short heuristic is allowed, just divide up the bins proportionally
    if (short_heuristic) {
      if (length_out == 1) {
        return(0)
      } else {
        prop_head <- linear_head_length /
          (linear_head_length + linear_tail_length)
        linear_head_pos <- min(round(prop_head * length_out), length_out - 1)
        linear_tail_pos <- linear_head_pos + 1
      }
    } else {
      stop(paste("total_width too short for rescaling",
                 "(<linear_head_length + linear_tail_length)"))
    }
  } else {
    # Get the position on which linear scaling should end and restart
    # If the head/tail is not an exact multiple of bin width, round up
    linear_head_pos <- ceiling(linear_head_length / bin_size)
    linear_tail_pos <- length_out - ceiling(linear_tail_length / bin_size) + 1
  }

  # If the transcript is too short such that the tail runs in to the
  # head, then just set the tail start position to be right after the head
  if (linear_tail_pos < linear_head_pos) {
    linear_tail_pos <- linear_head_pos + 1
  }

  # Generate a rescaled sequence that covers [0, 1] when the various
  # linear scalings are considered
  rescaled_pos <- numeric(length_out)
  rescaled_pos[1:linear_head_pos] <- seq(from = 0, to = head_percent,
                         length.out = linear_head_pos)
  rescaled_pos[linear_tail_pos:length_out] <-
    seq(from = 1 - tail_percent, to = 1, length.out = length_out - linear_tail_pos + 1)
  mid_pos_len <- linear_tail_pos - linear_head_pos - 1 + 2
  # Trim off extra start and end position of mid_pos that were added so that
  # first and last positions did not overlap the head and tail
  if (mid_pos_len > 2) {
    rescaled_pos[(linear_head_pos + 1):(linear_tail_pos - 1)] <-
      seq(from = head_percent, to = 1 - tail_percent,
          length.out = mid_pos_len)[c(-1, -mid_pos_len)]
  }
  return(rescaled_pos)
}

#' @title View transcript shape profile
#'
#' View transcript shape model natively or scaled to a specific gene length
#'
#' @param tsp A \link{transcript_shape_profile-class} object
#' @param gene_length A gene length, defaults to just plotting profile on
#' @param ... for extensibility
#' \code{[0,1]}
#' axis
#'
#' @name view
#' @rdname view
methods::setGeneric("view",
                    function(tsp, gene_length = NULL, ...) {
                      standardGeneric("view")
                    })

#' @rdname view
#' @export
methods::setMethod("view",
                   signature(tsp = "transcript_shape_profile"),
                   function(tsp, gene_length = NULL, ...) {
                     if (!is.null(gene_length)) {
                       pos <- rescale_fixed_width_position(
                         bin_size = tsp@bin_size,
                         length_out = round(gene_length / tsp@bin_size),
                         linear_head_length = tsp@linear_head_length,
                         linear_tail_length = tsp@linear_tail_length,
                         head_percent = tsp@head_percent,
                         tail_percent = tsp@tail_percent
                      )
                     } else {
                       grid_points <- 1e3 # default resolution for plot
                       pos <- seq(0, 1, length.out = grid_points)
                     }
                     # Predict and rescale shape profile
                     pred <- data.table::data.table(pos = pos)
                     pred[, y := stats::predict(tsp@shape_profile, pred)]
                     pred[, y := y * tsp@shape_scale_factor]

                     if (!is.null(gene_length)) {
                       pred[, pos := seq(0, by = tsp@bin_size,
                                         length.out = round(gene_length / tsp@bin_size))]
                       xlabel <- "Location (kb)"
                       pred[pos <= tsp@linear_head_length, region := "head"]
                       pred[pos >= gene_length - tsp@linear_tail_length,
                            region := "tail"]
                       pred[is.na(region), region := "mid"]
                       pred[, pos := pos / 1e3]
                     } else {
                       xlabel <- "Relative location"
                       pred[pos <= tsp@head_percent, region := "head"]
                       pred[pos >= 1 - tsp@tail_percent, region := "tail"]
                       pred[is.na(region), region := "mid"]
                     }

                     bin_label <- paste("bin width:", tsp@bin_size)

                     g <- ggplot2::ggplot(data = pred,
                                          ggplot2::aes(x = pos, y = y)) +
                       ggplot2::geom_line(ggplot2::aes(color = region)) +
                       ggplot2::theme_bw() +
                       ggplot2::ylim(0, max(pred$y)) +
                       ggplot2::xlab(xlabel) +
                       ggplot2::ylab("Relative height") +
                       ggplot2::geom_label(
                         mapping = ggplot2::aes(x = 0.03 * max(pred$pos),
                                      y = 0.2, label = bin_label),
                         hjust = 0, vjust = 0,
                         size = 4) +
                       ggplot2::scale_color_discrete(
                         guide = ggplot2::guide_legend(
                         direction = "horizontal"
                       )) +
                       ggplot2::theme(
                         legend.position = c(0.5, 0.2 / max(pred$y)),
                         legend.justification = c("center", "bottom"),
                         legend.box.just = "center",
                         text = ggplot2::element_text(size = 16)
                       )
                     return(g)
                   }
)

#' @title Change shape of transcript models
#'
#' Create a \code{\link{transcript_quantifier-class}} that uses the shape
#' profile from a \code{\link{transcript_shape_profile-class}}
#'
#' @param tq A \code{\link{transcript_quantifier-class}} object
#' @param tsp A \link{transcript_shape_profile-class} object
#'
#' @name apply_shape_profile
#' @rdname apply_shape_profile
#' @export
methods::setGeneric("apply_shape_profile",
                    function(tq, tsp) {
                      standardGeneric("apply_shape_profile")
                    })

#' @rdname apply_shape_profile
methods::setMethod("apply_shape_profile",
  signature(tq = "transcript_quantifier",
            tsp = "transcript_shape_profile"),
  function(tq, tsp) {
    if (tq@bin_size != tsp@bin_size) {
      stop("bin_size for transcript_quantifier and transcript_shape_profile ",
           "are not equal")
    }

    message("Caching shape profile function ...")
    # Pre-compute loess values at fine grid resolution for fast querying
    pre_profile <- stats::predict(
      tsp@shape_profile,
      newdata = seq(0, 1, by = 1e-05)) * tsp@shape_scale_factor
    grid_length <- length(pre_profile)

    message("Calculating transcript positional statistics ...")
    # Store group starts
    bin_starts <- GenomicRanges::start(tq@bins)
    start_idx <- c(1, cumsum(S4Vectors::elementNROWS(bin_starts)) + 1)
    start_idx <- start_idx[-length(start_idx)]
    grp_start <- BiocGenerics::unlist(bin_starts)[start_idx]
    # Store strands
    grp_strand <- unlist(S4Vectors::runValue(GenomicRanges::strand(tq@bins)))

    # Store transcript info
    tx_info <- data.table::data.table(
      tx_start = GenomicRanges::start(tq@transcripts),
      tx_end = GenomicRanges::end(tq@transcripts),
      tx_strand = BiocGenerics::as.vector(GenomicRanges::strand(tq@transcripts)),
      tx_name = GenomicRanges::values(tq@transcripts)[, tq@column_identifiers[1]])

    tx_info <- merge(
      tx_info,
      data.table::as.data.table(tq@transcript_model_key),
      by = "tx_name"
    )
    data.table::setkey(tx_info, "group")

    # Compute expected per model transcript start and end bins
    # Average over transcripts when there are multiple transcripts in a single
    # model
    tx_info[, expected_bin_start := mean((tx_start - grp_start[group] + 1) /
                                    tq@bin_size), by = c("group", "model")]
    tx_info[, expected_bin_start := ceiling(expected_bin_start)]
    tx_info[, expected_bin_end := mean((tx_end - grp_start[group] + 1) /
                                           tq@bin_size), by = c("group", "model")]
    tx_info[, expected_bin_end := ceiling(expected_bin_end)]
    data.table::setkey(tx_info, "tx_name")

    message("Updating transcript models ...")
    # Get indicies of masked bins per transcript group
    tx_remodel <- mapply(function(models, strand, masks) {
      for (m in seq_len(ncol(models))) {
        non_z <- which(models[, m] != 0)
        if (length(non_z) != 0) {
          # Get transcript name
          tx_name <- colnames(models)[m]
          # Compute expected properties of transcript
          expected_range <- unlist(tx_info[tx_name, list(
            expected_bin_start, expected_bin_end)], use.names = FALSE)
          expected_len <- expected_range[2] - expected_range[1] + 1
          # Compute which bins to omit by finding which indicies inside the
          # expected range of the transcript are zero
          omit <- -which(models[expected_range[1]:expected_range[2], m] == 0)

          # Compute pre-computed value closest to input value
          in_val <- rescale_fixed_width_position(
            bin_size = tsp@bin_size,
            length_out = expected_len,
            linear_head_length = tsp@linear_head_length,
            linear_tail_length = tsp@linear_tail_length,
            head_percent = tsp@head_percent,
            tail_percent = tsp@tail_percent
          )
          # Convert x values to lookup bins in the cached loess calculations
          lookup_ind <- round(in_val * (grid_length - 1)) + 1
          # Remove in_val indicies that have been removed from the model, either
          # via rounding or masking by comparing with the observed vs. expected
          # start and end bins
          if (strand == "+") {
            if (length(omit) != 0) {
              models[non_z, m] <- pre_profile[lookup_ind][omit]
            } else {
              models[non_z, m] <- pre_profile[lookup_ind]
            }
          } else if (strand == "-") {
            if (length(omit) != 0) {
              models[non_z, m] <- pre_profile[rev(lookup_ind)][omit]
            } else {
              models[non_z, m] <- pre_profile[rev(lookup_ind)]
            }
          }
        }
      }
      return(models)
    }, models = tq@models, strand = grp_strand, masks = tq@masks,
    SIMPLIFY = FALSE)
    tq@models <- tx_remodel
    return(tq)
  })

#' @title Slims down loess object
#'
#' Removes parts of \code{loess} object that are not required for prediction
#' or object printing. Edited down from code in the "strip" package.
#'
#' @param object a loess object
#'
#' @name strip_loess
strip_loess <- function(object) {
  if (class(object) != "loess") {
    stop("must be loess object")
  }
  # Remove bits of object not needed for prediction
  object$y <- NULL
  object$model <- NULL
  object$fitted <- NULL
  object$residuals <- NULL
  object$weights <- NULL
  object$data <- NULL
  object$one.delta <- NULL # nolint
  object$two.delta <- NULL # nolint
  attr(object$terms,".Environment") <- NULL # nolint
  return(object)
}
