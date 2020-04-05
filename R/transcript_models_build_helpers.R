#' Create bins for multiple groups of overlapped transcripts
#'
#' Function for generating bins for mutiple groups of overlapped transcripts.
#'
#' @param transcript_groups A \code{\link[GenomicRanges]{GRangesList-class}}
#' object,each element in the list contains transcripts from one gene or several
#' genes which are overlapped.
#' @param bin_size An integer, used to tile the gene region. Default is 250bp.
#' @return A \code{\link[GenomicRanges]{GRangesList-class}} object, containing
#' the binned region.
#' @export
create_bins <- function(transcript_groups, bin_size = 250) {
  # check input class
  if (!methods::is(transcript_groups, "GRangesList")) {
    stop("transcript_groups is not a GRangesList object")
  }
  if (!methods::is(bin_size, "numeric") || bin_size < 0) {
    stop("bin_size is not a positive number")
  }

  # First reduce transcript groups to single range using large gapwidth to
  # ensure single range
  red_tx_grps <- unlist(GenomicRanges::reduce(transcript_groups, min.gapwidth = 1e9))
  # Pre-compute number of bins in each range
  bin_count <- ceiling(GenomicRanges::width(red_tx_grps) / bin_size)
  # Pre-build rle for each element
  chrom <-
    S4Vectors::Rle(as.vector(GenomicRanges::seqnames(red_tx_grps)), bin_count)
  strnd <-
    S4Vectors::Rle(as.vector(GenomicRanges::strand(red_tx_grps)), bin_count)
  # Create Granges as single GRanges
  tiles <- GenomicRanges::GRanges(
    chrom,
    IRanges::IRanges(start = unlist(
      seqv(from = GenomicRanges::start(red_tx_grps),
           to = GenomicRanges::end(red_tx_grps),
           by = bin_size)),
      width = bin_size),
    strand = strnd)
  # Return split GRangesList
  return(GenomicRanges::split(x = tiles, rep(names(transcript_groups), bin_count))
  )
}

#' Vectorized seq.default
#'
#' Version of the \code{seq} function that takes vectorized arguments
#'
#' @inheritParams base::seq
seqv <- Vectorize(seq.default, vectorize.args = c("from", "to"))

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

#' Create transcript models
#'
#' Function for generating binned transcript models given a set of bins
#'
#' @param transcripts A \code{\link[GenomicRanges]{GRanges-class}}
#' object
#' @param bins A \code{\link[GenomicRanges]{GRangesList-class}} object,
#' each element in the list contains the bins for each corresponding group of
#' transcripts
#' @param bin_size width of bins
#' @param transcript_name_column A string that indicates which column in the
#' GRanges object contain the transcript names
#' @return A list of matrices where each row in the matrix corresponds to a
#' bin and each column is a transcript
#' @export
create_transcript_models <- function(transcripts, bins, bin_size,
                                     transcript_name_column) {
  # check input class
  if (!methods::is(transcripts, "GRanges")) {
    stop("transcripts is not a GRanges object")
  }
  # Compute percent overlap of each transcript per bin and cast to matrix
  # First find all pairwise overlaps
  ovr <- IRanges::findOverlapPairs(transcripts, bins)
  # Compute bin percent overlap for all intersections
  ovr_intersect <- IRanges::pintersect(ovr, drop.nohit.ranges = FALSE)
  ovr_val <- GenomicRanges::width(ovr_intersect) / bin_size

  # Compute the dimensions of model matrix for each loci
  d <- S4Vectors::elementNROWS(ovr_val)
  # Pre-extract transcript names for ease of access
  transcript_names <-
    GenomicRanges::values(transcripts)[, transcript_name_column]
  # Reshape intersection output into matricies of the correct dimension
  tx_matrix_models <- lapply(S4Vectors::split(d, names(d))[names(bins)], function(x) {
      # Retrive rows in transcripts that each intersection pair came from
      lookup_ind <- which(names(ovr_val) == names(x[1]))
      # Place into matrix of correct dimension
      m <- matrix(data = BiocGenerics::unlist(ovr_val[lookup_ind]),
                  nrow = x[1],
                  ncol = length(x))
      # Set column names to correct transcript names
      colnames(m) <- transcript_names[lookup_ind]
      return(m)
    }
  )
  return(tx_matrix_models)
}

#' @title Group transcripts
#'
#' @description
#' Creates groups of transcripts on the same strand based on
#' proximity. The groups are constructed of connected transcripts
#' such that a pair of transcripts are considered connected if
#' they are within a given distance \emph{d} of each other. The
#' group then consists of all transcripts that are connected to
#' at least one other member.
#'
#' @param transcript_granges \code{\link[GenomicRanges]{GRanges-class}}
#' object
#' @param distance the distance within which two transcripts are
#' considered connected
#'
#' @return A \code{\link[GenomicRanges]{GenomicRangesList-class}}
#' object
#'
#' @rdname group_transcripts
#' @export
group_transcripts <- function(transcript_granges, distance = 0) {
  # Check that txdb is an actual transcript database object
  if (!methods::is(transcript_granges, "GRanges")) {
    stop("transcript is not a TxDb object")
  }
  # Ensure everything is sorted, otherwise later assignment assumptions fail
  transcript_granges <- GenomeInfoDb::sortSeqlevels(transcript_granges)
  transcript_granges <- GenomicRanges::sort(transcript_granges)
  # Assign a unique id to each row
  transcript_granges$unique_id <- seq_len(length(transcript_granges))
  # Create expanded granges
  tx_granges_expand <- GenomicRanges::resize(
    transcript_granges,
    width = GenomicRanges::width(transcript_granges) + floor(distance / 2),
    fix = "end")
  tx_granges_expand <- GenomicRanges::resize(
    tx_granges_expand,
    width = GenomicRanges::width(tx_granges_expand) + ceiling(distance / 2),
    fix = "start")
  tx_granges_expand <- GenomicRanges::trim(tx_granges_expand)
  # Reduce granges to get unions of overlapping ranges which will become the
  # transcript groups. The -/+ 1 bit is to assure that perfectly adjoining
  # groups are not merged
  GenomicRanges::end(tx_granges_expand) <- GenomicRanges::end(tx_granges_expand) - 1
  tx_reduce <- GenomicRanges::reduce(tx_granges_expand, ignore.strand = FALSE)
  GenomicRanges::end(tx_reduce) <- GenomicRanges::end(tx_reduce) + 1
  # Get the group assignments
  group_assignment <- GenomicRanges::findOverlaps(query = transcript_granges,
                                                  subject = tx_reduce)
  # Create group vector string
  uid <- paste0(S4Vectors::subjectHits(group_assignment), "_",
                GenomicRanges::strand(transcript_granges))

  # Split the transcripting into their groups
  gr_groups <- GenomicRanges::split(transcript_granges, f = uid)
  return(gr_groups)
}

#' Reduce transcript models
#'
#' Combine identical transcript models and get non-redundent and reduced
#' transcript models
#'
#' @param transcript_models_ls A list of matrices where each row in the matrix
#' corresponds to a bin and each column is a transcript.
#' @param bin_operation Three different modes to deal with decimals in the
#' transript model (due to partial overlap of the first or last exon and bins).
#' Either "ceiling", "floor", or "round" (default: "round").
#'
#' @return  A list of matrices holding the reduced transcript models and a
#' dataframe holding each transcript belongs to which group and reduced model.
#' @export
reduce_transcript_models <-
  function(transcript_models_ls,
           bin_operation = c("round", "floor", "ceiling")) {
    # check transcript class
    tx_model_class <- unique(sapply(transcript_models_ls, class))
    if (length(tx_model_class) != 1 || tx_model_class != "matrix") {
      stop("invalid input for trascript models")
    }
    # check bin_operation
    if (length(bin_operation) > 1) bin_operation <- bin_operation[[1]]
    if (!bin_operation %in% c("ceiling", "floor", "round")) {
      stop("invalid arguement for bin_operation, it should be 'ceiling', 'floor', ",
           "or 'round'.")
    }
    # deal with decimals in transcript models
    integer_models <- lapply(transcript_models_ls, function(x, fun) {
      do.call(fun, list(x))
    }, fun = bin_operation)
    # compute reduced transcript groups
    tx_groups <- lapply(integer_models, group_trancript_models)
    # get reduced tx models
    get_reduced_tx_models <- function(tx_models, tx_groups) {
      tx_models[, unlist(lapply(tx_groups, function(x) x[[1]])), drop = FALSE]
    }
    reduced_tx_models <- mapply(get_reduced_tx_models,
                                integer_models,
                                tx_groups, SIMPLIFY = FALSE)
    # create group and model identifier for each transcript
    group_len <- lapply(tx_groups, lengths)
    tx_num <- sapply(group_len, sum)
    tx_group_model <- data.frame(tx_name = unlist(tx_groups),
                                 group = rep(seq_along(tx_num), tx_num),
                                 model = unlist(mapply(rep, lapply(group_len, seq_along),
                                                       group_len, SIMPLIFY = FALSE)))
    rownames(tx_group_model) <- NULL
    return(list(reduced_tx_models, tx_group_model))
    }

#' Group transcripts
#'
#' Group transripts if they have identical transcript model.
#'
#' @param tx_models A processed matrix without decimals in transcript model
#' @seealso \code{\link{reduce_transcript_models}},
#' @return  A list of grouped transcript names.

group_trancript_models <- function(tx_models) {
  pool <- colnames(tx_models)
  groups <- list()
  while (length(pool) > 1) {
    identical_models <- which(
      colSums(abs(tx_models[, pool[-1], drop = FALSE] - tx_models[, pool[1]])) == 0
    )
    groups <- append(groups, list(c(pool[1], names(identical_models))))
    pool <- pool[-c(1, identical_models + 1)]
  }
  if (length(pool) == 1) groups <- append(groups, pool)
  return(groups)
}
