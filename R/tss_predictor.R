#' @title Generate GRangesList of windows
#'
#' @description  Generates a set of windows for each element in a
#' \code{\link[GenomicRanges]{GRanges-class}} object. Uses the center of each
#' region as the 0th co-ordinate in generating windows.
#'
#' @param sites a \link[GenomicRanges]{GRanges-class} containing the sites of
#' interest. The central site of this range will be the focal site used to
#' generate bins
#' @param bins a single value or vector of length equal to bin_width indicating
#' the number of bins to generate, with the central bin centered on the
#' focal site, must be odd
#' @param bin_width a vector of one or more odd valued integers indicating the
#' width of the bins to be generated. If multiple value are specified, windows are
#' generated for all sizes specified
#'
#' @return A list of vectors with each one corresponding to one set of bins and
#' each element of a vector corresponding to a bin
#'
#' @name window_regions
#'
#' @export
window_regions <- function(sites, bins = 21L, bin_width = 51L) {
  ##** Begin checks **##
  if (any(bins %% 2 == 0)) {
    stop("bins must be odd-valued")
  }

  if (!class(sites) == "GRanges") {
    stop("sites must be a GRanges object")
  }

  if (any(!is.integer(bin_width)) && any(bin_width != as.integer(bin_width))) {
    stop("bin_width must be an integer")
  }

  if (any(bin_width %% 2 == 0)) {
    stop("all bin_width values must odd")
  }

  if (length(bins) > 1 && length(bins) != length(bin_width)) {
    stop("If multiple values are specified for bins there must be an equal ",
         "number of bin_width values specified")
  }

  ##** End checks **##

  # Make sure length of bins matches length of bin_width
  if (length(bins) != length(bin_width)) {
    bins <- rep(bins, length(bin_width))
  }
  # Compute focal sites
  focal_sites <- as.integer(round(
    (GenomicRanges::end(sites) + GenomicRanges::start(sites)) / 2))
  # Preallocate start vector
  bin_starts <- integer(sum(bins * length(sites)))
  # Pre-compute start adjustments and unlist
  start_adjust <- unlist(mapply(function(bins, bin_width) {
    as.integer(seq(from = - (bin_width - 1) / 2 - ((bins - 1) / 2) * bin_width,
                   by = bin_width, length.out = bins))
  }, bins = as.list(bins), bin_width = as.list(bin_width), SIMPLIFY = FALSE))
  # Generate vector of widths per site
  widths <- rep(bin_width, times = bins)
  # Pre-calculate the bin column indicies
  bin_idx <- unlist(lapply(as.list(bins), function(x) seq_len(x)))
  # Iterate over start sites
  idx <- 1
  stepsize <- sum(bins)
  for (i in seq_along(focal_sites)) {
    bin_starts[idx:(idx + stepsize - 1)] <-
      as.integer(focal_sites[i] + start_adjust)
    idx <- idx + stepsize
  }
  # Create GRanges object out of the starts
  gr <- GenomicRanges::GRanges(
    seqnames = rep(GenomicRanges::seqnames(sites), each = stepsize),
    ranges = IRanges::IRanges(start = bin_starts,
                              width = rep(widths, length(focal_sites))),
    strand = rep(GenomicRanges::strand(sites), each = stepsize),
    site = rep(factor(seq_along(focal_sites)), each = stepsize),
    bin_idx = rep(bin_idx, length(focal_sites))
  )
  gr <- gr[GenomicRanges::start(gr) > 0]
  # Split into GRangesList
  return(GenomicRanges::split(gr, gr$site))
}

#' @title Generate GRangesList of windows
#'
#' @description  Generates a set of windows for each transcript in a
#' \code{\link[GenomicRanges]{GRanges-class}} object. Uses the TSS of each
#' region as the 0th co-ordinate in generating windows.
#'
#' @param transcripts a \code{\link[GenomicRanges]{GRanges-class}} object
#' @param bigwig_plus the path to a bigwig for reads on the plus strand
#' @param bigwig_minus the path to a bigwig for reads on the minus strand
#' @param remove_incomplete remove regions where any window features go outside
#' the chromosomal range
#' @param remove_all_zero remove TSS that have only 0 values in their feature
#' vector
#' @inheritParams window_regions
#'
#' @return A list of vectors with each one corresponding to one set of bins and
#' each element of a vector corresponding to a bin
#'
#' @name collect_tss_features
#'
#' @export
collect_tss_features <- function(transcripts, bigwig_plus, bigwig_minus,
                                 bins = 21L, bin_width = 51L,
                                 remove_incomplete = TRUE,
                                 remove_all_zero = FALSE) {
  # Check that the same number of bigwig_plus and bigwig_minus files are
  # supplied
  if (length(bigwig_plus) != length(bigwig_minus)) {
    stop("equal numbers of bigwig_plus and biwig_minus files must be supplied")
  }

  # Check that all transcripts are stranded
  if (!all(S4Vectors::runValue(GenomicRanges::strand(transcripts)) %in% c("+", "-"))) {
    stop("all transcripts must have a +/- strand")
  }

  if (bins < 3) {
    stop("bins must be >= 3")
  }

  # Get TSS of transcripts
  tss <- GenomicRanges::promoters(transcripts, upstream = 0, downstream = 1)
  # Window sites
  windows <- window_regions(tss, bins = bins, bin_width = bin_width)
  names(windows) <- as.character(seq_along(windows))
  # Filter out incomplete feature vectors
  if (remove_incomplete) {
    expected_length <- sum(bins)
    remove_sites <- which(S4Vectors::elementNROWS(windows) != expected_length)
    if (length(remove_sites) > 0) {
      windows <- windows[-remove_sites]
      # Filter down TSS to only the sites that were retained
      tss <- tss[-remove_sites]
    }
    message("Removed ", length(remove_sites),
            " sites with incomplete feature vectors")
  } else {
    stop("Not able to handle ragged input vectors yet")
  }

  if (length(windows) == 0) {
    stop("No sites remaining")
  }

  # Collect sense and antisense counts
  # Reverse bin order such that if windows were indexed ..., -1, 0, 1, ... the
  # sense strand always runs in the positive direction and the anti-sense strand
  # always runs in the negative direction
  transcript_strand <- as.vector(GenomicRanges::strand(tss))
  is_sense_positive <- (transcript_strand == "+")
  sense_counts <- list()
  antisense_counts <- list()
  # Iterate over bigwig pairs
  message("Retrieving data ...")
  for (i in seq_along(bigwig_plus)) {
    sense_counts[[i]] <- c(
      summarize_bigwig(bigwig_plus[i], windows[is_sense_positive], "sum"),
      summarize_bigwig(bigwig_minus[i],
                       S4Vectors::revElements(windows[!is_sense_positive]),
                       "sum"))
    sense_counts[[i]] <-
      sense_counts[[i]][order(as.integer(names(sense_counts[[i]])))]

    antisense_counts[[i]] <- c(
      summarize_bigwig(bigwig_plus[i],
                       S4Vectors::revElements(windows[!is_sense_positive]),
                       "sum"),
      summarize_bigwig(bigwig_minus[i], windows[is_sense_positive], "sum"))
    antisense_counts[[i]] <-
      antisense_counts[[i]][order(as.integer(names(antisense_counts[[i]])))]
  }
  # Bind all sample matricies into a single matrix
  sense_counts <-
    abs(do.call("rbind", unlist(sense_counts, recursive = FALSE)))
  antisense_counts <-
    abs(do.call("rbind", unlist(antisense_counts, recursive = FALSE)))
  # Replicate TSS length to match
  tss <- rep(tss, length(bigwig_plus))

  message("Row scaling features ...")
  # Get row means
  mean_scale <- (matrixStats::rowMeans2(sense_counts) +
                   matrixStats::rowMeans2(antisense_counts)) / 2

  # Get row standard deviations
  sd_scale <- matrixStats::rowSds(cbind(sense_counts, antisense_counts))

  # Remove rows with mean zero as that implies all entries are 0
  if (remove_all_zero) {
    keep <- which(mean_scale > 0)
    if (length(keep) != length(mean_scale)) {
      message("Removing sites with all zero entries in their feature vector: ",
              length(keep), " of ", length(mean_scale), " remaining")
      sense_counts <- sense_counts[keep, ]
      antisense_counts <- antisense_counts[keep, ]
      mean_scale <- mean_scale[keep]
      sd_scale <- sd_scale[keep]
      tss <- tss[keep]
    } else {
      stop("No non-zero entries in feature matrix")
    }
  }

  # Scale feature matricies
  sense_counts <- (sense_counts - mean_scale) / sd_scale
  antisense_counts <- (antisense_counts - mean_scale) / sd_scale
  # Force all values in rows with no variance to be 0
  sense_counts[sd_scale == 0, ] <- 0
  antisense_counts[sd_scale == 0, ] <- 0

  # Place counts into a 3D array
  feat_arr <- array(0, dim = c(nrow(sense_counts), bins, 2))
  feat_arr[, , 1] <- sense_counts
  feat_arr[, , 2] <- antisense_counts

  attr(feat_arr, "scaled:center") <- mean_scale
  attr(feat_arr, "scaled:scale") <- sd_scale
  return(list(feature_array = feat_arr, tss = tss))
}

#' @title Build labeled data set
#'
#' @description Builds and labels feature vectors for both active and inactive
#' TSS.
#'
#' @param labeled_transcripts a \code{\link[GenomicRanges]{GRanges-class}} object which
#' has a logical valued metadata column \code{active_tss} indicating if the TSS of the
#' corresponding transcript is active (TRUE -> active)
#' @inheritParams collect_tss_features
#'
#' @return A list with two elements containing the feature arrays for the active
#' and inactive TSS respectively
#' @name labeled_feature_set
#'
#' @export
labeled_feature_set <- function(labeled_transcripts, bigwig_plus,
                                bigwig_minus, bins = 251L, bin_width = 21L) {

  # Drop inactive seqlevels to avoid spurious warnings
  GenomeInfoDb::seqlevels(labeled_transcripts) <-
    GenomeInfoDb::seqlevelsInUse(labeled_transcripts)

  # Check that active_tss column is present
  if (!"active_tss" %in% colnames(GenomicRanges::mcols(labeled_transcripts))) {
    stop("labeled_transcripts must contain column active_tss")
  }

  if (!is.logical(GenomicRanges::mcols(labeled_transcripts)$active_tss)) {
    stop("active_tss must be column of type logical")
  }

  # Subset transcript to only those not NA
  pre_filter <- length(labeled_transcripts)
  labeled_transcripts <-
    labeled_transcripts[!is.na(labeled_transcripts$active_tss)]
  post_filter <- length(labeled_transcripts)
  message("Removed NA transcripts: ", post_filter, " of ", pre_filter,
          " remaining")

  post_table <- table(labeled_transcripts$active_tss)
  if (post_table[1] == 0) {
    stop("No inactive transcripts remaining")
  } else if (post_table[2] == 0) {
    stop("No active transcripts remaining")
  }

  # Get TSS of transcripts
  tss <- GenomicRanges::promoters(labeled_transcripts,
                                  upstream = 0, downstream = 1)

  # Get features for the postive and negative sets of sites
  active_features <- collect_tss_features(tss[tss$active_tss], bigwig_plus,
                                          bigwig_minus, bins = bins,
                                          bin_width = bin_width)
  inactive_features <- collect_tss_features(tss[!tss$active_tss], bigwig_plus,
                                            bigwig_minus, bins = bins,
                                            bin_width = bin_width)
  return(list(active_features, inactive_features))
}

#' @title Generate training labels for TSS
#'
#' @description Labels active TSS from by generating an empirical distribution
#' of *-cap sequencing count data across all TSS and labeling the upper tail
#' as active and sites with no *-cap data as inactive.
#' @param transcripts transcripts a \code{\link[GenomicRanges]{GRanges-class}}
#' that holds all the transcript coordinates
#' @param cap_bigwig_plus PRO-cap/GRO-cap data from plus strand
#' @param cap_bigwig_minus PRO-cap/GRO-cap data from minus strand
#' @param ctrl_bigwig_plus PRO-seq/GRO-seq data from plus strand
#' @param ctrl_bigwig_minus PRO-seq/GRO-seq data from minus strand
#' @param ctrl_cdf_alpha foo
#' @param radius The radius around the TSS in which to collect counts
#' @return Return transcripts with an added active_tss column indicating if the
#' tss is active
#'
#' @include data_handling.R
#' @name label_transcripts
#'
#' @export
label_transcripts <- function(transcripts,
                              cap_bigwig_plus, cap_bigwig_minus,
                              ctrl_bigwig_plus, ctrl_bigwig_minus,
                              ctrl_cdf_alpha = 0.01, radius = 150) {

  # Drop inactive seqlevels to avoid spurious warnings
  GenomeInfoDb::seqlevels(transcripts) <-
    GenomeInfoDb::seqlevelsInUse(transcripts)

  # Get TSS of transcripts
  tss <- GenomicRanges::promoters(transcripts, upstream = radius,
                                  downstream = radius)
  # Default tss to being inactive
  tss$active_tss <- NA

  # Compute total reads in each bigwig file pair
  cap_reads <- total_coverage(cap_bigwig_plus) + total_coverage(cap_bigwig_minus)
  ctrl_reads <- total_coverage(ctrl_bigwig_plus) + total_coverage(ctrl_bigwig_minus)
  # Convert to units of millions
  cap_rpm <- cap_reads / 1e6
  ctrl_rpm <- ctrl_reads / 1e6

  # Get *-cap counts around promoters
  tss_plus <- GenomicRanges::strand(tss) == "+"
  cap_p <-
    abs(summarize_bigwig(cap_bigwig_plus, bins = tss[tss_plus])) / cap_rpm
  cap_m <-
    abs(summarize_bigwig(cap_bigwig_minus, bins = tss[!tss_plus])) / cap_rpm
  ctrl_p <-
    abs(summarize_bigwig(ctrl_bigwig_plus, bins = tss[tss_plus])) / ctrl_rpm
  ctrl_m <-
    abs(summarize_bigwig(ctrl_bigwig_minus, bins = tss[!tss_plus])) / ctrl_rpm

  # Create an empirical cdf
  ctrl_ecdf <- stats::ecdf(c(ctrl_p, ctrl_m))

  # Label the positive and negative training examples based on the tail
  # probability of the values observed within radius of the tss under the
  # control empirical CDF
  tss[tss_plus][1 - ctrl_ecdf(cap_p) <= ctrl_cdf_alpha]$active_tss <- TRUE
  tss[!tss_plus][1 - ctrl_ecdf(cap_m) <= ctrl_cdf_alpha]$active_tss <- TRUE
  tss[tss_plus][cap_p == 0]$active_tss <- FALSE
  tss[!tss_plus][cap_m == 0]$active_tss <- FALSE

  # Update transcripts with tss info
  transcripts$active_tss <- tss$active_tss

  # Return transcripts GRanges with active_tss column indicating if the tss is
  # active
  return(transcripts)
}

#' @title Train TSS activity predictor
#'
#' @description A function that trains and returns model to predict TSS
#' activity using the same convolutional architecture as the pre-trained
#' predictor that ships with this package. Just here in case a user wants to
#' try training the model on alternative data/labels.
#' @param train_features a 3-dimensional numeric array where the dimensions are
#' [transcript, bin, strand]. (IMPORTANT: bin must be >= 25)
#' @param train_labels a 0/1 encoded vector of labels with 1 per transcript
#' @param train logical. Whether to train model
#'
#' @name tss_predictor
#' @importFrom keras %>%
#'
#' @return a trained keras model
#' @export
tss_predictor <- function(train_features, train_labels, train = TRUE) {
  # Checks
  if (length(dim(train_features)) != 3) {
    stop("train_features does not have 3 dimensions")
  }
  if (length(train_labels) != dim(train_features)[1]) {
    stop("unequal number of transcripts and train_labels")
  }
  # Setup a callback for early stopping
  cb <- keras::callback_early_stopping(monitor = "val_loss", patience = 20,
                                       restore_best_weights = TRUE, verbose = 1)

  # Input dimensionality
  bins <- dim(train_features)[2]
  channels <- dim(train_features)[3]

  if (bins < 25) {
    stop("Features must have at least 25 bins")
  }

  # Initalize bias
  b_0 <- log(sum(train_labels) / (length(train_labels) - sum(train_labels)))
  b_0 <- keras::initializer_constant(value = b_0)

  # Create sequential model
  model <- keras::keras_model_sequential(name = "conv_pool_conv")
  model %>%
    keras::layer_conv_1d(filters = 10, kernel_size = 15,
                  padding = "same", input_shape = c(bins, channels)) %>%
    keras::layer_activation("relu") %>%
    keras::layer_max_pooling_1d(pool_size = 5) %>%
    # another 1-D convolution layer
    keras::layer_conv_1d(filters = 10, kernel_size = 10,
                         padding = "same") %>%
    # Defining a Pooling layer which reduces the dimensions of the
    # features map and reduces the computational complexity of the model
    keras::layer_max_pooling_1d(pool_size = 5) %>%
    #dropout layer to avoid overfitting
    keras::layer_dropout(0.25) %>%
    keras::layer_flatten() %>%
    keras::layer_dense(units = 15, activation = "relu") %>%
    keras::layer_dense(units = 1, activation = "sigmoid", bias_initializer = b_0)
  # compile model
  model %>% keras::compile(
    optimizer = "adam",
    loss = "binary_crossentropy",
    metrics =  "accuracy"
  )
  # fit model
  if (train) {
    model %>% keras::fit(train_features, train_labels, epochs = 300, verbose = 3,
                         validation_split = 0.1, batch_size = 300, callbacks = cb)
  }
  # Return fitted model
  return(model)
}

#' @title Return list of inactive TSS
#'
#' @description Predicts inactive transcripts based on GRO/PRO-seq data using a p
#' re-trained convolutional neural net. For stability reasones, this function enforces
#' the use of the CPU if a TF session has not already been started. Given the small size
#' of the network this shouldn not impact performance significantly.
#' @inheritParams add_data
#'
#' @return a vector of inactive transcript identifiers
#'
#' @name predict_inactive_transcripts
#' @include data_handling.R
#' @export
predict_inactive_transcripts <- function(tq, bigwig_plus, bigwig_minus) {
  # If TF session not already started this will force usage of CPU
  if (Sys.info()["sysname"] == "Windows") {
    Sys.setenv(CUDA_VISIBLE_DEVICES = "-1")
  } else {
    Sys.setenv(CUDA_DEVICE_ORDER = "PCI_BUS_ID")
  }
  Sys.setenv(TF_CPP_MIN_LOG_LEVEL = 3)
  model_id <- "conv_pool_conv_14e9c7d5"
  # The paths to relevant bigwig files
  ml_model <- system.file("extdata", model_id,
                     package = "tuSelecter2")
  # Get input dimensions and use it retrieve correct inputs
  config <- yaml::read_yaml(paste0(ml_model, "_input_info.yaml"))
  # Get features
  message("Collecting features ...")
  features <- collect_tss_features(transcripts = tq@transcripts,
                                   bigwig_plus, bigwig_minus,
                                   bins = config$bins,
                                   bin_width = config$bin_width)
  # load model
  model <- keras::load_model_hdf5(ml_model, compile = TRUE)
  # Predict classes
  y_pred <- keras::predict_classes(model, features$feature_array)
  inactive_tss <- S4Vectors::elementMetadata(
    features$tss[which(y_pred == 0)])[, tq@column_identifiers[1]]
  return(inactive_tss)
}
