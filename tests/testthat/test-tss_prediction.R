# Load data
suppressMessages({
  txdb_path_ds <- system.file("extdata", "test_double_strand.txdb",
                              package = "tuSelecter2")
  txdb_ds <- AnnotationDbi::loadDb(file = txdb_path_ds)
  gr_ds <- GenomicFeatures::transcripts(txdb_ds, c("tx_name", "gene_id"))
  gr_ds$gene_id <- unlist(gr_ds$gene_id)
})

# Create transcript_quantifier object
tq <- suppressMessages(transcript_quantifier(gr_ds, bin_size = 50,
                                             transcript_name_column = "tx_name",
                                             mask_start_bins = c(5, 5)))

# The paths to relevant bigwig files
bwp <- system.file("extdata", "test_double_strand_plus.bw",
                   package = "tuSelecter2")
bwm <- system.file("extdata", "test_double_strand_minus.bw",
                   package = "tuSelecter2")
bw_multi <- system.file("extdata/test_multichrom.bw", package = "tuSelecter2")

test_that("windowing regions", {
  tmp <- GenomicRanges::GRanges(seqnames = c(1, "Z"),
                                IRanges::IRanges(start = c(500, 750), width = 1))
  # Create test windows
  test_win_1 <- window_regions(tmp, bins = 5, bin_width = 11)
  test_win_2 <- window_regions(tmp, bins = 5, bin_width = c(11, 13))
  test_win_3 <- window_regions(tmp, bins = c(5, 11), bin_width = c(11, 13))
  # Check dimensionality of outputs
  expect_length(test_win_1, 2)
  expect_length(test_win_2, 2)
  expect_length(test_win_3, 2)
  expect_length(test_win_1[[1]], 5)
  expect_length(test_win_2[[2]], 10)
  expect_length(test_win_3[[1]], 16)
  # check window width
  expect_equal(GenomicRanges::width(test_win_1[[1]]), rep(11, 5))
  expect_equal(GenomicRanges::width(test_win_3[[1]]), rep(c(11, 13), c(5, 11)))
  # Check that center of windows is correct
  core <- BiocGenerics::unlist(GenomicRanges::resize(GenomicRanges::reduce(test_win_1),
                                             width = 1,
                                             fix = "center"))
  expect_equal(GenomicRanges::start(core), GenomicRanges::start(tmp))
  # check for errors
  expect_error(window_regions(tmp, bins = c(11, 13)))
  expect_error(window_regions(tmp, bins = 10), "bins must be odd-valued")
  expect_error(window_regions(tmp, bin_width = 10), "all bin_width values must odd")
})


test_that("feature collection", {
  tmp <- GenomicRanges::GRanges(seqnames = c(1, "Z"),
                                IRanges::IRanges(start = c(500, 750), width = 1),
                                strand = "+")
  feat <- suppressMessages(collect_tss_features(tmp, bigwig_plus = bw_multi,
                                                bigwig_minus = bw_multi,
                                                bins = 5, bin_width = 11))
  expect_equal(dim(feat$feature_array), c(2, 5, 2))
  # test feature where the first element goes out of range
  feat_oor <- suppressMessages(collect_tss_features(tmp, bigwig_plus = bw_multi,
                                                    bigwig_minus = bw_multi,
                                                    bins = 21, bin_width = 51))
  expect_equal(dim(feat_oor$feature_array), c(1, 21, 2))
  expect_equivalent(tmp[1], feat_oor$tss)
  # Test for error with incomplete vector
  expect_error(suppressMessages(collect_tss_features(tmp, bigwig_plus = bw_multi,
                                                    bigwig_minus = bw_multi,
                                                    bins = 21, bin_width = 51,
                                                    remove_incomplete = FALSE)))
  # No features remaining
  expect_error(suppressMessages(collect_tss_features(tmp, bigwig_plus = bw_multi,
                                                     bigwig_minus = bw_multi,
                                                     bins = 21, bin_width = 501)),
               "No sites remaining")
})

test_that("labeled feature collection", {
  tmp <- GenomicRanges::GRanges(seqnames = c(1, "Z"),
                                IRanges::IRanges(start = c(500, 750), width = 1),
                                strand = "+",
                                active_tss = c(TRUE, FALSE))
  tmp_fail <- GenomicRanges::GRanges(seqnames = c(1, "Z"),
                                IRanges::IRanges(start = c(500, 750), width = 1),
                                strand = "+")
  out <- suppressMessages(labeled_feature_set(tmp, bigwig_plus = bw_multi,
                                        bigwig_minus = bw_multi,
                                        bins = 11, bin_width = 7))
  expect_length(out, 2)
  expect_error(suppressMessages(labeled_feature_set(tmp_fail,
                                                    bigwig_plus = bw_multi,
                                                    bigwig_minus = bw_multi,
                                                    bins = 11, bin_width = 7)))
  expect_error(suppressMessages(labeled_feature_set(tmp[1], bigwig_plus = bw_multi,
                                                    bigwig_minus = bw_multi,
                                                    bins = 11, bin_width = 7)))
})


test_that("train and predict inactive transcripts", {
  skip_if_not(tensorflow::tf_config()$available)
  # Ensure that network tests with CPU, not GPU
  if (Sys.info()["sysname"] == "Windows") {
    Sys.setenv(CUDA_VISIBLE_DEVICES = "-1")
  } else {
    Sys.setenv(CUDA_DEVICE_ORDER = "PCI_BUS_ID")
  }
  tx <- tq@transcripts
  tx$active_tss <- c(T, F, F, T)
  feat <- suppressMessages(collect_tss_features(tx, bwp, bwm, bins = 25, bin_width = 1))
  mod <- tss_predictor(train_features = feat$feature_array,
                      train_labels = as.integer(feat$tss$active_tss), FALSE)
  expect_s3_class(mod, "keras.engine.sequential.Sequential")
  # Keep tensorflow quiet
  Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3")
  sink(file = nullfile())
  pred <- predict_inactive_transcripts(tq, bwp, bwm)
  sink()
  # check for obviously inactive tss, the model really should get these right
  expect_equal(sort(pred), c("t1.2", "t2.1"))
})
