context("Test transcriptModel construction")

# Load in test txdb
txdb_path_ss <- system.file("extdata", "test_single_strand.txdb",
                         package = "tuSelecter2")
txdb_ss <- AnnotationDbi::loadDb(file = txdb_path_ss)
txdb_path_ds <- system.file("extdata", "test_double_strand.txdb",
                         package = "tuSelecter2")
txdb_ds <- AnnotationDbi::loadDb(file = txdb_path_ds)

## Convert txdb to granges
gr_ss <- GenomicFeatures::transcripts(txdb_ss)
gr_ds <- GenomicFeatures::transcripts(txdb_ds)

test_that("Transcripts group correctly (single strand)", {
  tx_grp <- group_transcripts(gr_ss)
  # check that two groups are produced
  expect_equal(length(tx_grp), 2)
  # Check that there are the correct number of transcripts per group
  expect_equal(unlist(lapply(tx_grp, length), use.names = FALSE), c(4, 2))
  # Check that every member of the group overlaps at least one other member
  expect_equal({
      over_df <- as.data.frame(GenomicRanges::findOverlaps(tx_grp[[1]]))
      over_df <- over_df[over_df$queryHits != over_df$subjectHits, ]
      length(unique(over_df$queryHits))
    }, 4)
  # Create new groups with a non-zero distance option
  tx_grp_expand_1 <- group_transcripts(gr_ss, distance = 1000)
  tx_grp_expand_2 <- group_transcripts(gr_ss, distance = 1001)
  # Check that the correct number of groups are formed
  expect_equal(length(tx_grp_expand_1), 2)
  expect_equal(length(tx_grp_expand_2), 1)
})

test_that("Transcripts group correctly (double strand)", {
  gr_ds_mod <- gr_ds
  GenomicRanges::end(gr_ds_mod)[4] <- 5300
  tx_grp <- group_transcripts(gr_ds_mod)
  # check that two groups are produced
  expect_equal(length(tx_grp), 2)
  # Check that there are the correct number of transcripts per group
  expect_equal(unlist(lapply(tx_grp, length), use.names = FALSE), c(2, 2))
})

test_that("Transcripts are binned correctly", {
  tx_grp <- group_transcripts(gr_ss)
  # Check correct output length
  expect_equal(length(create_bins(tx_grp)), 2)
  # Check for the correct number of bins with various bin sizes
  expect_equal(length(create_bins(tx_grp[1], bin_size = 50)[[1]]),
               ceiling((5999 - 1000) / 50))
  expect_equal(length(create_bins(tx_grp[1], bin_size = 17)[[1]]),
               ceiling((5999 - 1000) / 17))
  # Check that the bins are the correct widths
  expect_true(
    all(S4Vectors::width(create_bins(tx_grp, bin_size = 17)[[1]]) == 17))
  # Check error catching
  expect_error(create_bins(tx_grp[[1]]),
               "transcript_groups is not a GRangesList object")
  expect_error(create_bins(tx_grp, bin_size = -1),
    "bin_size is not a positive number")
})

test_that("Transcript model generation", {
  # Some simple test cases
  test_tx <- GenomicRanges::GRangesList(
    GenomicRanges::GRanges(c(1, 1), IRanges::IRanges(c(1, 26), c(100, 100))),
    GenomicRanges::GRanges(c(2, 2), IRanges::IRanges(c(301, 351), c(500, 475)))
  )
  # Create bins at different sizes
  tx_bins_25 <- create_bins(test_tx, bin_size = 25)
  tx_bins_50 <- create_bins(test_tx, bin_size = 50)
  # Create models
  tx_models_25 <- create_transcript_models(test_tx, tx_bins_25)
  tx_models_25_short <- create_transcript_models(test_tx[1], tx_bins_25[1])
  tx_models_50 <- create_transcript_models(test_tx, tx_bins_50)
  # The true models
  true_models_25 <- list(
    t(matrix(c(1, 1, 1, 1,
             0, 1, 1, 1),
           byrow = T, ncol = 4)),
    t(matrix(c(1, 1, 1, 1, 1, 1, 1, 1,
             0, 0, 1, 1, 1, 1, 1, 0),
           byrow = T, ncol = 8))
  )
  true_models_50 <- list(
    t(matrix(c(1, 1,
             0.5, 1),
           byrow = T, ncol = 2)),
    t(matrix(c(1, 1, 1, 1,
             0, 1, 1, 0.5),
           byrow = T, ncol = 4))
  )
  # Test that the output is correct
  tmp25 <- mapply(expect_equivalent,
                 object = tx_models_25, expected = true_models_25)
  tmp25_short <- mapply(expect_equivalent,
                  object = tx_models_25_short, expected = true_models_25[1])
  tmp50 <- mapply(expect_equivalent,
                  object = tx_models_50, expected = true_models_50)
  # Test error catching
  expect_error(create_transcript_models(test_tx[[1]], tx_bins_25))
  expect_error(create_transcript_models(test_tx, tx_bins_25[[1]]))
  expect_error(create_transcript_models(test_tx, tx_bins_25[1]))
})
