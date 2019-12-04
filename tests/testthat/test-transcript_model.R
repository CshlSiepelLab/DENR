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
  # Check error catching
  expect_error(create_bins(tx_grp[[1]]), "gr_ls is not a GRangesList object")
  expect_error(create_bins(tx_grp, bin_size = -1),
    "bin_size is not a positive number")
})
