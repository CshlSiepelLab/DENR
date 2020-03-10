# Load in test txdb
txdb_path_ss <- system.file("extdata", "test_single_strand.txdb",
                         package = "tuSelecter2")
txdb_ss <- AnnotationDbi::loadDb(file = txdb_path_ss)
txdb_path_ds <- system.file("extdata", "test_double_strand.txdb",
                         package = "tuSelecter2")
txdb_ds <- AnnotationDbi::loadDb(file = txdb_path_ds)

## Convert txdb to granges
gr_ss <- GenomicFeatures::transcripts(txdb_ss, c("tx_name", "gene_id"))
gr_ds <- GenomicFeatures::transcripts(txdb_ds, c("tx_name", "gene_id"))
gr_ss$gene_id <- as.character(gr_ss$gene_id)
gr_ss$tx_name <- as.character(gr_ss$tx_name)

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
  # Additional test case for granges where there are three transcripts
  # where two non-overlapping transcripts are nested within a third
  gr_3 <- GenomicRanges::GRanges(1,
                                 IRanges::IRanges(
                                   c(1, 1, 500),
                                   c(1000, 200, 700)
                                 ))
  expect_equal(length(group_transcripts(gr_3)), 1)
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
               ceiling((5999 - 1000 + 1) / 50))
  expect_equal(length(create_bins(tx_grp[1], bin_size = 17)[[1]]),
               ceiling((5999 - 1000 + 1) / 17))
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
    grp1 = GenomicRanges::GRanges(c(1, 1),
                           IRanges::IRanges(c(1, 26), c(100, 100)),
                           tx_name = c("t1", "t2")),
    grp2 = GenomicRanges::GRanges(c(2, 2),
                           IRanges::IRanges(c(301, 351), c(500, 475)),
                           tx_name = c("t3", "t4"))
  )
  # Create bins at different sizes
  tx_bins_25 <- create_bins(test_tx, bin_size = 25)
  tx_bins_50 <- create_bins(test_tx, bin_size = 50)
  # Create models
  tx_models_25 <- create_transcript_models(test_tx, tx_bins_25, "tx_name")
  tx_models_25_short <- create_transcript_models(test_tx[1], tx_bins_25[1],
                                                 "tx_name")
  tx_models_50 <- create_transcript_models(test_tx, tx_bins_50, "tx_name")
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
  expect_error(create_transcript_models(test_tx[[1]], tx_bins_25, "tx_name"))
  expect_error(create_transcript_models(test_tx, tx_bins_25[[1]], "tx_name"))
  expect_error(create_transcript_models(test_tx, tx_bins_25[1], "tx_name"))
})

test_that("Transcript mask construction", {
  # Some simple test cases
  test_tx <- GenomicRanges::GRangesList(
    grp1 = GenomicRanges::GRanges(c(2, 2),
                                  IRanges::IRanges(c(301, 356), c(500, 475)),
                                  tx_name = c("t1", "t2"))
  )
  # Create bins at different sizes
  tx_bins_25 <- create_bins(test_tx, bin_size = 25)
  # Create models
  tx_models_25 <- create_transcript_models(test_tx,
                                           tx_bins_25, "tx_name")
  # Create masks
  mask_0_0_plus <- create_model_masks(transcript_models = tx_models_25,
                                      strand = "+")
  mask_1_0_plus <- create_model_masks(transcript_models = tx_models_25,
                                      strand = "+", c(1, 0), c(0, 0))
  mask_0_1_plus <- create_model_masks(transcript_models = tx_models_25,
                                      strand = "+", c(0, 0), c(1, 0))
  mask_1_1_plus <- create_model_masks(transcript_models = tx_models_25,
                                      strand = "+", c(1, 0), c(1, 0))
  mask_1_2_3_4_plus <- create_model_masks(transcript_models = tx_models_25,
                                          strand = "+", c(1, 2), c(3, 4))
  mask_0_0_minus <- create_model_masks(transcript_models = tx_models_25,
                                      strand = "-")
  mask_1_0_minus <- create_model_masks(transcript_models = tx_models_25,
                                      strand = "-", c(1, 0), c(0, 0))
  mask_0_1_minus <- create_model_masks(transcript_models = tx_models_25,
                                      strand = "-", c(0, 0), c(1, 0))
  mask_1_1_minus <- create_model_masks(transcript_models = tx_models_25,
                                      strand = "-", c(1, 0), c(1, 0))
  mask_1_2_3_4_minus <- create_model_masks(transcript_models = tx_models_25,
                                          strand = "-", c(1, 2), c(3, 4))

  # Test that masks are correct
  expect_equivalent(mask_0_0_plus, list(integer(0)))
  expect_equivalent(mask_1_0_plus, list(c(1, 3)))
  expect_equivalent(mask_0_1_plus, list(c(7, 8)))
  expect_equivalent(mask_1_1_plus, list(c(1, 3, 7, 8)))
  expect_equivalent(mask_1_2_3_4_plus, list(c(1:3, 5:8)))
  expect_equivalent(mask_0_0_minus, list(integer(0)))
  expect_equivalent(mask_1_0_minus, list(c(7, 8)))
  expect_equivalent(mask_0_1_minus, list(c(1, 3)))
  expect_equivalent(mask_1_1_minus, list(c(1, 3, 7, 8)))
  expect_equivalent(mask_1_2_3_4_minus, list(c(1:5, 7, 8)))

  # Test errors
  expect_error(create_model_masks(tx_models_25, "+", c(-1, 0), c(0, 0)),
               "the length of input vector")
  expect_error(create_model_masks(tx_models_25[[1]], "+"),
               "transcript models must be a list of matricies")
  expect_error(create_model_masks(tx_models_25, c("+", "+")),
               "strand and transcript_models must be the same length")
  expect_error(create_model_masks(tx_models_25, c("*")),
               "strand must be a vector containing only '\\+','\\-'")
})

test_that("Transcript masking", {
  # Some simple test cases
  test_tx <- GenomicRanges::GRangesList(
    grp1 = GenomicRanges::GRanges(c(2, 2),
                                  IRanges::IRanges(c(301, 356), c(500, 475)),
                           tx_name = c("t1", "t2"))
  )
  # Create bins at different sizes
  tx_bins_25 <- create_bins(test_tx, bin_size = 25)
  # Create models
  tx_models_25 <- create_transcript_models(test_tx, tx_bins_25, "tx_name")
  # Create masks for plus strand
  mask_1_0_plus <- create_model_masks(transcript_models = tx_models_25,
                                      strand = "+", c(1, 0), c(0, 0))
  masked_1_0_p <- mask_transcripts(transcript_models = tx_models_25,
                                   masks = mask_1_0_plus)
  # The true masked model for plus strand
  true_masked_model_plus <- list(
    t(matrix(c(0, 1, 0, 1, 1, 1, 1, 1,
               0, 0, 0, 1, 1, 1, 1, 0),
             byrow = T, ncol = 8))
  )
  # Create masks for minus strand
  mask_1_2_3_4_minus <- create_model_masks(transcript_models = tx_models_25,
                                           strand = "-", c(1, 2), c(3, 4))
  masked_1_2_3_4_m <- mask_transcripts(transcript_models = tx_models_25,
                                   masks = mask_1_2_3_4_minus)
  # The true masked model for minus strand
  true_masked_model_minus <- list(
      t(matrix(c(0, 0, 0, 0, 0, 1, 0, 0,
                 0, 0, 0, 0, 0, 1, 0, 0),
               byrow = T, ncol = 8))
  )
  # Test for correct output
  expect_equivalent(masked_1_0_p, true_masked_model_plus)
  expect_equivalent(masked_1_2_3_4_m, true_masked_model_minus)

  # Test errors
  expect_error(mask_transcripts(rep(tx_models_25, 2), mask_1_0_plus),
               "masks and transcript_models must be the same length")
  expect_error(mask_transcripts(tx_models_25[[1]], mask_1_0_plus),
               "transcript models must be a list of matricies")
  })

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

test_that("transcript_quantifier object construction", {
  # Check that models can be constructed correctly
  expect_s4_class(transcript_quantifier(transcripts = gr_ss,
                                        transcript_name_column = "tx_name",
                                        gene_name_column = "gene_id"),
                  "transcript_quantifier")
  # check for errors
  gr_ss_err <- gr_ss
  gr_ss_err$tx_name <- as.list(gr_ss_err$tx_name)
  gr_ss_err2 <- gr_ss
  gr_ss_err2$gene_id <- as.list(gr_ss_err$gene_id)
  expect_error(transcript_quantifier(transcripts = gr_ss_err,
                                        transcript_name_column = "tx_name",
                                        gene_name_column = "gene_id"),
                  "must be of class 'character'")
  expect_error(transcript_quantifier(transcripts = gr_ss_err2,
                                     transcript_name_column = "tx_name",
                                     gene_name_column = "gene_id"),
               "must be of class 'character'")
  expect_error(transcript_quantifier(transcripts = gr_ss,
                                transcript_name_column = "tx_name",
                                bin_size = -1))
  expect_error(transcript_quantifier(transcripts = gr_ss,
                                transcript_name_column = "tx_name",
                                distance = -1))
  expect_error(transcript_quantifier(transcripts = gr_ss,
                                transcript_name_column = "foo"),
               "transcripts does not have a column matching foo")
  expect_error(transcript_quantifier(transcripts = gr_ss,
                                     transcript_name_column = "tx_name",
                                     gene_name_column = "foo"),
               "transcripts does not have a column matching foo")
  # Try some illegal object modifications
  tq <- transcript_quantifier(transcripts = gr_ss,
                   transcript_name_column = "tx_name")
})
