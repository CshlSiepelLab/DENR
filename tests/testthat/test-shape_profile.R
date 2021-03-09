# Load in test txdb
txdb_path <- system.file("extdata", "test_transcript_filter.txdb",
                            package = "DENR")
txdb <- AnnotationDbi::loadDb(file = txdb_path)
gr <- GenomicFeatures::transcripts(txdb, c("tx_name", "gene_id"))

# Load in test bigwigs
bw_200 <- system.file("extdata", "test_transcript_filter_301_500.bw",
                         package = "DENR")
bw_1000 <- system.file("extdata", "test_transcript_filter_1_1000.bw",
                      package = "DENR")

# test inputs
test_that("covered_bases returns correct values", {
  test1 <- GenomicRanges::GRanges("1", IRanges::IRanges(c(1, 7, 7), c(10, 9, 18)),
                   strand = c("+", "+", "-"))
  test2 <- GenomicRanges::GRanges("1", IRanges::IRanges(5, 15), strand = "+")
  test3 <- GenomicRanges::GRanges("1", IRanges::IRanges(5, 15), strand = "-")
  test4 <- GenomicRanges::GRanges(c(1, 1, 2),
                                  IRanges::IRanges(c(1, 1, 1),
                                                   c(10, 10, 10)),
                                  strand = c("+", "-", "+"))
  test5 <- GenomicRanges::GRanges(c(1, 1, 2, 2),
                                  IRanges::IRanges(c(5, 8, 1, 100),
                                                   c(10, 10, 7, 200)),
                                  strand = c("+", "-", "+", "+"))
  expect_equal(covered_bases(test1, test2), c(6, 3, 0))
  expect_equal(covered_bases(test1, test3), c(0, 0, 9))
  expect_equal(covered_bases(test4, test5), c(6, 3, 7))
})

test_that("model_agreement returns correct values", {
  tx_gr <- group_transcripts(gr)
  m_200 <- model_agreement(tx_gr, bw_200, bw_200)
  m_1000 <- model_agreement(tx_gr, bw_1000, bw_1000)
  expect_equal(m_200$percent_transcribed, c(0.2, 1))
  expect_equal(m_200$percent_match, c(0.2, 1))
  expect_equal(m_1000$percent_transcribed, c(1, 1))
  expect_equal(m_1000$percent_match, c(1, 0.2))
})

test_that("rescale_fixed_width_position returns correct values", {
  ref_1 <- c(0, 0.2, seq(0.2, 0.8, length.out = 7)[c(-1, -7)], 0.8, 0.9, 0.1)
  test_1 <- rescale_fixed_width_position(bin_size = 1, length_out = 10,
                                         linear_head_length = 2,
                                         linear_tail_length = 3,
                                         head_percent = 0.2,
                                         tail_percent = 0.2,
                                         short_heuristic = TRUE)
  ref_2 <- c(0, 0.2, 0.8, 1)
  test_2 <- rescale_fixed_width_position(bin_size = 1, length_out = 4,
                                         linear_head_length = 2,
                                         linear_tail_length = 3,
                                         head_percent = 0.2,
                                         tail_percent = 0.2,
                                         short_heuristic = TRUE)
  expect_equal(test_2, ref_2)
  expect_error({
    rescale_fixed_width_position(bin_size = 1, length_out = 4,
                                 linear_head_length = 2,
                                 linear_tail_length = 3,
                                 head_percent = 0.2,
                                 tail_percent = 0.2,
                                 short_heuristic = FALSE)
  })
})

test_that("construction of transcript_shape_profile object and plotting", {
  tsp <- suppressMessages(
    transcript_shape_profile(transcripts = gr,
                                  bigwig_plus = bw_1000, bigwig_minus = bw_1000,
                                  bin_size = 10, min_transcript_length = 100,
                                  linear_head_length = 50,
                                  linear_tail_length = 50,
                                  span = 0.4)
  )
  expect_s4_class(tsp, "transcript_shape_profile")
  expect_equal(GenomicRanges::width(tsp@transcripts), 1e3)
  # check for error when no transcripts remaining
  expect_error(suppressMessages(transcript_shape_profile(transcripts = gr,
                             bigwig_plus = bw_1000, bigwig_minus = bw_1000,
                             bin_size = 10, min_transcript_length = 1e4,
                             linear_head_length = 50,
                             linear_tail_length = 50,
                             span = 0.4)))
  # check for error when bin too large
  expect_error(suppressMessages(
    transcript_shape_profile(transcripts = gr,
                             bigwig_plus = bw_1000, bigwig_minus = bw_1000,
                             bin_size = 100, min_transcript_length = 100,
                             linear_head_length = 50,
                             linear_tail_length = 50,
                             span = 0.4)
  ))
  # check plotting
  expect_s3_class(view(tsp), "ggplot")
  expect_s3_class(view(tsp, gene_length = 200), "ggplot")
})

test_that("transcript_shape_profile reshapes", {
  tsp <- suppressMessages(transcript_shape_profile(transcripts = gr,
                             bigwig_plus = bw_1000, bigwig_minus = bw_1000,
                             bin_size = 10, min_transcript_length = 100,
                             linear_head_length = 50,
                             linear_tail_length = 50,
                             span = 0.4))
  # Modify tsp manually for testing purposes
  tsp@shape_scale_factor <- 0.5
  # Create tq object
  tq <- transcript_quantifier(gr, "tx_name", bin_size = 10)
  # Modify tq object
  tq_shape <- suppressMessages(apply_shape_profile(tq, tsp))
  # Test for 1/2 model equality
  expect_equivalent(tq@models[[1]] * 0.5, tq_shape@models[[1]])
  # Same test as before but with masking
  tq <- transcript_quantifier(gr, "tx_name", bin_size = 10,
                              mask_start_bins = c(5, 5))
  tq_shape <- suppressMessages(apply_shape_profile(tq, tsp))
  # Test for 1/2 model equality
  expect_equivalent(tq@models[[1]] * 0.5, tq_shape@models[[1]])
  # Same test as before but on minus strand
  GenomicRanges::strand(gr) <- "-"
  tq <- transcript_quantifier(gr, "tx_name", bin_size = 10,
                              mask_start_bins = c(5, 5))
  tq_shape <- suppressMessages(apply_shape_profile(tq, tsp))
  # Test for 1/2 model equality
  expect_equivalent(tq@models[[1]] * 0.5, tq_shape@models[[1]])
})
