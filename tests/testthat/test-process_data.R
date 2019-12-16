# Load in test txdb
txdb_path_ss <- system.file("extdata", "test_single_strand.txdb",
                            package = "tuSelecter2")
txdb_ss <- AnnotationDbi::loadDb(file = txdb_path_ss)

## Convert txdb to granges
gr_ss <- GenomicFeatures::transcripts(txdb_ss)

# Import bigwigs
bw_ss <- system.file("extdata/test_single_strand.bw", package = "tuSelecter2")

# Create bins
tx_bins_25 <- create_bins(transcript_groups = group_transcripts(gr_ss),
                          bin_size = 25)

test_that("assert_chromosome_exists", {
  expect_silent(assert_chromosome_exists(1, bigwig_file = bw_ss))
  expect_silent(assert_chromosome_exists("1", bigwig_file = bw_ss))
  expect_error(assert_chromosome_exists(2, bigwig_file = bw_ss))
})

test_that("bigwig summaries have the correct type/dimensions", {
  expect_type(summarize_bigwig(bigwig_file = bw_ss, bins = tx_bins_25,
                   summary_operation = "sum"), "list")
  expect_type(summarize_bigwig(bigwig_file = bw_ss, bins = tx_bins_25,
                   summary_operation = "mean"), "list")
  expect_type(summarize_bigwig(bigwig_file = bw_ss, bins = tx_bins_25[[1]],
                   summary_operation = "median"), "list")
  expect_type(summarize_bigwig(bigwig_file = bw_ss, bins = tx_bins_25[[1]],
                               summary_operation = "min"), "list")
  expect_length(summarize_bigwig(bigwig_file = bw_ss, bins = tx_bins_25,
                               summary_operation = "sum"), 2)
  expect_length(summarize_bigwig(bigwig_file = bw_ss, bins = tx_bins_25[[1]],
                               summary_operation = "median"), 1)
})

test_that("summarize_bigwig parallel implementations", {
  expect_type(summarize_bigwig(bigwig_file = bw_ss, bins = tx_bins_25,
                               summary_operation = "sum", threads = 2), "list")
})

test_that("summarize_bigwig error catching", {
  expect_error(summarize_bigwig(bigwig_file = "foo.bw", bins = tx_bins_25,
                               summary_operation = "sum"),
               "bigwig_file path does not exist")
  expect_error(summarize_bigwig(bigwig_file = bw_ss, bins = tx_bins_25,
                                 summary_operation = "foo"))
  # Create GRanges with multiple chromosomes
  tx_seqnames_mix <- tx_bins_25[[1]]
  GenomeInfoDb::seqlevels(tx_seqnames_mix) <- c("1", "2")
  GenomicRanges::seqnames(tx_seqnames_mix) <- rep(c(1, 2), each = 100)
  expect_error(summarize_bigwig(bigwig_file = bw_ss, bins = tx_seqnames_mix),
               "One or more groups of contained multiple chromosomes")
})

# A simple two chromosome (1, Z) bigwig where [1:1000] = 1 and [1:500] = 1
# respectively
multichrom_bw <- system.file("extdata", "test_multichrom.bw",
                             package = "tuSelecter2")
multichrom_gr <- GenomicRanges::GRanges(c(1, "Z"),
                                        IRanges::IRanges(1, 1e3))
multichrom_gr_ls <- GenomicRanges::split(multichrom_gr,
                                      GenomicRanges::seqnames(multichrom_gr))

test_that("multichrom bigwig sum correct values", {
  true_values <- as.list(c(1e3, 500))
  expect_equivalent(summarize_bigwig(bigwig_file = multichrom_bw,
                                     bins = multichrom_gr_ls,
                                     summary_operation = "sum"),
                    true_values
  )
})
