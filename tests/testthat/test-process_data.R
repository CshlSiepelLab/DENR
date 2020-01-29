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
  expect_error(assert_chromosome_exists(c(1, 2), bigwig_file = bw_ss))
})

test_that("bigwig summaries have the correct type/dimensions", {
  expect_type(summarize_bigwig(bigwig_file = bw_ss, bins = tx_bins_25,
                   summary_operation = "sum"), "list")
  expect_type(summarize_bigwig(bigwig_file = bw_ss, bins = tx_bins_25,
                   summary_operation = "mean"), "list")
  expect_type(summarize_bigwig(bigwig_file = bw_ss, bins = tx_bins_25[[1]],
                   summary_operation = "median"), "double")
  expect_type(summarize_bigwig(bigwig_file = bw_ss, bins = tx_bins_25[[1]],
                               summary_operation = "min"), "double")
  expect_length(summarize_bigwig(bigwig_file = bw_ss, bins = tx_bins_25,
                               summary_operation = "sum"), 2)
  expect_length(summarize_bigwig(bigwig_file = bw_ss, bins = tx_bins_25[[1]],
                               summary_operation = "median"),
                length(tx_bins_25[[1]]))
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

# A test for adding data to the transcript_quantifier object
bw_plus <- system.file("extdata",
                       "test_double_strand_plus.bw", package = "tuSelecter2")
bw_minus <- system.file("extdata",
                        "test_double_strand_minus.bw", package = "tuSelecter2")

txdb_path_ds <- system.file("extdata", "test_double_strand.txdb",
                            package = "tuSelecter2")
txdb_ds <- AnnotationDbi::loadDb(file = txdb_path_ds)
gr_ds <-
    GenomicFeatures::transcripts(txdb_ds, c("tx_name", "gene_id"))

tq <- transcript_quantifier(gr_ds,
                            bin_size = 50,
                            transcript_name_column = "tx_name")

test_that("counts are added correctly and in the same order as bins", {
    tq_added_data <-
        add_data(tq, bigwig_plus = bw_plus, bigwig_minus = bw_minus)
    expect_equal(names(tq_added_data@counts), names(tq@bins))
})
