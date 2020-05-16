# Load in test txdb
txdb_path_ss <- system.file("extdata", "test_single_strand.txdb",
                            package = "tuSelecter2")
txdb_ss <- AnnotationDbi::loadDb(file = txdb_path_ss)

## Convert txdb to granges
gr_ss <- GenomicFeatures::transcripts(txdb_ss)

# Import bigwigs
bw_ss <- system.file("extdata/test_single_strand.bw", package = "tuSelecter2")
bw_multi <- system.file("extdata/test_multichrom.bw", package = "tuSelecter2")

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

test_that("total_coverage counted correctly", {
  expect_equal(total_coverage(bw_multi), 1500)
  expect_error(total_coverage(txdb_path_ss), "txdb is an unsupported file type")
  expect_error(total_coverage("foo"), "file does not exist")
})

test_that("upstream polymerase ratio calculations", {
  # Load in test txdb
  txdb_path <- system.file("extdata", "test_upr.txdb",
                           package = "tuSelecter2")
  txdb <- AnnotationDbi::loadDb(file = txdb_path)
  gr <- GenomicFeatures::transcripts(txdb, c("tx_name", "gene_id"))

  # Load in test bigwigs
  bw_plus <- system.file("extdata", "test_upr_plus.bw",
                        package = "tuSelecter2")
  bw_minus <- system.file("extdata", "test_upr_minus.bw",
                         package = "tuSelecter2")

  tq <- transcript_quantifier(gr, bin_size = 50,
                              transcript_name_column = "tx_name")
  upr <- upstream_polymerase_ratio(tq, bigwig_plus = bw_plus, bigwig_minus = bw_minus,
                                   up_width = 500, up_shift = 0, body_width = 500,
                                   body_shift = 0)

  expected_upr <- c(log2(1 + 1e-3) - log2(1e-3),
                    log2(1 + 1e-3) - log2(1 + 1e-3),
                    log2(2 + 1e-3) - log2(1e-3))
  expect_equivalent(upr, expected_upr)
})

test_that("bigwig seqinfo applied correctly", {
  # Set up GRanges with different seqstyles
  gr_ucsc <- gr_ss
  GenomeInfoDb::seqlevelsStyle(gr_ucsc) <- "UCSC"
  gr_ensembl <- gr_ss
  GenomeInfoDb::seqlevelsStyle(gr_ensembl) <- "Ensembl"
  bw_seqstyle <- GenomeInfoDb::seqlevelsStyle(rtracklayer::BigWigFile(bw_ss))[1]
  bw_seqinfo <- GenomeInfoDb::seqinfo(rtracklayer::BigWigFile(bw_ss))

  # Check that seqinfo is equivalent under varying seqstyles
  expect_equivalent(
    seqinfo(tuSelecter2:::apply_bigwig_seqinfo(x = gr_ucsc, bigwig_file = bw_ss)),
    bw_seqinfo)
  expect_equivalent(
    seqinfo(tuSelecter2:::apply_bigwig_seqinfo(x = gr_ensembl, bigwig_file = bw_ss)),
    bw_seqinfo)
})
