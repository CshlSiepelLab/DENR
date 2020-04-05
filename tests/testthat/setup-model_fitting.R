# Load data
suppressMessages({
  txdb_path_ds <- system.file("extdata", "test_double_strand.txdb",
                              package = "tuSelecter2")
  txdb_ds <- AnnotationDbi::loadDb(file = txdb_path_ds)
  gr_ds <- GenomicFeatures::transcripts(txdb_ds, c("tx_name", "gene_id"))
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
# Add the data
tq <- add_data(tq = tq,
               bigwig_plus = bwp,
               bigwig_minus = bwm,
               summary_operation = "mean")
