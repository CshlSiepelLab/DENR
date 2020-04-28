# Load data
suppressMessages({
  txdb_path_ds <- system.file("extdata", "test_double_strand.txdb",
                              package = "tuSelecter2")
  txdb_ds <- AnnotationDbi::loadDb(file = txdb_path_ds)
  gr_ds <- GenomicFeatures::transcripts(txdb_ds, c("tx_name", "gene_id"))
  gr_ds$gene_id <- as.character(gr_ds$gene_id)
})

# Create transcript_quantifier object
tq <- transcript_quantifier(gr_ds, bin_size = 50,
                            transcript_name_column = "tx_name",
                            gene_name_column = "gene_id",
                            mask_start_bins = c(5, 5))

# The paths to relevant bigwig files
bwp <- system.file("extdata", "test_double_strand_plus.bw",
                   package = "tuSelecter2")
bwm <- system.file("extdata", "test_double_strand_minus.bw",
                   package = "tuSelecter2")
# Add the data
tq <- add_data(tq = tq,
               bigwig_plus = bwp,
               bigwig_minus = bwm)

test_that("get_transcripts", {
  expect_equivalent(tuSelecter2:::get_transcripts(tq = tq,
                                gene_name = "g1"),
                    gr_ds[gr_ds$tx_name %in% c("t1.1", "t1.2", "t2.1")])
  expect_equivalent(tuSelecter2:::get_transcripts(tq = tq,
                                                  chrom = 1,
                                                  start = 1000,
                                                  end = 10000),
                    gr_ds)
  expect_equivalent(tuSelecter2:::get_transcripts(tq = tq,
                                                  chrom = 1,
                                                  start = 1000,
                                                  end = 10000,
                                                  strand = "+"),
                    gr_ds[GenomicRanges::strand(gr_ds) == "+"])
  expect_equivalent(tuSelecter2:::get_transcripts(tq = tq,
                                                  chrom = 1,
                                                  start = 1000,
                                                  end = 10000,
                                                  strand = "-"),
                    gr_ds[GenomicRanges::strand(gr_ds) == "-"])
  expect_setequal(tuSelecter2:::get_transcripts(tq = tq,
                                                  chrom = 1,
                                                  start = 7000,
                                                  end = 10000)$tx_name,
                    c("t1.1", "t1.2"))
})

test_that("get_masks", {
  expect_equivalent(tuSelecter2:::get_masks(tq, "t1.1"),
                    tq@bins[[2]][tq@masks[[2]]])
  expect_equivalent(tuSelecter2:::get_masks(tq, c("t2.1", "t1.2")),
                    c(tq@bins[[1]][tq@masks[[1]]],
                      (tq@bins[[2]][tq@masks[[2]]])))
})

test_that("get_data", {
  test_dat <- tuSelecter2:::get_data(tq, chrom = 1, start = 1,
                                           end = 6500)
  # Seperately compute the real bins and counts
  test_granges <- GenomicRanges::GRanges(1, IRanges::IRanges(1, 6500))
  hits <- lapply(tq@bins, function(x)
    GenomicRanges::findOverlaps(x, test_granges))
  test_counts <- mapply(function(count, hit) {
   return(count[S4Vectors::queryHits(hit)])
  }, tq@counts, hits)
  # compare the two outputs
  expect_equivalent(tq@counts$`1_-:1`, test_counts$`1_-:1`)
  expect_equivalent(tq@counts$`1_+:1`[1:27], test_counts$`1_+:1`)
})

test_that("set_datatrack_ylim", {
    dtrack_plus <- Gviz::DataTrack(rtracklayer::import(bwp), strand = "+")
    dtrack_plus <- set_datatrack_ylim(dtrack_plus, ylim = c(0, 10))
    dtrack_minus <- Gviz::DataTrack(rtracklayer::import(bwm), strand = "-")
    dtrack_minus <- set_datatrack_ylim(dtrack_minus, ylim = c(0, 10))
    # compare the ylim and see if it has been set
    expect_equivalent(as.vector(Gviz::displayPars(dtrack_plus)$ylim), c(0, 10))
    expect_equivalent(as.vector(Gviz::displayPars(dtrack_minus)$ylim), c(10, 0))
})

test_that("plot_model", {
  tmp <- tempfile()
  grDevices::pdf(tmp)
  expect_silent({
    pt1 <- plot_model(tq, gene_name = c("g1", "g2"))
  })
  tq_fitted <- fit(tq)
  expect_silent({
    pt2 <- plot_model(tq_fitted, gene_name = c("g1", "g2"))
  })
  expect_silent({
    pt3 <- plot_model(tq_fitted, gene_name = "g1",
                      bigwig_plus = bwp, bigwig_minus = bwm)
  })
  # Test errors
  tq_err <- tq
  tq_err@column_identifiers[2] <- NA
  expect_error(plot_model(tq, gene_name = "foo"),
               "Gene foo does not exist")
  expect_error(plot_model(tq_err, gene_name = "g2"),
               "transcript_quantifier object not built with gene ids")
  expect_error(plot_model(tq, chrom = 1, start = 1, end = 1000),
               "Positional query does not intersect any transcripts")
  expect_error(plot_model(tq, chrom = c(1, 2), start = 1, end = 1000),
               "Only one value may be specified for each positional option")
  expect_error(plot_model(tq, gene_name = "g1",
                          chrom = 1, start = 1, end = 1000),
               "Only gene name OR positional information can be specified")
  expect_error(plot_model(tq), "Incorrect positional specification")
  expect_error(plot_model(tq, chrom = 1),
               "Incorrect positional specification")
  expect_error(plot_model(tq, chrom = 1, start = 1),
               "Incorrect positional specification")
  expect_error(plot_model(tq_fitted, gene_name = "g1",
                          bigwig_plus = "A"),
               "File A does not exist")
  grDevices::dev.off()
  unlink(tmp, TRUE)
})
