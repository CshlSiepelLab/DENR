# Load data
suppressMessages({
  txdb_path_ds <- system.file("extdata", "test_double_strand.txdb",
                              package = "DENR")
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
                   package = "DENR")
bwm <- system.file("extdata", "test_double_strand_minus.bw",
                   package = "DENR")
# Add the data
tq <- add_data(tq = tq,
               bigwig_plus = bwp,
               bigwig_minus = bwm)

# Create transcript_quantifier object with gene annotations, also double all transcripts
# to ensure that duplicate models are being summed over correctly
gr_ds_doubled <- rep(gr_ds, 2)
gr_ds_doubled$tx_name <- paste(gr_ds_doubled$tx_name, seq_along(gr_ds_doubled),
                               sep = ".")
tq_gene <- suppressMessages(transcript_quantifier(gr_ds_doubled, bin_size = 50,
                                             transcript_name_column = "tx_name",
                                             gene_name_column = "gene_id",
                                             mask_start_bins = c(5, 5)))
# Add the data
tq <- add_data(tq = tq,
               bigwig_plus = bwp,
               bigwig_minus = bwm)

# Add the data
tq_gene <- add_data(tq = tq_gene,
               bigwig_plus = bwp,
               bigwig_minus = bwm)

num_grad <- function(x, models, data, transform, eps = 1e-6) {
  grad <- numeric(length(x))
  for (i in seq_along(x)) {
    xp <- x
    xm <- x
    xp[i] <- x[i] + eps
    xm[i] <- x[i] - eps
    ssp <- sum_squares_lasso(xp, models, data, transform = transform)
    ssm <- sum_squares_lasso(xm, models, data, transform = transform)
    grad[i] <- (ssp - ssm) / (xp[i] - xm[i])
  }
  return(grad)
}
