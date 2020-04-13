test_that("Model fitting", {
  # Check correct return class indicating fitting success
  expect_error({
    fit(tq, lambda = -1)
  })
  expect_s4_class({
    tq_fitted <<- fit(tq)
  }, "transcript_quantifier")
  # Check that return output has correct dimensionality
  expect_equal(c(length(tq_fitted@model_abundance),
                 length(tq_fitted@model_abundance[[1]]),
                 length(tq_fitted@model_abundance[[2]])),
               c(2, 2, 2))
  tq_zero <- fit(tq, inactive_transcripts = "t1.1")
  expect_equal(tq@model_abundance[[1]][1], 0)
})

test_that("Gradient estimation", {
  tq_gr <- tq
  tq_gr@model_abundance[[1]] <- c(0.5, 0.5)
  gr_ana_log <- sum_squares_grad(x = tq@model_abundance[[1]], models = tq@models[[1]],
                data = tq@counts[[1]], transform = "log")
  gr_num_log <- num_grad(x = tq@model_abundance[[1]], models = tq@models[[1]],
            data = tq@counts[[1]], transform = "log")
  gr_ana_id <- sum_squares_grad(x = tq@model_abundance[[1]], models = tq@models[[1]],
                        data = tq@counts[[1]], transform = "identity")
  gr_num_id <- num_grad(x = tq@model_abundance[[1]], models = tq@models[[1]],
                         data = tq@counts[[1]], transform = "identity")
  expect_equal(gr_num_id, gr_ana_id)
  expect_equal(gr_num_id, gr_ana_id)
})

# Check that abundance table output is correct
a_tab <- transcript_abundance(tq_fitted)
lookup <- tq_fitted@transcript_model_key[tq_fitted@transcript_model_key$tx_name
                                         == "t2.1", ]

test_that("transcript abundance table", {
  expect_equivalent(as.character(a_tab$transcript_name),
                    tq@transcripts$tx_name)
  expect_equivalent(a_tab[a_tab$transcript_name == "t2.1", ]$abundance,
                    tq_fitted@model_abundance[[lookup$group]][lookup$model])
})

test_that("gene abundance table", {
  tq_gene_fitted <- fit(tq_gene)
  t_tab <- transcript_abundance(tq_gene_fitted)
  g_tab <- gene_abundance(tq_gene_fitted)
  half_t <- t_tab[, .(abundance = sum(abundance) / 2), by = "gene_name"]
  expect_equivalent(g_tab, half_t)
  expect_error(gene_abundance(tq_fitted))
})
