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
})

# Check that abundance table output is correct
a_tab <- abundance_table(tq_fitted)
lookup <- tq_fitted@transcript_model_key[tq_fitted@transcript_model_key$tx_name
                                         == "t2.1", ]

test_that("Abundance table", {
  expect_equivalent(as.character(a_tab$transcript_name),
                    tq@transcripts$tx_name)
  expect_equivalent(a_tab[a_tab$transcript_name == "t2.1", ]$abundance,
                    tq_fitted@model_abundance[[lookup$group]][lookup$model])
})
