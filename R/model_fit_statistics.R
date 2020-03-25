#' predict_transcripts
#'
#' Predicts count data profile per transcript for a single loci
#'
#' @inheritParams predict_locus
#' @param abundance a vector of abundance values, one per transcript
#'
#' @return an matrix of predicted abundances per site (rows) and per model
#' (columns)
predict_transcripts <- function(models, abundance) {
  out <- models
  if (ncol(models) != length(abundance)) {
    stop("Different number of abundances and models")
  }
  for (a in seq_along(abundance)) {
    out[, a] <- out[, a] * abundance[a]
  }
  return(out)
}

#' tx_variance_explained
#'
#' Calculates fraction of variance in the data explained per
#' transcript/combination of \code{n} transcripts at a locus.
#'
#' @inheritParams predict_locus
#' @param abundance a vector of abundance values, one per model
#' @param data a vector of count data
#' @param masks an integer vector of positions to be masked
#' @param n the number of model with non-zero abundance that can be used to
#' explain variance. If fewer models available it will use all available models.
#' Note that some models correspond to multiple transcripts due to high
#' similarity in co-ordinates.
#'
#' @return a vector of the variance explained per model where the names indicate
#' the models each entry corresponds to. When there are multiple models being
#' used they are seperated by a "_" delimiter
tx_variance_explained <- function(models, abundance, data, masks, n = 1) {
  if (!is.numeric(abundance)) {
    stop("abundance must be a numeric vector")
  }
  if (length(abundance) != ncol(models)) {
    stop("Unequal length of abundance vector and number of models")
  }
  data[masks] <- 0
  # Get all the elements which are non-zero in the abundance vector
  non_zero_models <- which(abundance > 0)
  if (length(non_zero_models) == 0) {
    loc_pred <- NULL
  } else {
    if (n > length(non_zero_models)) {
      # If trying to sample more non-zero models than there are decrease n to
      # max viable number to sample (resulting in a single sample)
      n <- length(non_zero_models)
    }
    # generate all combinations of that abundance vector with n non-zero
    # elements. Use negative numbers so to avoid combn treating single positive
    # integers as sequence lengths
    a_mat <- utils::combn(x = -non_zero_models,
                   m = n,
                   FUN = function(x, a) {
                     if (length(x) != 0) {
                       a[setdiff(non_zero_models, abs(x))] <- 0
                     }
                     return(a)
                   }, a = abundance, simplify = TRUE)
    a_mat <- as.matrix(a_mat)
    # Compute locus level variation explained
    loc_pred <- locus_variance_explained(models, a_mat, data, masks)
    names(loc_pred) <-
      apply(a_mat, 2, function(x) paste(which(x > 0), collapse = "_"))
  }
  return(loc_pred)
}

#' locus_variance_explained
#'
#' Calculates fraction of variance in the data explained jointly by all models
#' at a locus
#'
#' @inheritParams predict_locus
#' @param abundance a matrix or vector of abundance values. When it is a
#' matrix each column is a set of abundance values for all models at that loci.
#' This is a more efficient way to calculate the variance explained at the locus
#' for multiple different combinations of abundance values. If it is a vector it
#' must have length equal to the number of models at the loci.
#' @inheritParams tx_variance_explained
#'
#' @include model_fitting.R
#' @return a vector of the total variance explained at the locus by all models
#' with one entry corresponding to each column of the aundance matrix or a
#' single value if abundance is a vector
locus_variance_explained <- function(models, abundance, data, masks) {
  if (is.matrix(abundance)) {
    if (nrow(abundance) != ncol(models)) {
      stop("Unequal number of rows of abundance matrix and number of models")
    }
  } else if (is.numeric(abundance)) {
    if (length(abundance) != ncol(models)) {
      stop("Unequal length of abundance vector and number of models")
    }
  }
  data[masks] <- 0
  locus_pred <- predict_locus(models, abundance)
  residual_model_variance <- apply(locus_pred, 2, function(x) {
    sum((x - data)^2)
    }
  )
  null_model_variance <- sum(data^2)
  variance_explained <- null_model_variance - residual_model_variance
  return(variance_explained)
}

#' locus_statistics
#'
#' Calculates a variety of statistics at the locus and transcript level
#' describing how well the model explains the variance in the data
#'
#' @inheritParams predict_locus
#' @inheritParams fit
#' @inheritParams tx_variance_explained
#'
#' @return A table describing the variance at each locus and the variance
#' explained by the model(s)
#' @export
locus_statistics <- function(transcript_quantifier, n = 1L) {
  if (as.integer(n) != n || n < 1) {
    stop("n must be an integer greater than 1")
  }

  tq <- transcript_quantifier
  # Compute total variance at a locus
  total_var <- unlist(lapply(tq@counts, function(x) sum(x^2)))

  # Compute number of models per locus
  model_count <- unlist(lapply(tq@model_abundance, length))

  # Compute number of models per locus
  nonzero_model_count <- unlist(lapply(tq@model_abundance, function(x) {
    sum(x > 0)
  }))

  # Marginal variance explained by each n transcript model(s)
  tx_var <- mapply(tx_variance_explained,
                models = tq@models,
                abundance = tq@model_abundance,
                data = tq@counts,
                masks = tq@masks,
                n = n)
  tx_unlisted <- unlist(tx_var, use.names = TRUE)
  # Extract the ids where the second column is the group and the third is the
  # models
  ids <- stringr::str_match(names(tx_unlisted), "(\\d+_[+-]).([\\d_]+)")

  # Variance explained by all variance models jointly per locus
  locus_var <-  mapply(locus_variance_explained,
                                 tq@models,
                                 tq@model_abundance,
                                 tq@counts,
                                 tq@masks)

  # Place variance explained in table
  group <- as.integer(factor(ids[, 2], levels = names(tq@models)))
  tx_var_table <- data.table::data.table(
    group,
    total_models = model_count[group],
    nonzero_models = nonzero_model_count[group],
    models = ids[, 3],
    variance_explained_model = tx_unlisted,
    variance_explained_all = locus_var[group],
    total_locus_variance = total_var[group])
  return(tx_var_table)
}
