#' @title Mask data
#'
#' @description masks positions in data vector at indices specified in masks
#' vector by replacing them with zeros.
#'
#' @param data a data vector
#' @param masks a vector of masked indicies
#'
#' @name mask_data
#' @rdname mask_data
mask_data <- function(data, masks) {
  data[masks] <- 0
  return(data)
}

#' @title Predicts signal at locus level
#'
#' @description outputs prediction for read counts across entire locus
#'
#' @param models the matrix of transcript models
#' @param abundance the transcript model abundances (can be either a vector or
#' a matrix where each column is a vecotor of abundances)
#'
#' @name predict_locus
#' @rdname predict_locus
predict_locus <- function(models, abundance) {
  return(models %*% abundance)
}

#' @title Sum of squares
#'
#' @description Computes sum-of-squares for a
#'
#' @param x the transcript model abundances
#' @param models the matrix of transcript models
#' @param data a vector of the observed data
#' @param lambda the weight of the lasso penalty
#' @param transform how to transform the data (either "identity" or "log")
#'
#' @name sum_squares_lasso
#' @rdname sum_squares_lasso
sum_squares_lasso <- function(x, models, data, lambda = 0,
                              transform) {
  # Normalize the sum-of-squares by the length of the region so the
  # lasso penalty is comparable across regions of different length
  locus_pred <- predict_locus(models, x)
  if (transform == "identity") {
      sum_sq <- sum((data - locus_pred)^2)
  } else if (transform == "log") {
      sum_sq <- sum((log(data + 1e-3) - log(locus_pred + 1e-3))^2)
  } else {
    stop("Invalid transform specified")
  }
  return(sum_sq)
}

#' @title Fit model
#'
#' @description Estimates trancript abundances for a given
#' \code{\link{transcript_quantifier-class}} object under a lasso penalty
#'
#' @param verbose if TRUE shows progress bar for fitting (default: FALSE)
#' @inheritParams add_data
#' @inheritParams sum_squares_lasso
#'
#' @include transcript_quantifier-class.R
#' @name fit
#' @rdname fit
#'
#' @export
methods::setGeneric("fit",
                    function(tq, lambda = 0,
                             transform = "log", verbose = FALSE) {
                      standardGeneric("fit")
                    })

#' @rdname fit
methods::setMethod("fit",
  signature(tq = "transcript_quantifier"),
  function(tq, lambda = 0, transform = "log",
           verbose = FALSE) {
    if (lambda < 0) {
      stop("lambda must be positive")
    }

    if (length(tq@counts) == 0) {
      stop("No count data has been entered for this model")
    }

    if (length(tq@counts) != length(tq@models)) {
      stop("There is data for a different number of loci than there are models")
    }

    if (!is.logical(verbose)) {
      stop("verbose most be either TRUE or FALSE")
    }

    # Handle transform option
    transform_opts <- c("identity", "log")
    if (length(transform) > 1) {
      warning("Multiple transform options specified, using the first")
      transform <- transform[1]
    }
    if (!transform %in% transform_opts) {
      stop(paste("transform must be one of these:",
           paste(transform_opts, collapse = ", ")))
    }

    # Zip together models, counts, and abundances
    sufficient_values <- mapply(function(x, y, z, za) {
      list(abundance = x, models = y, counts = z, masks = za)
    }, tq@model_abundance, tq@models, tq@counts, tq@masks, SIMPLIFY = FALSE)

    estim <- list()
    if (verbose) {
      message("Estimating abundance ...")
      pb <- utils::txtProgressBar(min = 1, max = length(sufficient_values), style = 3)
    }

    for (i in seq_along(sufficient_values)) {
      sv <- sufficient_values[[i]]
      opt_result <- stats::optim(sv$abundance, fn = sum_squares_lasso,
                                 models = sv$models,
                                 data = mask_data(sv$counts, sv$masks),
                                 lambda = lambda,
                                 transform = transform,
                                 lower = rep(0, length(sv$abundance)),
                                 method = "L-BFGS-B")
      estim[[i]] <- opt_result$par
      names(estim[[i]]) <- colnames(sv$models)
      if (verbose) {
        utils::setTxtProgressBar(pb, i)
      }
    }

    tq@model_abundance <- estim
    return(tq)
  }
)
