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
  if (transform == "identity") {
      ss <- mean((data - (models %*% x))^2) + lambda * sum(abs(x))
  } else if (transform == "log") {
      ss <- mean((log(data + 1e-3) - log((models %*% x) + 1e-3))^2) +
        lambda * sum(abs(x))
  } else {
    stop("Invalid transform specified")
  }
  return(ss)
}

#' @title Fit model
#'
#' @description Estimates trancript abundances for a given
#' \code{\link{transcript_quantifier-class}} object under a lasso penalty
#'
#' @param threads number of threads to use in model fitting
#' @inheritParams add_data
#' @inheritParams sum_squares_lasso
#'
#' @include transcript_quantifier-class.R
#' @name fit
#' @rdname fit
#'
#' @export
methods::setGeneric("fit",
                    function(transcript_quantifier, lambda = 0, threads = 1,
                             transform = "log") {
                      standardGeneric("fit")
                    })

#' @rdname fit
methods::setMethod("fit",
  signature(transcript_quantifier = "transcript_quantifier"),
  function(transcript_quantifier, lambda = 0, threads = 1,
           transform = "log") {
    if (lambda < 0) {
      stop("lambda must be positive")
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

    # Alias for ease of use
    tq <- transcript_quantifier
    # Zip together models, counts, and abundances
    sufficient_values <- mapply(function(x, y, z, za) {
      list(abundance = x, models = y, counts = z, masks = za)
    }, tq@model_abundance, tq@models, tq@counts, tq@masks, SIMPLIFY = FALSE)
    # Create iterator over list of abundances
    sv_iter <- iterators::iter(sufficient_values)

    if (threads > 1) {
      # Create parallel cluster and register it with the foreach backend
      cluster <- snow::makeCluster(threads)
      doSNOW::registerDoSNOW(cluster)
      `%doloop%` <- foreach::`%dopar%`
    } else {
      `%doloop%` <- foreach::`%do%`
    }

    estim <- foreach::foreach(sv = sv_iter) %doloop% {
      opt_result <- stats::optim(sv$abundance, fn = sum_squares_lasso,
                          models = sv$models,
                          data = mask_data(sv$counts, sv$masks),
                          lambda = lambda,
                          transform = transform,
                          lower = rep(0, length(sv$abundance)),
                          method = "L-BFGS-B")
      abundance_est <- opt_result$par
      names(abundance_est) <- colnames(sv$models)
      return(abundance_est)
    }

    if (threads > 1) {
      snow::stopCluster(cluster)
    }

    tq@model_abundance <- estim
    return(tq)
  }
)

## Appease R CMD check
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("sv"))
}
