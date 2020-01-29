#' @title Sum of squares
#'
#' @description Computes sum-of-squares for a
#'
#' @param x the transcript model abundances
#' @param models the matrix of transcript models
#' @param data a vector of the observed data
#' @param lambda the weight of the lasso penalty
#'
#' @name sum_squares
#' @rdname sum_squares
sum_squares_lasso <- function(x, models, data, lambda = 0) {
  # Normalize the sum-of-squares by the length of the region so the
  # lasso penalty is comparable across regions of different length
  return(
    mean((data - (models %*% x))^2) +
      lambda * sum(abs(x))
  )
}

#' @title Fit model
#'
#' @description Estimates trancript abundances for a given
#' \code{\link{transcript_quantifier-class}} object under a lasso penalty
#'
#' @param threads number of threads to use in model fitting
#' @inheritParams add_data
#' @inheritParams sum_squares
#'
#' @include transcript_quantifier-class.R
#' @name fit
#' @rdname fit
#'
#' @export
methods::setGeneric("fit",
                    function(transcript_quantifier, lambda = 0, threads = 1) {
                      standardGeneric("fit")
                    })

#' @rdname fit
methods::setMethod("fit",
  signature(transcript_quantifier = "transcript_quantifier"),
  function(transcript_quantifier, lambda = 0, threads = 1) {
    if (lambda < 0) {
      stop("lambda must be positive")
    }
    # Alias for ease of use
    tq <- transcript_quantifier
    # Zip together models, counts, and abundances
    sufficient_values <- mapply(function(x, y, z) {
      list(abundance = x, models = y, counts = z)
    }, tq@model_abundance, tq@models, tq@counts, SIMPLIFY = FALSE)
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
                          data = sv$counts,
                          lambda = lambda,
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
