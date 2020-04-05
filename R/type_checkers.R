#' is_matrix_list
#'
#' Checks that all elements of list are a matrix. Will return FALSE for empty
#' list.
#' @param l a list
#'
#' @return boolean
is_matrix_list <- function(l) {
  out <- all(unlist(lapply(l, is.matrix))) && length(l) > 0
  return(out)
}

#' is_vector_list
#'
#' Checks that all elements of list are a matrix. Will return FALSE for empty
#' list.
#' @param l a list
#'
#' @return boolean
is_vector_list <- function(l) {
  out <- all(unlist(lapply(l, is.vector))) && length(l) > 0
  return(out)
}

#' is_strand_vector
#'
#' Checks that all elements are valid strand values ("+","-","*"). Returns FALSE
#' for empty vector.
#' @param x a vector
#' @param allow_star Is "*" allowed in the vector, default TRUE.
#'
#' @return boolean
is_strand_vector <- function(x, allow_star = TRUE) {
    step_one <- (is.vector(x) && length(x) > 0)
    if (allow_star == TRUE) {
        return(step_one && all(as.character(x) %in% c("+", "-", "*")))
    } else {
        return(step_one && all(as.character(x) %in% c("+", "-")))
    }
}

#' matrix_list_dim_equal
#'
#' Checks that all elements of a pair of lists have the same dimensionality
#' @param l1 a list
#' @param l2 a list
#'
#' @return boolean
matrix_list_dim_equal <- function(l1, l2) {
  if (!is_matrix_list(l1) || !is_matrix_list(l2)) {
    stop("l1 and l2 must both be lists of matrices")
  }
  if (length(l1) != length(l2)) {
    stop("lists must be of equal length")
  }
  l1_dim <- lapply(l1, dim)
  l2_dim <- lapply(l2, dim)
  all_dim_equal <- all(mapply(identical, l1_dim, l2_dim, SIMPLIFY = T))
  return(all_dim_equal)
}

#' is_length_two_vector
#'
#' Checks that all elements in the vector is integer and the vector length is 2.
#' @param x a vector
#'
#' @return boolean
is_length_two_vector <- function(x) {
    out <- all(!is.na(x)) && length(x) == 2 && all(x >= 0) && all(x - floor(x) == 0)
    return(out)
}
