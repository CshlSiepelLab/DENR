context("Test type checkers")


test_that("Matrix list checker", {
  ml <- list(matrix(1))
  expect_true(is_matrix_list(ml))
  expect_false(is_matrix_list(unlist(ml)))
  expect_false(is_matrix_list(c(ml, numeric(1))))
  expect_false(is_matrix_list(list()))
})

test_that("Strand vector checker", {
  sv <- c("+", "-", "*")
  expect_true(is_strand_vector(sv))
  expect_false(is_strand_vector(c(sv, "%")))
  expect_false(is_strand_vector(character(0)))
})

test_that("Matrix list dimensions checker", {
  m_1_1 <- matrix(1)
  m_2_1 <- matrix(c(1, 2))
  v_1 <- integer(1)
  expect_true(matrix_list_dim_equal(list(m_1_1, m_2_1),
              list(m_1_1, m_2_1)))
  expect_false(matrix_list_dim_equal(list(m_2_1, m_2_1),
                                     list(m_1_1, m_2_1)))
  expect_error(matrix_list_dim_equal(m_1_1, list(m_1_1)),
               "l1 and l2 must both be lists of matrices")
  expect_error(matrix_list_dim_equal(list(m_2_1),
                                     list(m_1_1, m_2_1)),
               "lists must be of equal length")
})
