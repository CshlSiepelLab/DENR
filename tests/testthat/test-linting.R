testthat::test_that("Style should be lint-free", {
    testthat::skip_if_not(
        requireNamespace("lintr", quietly = TRUE),
        message = "Package lintr must be installed!"
    )
    # Create list of excluded files
    excluded_files <- c(
        list.files(file.path(lint_path, "data"),
                   recursive = TRUE, full.names = TRUE),
        list.files(file.path(lint_path, "doc"),
                   recursive = TRUE, full.names = TRUE),
        list.files(file.path(lint_path, "inst/doc"),
                   recursive = TRUE, full.names = TRUE),
        list.files(file.path(lint_path, "man"),
                   recursive = TRUE, full.names = TRUE),
        list.files(file.path(lint_path, "vignettes"),
                   recursive = TRUE, full.names = TRUE)
    )
    my_linters <- lintr::with_defaults(cyclocomp_linter =
                                         lintr::cyclocomp_linter(30),
                                       line_length_linter =
                                         lintr::line_length_linter(89)
                                       )
    # Check for lint
    lintr::expect_lint_free(
        path = lint_path,
        relative_path = TRUE,
        exclusions = excluded_files,
        linters = my_linters
    )
})
