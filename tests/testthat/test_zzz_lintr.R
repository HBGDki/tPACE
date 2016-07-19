
# https://github.com/jimhester/lintr
if (FALSE && requireNamespace("lintr", quietly = TRUE)) {
  context("lints")
  test_that("Package Style", {
    lintr::expect_lint_free(cache = TRUE)
  })
}
