test_that("unpenalized_estimate runs on toy data", {
  set.seed(123)
  X <- matrix(rnorm(3 * 80), ncol = 3)  # 80 rows × 3 cols
  N <- ncol(X); p <- 2; r <- 1
  Phi_block <- diag(1, N)               # simple structure
  
  fit <- unpenalized_estimate(
    data = X, p = p, r = r,
    zero_mean = TRUE,
    phi_f_structure = Phi_block,
    control = list(maxit = 10)
  )
  
  expect_s3_class(fit, "tvvar_fit")
  expect_named(fit$ic, c("AIC","AICc","BIC"))
})