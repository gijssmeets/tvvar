test_that("opti_fct runs and returns scalar for purpose='optim'", {
  set.seed(1)
  # tiny toy data
  VAR.data <- matrix(rnorm(100), ncol = 2)
  n_var <- ncol(VAR.data); p <- 1L; r <- 1L
  
  # minimal masks / shapes
  Phi.f.mask <- matrix(0, nrow = n_var, ncol = n_var * p + 1L)
  
  # params: length must match your unpacking
  # 2 (A,B) + r (phi_r) + sum(Phi.f!=0)=0 + Nstar + (n_var^2*p) for zero.mean
  Nstar <- n_var * (n_var + 1L) / 2L
  pars <- c(rnorm(2), rnorm(r), rnorm(Nstar), rnorm(n_var^2 * p))
  
  val <- opti_fct(
    params = pars,
    VAR.data = VAR.data,
    Phi.f.mask = Phi.f.mask,
    zero.mean = TRUE,
    Smooth = FALSE,
    purpose = "optim",
    r = r, p = p, n_var = n_var
  )
  
  expect_true(is.numeric(val) && length(val) == 1L && is.finite(val))
})

test_that("opti_fct returns named list for purpose='eval'", {
  set.seed(1)
  VAR.data <- matrix(rnorm(100), ncol = 2)
  n_var <- ncol(VAR.data); p <- 1L; r <- 1L
  Phi.f.mask <- matrix(0, nrow = n_var, ncol = n_var * p + 1L)
  Nstar <- n_var * (n_var + 1L) / 2L
  pars <- c(rnorm(2), rnorm(r), rnorm(Nstar), rnorm(n_var^2 * p))
  
  res <- opti_fct(
    params = pars,
    VAR.data = VAR.data,
    Phi.f.mask = Phi.f.mask,
    zero.mean = TRUE,
    Smooth = FALSE,
    purpose = "eval",
    r = r, p = p, n_var = n_var
  )
  
  expect_type(res, "list")
  expect_true(all(c("average.L","full.L","Phi.c","array.filtered.H","momega","Nstar") %in% names(res)))
})