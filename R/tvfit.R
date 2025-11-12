# ---- small helper ----
`%||%` <- function(x, y) if (!is.null(x)) x else y

# ---- shared setup ----
.tvvar_shared_setup <- function(data, p, r, zero.mean, factor.structure) {
  Y <- as.matrix(data)
  N <- ncol(Y); T.fin <- nrow(Y)
  phi_arr      <- .normalize_phi_structure(factor.structure, N, r)
  Phi.f.array <- make.Phi.f(structure = phi_arr, lag.order = p)
  cfg <- list(dim.VAR = N, lag.order = p, number.factors = r, zero.mean = isTRUE(zero.mean))
  list(Y = Y, N = N, T.fin = T.fin, Phi.f.array = Phi.f.array, cfg = cfg)
}

# ---- shared initialization ----
.tvvar_build_par_init <- function(Y, p, r, zero.mean, Phi.f.array,
                                  init = c("default", "random", "custom"),
                                  init_list = NULL) {
  init <- match.arg(init)
  N <- ncol(Y)
  logit <- function(x) log(x/(1 - x))
  
  # Defaults
  A0 <- 0.10; B0 <- 0.80
  phi_r0 <- if (r > 0) rep(0.95, r) else numeric(0L)
  Omega0 <- matrix(0.2, N, N); diag(Omega0) <- 0.3
  
  # Phi_c0 via static VAR (Gorgi-style)
  Y_named <- as.data.frame(Y)
  colnames(Y_named) <- paste0("y", seq_len(ncol(Y_named)))
  B.t <- matrix(0, nrow = N, ncol = (N * p + 1))
  for (j in 1:N) {
    B.t[j, ] <- stats::coefficients(vars::VAR(Y_named, p = p))[[j]][, 1]
  }
  A.t <- B.t[, -(N * p + 1), drop = FALSE]
  if (isFALSE(zero.mean)) A.t <- cbind(B.t[, (N * p + 1), drop = FALSE], A.t)
  Phi_c0 <- A.t
  
  # mask & free params
  Phi.f.as.matrix <- matrix(Phi.f.array, nrow = N)
  Phi.f.mask  <- (Phi.f.as.matrix != 0)
  n_phi_free  <- sum(Phi.f.mask)
  
  # --- handle init type ---
  if (init == "random") {
    set.seed(NULL)
    A0 <- runif(1, 0.05, 0.95)
    B0 <- runif(1, 0.05, 0.95)
    phi_r0 <- if (r > 0) runif(r, 0.5, 0.98) else numeric(0L)
    L <- matrix(rnorm(N^2, sd = 0.2), N, N)
    diag(L) <- abs(rnorm(N, mean = 0.6, sd = 0.1))
    Omega0 <- tcrossprod(L)
    s <- mean(diag(Omega0)); if (is.finite(s) && s > 0) Omega0 <- Omega0 * (0.3 / s)
    Phi_c0 <- Phi_c0 + matrix(rnorm(length(Phi_c0), sd = 0.05), nrow = nrow(Phi_c0))
    Phi_f0 <- rnorm(n_phi_free, mean = 0, sd = 0.05)
  } else if (init == "custom") {
    if (is.null(init_list)) stop("Provide init_list when init='custom'.")
    A0 <- init_list$A %||% A0
    B0 <- init_list$B %||% B0
    phi_r0 <- init_list$phi_r %||% phi_r0
    Omega0 <- init_list$Omega %||% Omega0
    Phi_c0 <- init_list$Phi_c %||% Phi_c0
    Phi_f0 <- init_list$Phi_f %||% rep(0.01, n_phi_free)
  } else {
    Phi_f0 <- rep(0.01, n_phi_free)
  }
  
  par_init <- c(
    if (length(phi_r0)) c(logit(A0), logit(B0), logit(phi_r0)) else c(logit(A0), logit(B0)),
    Phi_f0,
    vech(Omega0),
    as.vector(Phi_c0)
  )
  
  list(
    par_init   = par_init,
    Phi_c0     = Phi_c0,
    Phi_f_mask = Phi.f.mask,
    n_phi_free = n_phi_free,
    init_mode  = init
  )
}

#' Time-varying VAR estimation (ML or EM)
#'
#' @description
#' Fits a time-varying VAR (TV–VAR) using either direct maximum likelihood (ML)
#' or an EM routine, without sparsity penalties. Produces an object compatible
#' with downstream functions such as [tvirf()] and [tvpred()].
#'
#' @param data Numeric matrix \code{T x N} (rows = time, cols = variables).
#' @param p Integer VAR lag order.
#' @param r Integer number of latent factors (state dimension).
#' @param zero.mean Logical; if \code{TRUE}, intercepts in \eqn{\Phi^c} are fixed at 0.
#' @param factor_structure Specification of the factor loading structure Φᶠ.
#'   See \code{\link{tvpenfit}} for details on accepted forms and identification restrictions.
#' @param method One of \code{"ML"} (direct BFGS on likelihood) or \code{"EM"}.
#' @param control List of optimizer controls passed to ML (e.g. \code{list(maxit=2000, trace=0)}).
#' @param em_control List of EM controls
#'   (e.g. \code{list(max_iter=200, tol=1e-3, trace=TRUE)}).
#' @param init Initialization mode: \code{"default"}, \code{"random"}, or \code{"custom"}.
#'   See Details.
#' @param init_list Named list of custom initial values used when \code{init="custom"}.
#'   Recognized elements: \code{A}, \code{B}, \code{phi_r}, \code{Omega} (\eqn{N \times N}),
#'   \code{Phi_c} (\eqn{N \times (1+Np)} or \eqn{N \times Np} if \code{zero.mean=TRUE}),
#'   and \code{Phi_f} (vector for the free \eqn{\Phi^f} entries in mask order).
#'
#' @details
#' Initialization:
#' \itemize{
#'   \item \strong{default}: \eqn{A=0.10}, \eqn{B=0.80}, \eqn{\phi_r=0.95}, \eqn{\Omega} with
#'     0.3 on diagonal and 0.2 off-diagonal, and small 0.01 values for free \eqn{\Phi^f}.
#'     \eqn{\Phi^c} is initialized via a static VAR (Gorgi-style), respecting \code{zero.mean}.
#'   \item \strong{random}: draws reasonable random starting values (bounded and PD).
#'   \item \strong{custom}: take values from \code{init_list} (fallback to defaults if a field is missing).
#' }
#'
#' The returned object has class \code{"tvfit"} and contains:
#' \itemize{
#'   \item \code{$estimate}: list with \code{A}, \code{B}, \code{phi_r}, \code{Phi_f}, \code{Phi_c}, \code{Omega}.
#'   \item \code{$lik}: average and summed log-likelihood.
#'   \item \code{$ic}: AIC, AICc, BIC (computed in the usual unpenalized way).
#'   \item \code{$eval}: last evaluation bundle (e.g., filtered states).
#'   \item \code{$optim}: optimizer details.
#'   \item \code{$vcov}: asymptotic covariance matrix (if computed in your ML core).
#'   \item \code{$theta}: final parameter vector (unconstrained scale).
#'   \item \code{$meta}: dims, method, free pattern, and \code{time} (runtime in seconds).
#' }
#'
#' @return A \code{tvfit} list (same structure as from \code{tvpenfit()}).
#' @seealso [tvpenfit()], [tvirf()], [tvpred()]
#' @export
#'
#' @examples
#' \dontrun{
#' # ML with default init
#' fit_ml <- tvfit(
#'   data = simdata$Y, p = 1, r = 1, zero.mean = TRUE,
#'   factor.structure = matrix(1, 2, 2), method = "ML"
#' )
#'
#' # EM with random init
#' fit_em <- tvfit(
#'   data = simdata$Y, p = 1, r = 1, zero.mean = TRUE,
#'   factor.structure = matrix(1, 2, 2), method = "EM",
#'   init = "random"
#' )
#'
#' # ML with custom init for A and B only (others fall back to defaults)
#' fit_ml_custom <- tvfit(
#'   data = simdata$Y, p = 1, r = 1, zero.mean = TRUE,
#'   factor.structure = matrix(1, 2, 2), method = "ML",
#'   init = "custom", init_list = list(A = 0.2, B = 0.7)
#' )
#' }
tvfit <- function(data, p = 1, r = 1, zero.mean = TRUE,
                  factor.structure = NULL,
                  method = c("ML","EM"),
                  control = list(maxit = 2000, trace = 0),
                  em_control = list(max_iter = 200, tol = 1e-3, trace = TRUE),
                  init = c("default", "random", "custom"),
                  init_list = NULL) {
  method <- match.arg(method)
  init <- match.arg(init)
  
  
  if (is.null(factor.structure)) {
    factor.structure <- matrix(1, nrow = dim(data), ncol = dim(data))
    #TODO: Zero-mean = FALSE option
  }
  
  ss <- .tvvar_shared_setup(data, p, r, zero.mean, factor.structure)
  initvals <- .tvvar_build_par_init(ss$Y, p, r, zero.mean, ss$Phi.f.array, init, init_list)
  
  start_time <- Sys.time()
  fit <- if (identical(method, "ML")) {
    .unpen_ml(ss, initvals, control)
  } else {
    .unpen_em(ss, initvals, em_control)
  }
  fit$meta$time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  fit$meta$init <- initvals$init_mode
  fit$meta$data <- data
  fit
}

# ---- ML ----
.unpen_ml <- function(ss, init, control) {
  Y <- ss$Y; cfg <- ss$cfg; Phi.f.array <- ss$Phi.f.array
  p <- cfg$lag.order; r <- cfg$number.factors; N <- cfg$dim.VAR
  
  opt <- optim(
    par        = init$par_init,
    fn         = opti.fct,
    par_fixed  = NaN,
    VAR.data   = Y,
    Phi.f.array= Phi.f.array,
    cfg        = cfg,
    Smooth     = FALSE,
    purpose    = "optim",
    method     = "BFGS",
    control    = control
  )
  
  ev <- opti.fct(
    par_free    = opt$par,
    par_fixed   = NaN,
    VAR.data    = Y,
    Phi.f.array = Phi.f.array,
    cfg         = cfg,
    Smooth      = TRUE,
    purpose     = "eval"
  )
  
  p_out <- opt$par
  my.hess <- optimHess(
    p_out, hessian.fct.untr,
    VAR.data = Y, Phi.f.array = Phi.f.array,
    zero.mean = cfg$zero.mean, cfg = cfg
  )
  cov.matrix <- tryCatch(solve(my.hess * ss$T.fin), error = function(e) diag(1e-8, length(p_out)))
  
  # unpack
  inv_logit <- plogis
  idx <- 1L
  A <- inv_logit(opt$par[idx]); idx <- idx + 1L
  B <- inv_logit(opt$par[idx]); idx <- idx + 1L
  phi_r <- if (r > 0) inv_logit(opt$par[idx:(idx + r - 1L)]) else numeric(0L)
  idx <- idx + length(phi_r)
  
  Phi.f.as.matrix <- matrix(Phi.f.array, nrow = N)
  n_phi_free <- sum(Phi.f.as.matrix != 0)
  Phi.f.est <- matrix(0, nrow = nrow(Phi.f.as.matrix), ncol = ncol(Phi.f.as.matrix))
  if (n_phi_free > 0) {
    Phi.f.est[Phi.f.as.matrix != 0] <- opt$par[idx:(idx + n_phi_free - 1L)]
    idx <- idx + n_phi_free
  }
  Phi.f.est <- array(Phi.f.est, dim = dim(Phi.f.array))
  
  Nstar <- N * (N + 1L) / 2L
  L_vec <- opt$par[idx:(idx + Nstar - 1L)]; idx <- idx + Nstar
  Lmat  <- matrix(D.matrix(N) %*% L_vec, ncol = N)
  Lmat[upper.tri(Lmat)] <- 0
  Omega <- Lmat %*% t(Lmat)
  
  need <- if (isTRUE(cfg$zero.mean)) N^2 * p else (N^2 * p + N)
  Phi_c_vec <- opt$par[idx:(idx + need - 1L)]
  Phi_c <- if (isTRUE(cfg$zero.mean)) cbind(0, matrix(Phi_c_vec, nrow = N)) else matrix(Phi_c_vec, nrow = N)
  
  K <- length(c(A, B, phi_r)) + n_phi_free + Nstar + length(Phi_c_vec)
  n_eff <- (ss$T.fin - p) * N
  AIC  <-  2 * K + 2 * ev$average.L * (ss$T.fin - p)
  AICc <-  AIC + (2*K^2 + 2*K) / (max(1, n_eff - K - 1))
  BIC  <-  log(n_eff) * K + 2 * ev$average.L * (ss$T.fin - p)
  
  out <- list(
    estimate = list(A = A, B = B, phi_r = phi_r,
                    Phi_f = Phi.f.est, Phi_c = Phi_c, Omega = Omega),
    lik = list(avg = ev$average.L, sum = ev$full.L),
    ic  = list(AIC = AIC, AICc = AICc, BIC = BIC),
    eval = ev,
    optim = list(par = opt$par, value = opt$value, convergence = opt$convergence,
                 counts = opt$counts, message = opt$message),
    vcov = cov.matrix,
    theta = p_out,
    meta = list(N = N, p = p, r = r, zero.mean = isTRUE(cfg$zero.mean),
                nobs = ss$T.fin, method = "unpenalized", Phi.f.array = Phi.f.array)
  )
  class(out) <- "tvfit"
  out
}


.unpen_em <- function(ss, init, em_control) {
  Y <- ss$Y; cfg <- ss$cfg; Phi.f.array <- ss$Phi.f.array
  p <- cfg$lag.order; r <- cfg$number.factors; N <- cfg$dim.VAR
  
  tol       <- em_control$tol %||% 1e-3
  max_iter  <- em_control$max_iter %||% 200
  do_trace  <- isTRUE(em_control$trace)
  
  # Bounds for (A,B,phi_r)
  lower_bounds <- rep(-Inf, length(init$par_init))
  upper_bounds <- rep( Inf, length(init$par_init))
  n_head <- 2 + r
  if (n_head > 0) { lower_bounds[1:n_head] <- -10; upper_bounds[1:n_head] <- 10 }
  
  params   <- init$par_init
  obj_prev <- Inf
  iter     <- 1L
  repeat {
    # --- E-step ---
    invisible(opti.fct(
      par_free    = params,
      par_fixed   = NaN,
      VAR.data    = Y,
      Phi.f.array = Phi.f.array,
      cfg         = cfg,
      Smooth      = TRUE,
      purpose     = "eval"
    ))
    
    # --- M-step ---
    opt <- optim(
      par        = params,
      fn         = opti.fct,
      par_fixed  = NaN,
      VAR.data   = Y,
      Phi.f.array= Phi.f.array,
      cfg        = cfg,
      Smooth     = FALSE,
      purpose    = "optim",
      method     = "L-BFGS-B",
      lower      = lower_bounds,
      upper      = upper_bounds,
      control    = list(maxit = 500, trace = 0)
    )
    params <- opt$par
    obj    <- opt$value
    
    # --- convergence check ---
    rel <- 2 * abs(obj - obj_prev) / (abs(obj) + abs(obj_prev))
    if (do_trace) message(sprintf("[EM %03d] obj=%.6f  rel=%.3e", iter, obj, rel))
    if (!is.finite(rel) || rel < tol || iter >= max_iter) break
    obj_prev <- obj; iter <- iter + 1L
  }
  
  # --- final evaluation (with smoother) ---
  ev <- opti.fct(
    par_free    = params,
    par_fixed   = NaN,
    VAR.data    = Y,
    Phi.f.array = Phi.f.array,
    cfg         = cfg,
    Smooth      = TRUE,
    purpose     = "eval"
  )
  
  # --- Hessian for vcov ---
  my.hess <- optimHess(
    params, hessian.fct.untr,
    VAR.data = Y, Phi.f.array = Phi.f.array,
    zero.mean = isTRUE(cfg$zero.mean), cfg = cfg
  )
  cov.matrix <- tryCatch(solve(my.hess * ss$T.fin), error = function(e) diag(1e-8, length(params)))
  
  # --- unpack parameters ---
  inv_logit <- plogis
  idx <- 1L
  A  <- inv_logit(params[idx]); idx <- idx + 1L
  B  <- inv_logit(params[idx]); idx <- idx + 1L
  phi_r <- if (r > 0) { x <- inv_logit(params[idx:(idx + r - 1L)]); idx <- idx + r; x } else numeric(0L)
  
  Phi.f.as.matrix <- matrix(Phi.f.array, nrow = N)
  n_phi_free <- sum(Phi.f.as.matrix != 0)
  
  Phi.f.est <- matrix(0, nrow = nrow(Phi.f.as.matrix), ncol = ncol(Phi.f.as.matrix))
  if (n_phi_free > 0) {
    Phi.f.est[Phi.f.as.matrix != 0] <- params[idx:(idx + n_phi_free - 1L)]
    idx <- idx + n_phi_free
  }
  Phi.f.est <- array(Phi.f.est, dim = dim(Phi.f.array))
  
  Nstar <- N * (N + 1L) / 2L
  L_vec <- params[idx:(idx + Nstar - 1L)]; idx <- idx + Nstar
  Lmat  <- matrix(D.matrix(N) %*% L_vec, ncol = N)
  Lmat[upper.tri(Lmat)] <- 0
  Omega <- Lmat %*% t(Lmat)
  
  need <- if (isTRUE(cfg$zero.mean)) N^2 * p else (N^2 * p + N)
  Phi_c_vec <- params[idx:(idx + need - 1L)]
  Phi_c <- if (isTRUE(cfg$zero.mean)) cbind(0, matrix(Phi_c_vec, nrow = N)) else matrix(Phi_c_vec, nrow = N)
  
  # --- Information criteria ---
  K <- length(c(A, B, phi_r)) + n_phi_free + Nstar + length(Phi_c_vec)
  n_eff <- (ss$T.fin - p) * N
  AIC  <-  2 * K + 2 * ev$average.L * (ss$T.fin - p)
  AICc <-  AIC + (2*K^2 + 2*K) / (max(1, n_eff - K - 1))
  BIC  <-  log(n_eff) * K + 2 * ev$average.L * (ss$T.fin - p)
  
  # --- final output (same as ML) ---
  out <- list(
    estimate = list(A = A, B = B, phi_r = phi_r,
                    Phi_f = Phi.f.est, Phi_c = Phi_c, Omega = Omega),
    lik = list(avg = ev$average.L, sum = ev$full.L),
    ic  = list(AIC = AIC, AICc = AICc, BIC = BIC),
    eval = ev,
    optim = list(par = params, value = obj, convergence = NA,
                 counts = NA, message = sprintf("EM finished in %d iterations", iter)),
    vcov  = cov.matrix,
    theta = params,
    meta = list(N = N, p = p, r = r, zero.mean = isTRUE(cfg$zero.mean),
                nobs = ss$T.fin, method = "EM", Phi.f.array = Phi.f.array)
  )
  class(out) <- "tvfit"
  out
}

#' Extract all coefficients from a tvfit model
#'
#' @description
#' Returns all estimated coefficients from a fitted \code{tvfit} object as one
#' consolidated table. Unlike \code{summary.tvfit()}, this method does not print
#' block sections — it merges all parameter groups into a single data frame.
#'
#' @param object A fitted object of class \code{"tvfit"}.
#' @param digits Number of digits to round for display (default 6).
#' @param ... Ignored; for S3 compatibility.
#'
#' @return A data frame with columns:
#' \itemize{
#'   \item \code{parameter} — parameter name
#'   \item \code{estimate} — estimated value (transformed for scalar parameters)
#'   \item \code{std.error} — standard error (if available)
#'   \item \code{z.value} — z-statistic (if applicable)
#'   \item \code{p.value} — p-value (if applicable)
#' }
#' @export
coef.tvfit <- function(object, digits = 6, ...) {
  stopifnot(inherits(object, "tvfit"))
  
  theta <- object$theta %||% object$optim$par
  V     <- object$vcov
  if (is.list(V) && !is.null(V$untransformed)) V <- V$untransformed
  
  map  <- .tvvar_param_index(object)
  npar <- length(map$names)
  
  # Ensure theta length matches
  if (!length(theta) || length(theta) != npar) {
    theta <- theta[seq_len(min(length(theta), npar))]
    if (length(theta) < npar)
      theta <- c(theta, rep(NA_real_, npar - length(theta)))
  }
  
  have_V <- is.matrix(V) && all(dim(V) == c(length(theta), length(theta)))
  se_raw <- if (have_V) sqrt(pmax(diag(V), 0)) else rep(NA_real_, length(theta))
  
  est_nat <- theta
  se_nat  <- se_raw
  
  # Transform head scalars (A, B, phi_r) to (0,1) with delta-method SEs
  idx <- map$idx
  if (length(idx$head)) {
    inv_logit <- stats::plogis
    h <- idx$head
    p <- inv_logit(theta[h])
    est_nat[h] <- p
    if (have_V) {
      grad <- p * (1 - p)
      se_nat[h] <- se_raw[h] * grad
    }
  }
  
  z   <- if (all(!is.na(se_nat))) est_nat / se_nat else rep(NA_real_, length(theta))
  pvl <- if (all(!is.na(z))) 2 * stats::pnorm(-abs(z)) else rep(NA_real_, length(theta))
  
  out <- data.frame(
    parameter = map$names,
    estimate  = round(est_nat, digits),
    std.error = round(se_nat, digits),
    z.value   = round(z, digits),
    p.value   = round(pvl, digits),
    stringsAsFactors = FALSE
  )
  
  out
}
