#' Penalized time-varying VAR estimation (ECM)
#'
#' Fits a penalized time-varying VAR (TV-VAR) model using the ECM algorithm.
#' Optionally uses adaptive weights for shrinkage. The estimated object can be
#' passed to functions such as `tvirf()` and `tvpred()` for IRFs and forecasting.
#'
#' @param data Numeric matrix (T x N)
#' @param p Integer VAR lag order
#' @param r Integer number of factors
#' @param zero.mean Logical; if TRUE, intercepts fixed at 0
#' @param lambda_penalty Numeric scalar; L1 penalty level for Phi.f
#' @param penalty_type "adaptive" or "regular"
#' @param factor.structure Structure specification for Φᶠ (the factor loading component).
#'
#'   Defines which autoregressive coefficients are allowed to vary with each
#'   latent factor and which are fixed to zero. Each entry of
#'   \code{factor.structure} indicates whether a coefficient is estimated (\code{1})
#'   or fixed to zero (\code{0}). The input can be:
#'   \itemize{
#'     \item a single matrix of dimension N×N (applied to all factors), or
#'     \item a list or 3D array of r matrices, one per factor.
#'   }
#'
#'   All coefficients corresponding to entries equal to \code{1} are *eligible*
#'   for penalization, depending on the chosen penalty type and strength.
#'   Coefficients fixed to \code{0} are excluded from estimation entirely.
#'
#'   The factors must be ordered such that the identification condition holds:
#'   the first r×r block of the factor-loading matrix (Λ₀) is lower-triangular
#'   with strictly positive diagonal elements. These diagonal elements are
#'   internally forced to remain nonzero and *not penalized* to ensure identification.
#'
#'   In practice, this means that the first variable primarily loads on the first
#'   factor, the second on the second factor, and so on. If the user provides a
#'   custom factor structure, they should ensure this ordering and nonzero
#'   diagonal pattern before estimation.
#'
#'   If no structure is provided, \code{tvpenfit()} uses the default
#'   \code{\link{init.factor.structure}()} pattern, which satisfies these
#'   identification requirements.
#'
#' @details
#'   In the penalized model, identification is maintained by fixing the diagonal
#'   of Λ₀ as nonzero and unpenalized, while all remaining free coefficients
#'   (entries equal to 1 in \code{factor.structure}) can be subject to L1 or
#'   adaptive L1 regularization. The user must ensure factor ordering consistent
#'   with this restriction.
#'
#' @seealso
#'   \code{\link{tvfit}} for the unpenalized version, which automatically checks
#'   whether a valid lower-triangular identification can be achieved up to
#'   permutation.
#' @param init "default", "random", or "custom"
#' @param init_list Named list for custom init (A, B, phi_r, Omega, Phi_c, Phi_f)
#'
#' @return Object of class `tvfit`.
#' @seealso [tvfit()], [tvirf()], [tvpred()]
#' @export
tvpenfit <- function(data,
                     p = 1,
                     r = 1,
                     zero.mean = TRUE,
                     lambda_penalty = 0.01,
                     penalty_type = c("adaptive", "regular"),
                     factor.structure = NULL,
                     init = c("default", "random", "custom"),
                     init_list = NULL) {
  start_time   <- Sys.time()
  penalty_type <- match.arg(penalty_type)
  init         <- match.arg(init)
  
  ## --- check id factor.structure ----
  if (!is.null(factor.structure)) {
    ok_id <- check_identification(factor.structure, dim.VAR = ncol(Y),
                                  number.factors = r, lag.order = p)
    if (!ok_id) stop("The specified factor.structure is not identified. 
                   Ensure Λ₀ is lower-triangular with positive diagonal.")
  }
  
  ## ---- data & dims ----
  Y <- as.matrix(data)
  if (!is.numeric(Y) || length(dim(Y)) != 2L)
    stop("`data` must be a numeric T x N matrix.")
  T.fin <- nrow(Y); N <- ncol(Y)
  
  ## ---- structure setup (aligned with unpen ML/EM) ----
  # default structure: all ones per factor if nothing provided
  base_struct <- if (is.null(factor.structure)) matrix(1, N, N) else factor.structure
  phi_arr <- .normalize_phi_structure(base_struct, N = N, r = r)
  Phi.f.array <- make.Phi.f(structure = phi_arr, lag.order = p)
  
  # mask for free Phi.f entries (legacy flattening)
  Phi.f.as.mat <- matrix(Phi.f.array, nrow = N)
  Phi.f.mask   <- (Phi.f.as.mat != 0)
  n_phi_free   <- sum(Phi.f.mask)
  
  ## ---- cfg (consistent with unpen ML/EM) ----
  cfg <- list(
    dim.VAR = N,
    lag.order = p,
    number.factors = r,
    zero.mean = isTRUE(zero.mean),
    T = T.fin,
    T.here = T.fin - p,
    Phi.f.array = Phi.f.array,
    Phi.f.array.mat.structure = Phi.f.as.mat,
    data = Y,
    penalty_type = penalty_type
  )
  
  ## ---- initialization (shared) ----
  if (!is.null(init_list) && !is.list(init_list)) init_list <- as.list(init_list)
  init_res <- .tvvar_build_par_init(
    Y = Y, p = p, r = r, zero.mean = zero.mean,
    Phi.f.array = Phi.f.array,
    init = init, init_list = init_list
  )
  params0     <- init_res$par_init
  n_phi_free0 <- init_res$n_phi_free  # for sanity (should equal n_phi_free)
  
  ## ---- adaptive or regular weights for Phi.f ----
  if (identical(penalty_type, "adaptive")) {
    fit_ml <- tvfit(
      data = Y, p = p, r = r, zero.mean = zero.mean,
      factor.structure = matrix(1, nrow = N, ncol = N),
      method = "ML", init = "default"
    )
    # extract ML Phi.f on the free pattern
    Phi.f_ml_free <- matrix(fit_ml$estimate$Phi_f, nrow = N)[Phi.f.mask]
    w.Phi.f <- (1 / (abs(Phi.f_ml_free) + 1e-6)) / 2
  } else {
    w.Phi.f <- 1
  }
  
  # do not penalize the diagonals of Lambda for identification purposes (see theory thesis)
  diag_idx <- calc_indices_of_lambda_diagonals(N, p, r)
  idx <- diag_idx$without_zeros
  idx <- idx[idx >= 1 & idx <= length(w.Phi.f)]
  w.Phi.f[idx] <- 0
  
  ## ---- run penalized ECM (your routine) ----
  fit_pen <- em_algorithm(
    params      = params0,
    data        = Y,
    Phi.f.array = Phi.f.array,
    cfg         = cfg,
    conditional = TRUE,
    lambda      = c(0, lambda_penalty),  # only penalize Phi.f
    weights.c   = 1,
    weights.f   = w.Phi.f
  )
  
  ## ---- one evaluation for likelihoods ----
  ev <- opti.fct(
    par_free    = fit_pen$par,
    par_fixed   = NaN,
    VAR.data    = Y,
    Phi.f.array = Phi.f.array,
    cfg         = cfg,
    Smooth      = TRUE,
    purpose     = "eval"
  )
  
  ## ---- unpack estimates (mirror unpenalized) ----
  inv_logit <- stats::plogis
  params <- fit_pen$par
  idx <- 1L
  A  <- inv_logit(params[idx]); idx <- idx + 1L
  B  <- inv_logit(params[idx]); idx <- idx + 1L
  phi_r <- if (r > 0) inv_logit(params[idx:(idx + r - 1L)]) else numeric(0L)
  idx <- idx + length(phi_r)
  
  Phi.f.est.mat <- matrix(0, nrow = nrow(Phi.f.as.mat), ncol = ncol(Phi.f.as.mat))
  if (n_phi_free > 0) {
    Phi.f.est.mat[Phi.f.mask] <- params[idx:(idx + n_phi_free - 1L)]
    idx <- idx + n_phi_free
  }
  Phi.f.est <- array(Phi.f.est.mat, dim = dim(Phi.f.array))
  
  Nstar <- N * (N + 1L) / 2L
  L_vec <- params[idx:(idx + Nstar - 1L)]; idx <- idx + Nstar
  Lmat  <- matrix(D.matrix(N) %*% L_vec, ncol = N)
  Lmat[upper.tri(Lmat)] <- 0
  Omega <- Lmat %*% t(Lmat)
  
  need <- if (isTRUE(zero.mean)) N^2 * p else (N^2 * p + N)
  Phi_c_vec <- params[idx:(idx + need - 1L)]
  Phi_c <- if (isTRUE(zero.mean))
    cbind(0, matrix(Phi_c_vec, nrow = N)) else matrix(Phi_c_vec, nrow = N)
  
  ## ---- ICs (naive under penalty) ----
  K <- length(c(A, B, phi_r)) + n_phi_free + Nstar + length(Phi_c_vec)
  n_eff <- (T.fin - p) * N
  AIC  <-  2 * K + 2 * ev$average.L * (T.fin - p)
  AICc <-  AIC + (2 * K^2 + 2 * K) / (max(1, n_eff - K - 1))
  BIC  <-  log(n_eff) * K + 2 * ev$average.L * (T.fin - p)
  
  end_time <- Sys.time()
  runtime_sec <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  ## ---- output (same structure as ML/EM) ----
  out <- list(
    estimate = list(A = A, B = B, phi_r = phi_r,
                    Phi_f = Phi.f.est, Phi_c = Phi_c, Omega = Omega),
    lik = list(avg = ev$average.L, sum = ev$full.L),
    ic  = list(AIC = AIC, AICc = AICc, BIC = BIC),
    eval = ev,
    optim = list(par = params, value = fit_pen$obj, convergence = NA,
                 counts = NA, message = sprintf("ECM-penalized (lambda=%.4g, %s)", lambda_penalty, penalty_type)),
    vcov  = matrix(NA_real_, length(params), length(params)),  # placeholder under L1
    theta = params,
    meta  = list(N = N, p = p, r = r, zero.mean = isTRUE(zero.mean),
                 nobs = T.fin, method = "penalized",
                 Phi.f.array = Phi.f.array,
                 lambda = lambda_penalty, weights = penalty_type,
                 init = init, time = runtime_sec, data = data)
  )
  class(out) <- "tvfit"
  out
}