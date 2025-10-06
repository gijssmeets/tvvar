#' Penalized estimation (ECM) with optional adaptive weights
#'
#' @param data            numeric matrix (T x N)
#' @param p               integer VAR lag order
#' @param r               integer number of factors
#' @param zero_mean       logical; if TRUE, intercepts fixed at 0
#' @param lambda_penalty  numeric scalar; L1 penalty level for Phi.f
#' @param penalty_type    "adaptive" or "regular"
#' @param Phi.f.structure optional free-pattern for Phi^f (3D array [N,N,r], list of r N x N, or N x (N*r) matrix)
#' @param Omega_init      optional N x N PD matrix for initialization (used by your init code if applicable)
#'
#' @return tvvar_fit list (same fields as unpenalized_estimate output)
#' @export
penalized_estimate <- function(data,
                               p,
                               r,
                               zero_mean = TRUE,
                               lambda_penalty = 0.01,
                               penalty_type = c("adaptive", "regular"),
                               Phi.f.structure = NULL,
                               Omega_init = NULL) {
  penalty_type <- match.arg(penalty_type)
  
  # ---- data & dims ----
  Y <- as.matrix(data)
  if (!is.numeric(Y) || length(dim(Y)) != 2L) stop("`data` must be a numeric T x N matrix.")
  T.fin <- nrow(Y); N <- ncol(Y)
  
  Phi.f.str <- array(1, dim = c(N,N,r))
  Phi.f.array <- make.Phi.f(structure = Phi.f.str, lag.order = p)
  Phi.f.as.mat <- matrix(Phi.f.array, nrow = N)      # same flattening as legacy code
  Phi.f.mask   <- (Phi.f.as.mat != 0)
  n_phi_free   <- sum(Phi.f.mask)

  # ---- minimal cfg for legacy code ----
  cfg <- list(
    dim.VAR = N,
    lag.order = p,
    number.factors = r,
    zero.mean = isTRUE(zero_mean),
    T = T.fin,
    T.here = T.fin - p,
    Phi.f.array = Phi.f.array,
    Phi.f.array.mat.structure = Phi.f.as.mat,
    data = Y,
    penalty_type = penalty_type
  )
  if (!is.null(Omega_init)) cfg$omega <- Omega_init
  
  # ---- initial parameters (your existing initializer) ----
  params0 <- initialize_parameters(Y, Phi.f.as.mat, cfg)
  
  # ---- choose weights for Phi.f ----
  if (identical(penalty_type, "adaptive")) {
    # *** ADAPTIVE WEIGHTS FROM UNPENALIZED ML ***
    fit_ml <- unpenalized_estimate(
      data = Y, p = p, r = r, zero_mean = zero_mean,
      phi_f_structure = maxtrix(1, nrow=N,ncol=N), method = "ML"
    )
    # extract free Phi.f entries (by mask order)
    Phi.f_ml_free <- matrix(fit_ml$estimate$Phi_f, nrow = N)[Phi.f.mask]
    w.Phi.f <- (1 / (abs(Phi.f_ml_free) + 1e-6)) / 2
    
    # optional: zero-out diagonal-group weights if that's your convention
    diag_idx <- calc_indices_of_lambda_diagonals(N, p, r)
    if (!is.null(diag_idx$without_zeros)) {
      idx <- diag_idx$without_zeros
      idx <- idx[idx >= 1 & idx <= length(w.Phi.f)]
      if (length(idx)) w.Phi.f[idx] <- 0
    }
  } else {
    # constant weights → plain lasso
    w.Phi.f <- 1
  }
  
  # ---- run penalized ECM (uses your em_algorithm and APG inside Phi.f C-step) ----
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
  
  # ---- evaluate once to collect likelihood etc. ----
  ev <- opti.fct(
    par_free    = fit_pen$par,
    par_fixed   = NaN,
    VAR.data    = Y,
    Phi.f.array = Phi.f.array,
    cfg         = cfg,
    Smooth      = FALSE,
    purpose     = "eval"
  )
  
  # ---- unpack (mirror unpenalized) ----
  inv_logit <- stats::plogis
  params <- fit_pen$par
  idx <- 1L
  A  <- inv_logit(params[idx]); idx <- idx + 1L
  B  <- inv_logit(params[idx]); idx <- idx + 1L
  phi_r <- if (r > 0) { x <- inv_logit(params[idx:(idx + r - 1L)]); idx <- idx + r; x } else numeric(0L)
  
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
  
  need <- if (isTRUE(zero_mean)) N^2 * p else (N^2 * p + N)
  Phi_c_vec <- params[idx:(idx + need - 1L)]
  Phi_c <- if (isTRUE(zero_mean)) cbind(0, matrix(Phi_c_vec, nrow = N)) else matrix(Phi_c_vec, nrow = N)
  
  # ---- ICs (same fields as unpenalized; naive under penalty) ----
  K <- length(c(A, B, phi_r)) + n_phi_free + Nstar + length(Phi_c_vec)
  n_eff <- (T.fin - p) * N
  AIC  <-  2 * K + 2 * ev$average.L * (T.fin - p)
  AICc <-  AIC + (2*K^2 + 2*K) / (max(1, n_eff - K - 1))
  BIC  <-  log(n_eff) * K + 2 * ev$average.L * (T.fin - p)
  
  out <- list(
    estimate = list(A = A, B = B, phi_r = phi_r,
                    Phi_f = Phi.f.est, Phi_c = Phi_c, Omega = Omega),
    lik = list(avg = ev$average.L, sum = ev$full.L),
    ic  = list(AIC = AIC, AICc = AICc, BIC = BIC),
    eval = ev,
    optim = list(par = fit_pen$par, value = fit_pen$obj, convergence = NA,
                 counts = NA, message = sprintf("ECM-penalized (lambda=%.4g, weights=%s)", lambda_penalty, penalty_type)),
    meta = list(N = N, p = p, r = r, zero_mean = isTRUE(zero_mean),
                nobs = T.fin, method = "ECM-penalized",
                lambda = c(0, lambda_penalty), weights = penalty_type)
  )
  class(out) <- "tvvar_fit"
  out
}