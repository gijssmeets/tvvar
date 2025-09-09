#' @keywords internal
#' @noRd
Q.matrix.fun <- function(r, N, p) {
  aid.vec <- matrix(1:r, ncol = 1)
  aid.mat <- matrix(aid.vec %x% diag(N * p + 1))
  Qmat <- matrix(0, nrow = r * (N * p + 1)^2, ncol = r)
  for (i in 1:r) {
    Qmat[which(aid.mat == i), i] <- 1
  }
  Qmat
}

#' @keywords internal
#' @noRd
Phi.f <- function(structure, lag.order) {
  dim.VAR <- dim(structure)[1]
  number.factors <- dim(structure)[3]
  out <- array(0, c(dim.VAR, (dim.VAR * lag.order + 1), number.factors))
  for (j in 1:number.factors) {
    out[, , j] <- cbind(rep(0, dim.VAR),
                        matrix(rep(structure[, , j], lag.order), nrow = dim.VAR))
  }
  out
}

#' Build state-space inputs (internal)
#' @keywords internal
#' @noRd
build.statespace.form <- function(VAR.data, lag.order, number.factors, Phi.f, Phi.c) {
  dim.VAR <- ncol(VAR.data)
  T.fin <- nrow(VAR.data)
  Y <- t(VAR.data)
  
  # C++ helper builds lagged regressors
  Y.minus1 <- create_Y_minus1(Y, lag.order, T.fin, dim.VAR)
  
  # drop initial lags from Y to align
  Y <- Y[, -(1:lag.order), drop = FALSE]
  
  # dependent variable
  ytilde <- Y - Phi.c %*% Y.minus1
  
  # selection and Z array (C++)
  Q <- Q.matrix.fun(r = number.factors, N = dim.VAR, p = lag.order)
  Z <- create_Z(Y.minus1, Phi.f, Q, T.fin, dim.VAR, number.factors, lag.order)
  
  list(Z = Z, ytilde = ytilde)
}







#' Unpenalized estimation (ML) for TV-VAR — minimal wrapper
#' @param data numeric matrix (T x N), already transformed (demeaned if zero_mean=TRUE).
#' @param p integer, VAR lag order.
#' @param r integer, number of factors.
#' @param zero_mean logical; TRUE fixes intercept to 0.
#' @param phi_f_structure 3D array [N,N,r], or list of r (N x N) matrices, or block matrix N x (N*r).
#' @param control list for optim (e.g., list(maxit=2000, trace=0)).
#' @return list with estimates, likelihood, ICs, eval, optim info.
#' @export
unpenalized_estimate <- function(data, p, r, zero_mean = TRUE,
                                 phi_f_structure,
                                 control = list(maxit = 2000, trace = 0)) {
  stopifnot(is.matrix(data) || is.data.frame(data))
  Y <- as.matrix(data)
  n_var <- ncol(Y); T.fin <- nrow(Y)
  
  # normalize phi_f_structure -> 3D array [N,N,r]
  phi_arr <- .normalize_phi_structure(phi_f_structure, n_var, r)
  
  # build Phi^f array [N, N*p+1, r] and mask
  Phi.f.array <- Phi.f(structure = phi_arr, lag.order = p)
  Phi.f.mat  <- matrix(Phi.f.array, nrow = n_var)
  Phi.f.mask <- matrix(Phi.f.array, nrow = n_var); Phi.f.mask[Phi.f.mask != 0] <- 1L
  n_phi_free <- sum(Phi.f.mask != 0)
  
  # quick starts: (A,B,phi_r) in (0,1), Omega small diag, Phi_c by OLS
  logit <- function(x) log(x/(1 - x))
  A0 <- 0.10; B0 <- 0.80; phi_r0 <- rep(0.95, r)
  Omega0 <- diag(0.1, n_var)
  
  # Phi_c by static VAR (fast OLS)
  Ydep <- Y[(p + 1):T.fin, , drop = FALSE]
  Xlag <- do.call(cbind, lapply(1:p, function(L) Y[(p + 1 - L):(T.fin - L), , drop = FALSE]))
  X <- cbind(1, Xlag)
  Bhat <- vapply(seq_len(n_var), function(j) lm.fit(X, Ydep[, j])$coefficients, numeric(ncol(X)))
  Phi_c0 <- t(Bhat); if (isTRUE(zero_mean)) Phi_c0[, 1] <- 0
  
  # pack params (unconstrained)
  params0 <- c(logit(A0), logit(B0), logit(phi_r0),
               rep(0.01, n_phi_free),
               vech(Omega0),
               as.vector(Phi_c0))
  
  # ML optimize
  opt <- optim(
    par = params0,
    fn  = opti_fct,
    purpose   = "optim",
    VAR.data  = Y,
    Phi.f.mask = Phi.f.mask,
    zero.mean  = isTRUE(zero_mean),
    Smooth     = FALSE,
    r = r, p = p, n_var = n_var,
    method = "BFGS",
    hessian = FALSE,
    control = control
  )
  
  # evaluate at optimum
  ev <- opti_fct(
    params   = opt$par,
    VAR.data = Y,
    Phi.f.mask = Phi.f.mask,
    zero.mean  = isTRUE(zero_mean),
    Smooth     = FALSE,
    purpose = "eval",
    r = r, p = p, n_var = n_var
  )
  
  # transform back (mirror slicing in opti_fct)
  inv_logit <- stats::plogis
  idx <- 1L
  A  <- inv_logit(opt$par[idx]); idx <- idx + 1L
  B  <- inv_logit(opt$par[idx]); idx <- idx + 1L
  phi_r <- inv_logit(opt$par[idx:(idx + r - 1L)]); idx <- idx + r
  
  Phi.f.est <- matrix(0, nrow = nrow(Phi.f.mat), ncol = ncol(Phi.f.mat))
  if (n_phi_free > 0) {
    Phi.f.est[Phi.f.mask != 0] <- opt$par[idx:(idx + n_phi_free - 1L)]
    idx <- idx + n_phi_free
  }
  Phi.f.est <- array(Phi.f.est, dim = dim(Phi.f.array))
  
  Nstar <- n_var * (n_var + 1L) / 2L
  L_vec <- opt$par[idx:(idx + Nstar - 1L)]; idx <- idx + Nstar
  Lmat  <- matrix(D.matrix(n_var) %*% L_vec, ncol = n_var)
  Lmat[upper.tri(Lmat)] <- 0
  Omega <- Lmat %*% t(Lmat)
  
  need <- if (isTRUE(zero_mean)) n_var^2 * p else (n_var^2 * p + n_var)
  Phi_c_vec <- opt$par[idx:(idx + need - 1L)]
  Phi_c <- if (isTRUE(zero_mean)) cbind(0, matrix(Phi_c_vec, nrow = n_var)) else matrix(Phi_c_vec, nrow = n_var)
  
  # ICs
  K <- length(c(A, B, phi_r)) + n_phi_free + Nstar + length(Phi_c_vec)
  n_eff <- (T.fin - p) * n_var
  AIC  <-  2 * K + 2 * ev$average.L * (T.fin - p)
  AICc <-  AIC + (2*K^2 + 2*K) / (n_eff - K - 1)
  BIC  <-  log(n_eff) * K + 2 * ev$average.L * (T.fin - p)
  
  out <- list(
    estimate = list(A = A, B = B, phi_r = phi_r,
                    Phi_f = Phi.f.est, Phi_c = Phi_c, Omega = Omega),
    lik = list(avg = ev$average.L, sum = ev$full.L),
    ic  = list(AIC = AIC, AICc = AICc, BIC = BIC),
    eval = ev,
    optim = list(par = opt$par, value = opt$value, convergence = opt$convergence,
                 counts = opt$counts, message = opt$message),
    meta = list(N = n_var, p = p, r = r, zero_mean = isTRUE(zero_mean), nobs = T.fin, method = "ML")
  )
  class(out) <- "tvvar_fit"
  out
}

# internal: accept array/list/block matrix; enforce r if provided
#' @keywords internal
#' @noRd
.normalize_phi_structure <- function(phi, n_var, r = NULL) {
  if (is.array(phi) && length(dim(phi)) == 3L) {
    stopifnot(dim(phi)[1L] == n_var, dim(phi)[2L] == n_var)
    if (!is.null(r)) stopifnot(dim(phi)[3L] == r)
    return(phi)
  }
  if (is.list(phi)) {
    r0 <- length(phi); if (!is.null(r)) stopifnot(r0 == r) else r <- r0
    stopifnot(all(vapply(phi, function(m) all(dim(m) == c(n_var, n_var)), logical(1))))
    arr <- array(0, c(n_var, n_var, r))
    for (j in seq_len(r)) arr[, , j] <- phi[[j]]
    return(arr)
  }
  if (is.matrix(phi)) {
    stopifnot(nrow(phi) == n_var, ncol(phi) %% n_var == 0L)
    r0 <- ncol(phi) %/% n_var; if (!is.null(r)) stopifnot(r0 == r) else r <- r0
    arr <- array(0, c(n_var, n_var, r))
    for (j in seq_len(r)) {
      cols <- ((j - 1L) * n_var + 1L):(j * n_var)
      arr[, , j] <- phi[, cols, drop = FALSE]
    }
    return(arr)
  }
  stop("`phi_f_structure` must be a 3D array [N,N,r], a list of N x N matrices, or an N x (N*r) block matrix.")
}






#' Core optimizer for TV-VAR likelihood / evaluation
#'
#' Runs the Kalman filter (and optional smoother) and returns the
#' negative average log-likelihood (for optimization) or detailed
#' state outputs (for evaluation).
#'
#' @details Parameters are unconstrained and mapped internally
#'   (e.g., via `plogis`). Nonzero entries of `Phi.f` are filled
#'   using `Phi.f.mask`.
#'
#' @param params Numeric vector of unconstrained parameters (see Details).
#' @param VAR.data Matrix (T x N).
#' @param Phi.f.mask Matrix or array used to indicate nonzero entries of Phi.f
#'   (same shape as Phi.f when reshaped); values != 0 are filled from `params`.
#' @param zero.mean Logical; if TRUE, intercept is fixed at 0.
#' @param Smooth Logical; if TRUE, run Kalman smoother.
#' @param purpose "optim" (return negative avg loglik) or "eval" (return list).
#' @param r Integer, number of factors.
#' @param p Integer, VAR lag order.
#' @param n_var Integer, number of VAR variables (ncol(VAR.data)).
#'
#' @return If purpose="optim": scalar negative average loglik.
#'   If purpose="eval": list with likelihoods, states, and components.
#' @export
#' @importFrom stats plogis


opti_fct <- function(params, VAR.data, Phi.f.mask,
                     zero.mean = TRUE, Smooth = FALSE,
                     purpose = c("optim", "eval"),
                     r, p, n_var) {
  
  purpose <- match.arg(purpose)
  stopifnot(is.matrix(VAR.data), is.numeric(params), length(r) == 1, length(p) == 1, length(n_var) == 1)
  
  inv_logit <- function(x) plogis(x)
  
  ## ---- 1) unpack parameters -------------------------------------------------
  idx <- 1L
  A <- inv_logit(params[idx]); idx <- idx + 1L
  B <- inv_logit(params[idx]); idx <- idx + 1L
  
  # phi_r (pphi) on (0,1); diagonal unless r==1 (1x1)
  if (r > 0) {
    phi_r <- inv_logit(params[idx:(idx + r - 1L)]); idx <- idx + r
    pphi <- if (r == 1L) matrix(phi_r, 1, 1) else diag(phi_r, nrow = r)
    Q <- diag(r) - pphi %*% t(pphi)
    starta <- numeric(r); startP <- diag(r)
  } else {
    # r == 0: degenerate factor (keep 1x1 placeholders but we won't smooth)
    pphi <- matrix(0, 1, 1)
    Q <- diag(1) - pphi %*% t(pphi)
    starta <- numeric(1); startP <- diag(1)
  }
  
  # Fill Phi.f from masked positions
  Phi.f.as.mat <- matrix(Phi.f.mask, nrow = n_var)
  k_phi <- sum(Phi.f.as.mat != 0)
  stopifnot(idx + k_phi - 1L <= length(params))
  Phi.f <- matrix(0, nrow = nrow(Phi.f.as.mat), ncol = ncol(Phi.f.as.mat))
  Phi.f[Phi.f.as.mat != 0] <- params[idx:(idx + k_phi - 1L)]
  idx <- idx + k_phi
  
  # Omega (momega) via lower-tri param vector (N*(N+1)/2)
  Nstar <- n_var * (n_var + 1L) / 2L
  stopifnot(idx + Nstar - 1L <= length(params))
  # D.matrix should map vech to lower-tri; use your Rcpp D.matrix or R equivalent
  L_vec <- params[idx:(idx + Nstar - 1L)]; idx <- idx + Nstar
  m1 <- matrix(D.matrix(n_var) %*% L_vec, ncol = n_var)  # D.matrix provided by C++ or R
  m1[upper.tri(m1)] <- 0
  momega <- m1 %*% t(m1)
  firstH <- momega / (1 - A - B)
  covi <- 1
  
  # Phi.c (deterministic part)
  if (zero.mean) {
    need <- n_var^2 * p
    stopifnot(idx + need - 1L <= length(params))
    Phi.c <- cbind(rep(0, n_var),
                   matrix(params[idx:(idx + need - 1L)], nrow = n_var))
    idx <- idx + need
  } else {
    need <- n_var^2 * p + n_var
    stopifnot(idx + need - 1L <= length(params))
    Phi.c <- matrix(params[idx:(idx + need - 1L)], nrow = n_var)
    idx <- idx + need
  }
  
  ## ---- 2) build state-space inputs -----------------------------------------
  nf <- if (r == 0L) 1L else r
  ss <- build.statespace.form(VAR.data = VAR.data, lag.order = p,
                              number.factors = nf, Phi.f = Phi.f, Phi.c = Phi.c)
  ytilde <- ss$ytilde
  Z <- ss$Z
  T.here <- ncol(ytilde) - 1L
  
  ## ---- 3) main filter / likelihood (C++ backend) ----------------------------
  out <- my_loop_main(ytilde, Z, startP, covi, firstH, momega, pphi, A, B, Q)
  
  if (purpose == "optim") {
    return(-out$avlik)
  }
  
  ## ---- 4) (optional) RTS smoother ------------------------------------------
  res <- list(
    average.L = -out$avlik,
    full.L    = -out$sumlik,
    Phi.c     = Phi.c,
    array.filtered.H = out$aH,
    momega = vech(momega),   # assumes you provide vech()
    Nstar  = Nstar
  )
  
  last.obs <- ncol(out$pstate)
  if (r > 0) {
    pred_a <- matrix(out$pstate[, -last.obs, drop = FALSE], ncol = T.here)
    pred_P <- out$pstatevariance[, , -last.obs, drop = FALSE]
    filt_a <- out$fstate[, -last.obs, drop = FALSE]
    filt_P <- out$fstatevariance[, , -last.obs, drop = FALSE]
    L.arr  <- out$L[, , -last.obs, drop = FALSE]
    F.arr  <- out$F[, , -last.obs, drop = FALSE]
    v.mat  <- out$vmat[, -last.obs, drop = FALSE]
    pred.err.var <- F.arr
    
    res$predicted.state            <- pred_a
    res$predicted.state.variance   <- pred_P
    res$filtered.state             <- filt_a
    res$filtered.state.variance    <- filt_P
    res$Z                          <- Z
    res$v.mat                      <- v.mat
    res$pred.err.var               <- pred.err.var
    
    if (isTRUE(Smooth)) {
      rts_r <- matrix(0, nrow = r, ncol = T.here)
      rts_N <- array(0, dim = c(r, r, T.here))
      alpha_hat <- matrix(0, nrow = nrow(pred_a), ncol = ncol(pred_a))
      V_arr <- array(0, dim = c(r, r, T.here))
      smooth_err <- matrix(0, nrow = nrow(ytilde), ncol = T.here)
      
      for (t in T.here:1) {
        F_inv <- solve(F.arr[, , t])
        L_t   <- L.arr[, , t]
        v_t   <- v.mat[, t, drop = FALSE]
        a_t   <- pred_a[, t, drop = FALSE]
        P_t   <- pred_P[, , t, drop = FALSE]
        Z_t   <- Z[, , t, drop = FALSE]
        
        r_tm1 <- t(Z_t) %*% F_inv %*% v_t + t(L_t) %*% rts_r[, t, drop = FALSE]
        N_tm1 <- t(Z_t) %*% F_inv %*% Z_t + t(L_t) %*% rts_N[, , t] %*% L_t
        
        alpha_hat[, t] <- a_t + P_t %*% r_tm1
        V_arr[, , t]   <- P_t - P_t %*% N_tm1 %*% P_t
        
        if (t > 1) {
          rts_r[, t - 1] <- r_tm1
          rts_N[, , t - 1] <- N_tm1
        }
        smooth_err[, t] <- ytilde[, t] - Z_t %*% alpha_hat[, t, drop = FALSE]
      }
      
      res$smoothed.state            <- alpha_hat
      res$smoothed.state.variance   <- V_arr
      res$smooth.error              <- smooth_err
    }
  }
  
  res
}