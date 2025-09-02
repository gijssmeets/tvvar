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
  # D_matrix should map vech to lower-tri; use your Rcpp D_matrix or R equivalent
  L_vec <- params[idx:(idx + Nstar - 1L)]; idx <- idx + Nstar
  m1 <- matrix(D_matrix(n_var) %*% L_vec, ncol = n_var)  # D_matrix provided by C++ or R
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