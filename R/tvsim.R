#' Simulate data from a TV-VAR with scalar BEKK innovations
#'
#' @param T Number of time points (after burn-in).
#' @param N Number of variables.
#' @param p Lag order of the VAR.
#' @param r Number of factors.
#' @param A2 Scalar BEKK parameter A.
#' @param B2 Scalar BEKK parameter B.
#' @param pphi If r=1: scalar. If r>1: (rxr) matrix with pphi values on the diagonal.
#' @param omega Covariance matrix (N x N).
#' @param Phi.c Static VAR coefficients (N x (N*p)).
#' @param Phi.f Factor loadings matrix (N x (N*pr)).
#' @param burn_in Number of burn-in steps (default 500).
#' @return List with simulated data matrix Y (T x N), latent factors, shocks, and Phi.f structure.
#' @export
tvsim <- function(T = 200, N = 2, p = 1, r = 1,
                                A2 = 0.1, B2 = 0.75,
                                pphi = 0.95,
                                omega = matrix(c(0.3, 0.2, 0.2, 0.3), 2),
                                Phi.c = matrix(c(0.3, 0.15, 0, 0.3), 2),
                                Phi.f = matrix(c(0.2, 0, 0, 0.2), 2),
                                burn_in = 500) {
  stopifnot(is.matrix(omega), nrow(omega) == N, ncol(omega) == N)
  stopifnot(is.matrix(Phi.c), nrow(Phi.c) == N, ncol(Phi.c) == N * p)
  
  ## ---- coerce Phi.f to [N, N*p, r] ----
  if (length(dim(Phi.f)) == 3L) {
    stopifnot(dim(Phi.f)[1] == N, dim(Phi.f)[2] == N * p, dim(Phi.f)[3] == r)
    Phi_f_arr <- Phi.f
  } else if (is.matrix(Phi.f)) {
    if (ncol(Phi.f) == N * p && r == 1L) {
      Phi_f_arr <- array(Phi.f, dim = c(N, N * p, 1L))
    } else if (ncol(Phi.f) == N * p * r) {
      Phi_f_arr <- array(0, dim = c(N, N * p, r))
      for (k in seq_len(r)) {
        cols <- ((k - 1L) * N * p + 1L):(k * N * p)
        Phi_f_arr[,,k] <- Phi.f[, cols, drop = FALSE]
      }
    } else {
      stop("Phi.f must be N x (N*p) (when r=1) or N x (N*p*r) (when r>1).")
    }
  } else {
    stop("Phi.f must be a matrix or a 3D array.")
  }
  
  ## ---- set up factor AR(1) ----
  # pphi can be: scalar (r=1), length-r vector (diagonal), or r x r matrix
  if (r == 1L) {
    phi_mat <- matrix(as.numeric(pphi), 1, 1)
  } else if (is.matrix(pphi)) {
    stopifnot(all(dim(pphi) == c(r, r)))
    phi_mat <- pphi
  } else if (is.numeric(pphi) && length(pphi) == r) {
    phi_mat <- diag(pphi, r)
  } else {
    stop("`pphi` must be scalar (r=1), length-r numeric (diagonal), or r x r matrix.")
  }
  # innovation covariance so that Var(f_t) = I_r if |eig(phi)|<1:
  eta <- diag(r) - phi_mat %*% t(phi_mat)
  # small jitter if near-singular
  eta <- (eta + t(eta)) / 2
  if (min(eigen(eta, symmetric = TRUE, only.values = TRUE)$values) <= 1e-10) {
    eta <- eta + diag(1e-6, r)
  }
  
  ## ---- storage ----
  T.burn <- T + burn_in
  Y       <- matrix(0, T.burn, N)
  factors <- matrix(0, T.burn, r)
  u       <- matrix(0, T.burn, N)
  H.array <- array(0, dim = c(N, N, T.burn))
  
  ## ---- initialize H and first shock ----
  H.array[,,1] <- omega / (1 - A2 - B2)
  u[1, ] <- t(chol(H.array[,,1]) %*% rnorm(N))
  
  start_idx <- p + 1L
  
  for (i in start_idx:T.burn) {
    ## factor AR update: f_t = phi f_{t-1} + e_t, e_t ~ N(0, eta)
    if (r == 1L) {
      innov <- rnorm(1, sd = sqrt(eta[1, 1]))
    } else {
      innov <- MASS::mvrnorm(1, mu = rep(0, r), Sigma = eta)
    }

    factors[i, ] <- (phi_mat %*% factors[i - 1,]) + innov
    
    ## time-varying coefficient: Phi_t = Phi_c + sum_k Phi_f^{(k)} * f_{t,k}
    Phi_t <- Phi.c
    for (k in seq_len(r)) {
      Phi_t <- Phi_t + Phi_f_arr[,,k] * factors[i, k]
    }
    
    ## BEKK(1,1) recursion (scalar A2, B2)
    H.array[,,i] <- omega + B2 * H.array[,,i - 1] + A2 * (u[i - 1, ] %*% t(u[i - 1, ]))
    
    ## new structural shock u_t ~ N(0, H_t)
    u[i, ] <- t(chol(H.array[,,i])) %*% rnorm(N)
    
    ## VAR update with p lags (stacked as y_{t-1}, y_{t-2}, ..., y_{t-p})
    y_lags <- as.numeric(t(Y[(i - 1):(i - p), , drop = FALSE]))
    Y[i, ] <- Phi_t %*% matrix(y_lags, ncol = 1) + u[i, ]
  }
  
  ## ---- drop burn-in ----
  keep <- (burn_in + 1):T.burn
  list(
    Y       = Y[keep, , drop = FALSE],
    factors = factors[keep, , drop = FALSE],
    shocks  = u[keep, , drop = FALSE]
  )
}