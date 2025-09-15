#' Simulate data from a TV-VAR with scalar BEKK innovations
#'
#' @param T Number of time points (after burn-in).
#' @param N Number of variables.
#' @param p Lag order of the VAR.
#' @param r Number of factors.
#' @param A2 Scalar BEKK parameter A.
#' @param B2 Scalar BEKK parameter B.
#' @param pphi AR(1) coefficient for factors (scalar or matrix).
#' @param omega Covariance matrix (N x N).
#' @param Phi.c Static VAR coefficients (N x (N*p)).
#' @param Phi.f Factor loadings matrix (N x N).
#' @param burn_in Number of burn-in steps (default 500).
#' @return List with simulated data matrix Y (T x N), latent factors, shocks, and Phi.f structure.
#' @export
simulate_tvvar_data <- function(T = 200, N = 2, p = 1, r = 1,
                                A2 = 0.1, B2 = 0.75,
                                pphi = 0.95,
                                omega = matrix(c(0.3, 0.2, 0.2, 0.3), 2),
                                Phi.c = matrix(c(0.3, 0.15, 0, 0.3), 2),
                                Phi.f = matrix(c(0.2, 0, 0, 0.2), 2),
                                burn_in = 500) {
  T.burn <- T + burn_in
  Y <- matrix(0, nrow = T.burn, ncol = N)
  factors <- matrix(0, nrow = T.burn, ncol = r)
  u <- matrix(0, nrow = T.burn, ncol = N)
  
  # innovation variance for factors
  if (r == 1) {
    eta <- 1 - pphi^2
  } else {
    eta <- diag(r) - pphi %*% t(pphi)
  }
  
  H.array <- array(0, dim = c(N, N, T.burn))
  H.array[, , 1] <- omega / (1 - A2 - B2)
  u[1, ] <- t(chol(H.array[,,1]) %*% rnorm(N))
  
  start_idx <- p + 1
  
  for (i in start_idx:T.burn) {
    # factor AR update
    innov <- if (r == 1) {
      rnorm(1, sd = sqrt(eta))
    } else {
      MASS::mvrnorm(1, mu = rep(0, r), Sigma = eta)
    }
    factors[i, ] <- pphi %*% factors[i-1, ] + innov
    
    # time-varying coefficient
    Phi.t <- Phi.c + Phi.f * factors[i, ]
    
    # BEKK recursion
    H.array[,,i] <- omega + B2 * H.array[,,i-1] + A2 * u[i-1, ] %*% t(u[i-1, ])
    
    # new shock
    u[i, ] <- t(chol(H.array[,,i])) %*% rnorm(N)
    
    # VAR update
    if (i > p) {
      y_lags <- as.numeric(t(Y[(i - 1):(i - p), , drop = FALSE]))
      Y[i, ] <- Phi.t %*% matrix(y_lags, ncol = 1) + u[i, ]
    }
  }
  
  # drop burn-in
  Y <- Y[(burn_in+1):T.burn, , drop = FALSE]
  u <- u[(burn_in+1):T.burn, , drop = FALSE]
  factors <- factors[(burn_in+1):T.burn, , drop = FALSE]
  
  list(
    Y = Y,
    factors = factors,
    shocks = u
    )
}