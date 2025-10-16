
#' IRFs for penalized tvvar_fit (factor-uncertainty only)
#'
#' @param fit tvvar_fit from penalized_estimate() or unpenalized_estimate()
#' @param lags horizon (number of periods ahead)
#' @param T.max number of time points (from start) at which to compute IRFs
#' @param B Monte Carlo draws (factor simulations only)
#' @param shock_index index of shocked variable (1..N) [for convenience when using ek as unit vector]
#' @param shock_position numeric length-N shock vector (if given, overrides shock_index)
#' @param seed base RNG seed
#' @return list(IRF_lb, IRF_med, IRF_ub)
#' @export
irf_penalized <- function(fit,
                          lags = 10,
                          T.max = 5,
                          B = 500,
                          shock_index = 1,
                          shock_position = NULL,
                          seed = 1234) {
  stopifnot(inherits(fit, "tvvar_fit"))
  N <- fit$meta$N
  p <- fit$meta$p
  r <- fit$meta$r
  zero.mean <- isTRUE(fit$meta$zero.mean)
  
  # shock vector
  ek <- if (is.null(shock_position)) {
    v <- rep(0, N); v[shock_index] <- 1; v
  } else {
    as.numeric(shock_position)
  }
  if (length(ek) != N) stop("shock_position must be length N.")
  
  # point estimates (fixed across draws)
  Phi_c_est <- fit$estimate$Phi_c                       # N x (1+N*p)  (if zero.mean) or N x (N*p+1)
  Phi_f_est <- fit$estimate$Phi_f                       # N x (1+N*p) x r
  psi        <- as.numeric(fit$estimate$phi_r)          # length r on (0,1)
  if (length(dim(Phi_f_est)) != 3L) {
    stop("fit$estimate$Phi_f must be an array N x (1+N*p) x r.")
  }
  
  # Build Phi cube expected by irf_generator_cpp:
  #   slice 1: Phi^c (with intercept col),
  #   slices 2..(r+1): Phi^f_k (with intercept col)
  Phi <- array(0, dim = c(N, N*p + 1, r + 1))
  Phi[,,1] <- Phi_c_est
  for (k in 1:r) Phi[,,k+1] <- Phi_f_est[,,k]
  
  # Storage for bands
  IRF_lb  <- array(NA_real_, c(N, lags, T.max))
  IRF_ub  <- array(NA_real_, c(N, lags, T.max))
  IRF_med <- array(NA_real_, c(N, lags, T.max))
  
  # Quick stability check ingredients (no intercept)
  Phi_c_noint <- Phi_c_est[, -1, drop = FALSE]
  Phi_f_noint <- Phi_f_est[, -1, , drop = FALSE]  # N x (N*p) x r
  
  # Pull evaluation bundle (filtered states, etc.)
  ev <- fit$eval
  if (is.null(ev)) stop("fit$eval is NULL; call evaluation in your estimator.")
  
  # Loop over the first T.max time points
  for (t in 1:T.max) {
    # filtered state, covariance, and H_t at time t
    aT <- ev$filtered.state[, t, drop = TRUE]                   # r
    PT <- ev$filtered.state.variance[, , t, drop = FALSE]       # r x r
    PT <- matrix(PT, nrow = r, ncol = r)                        # ensure 2D
    HT <- ev$array.filtered.H[, , t, drop = FALSE]              # N x N
    HT <- matrix(HT, nrow = N, ncol = N)
    
    # stability (Lyapunov) — positive ⇒ unstable, we will NA those draws
    lyap <- tryCatch(
      check_lyapunov(Phi_c_noint, Phi_f_noint, psi, seed_ly = 12),
      error = function(e) Inf
    )
    
    # Monte Carlo over factor paths only — use irf_generator_cpp with fixed=FALSE
    # Note: that C++ averages across B; to get quantiles, do B calls with B=1.
    set.seed(seed + t)
    irf_draws <- array(NA_real_, c(N, lags, B))
    for (b in 1:B) {
      one <- irf_generator_cpp(
        Phi, psi, aT, PT, HT, ek,
        lags = lags,
        fixed = FALSE,
        B = 1,                # one factor-path draw per call
        seed_f = seed + t + 1000L * b
      )
      irf_draws[,,b] <- one$irf
    }
    
    # If unstable, drop the draws
    if (is.finite(lyap) && lyap > 0) {
      irf_draws[,] <- NA_real_
    }
    
    # Summarize to bands
    for (i in 1:N) {
      for (h in 1:lags) {
        vals <- irf_draws[i, h, ]
        IRF_lb[i, h, t]  <- stats::quantile(vals, 0.16, na.rm = TRUE)
        IRF_med[i, h, t] <- stats::quantile(vals, 0.50, na.rm = TRUE)
        IRF_ub[i, h, t]  <- stats::quantile(vals, 0.84, na.rm = TRUE)
      }
    }
  }
  
  list(IRF_lb = IRF_lb, IRF_med = IRF_med, IRF_ub = IRF_ub)
}




# impulse_response.R

#' Time-varying impulse responses (basic MC bands)
#'
#' Minimal IRF generator that mirrors the provided script.
#' Expects an unpenalized ML fit where \code{fit$theta$untransformed}
#' and \code{fit$vcov$untransformed} exist. Falls back to \code{fit$optim$par}
#' and a tiny diagonal covariance if missing.
#'
#' @param fit tvvar_fit from unpenalized_estimate(..., method = "ML")
#' @param lags horizon
#' @param T.max number of time points (from start) to compute
#' @param B number of Monte Carlo draws (set B=0 or 1 for plug-in only)
#' @param shock_index index of unit shock (1..N)
#' @param seed RNG seed base
#' @return list(IRF_lb, IRF_med, IRF_ub)
#' @importFrom MASS mvrnorm
#' @export
irf <- function(fit,
                lags = 10,
                T.max = 5,
                B = 500,
                shock_index = 1,
                shock_position = c(0,1)) {
  
  N <- fit$meta$N
  p <- fit$meta$p
  r <- fit$meta$r
  zero.mean <- fit$meta$zero.mean
  Phi.f.array <- fit$meta$Phi.f.array
  Phi.f.array.mat.structure <- matrix(Phi.f.array, nrow = N)
  Phi.f.array.mat.structure[Phi.f.array.mat.structure != 0] <- 1
  number.nonzero.phis <- sum(Phi.f.array.mat.structure != 0)
  
  VAR.data <- fit$meta$data

  IRF_lb=array(NA,c(N,lags, T.max))
  IRF_ub=array(NA,c(N,lags, T.max))
  IRF_med=array(NA,c(N,lags, T.max))

  for(t in 1:T.max){
    
    T.index=t #time point
    theta= fit$theta 
    V=fit$vcov 
    s.h= shock_index # see index
    ek= shock_position # shock position (here: third variable)

    irf = array(0,c(N,lags,B))
    Pi = array(0, c(N*p,N*p,B))
    ev <- numeric()
    lyapunov <- numeric()
    
    for(i in 1:B){
      
      s.h = s.h+1
      index=s.h*t
      set.seed(index)

      par_val = MASS::mvrnorm(n = 1, mu = theta, Sigma = V)
      
      #### obtain Phi,psi,aT,PT,HT evaluated at parameter value par_val  
      opti.eval <- opti.fct(
        par_free    = par_val,
        par_fixed   = NaN,
        VAR.data    = VAR.data,
        Phi.f.array = Phi.f.array,
        cfg         = list(
          zero.mean = zero.mean,
          dim.VAR = N,
          lag.order = p,
          number.factors = r
        ),
        Smooth      = FALSE,
        purpose     = "eval"
      )
      
      phi_est <- par_val[3:(3+r-1)]
      
      Phi_f_est <- Phi.f.array.mat.structure
      Phi_f_est[which(Phi.f.array.mat.structure !=0)]=par_val[(3+r):(3+r+number.nonzero.phis-1)]
      
      Phi_c_est <- opti.eval$Phi.c
      
      phi_est_tr <- exp(par_val[3:(3+r-1)])/(exp(par_val[3:(3+r-1)])+1)
      Phi_c_noint <- Phi_c_est[,-1]
      Phi_f_noint <- array(Phi_f_est, dim=c(N, (N*p+1), r))[,-1,,drop=FALSE]
      seed_ly = 12
    
      lyapunov[i] = check_lyapunov(Phi_c_noint, Phi_f_noint, phi_est_tr, seed_ly)
      
      Phi = array(cbind(Phi_c_est, Phi_f_est), c(N,N*p+1, (r+1)))
      psi=exp(phi_est)/(1+exp(phi_est))
      #psi=phi_est
      
      aT <- opti.eval$filtered.state[,T.index]
      
      # Changed this to make sure it is always a matrix
      PT <- opti.eval$filtered.state.variance[, , T.index, drop = FALSE]
      PT <- matrix(PT, nrow = r, ncol = r)
      
      HT <- opti.eval$array.filtered.H[,,T.index]
      
      irf.gen <- irf_generator_cpp(Phi, psi, aT, PT, HT, ek, lags, fixed = TRUE, B = B, seed_f = 1234)
      
      irf[,,i] = irf.gen$irf
      Pi[,,i] = irf.gen$PI
      ev[i] <- eigen(Pi[,,i])$values[1]
    }
    
    unstable <- which(lyapunov>0)
    irf[,,unstable]=NA
    
    for(i in 1:N){
      for(j in 1:lags){
        IRF_lb[i,j,t] = quantile(irf[i,j,], probs=0.16,na.rm=TRUE)
        IRF_ub[i,j,t] = quantile(irf[i,j,], probs=0.84,na.rm=TRUE)
        IRF_med[i,j,t] = quantile(irf[i,j,], probs=0.5,na.rm=TRUE)
      }
    }
    
  }
  
  return(list(IRF_lb=IRF_lb, IRF_med=IRF_med, IRF_ub=IRF_ub) )
}

#' Plot time-varying impulse responses
#'
#' @description
#' Quick plot for IRF objects returned by \code{impulse_response()}.
#' Draws the median response and a basic 68% band (16–84% quantiles).
#'
#' @param irf A list with components \code{IRF_med}, \code{IRF_lb}, \code{IRF_ub},
#'   each an array of shape \code{[N, lags, T.max]}.
#' @param t Integer time index in \code{1:T.max}. If \code{NULL}, uses the last time point.
#' @param j Integer series index in \code{1:N} to plot (row of the IRF arrays).
#' @param main_label Optional main title. If \code{NULL}, a default title is used.
#'
#' @return Invisibly returns \code{NULL}. Draws a base R plot.
#' @export
#'
#' @examples
#' \dontrun{
#' irf_ml <- impulse_response(simdata$Y, fit_ml, lags = 12, T.max = 5, B = 200)
#' plot_irf(irf_ml, t = 5, j = 1)
#' }
plot_irf <- function(irf, t = NULL, j = 1, main_label = NULL) {
  IRF_med <- irf$IRF_med
  IRF_lb  <- irf$IRF_lb
  IRF_ub  <- irf$IRF_ub
  
  if (is.null(IRF_med)) stop("irf$IRF_med is missing.")
  
  dims  <- dim(IRF_med)     # c(N, lags, T.used)
  N     <- dims[1]
  lags  <- dims[2]
  Tused <- dims[3]
  
  if (is.null(t)) t <- Tused
  if (!(j %in% seq_len(N))) stop(sprintf("`j` must be in 1..%d", N))
  if (!(t %in% seq_len(Tused))) stop(sprintf("`t` must be in 1..%d", Tused))
  
  h <- 0:(lags - 1)
  
  # cross-platform: don't force a new device; just plot
  op <- par(mgp = c(2, 0.8, 0), mar = 0.1 + c(3, 3, 3, 1))
  on.exit(par(op), add = TRUE)
  
  if (is.null(main_label)) main_label <- sprintf("IRF (series %d, t = %d)", j, t)
  
  plot(h, IRF_med[j, , t], type = "l", lwd = 2,
       col  = "steelblue",
       ylim = range(IRF_lb[j, , t], IRF_ub[j, , t], na.rm = TRUE),
       xlab = "Horizon", ylab = "", main = main_label)
  lines(h, IRF_ub[j, , t], lty = 2, col = "gray40", lwd = 2)
  lines(h, IRF_lb[j, , t], lty = 2, col = "gray40", lwd = 2)
  abline(h = 0, lwd = 1)
}
      