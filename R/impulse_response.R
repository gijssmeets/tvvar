#' Time-varying impulse responses (bands via MC)
#'
#' Works for both unpenalized ML fits and penalized fits.
#' For ML: identical to your original implementation.
#' For penalized: factors are simulated, static coefficients fixed.
#'
#' @param fit tvvar_fit
#' @param lags horizon
#' @param T.max number of time points (from start)
#' @param B number of Monte Carlo draws
#' @param shock_index index of unit shock (1..N)
#' @param shock_position numeric length-N shock vector (used by your original code)
#' @return list(IRF_lb, IRF_med, IRF_ub)
#' @importFrom MASS mvrnorm
#' @export
irf <- function(fit,
                lags = 10,
                T.max = 5,
                B = 500,
                shock_index = 1,
                shock_position = c(0,1)) {
  
  # --- Branch on estimator: keep original code for ML ---
  if (identical(fit$meta$method, "ML")) {
    
    ## -----------------------
    ## ORIGINAL ML CODE (verbatim)
    ## -----------------------
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
    
  } else {
    # --- Penalized branch: keep coefficients fixed at estimates; simulate factors only ---
    N <- fit$meta$N
    p <- fit$meta$p
    r <- fit$meta$r
    Phi_c_est <- fit$estimate$Phi_c          # N x (1+N*p)
    Phi_f_est <- fit$estimate$Phi_f          # N x (1+N*p) x r
    psi       <- as.numeric(fit$estimate$phi_r)  # length r in (0,1)
    if (length(dim(Phi_f_est)) != 3L) stop("fit$estimate$Phi_f must be an array [N x (1+N*p) x r].")
    
    # Build Phi cube: slice 1 = Phi^c, slices 2..r+1 = Phi^f_k
    Phi <- array(0, dim = c(N, N*p + 1, r + 1))
    Phi[,,1] <- Phi_c_est
    for (k in seq_len(r)) Phi[,,k+1] <- Phi_f_est[,,k]
    
    # shock vector (keep your original interface)
    ek <- shock_position
    if (length(ek) != N) stop("`shock_position` must be length N for penalized IRF.")
    
    IRF_lb  <- array(NA_real_, c(N, lags, T.max))
    IRF_ub  <- array(NA_real_, c(N, lags, T.max))
    IRF_med <- array(NA_real_, c(N, lags, T.max))
    
    ev <- fit$eval
    if (is.null(ev)) stop("fit$eval is NULL; evaluate model before IRFs.")
    
    # quick stability check uses no-intercept versions
    Phi_c_noint <- Phi_c_est[, -1, drop = FALSE]
    Phi_f_noint <- Phi_f_est[, -1, , drop = FALSE]  # N x (N*p) x r
    lyap_val <- tryCatch(
      check_lyapunov(Phi_c_noint, Phi_f_noint, psi, seed_ly = 12L),
      error = function(e) Inf
    )
    
    #message("IRF for penalized model: parameter uncertainty not propagated — bands reflect only factor variation.")
    set.seed(1234L)
    for (t in 1:T.max) {
      aT <- ev$filtered.state[, t]
      print(dim(ev$filtered.state))
      PT <- matrix(ev$filtered.state.variance[, , t, drop = FALSE], nrow = r, ncol = r)

      HT <- matrix(ev$array.filtered.H[,,t], nrow = N, ncol = N)

      
      # Draw factor paths only; average across B by taking quantiles of B=1 runs
      irf_draws <- array(NA_real_, c(N, lags, B))
      for (b in 1:B) {
        one <- irf_generator_cpp(
          Phi, psi, aT, PT, HT, ek,
          lags = lags,
          fixed = TRUE,   # <-- factor simulation inside C++
          B = 1,
          seed_f = 1000L * t + b
        )
        irf_draws[,,b] <- one$irf
      }
      
      if (is.finite(lyap_val) && lyap_val > 0) irf_draws[,] <- NA_real_
      
      for (i in 1:N) for (h in 1:lags) {
        z <- irf_draws[i, h, ]
        IRF_lb[i, h, t]  <- stats::quantile(z, 0.16, na.rm = TRUE)
        IRF_med[i, h, t] <- stats::quantile(z, 0.50, na.rm = TRUE)
        IRF_ub[i, h, t]  <- stats::quantile(z, 0.84, na.rm = TRUE)
      }
    }
    
    return(list(IRF_lb=IRF_lb, IRF_med=IRF_med, IRF_ub=IRF_ub))
  }
}

#' Plot IRFs (single or multiple shocks)
#' @param irf_obj result of irf()
#' @param t which time index (default 1)
#' @export
plot_irf <- function(irf_obj, t = 1) {
  # normalize to a list-of-IRFs
  if (!is.list(irf_obj) || is.null(attr(irf_obj, "tvvar_irf_multi"))) {
    irf_obj <- list(shock_1 = irf_obj)
  }
  shocks <- names(irf_obj)
  N <- irf_obj[[1]]$meta$N
  lags <- irf_obj[[1]]$meta$lags
  
  if (length(shocks) == 1L) {
    # 2 panels (responses of all N series to one shock)
    op <- par(mfrow = c(1, N), mgp = c(2, .8, 0), mar = .1 + c(3,3,3,1))
    on.exit(par(op), add = TRUE)
    k <- irf_obj[[1]]$meta$shock_var
    for (i in 1:N) {
      lb  <- irf_obj[[1]]$IRF_lb[i, , t]
      md  <- irf_obj[[1]]$IRF_med[i, , t]
      ub  <- irf_obj[[1]]$IRF_ub[i, , t]
      ylim <- range(c(lb, md, ub), na.rm = TRUE)
      plot(0:(lags-1), md, type = "l", xlab = "Horizon", ylab = "",
           lwd = 2, ylim = ylim, main = sprintf("Resp y%d to shock y%d (t=%d)", i, k, t))
      lines(0:(lags-1), lb, lty = "dashed")
      lines(0:(lags-1), ub, lty = "dashed")
      abline(h = 0, col = "grey40")
    }
  } else if (length(shocks) == 2L && N == 2L) {
    # 2x2 grid: row = responder, col = shocked series
    op <- par(mfrow = c(2, 2), mgp = c(2, .8, 0), mar = .1 + c(3,3,3,1))
    on.exit(par(op), add = TRUE)
    for (s in seq_along(shocks)) {
      k <- irf_obj[[s]]$meta$shock_var
      for (i in 1:2) {
        lb <- irf_obj[[s]]$IRF_lb[i, , t]
        md <- irf_obj[[s]]$IRF_med[i, , t]
        ub <- irf_obj[[s]]$IRF_ub[i, , t]
        ylim <- range(c(lb, md, ub), na.rm = TRUE)
        plot(0:(lags-1), md, type = "l", xlab = "Horizon", ylab = "",
             lwd = 2, ylim = ylim, main = sprintf("Resp y%d to shock y%d (t=%d)", i, k, t))
        lines(0:(lags-1), lb, lty = "dashed")
        lines(0:(lags-1), ub, lty = "dashed")
        abline(h = 0, col = "grey40")
      }
    }
  } else {
    stop("plot_irf() currently supports 1 shock (N panels) or 2 shocks with N=2 (2x2).")
  }
}
      