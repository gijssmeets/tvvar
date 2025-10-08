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
#' @export
impulse_response <- function(data,
                             fit,
                             lags = 10,
                             T.max = 10,
                             B = 500,
                             shock_index = 1,
                             shock_position = c(0,0,1),
                             seed = 1234) {
  
  N <- fit$meta$N
  p <- fit$meta$p
  r <- fit$meta$r
  zero.mean <- fit$meta$zero.mean
  Phi.f.array <- fit$meta$Phi.f.array
  Phi.f.array.mat.structure <- matrix(Phi.f.array, nrow = N)
  Phi.f.array.mat.structure[Phi.f.array.mat.structure != 0] <- 1
  number.nonzero.phis <- sum(Phi.f.array.mat.structure != 0)
  
  VAR.data <- data

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
      
      print(theta)
      print(V)
      par_val = mvrnorm(n = 1, mu = theta, Sigma = V)
      
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
      
      print(phi_est_tr)
      print(Phi_c_noint)
      print(Phi_f_noint)
      
      print(dim(phi_est_tr))
      print(dim(Phi_c_noint))
      print(dim(Phi_f_noint))
      
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
    
    
    print(t)
    
  }
  
  return(list(IRF_lb=IRF_lb, IRF_med=IRF_med, IRF_ub=IRF_ub) )
}


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
      