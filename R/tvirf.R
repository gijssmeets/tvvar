#' Time-varying impulse response function (TVIRF)
#'
#' Computes Monte Carlo impulse responses for a fitted time-varying VAR model.
#' For each specified horizon, it generates simulated responses to shocks at given
#' positions in the system. Uncertainty bands are based on parameter draws.
#'
#' @param fit A fitted tvvar model object (from \code{tvfit()} or \code{tvpenfit()}).
#' @param horizon Integer. Forecast horizon (number of steps ahead) over which to compute IRFs.
#'   The function includes responses up to and including this number of steps.
#' @param T.max Number of time points (from the start of the sample) at which to evaluate the IRFs.
#' @param B Number of Monte Carlo draws used for uncertainty estimation.
#' @param shock_index Integer index (1..N) of the shocked variable. Ignored if \code{shock_position} is supplied.
#' @param shock_position Optional numeric vector of length N specifying the exact shock magnitudes
#'   (non-zero entries indicate which series receive shocks).
#'
#' @return Either a single IRF object — a list with components
#'   \code{IRF_lb}, \code{IRF_med}, \code{IRF_ub}, and \code{meta} —
#'   or a named list of such objects (one per shocked series),
#'   with attribute \code{attr(*, "tvvar_irf_multi")} set to \code{TRUE}.
#'
#' @details
#' The function simulates the system’s response to a one-standard-deviation shock at the
#' indicated variable(s), conditional on the estimated time-varying parameters.
#' If multiple shocks are provided via \code{shock_position}, a joint IRF is computed.
#' The horizon is fully inclusive: for example, \code{horizon = 10} produces responses
#' at horizons 0 through 10.
#'
#' @importFrom MASS mvrnorm
#' @export
tvirf <- function(fit,
                  horizon = 10,
                  T.max = 5,
                  B = 500,
                  shock_index = 1,
                  shock_position = NULL) {
  
  stopifnot(is.list(fit), !is.null(fit$meta))
  N <- fit$meta$N
  p <- fit$meta$p
  r <- fit$meta$r
  
  if (is.null(shock_position)) {
    warning("`shock_position` not provided — defaulting to a unit shock in variable 1.")
    shock_position <- rep(0, N)
    shock_position[1] <- 1
  } else {
    shock_position <- as.numeric(shock_position)
    if (length(shock_position) != N)
      stop("`shock_position` must be a numeric vector of length N (", N, ").")
  }
  
  shocked_ids <- which(shock_position != 0)
  if (length(shocked_ids) == 0L)
    stop("No shocks requested: all entries of `shock_position` are zero.")
  
  # container for multi-shock output (only used if >1 shocks)
  out_multi <- setNames(vector("list", length(shocked_ids)),
                        paste0("shock_", shocked_ids))
  
  # ========= Unpenalized branch (your original code, per shock) =========
  if (identical(fit$meta$method, "unpenalized")) {
    
    Phi.f.array <- fit$meta$Phi.f.array
    Phi.f.array.mat.structure <- matrix(Phi.f.array, nrow = N)
    Phi.f.array.mat.structure[Phi.f.array.mat.structure != 0] <- 1
    number.nonzero.phis <- sum(Phi.f.array.mat.structure != 0)
    VAR.data <- fit$meta$data
    zero.mean <- fit$meta$zero.mean
    
    for (kk in seq_along(shocked_ids)) {
      k <- shocked_ids[kk]  # which series to shock
      ek <- rep(0, N); ek[k] <- 1
      
      IRF_lb  <- array(NA_real_, c(N, horizon, T.max))
      IRF_ub  <- array(NA_real_, c(N, horizon, T.max))
      IRF_med <- array(NA_real_, c(N, horizon, T.max))
      
      for (t in 1:T.max) {
        T.index <- t
        theta <- fit$theta
        V     <- fit$vcov
        s.h   <- k
        
        irf       <- array(0, c(N, horizon, B))
        Pi        <- array(0, c(N * p, N * p, B))
        ev        <- numeric()
        lyapunov  <- numeric()
        
        for (i in 1:B) {
          s.h   <- s.h + 1
          index <- s.h * t
          set.seed(index)
          
          par_val <- MASS::mvrnorm(n = 1, mu = theta, Sigma = V)
          
          # obtain Phi,psi,aT,PT,HT at par draw
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
          
          phi_est <- par_val[3:(3 + r - 1)]
          
          Phi_f_est <- Phi.f.array.mat.structure
          Phi_f_est[Phi.f.array.mat.structure != 0] <-
            par_val[(3 + r):(3 + r + number.nonzero.phis - 1)]
          
          Phi_c_est <- opti.eval$Phi.c
          
          phi_est_tr <- exp(par_val[3:(3 + r - 1)]) / (exp(par_val[3:(3 + r - 1)]) + 1)
          Phi_c_noint <- Phi_c_est[, -1]
          Phi_f_noint <- array(Phi_f_est, dim = c(N, (N * p + 1), r))[,-1,, drop = FALSE]
          seed_ly <- 12
          
          lyapunov[i] <- check_lyapunov(Phi_c_noint, Phi_f_noint, phi_est_tr, seed_ly)
          
          Phi <- array(cbind(Phi_c_est, Phi_f_est), c(N, N * p + 1, (r + 1)))
          psi <- exp(phi_est) / (1 + exp(phi_est))
          
          aT <- opti.eval$filtered.state[, T.index]
          
          PT <- opti.eval$filtered.state.variance[, , T.index, drop = FALSE]
          PT <- matrix(PT, nrow = r, ncol = r)
          
          HT <- opti.eval$array.filtered.H[, , T.index]
          
          irf.gen <- irf_generator_cpp(Phi, psi, aT, PT, HT, ek,
                                       horizon, fixed = TRUE, B = B, seed_f = 1234)
          
          irf[,,i] <- irf.gen$irf
          Pi[,,i]  <- irf.gen$PI
          ev[i]    <- eigen(Pi[,,i])$values[1]
        }
        
        unstable <- which(lyapunov > 0)
        if (length(unstable)) irf[,,unstable] <- NA_real_
        
        for (ii in 1:N) {
          for (jj in 1:horizon) {
            IRF_lb[ii, jj, t]  <- quantile(irf[ii, jj, ], probs = 0.16, na.rm = TRUE)
            IRF_ub[ii, jj, t]  <- quantile(irf[ii, jj, ], probs = 0.84, na.rm = TRUE)
            IRF_med[ii, jj, t] <- quantile(irf[ii, jj, ], probs = 0.50, na.rm = TRUE)
          }
        }
      }
      
      obj <- list(IRF_lb = IRF_lb, IRF_med = IRF_med, IRF_ub = IRF_ub,
                  meta   = list(N = N, horizon = horizon, T.max = T.max,
                                shock_var = k,
                                shock_position = shock_position,
                                method = "ML"))
      class(obj) <- "tvirf"
      
      if (length(shocked_ids) == 1L) {
        return(obj)                    # single-shock: return the object directly
      } else {
        out_multi[[kk]] <- obj         # multi: collect and return after loop
      }
    }
    
    # name columns by shocked variable if not already named
    if (is.null(names(out_multi)) || any(!nzchar(names(out_multi)))) {
      names(out_multi) <- paste0("shock_", shocked_ids)
    }
    
    # unify class for multi-return as well
    class(out_multi) <- "tvirf"
    return(out_multi)
  }
  
  # ========= Penalized branch (your current code, per shock) =========
  {
    
    
    Phi_c_est <- fit$estimate$Phi_c          # N x (1+N*p)
    Phi_f_est <- fit$estimate$Phi_f          # N x (1+N*p) x r
    psi       <- as.numeric(fit$estimate$phi_r)
    
    if (length(dim(Phi_f_est)) != 3L)
      stop("fit$estimate$Phi_f must be an array [N x (1+N*p) x r].")
    
    # quick stability check (no intercept)
    Phi_c_noint <- Phi_c_est[, -1, drop = FALSE]
    Phi_f_noint <- Phi_f_est[, -1, , drop = FALSE]
    lyap_val <- tryCatch(
      check_lyapunov(Phi_c_noint, Phi_f_noint, psi, seed_ly = 12L),
      error = function(e) Inf
    )
    
    ev <- fit$eval
    if (is.null(ev)) stop("fit$eval is NULL; evaluate model before IRFs.")
    
    warning(
      "Penalized IRFs do not include parameter uncertainty. ",
      "Reported variability (if any) reflects only factor simulations."
    )
    
    for (kk in seq_along(shocked_ids)) {
      k <- shocked_ids[kk]
      ek <- rep(0, N); ek[k] <- 1
      
      # Build Phi cube: slice 1 = Phi^c, slices 2..r+1 = Phi^f_k
      Phi <- array(0, dim = c(N, N * p + 1, r + 1))
      Phi[,,1] <- Phi_c_est
      for (j in seq_len(r)) Phi[,,j + 1] <- Phi_f_est[,,j]
      
      IRF_lb  <- array(NA_real_, c(N, horizon, T.max))
      IRF_ub  <- array(NA_real_, c(N, horizon, T.max))
      IRF_med <- array(NA_real_, c(N, horizon, T.max))
      
      set.seed(1234L)
      for (t in 1:T.max) {
        aT <- ev$filtered.state[, t]
        PT <- matrix(ev$filtered.state.variance[, , t, drop = FALSE], nrow = r, ncol = r)
        HT <- matrix(ev$array.filtered.H[, , t], nrow = N, ncol = N)
        
        irf_draws <- array(NA_real_, c(N, horizon, B))
        for (b in 1:B) {
          one <- irf_generator_cpp(
            Phi, psi, aT, PT, HT, ek,
            lags = horizon,
            fixed = TRUE,   # fixed coefficients; (factor sim would be fixed=FALSE if desired)
            B = 1,
            seed_f = 1000L * t + b
          )
          irf_draws[,,b] <- one$irf
        }
        
        if (is.finite(lyap_val) && lyap_val > 0) irf_draws[,] <- NA_real_
        
        for (ii in 1:N) for (jj in 1:horizon) {
          z <- irf_draws[ii, jj, ]
          IRF_lb[ii, jj, t]  <- stats::quantile(z, 0.16, na.rm = TRUE)
          IRF_med[ii, jj, t] <- stats::quantile(z, 0.50, na.rm = TRUE)
          IRF_ub[ii, jj, t]  <- stats::quantile(z, 0.84, na.rm = TRUE)
        }
      }
      
      obj <- list(IRF_lb = IRF_lb, IRF_med = IRF_med, IRF_ub = IRF_ub,
                  meta   = list(N = N, horizon = horizon, T.max = T.max,
                                shock_var = k,
                                shock_position = shock_position,
                                method = fit$meta$method))
      class(obj) <- "tvirf"
      
      # return single or collect multi
      if (length(shocked_ids) == 1L) return(obj)
      out_multi[[kk]] <- obj
    }
    
    attr(out_multi, "tvvar_irf_multi") <- TRUE
    class(out_multi) <- "tvirf"
    return(out_multi)
  }
}

#' Plot IRFs (single or multiple shocks)
#' @param irf_obj result of irf() (single IRF list or multi-shock list)
#' @param t which time index to plot (default 1)
#' @export
plot.tvirf <- function(irf_obj, t = 1) {
  if (is.list(irf_obj) && !is.null(irf_obj$IRF_lb)) {
    # single object -> wrap as one-element list
    nm <- paste0("shock_", irf_obj$meta$shock_var %||% 1L)
    irf_obj <- setNames(list(irf_obj), nm)
  } else if (is.list(irf_obj) && length(irf_obj) >= 1L && inherits(irf_obj[[1]], "tvirf")) {
    # multi already; ensure names
    if (is.null(names(irf_obj)) || any(!nzchar(names(irf_obj)))) {
      names(irf_obj) <- vapply(irf_obj, function(x)
        paste0("shock_", x$meta$shock_var %||% 1L), character(1))
    }
  } else {
    stop("plot.tvirf: input must be a 'tvirf' object or a list of 'tvirf' objects.")
  }
  
  # Basic dims from first entry
  first <- irf_obj[[1]]
  N     <- first$meta$N
  horizon  <- first$meta$horizon
  S     <- length(irf_obj)                # number of shocks to show (relevant only)
  H     <- 0:(horizon - 1)                   # horizons
  
  # Layout: rows = responders (1..N), cols = shocks (the ones you requested)
  op <- par(mfrow = c(N, S), mgp = c(2, 0.8, 0), mar = 0.1 + c(3, 3, 3, 1))
  on.exit(par(op), add = TRUE)
  
  # Iterate responders (rows)
  for (i in seq_len(N)) {
    # Iterate shocks (cols)
    for (s in seq_len(S)) {
      irf_s  <- irf_obj[[s]]
      kshock <- irf_s$meta$shock_var %||% s   # which series was shocked
      
      lb <- irf_s$IRF_lb[i, , t]
      md <- irf_s$IRF_med[i, , t]
      ub <- irf_s$IRF_ub[i, , t]
      
      # Y-limits per panel (robust to all-NA)
      vals <- c(lb, md, ub)
      if (all(is.na(vals))) {
        ylim <- c(-1, 1)
      } else {
        ylim <- range(vals, na.rm = TRUE)
        if (!all(is.finite(ylim))) ylim <- c(-1, 1)
        if (diff(ylim) == 0) ylim <- ylim + c(-1, 1) * max(1e-6, abs(ylim[1]) * 0.1)
      }
      
      # Title: response-by-shock
      main_lbl <- sprintf("Response y%d to shock y%d (t=%d)", i, kshock, t)
      
      plot(0:(horizon - 1), md, type = "l", lwd = 2, xlab = "Horizon", ylab = "",
           ylim = ylim, main = main_lbl)
      if (!all(is.na(lb))) lines(0:(horizon - 1), lb, lty = "dashed")
      if (!all(is.na(ub))) lines(0:(horizon - 1), ub, lty = "dashed")
      abline(h = 0, col = "grey40")
    }
  }
  
  invisible(NULL)
}