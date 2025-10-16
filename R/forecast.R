#' Multi-step forecasts for a tvvar_fit (unpenalized ML)
#'
#' Simulates recursive forecasts y_{T+1}, …, y_{T+h} using parameter draws
#' from the asymptotic normal N(theta, vcov) and simulated factor paths.
#' Volatility is propagated with a simple BEKK(1,1) recursion using simulated shocks.
#'
#' @param fit   tvvar_fit from unpenalized_estimate(..., method = "ML")
#' @param h     integer horizon (number of steps ahead)
#' @param B     number of Monte Carlo paths (set B=1 for plug-in path)
#' @param seed  optional integer seed (reproducibility)
#' @param intervals numeric probs to report (default 16%, 50%, 84%)
#' @param use_param_draws logical; if FALSE, uses plug-in theta only
#'
#' @return list with elements:
#'   \item{mean, median, lb, ub}{N x h matrices of forecast stats}
#'   \item{paths}{optional B x h x N array (only if B not too large)}
#' @importFrom MASS mvrnorm
#' @export
forecast_tvvar <- function(fit,
                           h = 8,
                           B = 500,
                           seed = NULL,
                           intervals = c(0.16, 0.50, 0.84),
                           use_param_draws = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  
  N   <- fit$meta$N
  p   <- fit$meta$p
  r   <- fit$meta$r
  zm  <- fit$meta$zero.mean
  Phi.f.array <- fit$meta$Phi.f.array        # N x (N*p+1) x r (mask-built)
  VAR.data    <- fit$meta$data
  T.fin <- nrow(VAR.data)
  
  # parameter center & covariance (as in your IRF)
  theta <- fit$theta
  V     <- fit$vcov
  if (is.null(theta) || is.null(V)) {
    # fallback (should not normally happen for ML)
    theta <- fit$optim$par
    V     <- diag(1e-6, length(theta))
  }
  
  # free-pattern bookkeeping (same as IRF)
  Phi.f.as.mat <- matrix(Phi.f.array, nrow = N)
  Phi.f.as.mat[Phi.f.as.mat != 0] <- 1
  n_phi_free <- sum(Phi.f.as.mat != 0)
  
  # helper: slice theta into (A,B,phi_r, Phi_f_free, Lvech, Phi_c_vec)
  inv_logit <- stats::plogis
  
  # storage for forecast draws
  Ypaths <- array(NA_real_, dim = c(B, h, N))
  
  # last p lags from the sample (stacked as y_{T}, y_{T-1}, ..., y_{T-p+1})
  y_hist <- VAR.data[(T.fin - p + 1):T.fin, , drop = FALSE]
  # companion-state lag vector in order: [y_{t-1}^{(1..N)}, y_{t-2}^{(1..N)}, ...]
  make_lagvec <- function(Ylag) as.numeric(t(Ylag[nrow(Ylag):(nrow(Ylag)-p+1), , drop=FALSE]))
  y_lags0 <- make_lagvec(rbind(y_hist[nrow(y_hist), , drop=FALSE], y_hist[1:(p-1), , drop=FALSE]))
  
  for (b in seq_len(B)) {
    # 1) draw parameters (or plug-in)
    par_b <- if (use_param_draws) MASS::mvrnorm(1, mu = theta, Sigma = V) else theta
    
    # 2) get evaluation pieces at (par_b)
    ev <- opti.fct(
      par_free    = par_b,
      par_fixed   = NaN,
      VAR.data    = VAR.data,
      Phi.f.array = Phi.f.array,
      cfg         = list(zero.mean = zm, dim.VAR = N, lag.order = p, number.factors = r),
      Smooth      = FALSE,
      purpose     = "eval"
    )
    Phi_c  <- ev$Phi.c                                # N x (1 + N*p)
    aT     <- if (r > 0) ev$filtered.state[, ncol(ev$filtered.state), drop=TRUE] else numeric(0)
    HT0    <- ev$array.filtered.H[, , ncol(ev$array.filtered.H)]
    # Rebuild A,B scalars (first two are unconstrained, map to (0,1))
    idx <- 1L
    A2  <- inv_logit(par_b[idx]); idx <- idx + 1L
    B2  <- inv_logit(par_b[idx]); idx <- idx + 1L
    psi <- if (r > 0) inv_logit(par_b[idx:(idx + r - 1L)]) else numeric(0L)
    idx <- idx + length(psi)
    
    # Phi_f free entries, rebuild Phi_f array (N x (N*p+1) x r)
    Phi_f_est <- Phi.f.as.mat
    if (n_phi_free > 0) {
      Phi_f_est[Phi.f.as.mat != 0] <- par_b[idx:(idx + n_phi_free - 1L)]
    }
    idx <- idx + n_phi_free
    Phi_f_arr <- array(Phi_f_est, dim = dim(Phi.f.array))   # includes intercept column; we’ll drop later
    
    # BEKK omega (approx via unconditional relation using HT0)
    # Unconditional: H ≈ omega / (1 - A2 - B2)  => omega ≈ (1 - A2 - B2) * H
    Omega_b <- drop((1 - A2 - B2)) * HT0
    
    # Drop intercept column for VAR propagation
    Phi_c_noint <- Phi_c[, -1, drop = FALSE]                  # N x (N*p)
    Phi_f_noint <- if (r > 0) Phi_f_arr[, -1, , drop = FALSE] else array(0, dim=c(N, N*p, 0))
    
    # initialize factor state and H
    f_t   <- if (r > 0) as.numeric(aT) else numeric(0)
    H_t   <- HT0
    y_lag <- y_lags0
    
    for (tstep in 1:h) {
      # factor evolve: f_{t+1} = psi * f_t + eta, eta ~ N(0, Q) with Q = I - diag(psi^2)
      if (r > 0) {
        Q <- diag(1 - psi^2, r)
        f_t <- psi * f_t + MASS::mvrnorm(1, mu = rep(0, r), Sigma = Q)
      }
      
      # time-varying VAR matrix at this step:
      Phi_t <- Phi_c_noint
      if (r > 0) {
        for (k in seq_len(r)) Phi_t <- Phi_t + Phi_f_noint[, , k] * f_t[k]
      }
      
      # shock ~ N(0, H_t), then update H_{t+1} via BEKK(1,1) with scalar A2,B2
      eps_t <- as.numeric(t(chol(H_t)) %*% rnorm(N))
      
      # y_{t+1} = Phi_t * [y_t, y_{t-1}, ..., y_{t-p+1}] + eps_t
      y_next <- matrix(Phi_t, nrow = N, ncol = N * p) %*% matrix(y_lag, ncol = 1)
      y_next <- as.numeric(y_next) + eps_t
      
      # store
      Ypaths[b, tstep, ] <- y_next
      
      # update lag vector: prepend y_next, drop oldest block
      if (p > 1) {
        y_lag <- c(y_next, y_lag[1:(N * (p - 1))])
      } else {
        y_lag <- y_next
      }
      
      # BEKK recursion for next H: H_{t+1} = Omega + B2 * H_t + A2 * eps_t eps_t'
      H_t <- Omega_b + B2 * H_t + A2 * (eps_t %*% t(eps_t))
    }
  }
  
  # collapse draws -> stats
  # (B x h x N) → for each variable n and horizon t, take quantiles across B
  qfun <- function(x, probs) stats::quantile(x, probs = probs, na.rm = TRUE)
  med  <- apply(Ypaths, c(2, 3), qfun, probs = 0.50)
  lb   <- apply(Ypaths, c(2, 3), qfun, probs = min(intervals))
  ub   <- apply(Ypaths, c(2, 3), qfun, probs = max(intervals))
  mu   <- apply(Ypaths, c(2, 3), mean, na.rm = TRUE)
  
  # return as N x h matrices (transpose: horizons in columns)
  list(
    median = t(med),
    lb     = t(lb),
    ub     = t(ub),
    mean   = t(mu)
    # You can also return paths = Ypaths if you want raw draws
  )
}


#' Plot forecasts with optional history and forecast boundary marker
#'
#' Plots forecasted means or medians (with CI if available), shows optional
#' recent history before the forecast, and adds a vertical line at the
#' forecast start. The last observed point is connected to the first forecast.
#'
#' @param fc list returned by `forecast_tvvar()`.
#' @param var integer index of variable (1..N).
#' @param y_hist optional numeric vector of observed values (same variable).
#' @param hist_n number of last historical observations to show (default 50).
#' @param main plot title.
#' @param ci_col shading color for confidence band.
#' @param line_col color for forecast line.
#' @param hist_col color for history line.
#' @param connect_col color for connecting line (default same as forecast line).
#' @param vline_col color for vertical split line (default grey40).
#' @param vline_lty linetype for split (default dotted, 2).
#' @param show_median logical; if TRUE, plot median if available.
#' @param ... passed to underlying plot functions.
#' @export
plot_forecast <- function(fc,
                          var = 1,
                          y_hist = NULL,
                          hist_n = 50,
                          main = NULL,
                          ci_col = "grey80",
                          line_col = "steelblue",
                          hist_col = "black",
                          connect_col = NULL,
                          vline_col = "grey40",
                          vline_lty = 2,
                          show_median = TRUE,
                          ...) {
  
  connect_col <- connect_col %||% line_col
  stopifnot(is.list(fc))
  
  center_mat <- if (show_median && !is.null(fc$median)) fc$median else fc$mean
  if (is.null(center_mat))
    stop("`fc` must have at least `mean` (N x h).")
  
  N <- nrow(center_mat); h <- ncol(center_mat)
  if (var > N) stop("`var` exceeds N.")
  center <- as.numeric(center_mat[var, ])
  lower  <- if (!is.null(fc$lb)) as.numeric(fc$lb[var, ]) else rep(NA_real_, h)
  upper  <- if (!is.null(fc$ub)) as.numeric(fc$ub[var, ]) else rep(NA_real_, h)
  
  # Optional history
  has_hist <- !is.null(y_hist)
  if (has_hist) {
    y_hist <- as.numeric(y_hist)
    hist_n <- min(hist_n, length(y_hist))
    y_hist_plot <- tail(y_hist, hist_n)
    x_hist <- seq.int(-hist_n + 1L, 0L)
    x_fc <- seq_len(h)
  } else {
    y_hist_plot <- NULL
    x_hist <- NULL
    x_fc <- seq_len(h)
  }
  
  # Try ggplot2
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    gg <- ggplot2::ggplot()
    
    if (has_hist) {
      df_hist <- data.frame(x = x_hist, y = y_hist_plot)
      gg <- gg + ggplot2::geom_line(data = df_hist, ggplot2::aes(x, y),
                                    color = hist_col, linewidth = 0.6)
    }
    
    df_fc <- data.frame(x = x_fc, center = center, lower = lower, upper = upper)
    
    if (!all(is.na(lower)) && !all(is.na(upper))) {
      gg <- gg + ggplot2::geom_ribbon(
        data = df_fc,
        ggplot2::aes(x = x, ymin = lower, ymax = upper),
        fill = ci_col, alpha = 0.5
      )
    }
    
    gg <- gg + ggplot2::geom_line(
      data = df_fc,
      ggplot2::aes(x = x, y = center),
      color = line_col, linewidth = 1.1
    )
    
    # Add connecting line (from last history to first forecast)
    if (has_hist) {
      gg <- gg +
        ggplot2::geom_segment(
          x = 0, xend = 1,
          y = tail(y_hist_plot, 1),
          yend = center[1],
          color = connect_col,
          linewidth = 1.0
        ) +
        ggplot2::geom_vline(xintercept = 0, linetype = vline_lty, color = vline_col)
    }
    
    gg <- gg +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::labs(
        title = main %||% sprintf("Forecast (series %d)", var),
        x = if (has_hist) "Time (0 = forecast start)" else "Horizon",
        y = if (show_median && !is.null(fc$median)) "Median" else "Mean"
      )
    
    return(gg)
  }
  
  # --- Base R fallback ---
  if (has_hist) {
    y_all <- c(y_hist_plot, center, lower, upper)
    y_rng <- range(y_all, na.rm = TRUE)
    
    plot(x_hist, y_hist_plot, type = "l", col = hist_col, lwd = 1,
         ylim = y_rng,
         xlab = "Time (0 = forecast start)",
         ylab = if (show_median && !is.null(fc$median)) "Median" else "Mean",
         main = main %||% sprintf("Forecast (series %d)", var),
         ...)
    
    if (!any(is.na(lower)) && !any(is.na(upper))) {
      polygon(c(x_fc, rev(x_fc)), c(lower, rev(upper)),
              col = adjustcolor(ci_col, alpha.f = 0.4), border = NA)
    }
    lines(x_fc, center, col = line_col, lwd = 2)
    
    # connecting line between last observed and first forecast
    segments(0, tail(y_hist_plot, 1), 1, center[1],
             col = connect_col, lwd = 2)
    
    abline(v = 0, lty = vline_lty, col = vline_col)
    
  } else {
    y_rng <- range(c(center, lower, upper), na.rm = TRUE)
    plot(x_fc, center, type = "l", col = line_col, lwd = 2,
         ylim = y_rng, xlab = "Horizon",
         ylab = if (show_median && !is.null(fc$median)) "Median" else "Mean",
         main = main %||% sprintf("Forecast (series %d)", var), ...)
    if (!any(is.na(lower)) && !any(is.na(upper))) {
      polygon(c(x_fc, rev(x_fc)), c(lower, rev(upper)),
              col = adjustcolor(ci_col, alpha.f = 0.4), border = NA)
      lines(x_fc, center, col = line_col, lwd = 2)
    }
  }
  
  invisible(NULL)
}