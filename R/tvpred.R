#' Time-varying VAR forecasting
#'
#' Generates multi-step forecasts for a fitted tvvar model.
#' For unpenalized (ML) fits, parameter uncertainty is included via draws from
#' the asymptotic covariance. For penalized fits, parameters are fixed, and
#' only factor and shock uncertainty is simulated.
#'
#' @param fit Fitted tvvar object from `tvfit()` (unpenalized) or `tvpenfit()` (penalized)
#' @param h Integer forecast horizon
#' @param B Number of Monte Carlo forecast paths
#' @param seed Optional random seed
#' @param intervals Quantiles to report for prediction intervals
#' @param use_param_draws Logical; ignored for penalized fits
#' @return Object of class `tvpred` containing median, mean, lower/upper bounds, and meta info.
#' @export
tvpred <- function(fit,
                   h = 8,
                   B = 500,
                   seed = NULL,
                   intervals = c(0.16, 0.50, 0.84),
                   use_param_draws = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  
  N   <- fit$meta$N
  p   <- fit$meta$p
  r   <- fit$meta$r
  VAR.data <- fit$meta$data
  T.fin <- nrow(VAR.data)
  
  is_penalized <- !identical(fit$meta$method, "unpenalized")
  
  if (is_penalized) {
    warning("Penalized forecast: parameter uncertainty ignored; only factor & shock variability simulated.")
    Phi_c_est <- fit$estimate$Phi_c
    Phi_f_est <- fit$estimate$Phi_f
    psi       <- as.numeric(fit$estimate$phi_r)
    ev        <- fit$eval
  } else {
    zm  <- fit$meta$zero.mean
    Phi.f.array <- fit$meta$Phi.f.array
    theta <- fit$theta
    V     <- fit$vcov
    if (is.null(theta) || is.null(V)) {
      theta <- fit$optim$par
      V     <- diag(1e-6, length(theta))
    }
    
    # bookkeeping for free Phi_f entries
    Phi.f.as.mat <- matrix(Phi.f.array, nrow = N)
    Phi.f.as.mat[Phi.f.as.mat != 0] <- 1
    n_phi_free <- sum(Phi.f.as.mat != 0)
  }
  
  inv_logit <- stats::plogis
  Ypaths <- array(NA_real_, dim = c(B, h, N))
  
  # last p lags
  y_hist <- VAR.data[(T.fin - p + 1):T.fin, , drop = FALSE]
  make_lagvec <- function(Ylag) as.numeric(t(Ylag[nrow(Ylag):(nrow(Ylag)-p+1), , drop = FALSE]))
  y_lags0 <- make_lagvec(rbind(y_hist[nrow(y_hist), , drop = FALSE], y_hist[1:(p-1), , drop = FALSE]))
  
  for (b in seq_len(B)) {
    
    if (!is_penalized) {
      # === Unpenalized (ML) branch ===
      par_b <- if (use_param_draws) MASS::mvrnorm(1, mu = theta, Sigma = V) else theta
      
      ev <- opti.fct(
        par_free    = par_b,
        par_fixed   = NaN,
        VAR.data    = VAR.data,
        Phi.f.array = Phi.f.array,
        cfg         = list(zero.mean = zm, dim.VAR = N, lag.order = p, number.factors = r),
        Smooth      = FALSE,
        purpose     = "eval"
      )
      
      Phi_c  <- ev$Phi.c
      HT0    <- ev$array.filtered.H[, , ncol(ev$array.filtered.H)]
      aT     <- if (r > 0) ev$filtered.state[, ncol(ev$filtered.state)] else numeric(0)
      
      # rebuild dynamic params
      idx <- 1L
      A2  <- inv_logit(par_b[idx]); idx <- idx + 1L
      B2  <- inv_logit(par_b[idx]); idx <- idx + 1L
      psi <- if (r > 0) inv_logit(par_b[idx:(idx + r - 1L)]) else numeric(0L)
      idx <- idx + length(psi)
      
      Phi_f_est <- Phi.f.as.mat
      if (n_phi_free > 0)
        Phi_f_est[Phi.f.as.mat != 0] <- par_b[idx:(idx + n_phi_free - 1L)]
      Phi_f_arr <- array(Phi_f_est, dim = dim(Phi.f.array))
      
      Omega_b <- drop((1 - A2 - B2)) * HT0
      Phi_c_noint <- Phi_c[, -1, drop = FALSE]
      Phi_f_noint <- if (r > 0) Phi_f_arr[, -1, , drop = FALSE] else array(0, dim = c(N, N*p, 0))
      
    } else {
      # === Penalized branch ===
      Phi_c_noint <- Phi_c_est[, -1, drop = FALSE]
      Phi_f_noint <- Phi_f_est[, -1, , drop = FALSE]
      HT0   <- ev$array.filtered.H[, , ncol(ev$array.filtered.H)]
      aT    <- if (r > 0) ev$filtered.state[, ncol(ev$filtered.state)] else numeric(0)
      A2 <- 0.2; B2 <- 0.75  # fixed scalar BEKK params (approx)
      Omega_b <- drop((1 - A2 - B2)) * HT0
    }
    
    # === Forecast recursion ===
    f_t <- if (r > 0) as.numeric(aT) else numeric(0)
    H_t <- HT0
    y_lag <- y_lags0
    
    for (tstep in 1:h) {
      if (r > 0) {
        Q <- diag(1 - psi^2, r)
        f_t <- psi * f_t + MASS::mvrnorm(1, mu = rep(0, r), Sigma = Q)
      }
      
      Phi_t <- Phi_c_noint
      if (r > 0)
        for (k in seq_len(r)) Phi_t <- Phi_t + Phi_f_noint[, , k] * f_t[k]
      
      eps_t <- as.numeric(t(chol(H_t)) %*% rnorm(N))
      y_next <- as.numeric(matrix(Phi_t, nrow = N, ncol = N*p) %*% matrix(y_lag, ncol = 1)) + eps_t
      
      Ypaths[b, tstep, ] <- y_next
      if (p > 1) y_lag <- c(y_next, y_lag[1:(N*(p - 1))]) else y_lag <- y_next
      H_t <- Omega_b + B2 * H_t + A2 * (eps_t %*% t(eps_t))
    }
  }
  
  # summarize forecast draws
  qfun <- function(x, probs) stats::quantile(x, probs = probs, na.rm = TRUE)
  med  <- apply(Ypaths, c(2, 3), qfun, probs = 0.50)
  lb   <- apply(Ypaths, c(2, 3), qfun, probs = min(intervals))
  ub   <- apply(Ypaths, c(2, 3), qfun, probs = max(intervals))
  mu   <- apply(Ypaths, c(2, 3), mean, na.rm = TRUE)
  
  fc <- list(
    median = t(med),
    lb     = if (is_penalized) NULL else t(lb),
    ub     = if (is_penalized) NULL else t(ub),
    mean   = t(mu),
    data   = VAR.data,
    meta   = list(N = N, h = h, B = B, method = if (is_penalized) "tvpred_pen" else "tvpred")
  )
  
  class(fc) <- c("tvpred", "tvvar_result")
  return(fc)
}

#' Plot forecasts with optional history and forecast boundary marker
#'
#' Plots forecasted means or medians (with CI if available), shows recent history,
#' and adds a vertical line at the forecast start. By default (`var = NULL`) it
#' plots all variables (one panel per series, vertically stacked).
#'
#' @param fc list returned by `tvpred()`. Must contain
#'   `median` or `mean` (N×h or h×N), optional `lb`/`ub`, and `data` (T×N).
#' @param var integer vector of variables to plot (1..N). Default `NULL` = all.
#' @param hist_n number of last historical observations to show (default 50).
#' @param main optional overall title (used in base mode).
#' @param ci_col shading color for confidence band.
#' @param line_col color for forecast line.
#' @param hist_col color for history line.
#' @param connect_col color for connecting line (default same as forecast line).
#' @param vline_col color for vertical split line (default grey40).
#' @param vline_lty linetype for split (default dotted, 2).
#' @param show_median logical; if TRUE and `fc$median` exists, plot medians; else means.
#' @param y_hist optional history override. If supplied, can be a vector (T) for
#'   a single series or a matrix (T×N) matching `fc$data`.
#' @param ... passed to underlying plot functions.
#' @export
plot.tvpred <- function(fc,
                          var = NULL,
                          hist_n = 50,
                          main = NULL,
                          ci_col = "grey80",
                          line_col = "steelblue",
                          hist_col = "black",
                          connect_col = NULL,
                          vline_col = "grey40",
                          vline_lty = 2,
                          show_median = TRUE,
                          y_hist = NULL,
                          ...) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
  connect_col <- connect_col %||% line_col
  stopifnot(is.list(fc))
  if (is.null(fc$data))
    stop("`fc$data` (T×N) was not found in `fc`. The forecast object must include the original data.")
  
  Y <- as.matrix(fc$data)           # T × N history
  Tn <- nrow(Y); N <- ncol(Y)
  
  # Choose center matrix
  center0 <- if (isTRUE(show_median) && !is.null(fc$median)) fc$median else fc$mean
  if (is.null(center0)) stop("`fc` must contain either `median` or `mean`.")
  
  # Normalize forecast matrices to N×h (rows = series, cols = horizons)
  to_Nxh <- function(M) {
    if (is.null(M)) return(NULL)
    M <- as.matrix(M)
    if (nrow(M) == N) {
      M
    } else if (ncol(M) == N) {
      t(M)
    } else {
      stop("Could not align forecast matrix with N = ", N, ". Expect N×h or h×N.")
    }
  }
  center <- to_Nxh(center0)
  lb     <- to_Nxh(fc$lb)
  ub     <- to_Nxh(fc$ub)
  h      <- ncol(center)
  
  # Pick series set
  if (is.null(var)) var <- seq_len(N)
  if (any(var < 1 | var > N)) stop("`var` must be indices in 1..N.")
  
  # Determine history source (override if provided)
  if (!is.null(y_hist)) {
    y_hist <- as.matrix(y_hist)
    if (nrow(y_hist) != Tn) stop("`y_hist` must have ", Tn, " rows to match `fc$data`.")
    if (ncol(y_hist) == 1L && length(var) == 1L) {
      Y_use <- matrix(NA_real_, nrow = Tn, ncol = N)
      Y_use[, var] <- y_hist[, 1]
    } else if (ncol(y_hist) == N) {
      Y_use <- y_hist
    } else {
      stop("`y_hist` must be a vector (T) for one series or a matrix (T×N).")
    }
  } else {
    Y_use <- Y
  }
  
  hist_n <- max(1L, min(hist_n, Tn))
  Hx <- seq_len(h)                     # 1..h
  
  # ggplot version (preferred)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    df_list <- vector("list", length(var))
    for (idx in seq_along(var)) {
      j <- var[idx]
      y_hist_j <- tail(Y_use[, j], hist_n)
      x_hist   <- seq.int(-hist_n + 1L, 0L)
      
      df_h <- data.frame(
        series = factor(paste0("y", j), levels = paste0("y", var)),
        x      = x_hist,
        value  = y_hist_j,
        part   = "history",
        lower  = NA_real_,
        upper  = NA_real_
      )
      
      df_f <- data.frame(
        series = factor(paste0("y", j), levels = paste0("y", var)),
        x      = Hx,
        value  = center[j, ],
        part   = "forecast",
        lower  = if (is.null(lb)) NA_real_ else lb[j, ],
        upper  = if (is.null(ub)) NA_real_ else ub[j, ]
      )
      
      df_list[[idx]] <- rbind(df_h, df_f)
    }
    df <- do.call(rbind, df_list)
    
    gg <- ggplot2::ggplot()
    
    # ribbons
    if (!is.null(lb) && !is.null(ub)) {
      gg <- gg +
        ggplot2::geom_ribbon(
          data = subset(df, part == "forecast"),
          ggplot2::aes(x = x, ymin = lower, ymax = upper),
          fill = ci_col, alpha = 0.45
        )
    }
    
    # forecast + history lines
    gg <- gg +
      ggplot2::geom_line(
        data = subset(df, part == "forecast"),
        ggplot2::aes(x = x, y = value),
        linewidth = 1.1, color = line_col
      ) +
      ggplot2::geom_line(
        data = subset(df, part == "history"),
        ggplot2::aes(x = x, y = value),
        linewidth = 0.7, color = hist_col
      )
    
    # connector per series
    connectors <- do.call(rbind, lapply(split(df, df$series), function(dsi) {
      dh <- subset(dsi, part == "history")
      dfc <- subset(dsi, part == "forecast")
      if (nrow(dh) == 0 || nrow(dfc) == 0) return(NULL)
      data.frame(
        series = unique(dsi$series),
        x     = 0, xend = 1,
        y     = dh$value[which.max(dh$x)],
        yend  = dfc$value[1]
      )
    }))
    
    if (!is.null(connectors) && nrow(connectors)) {
      gg <- gg + ggplot2::geom_segment(
        data = connectors,
        ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
        linewidth = 1.05, color = connect_col
      )
    }
    
    # vertical split + facets (VERTICAL)
    gg <- gg +
      ggplot2::geom_vline(xintercept = 0, linetype = vline_lty, color = vline_col) +
      ggplot2::facet_wrap(~ series, ncol = 1, scales = "free_y") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(
        title = main %||% "Forecasts",
        x = "Time (0 = forecast start)",
        y = if (isTRUE(show_median) && !is.null(fc$median)) "Median" else "Mean"
      )
    
    return(gg)
  }
  
  ## --- Base R fallback (stacked vertically) ---
  K <- length(var)
  op <- par(mfrow = c(K, 1), mgp = c(2, 0.8, 0), mar = 0.1 + c(3,3,2,1))
  on.exit(par(op), add = TRUE)
  
  for (j in var) {
    y_hist_j <- tail(Y_use[, j], hist_n)
    x_hist   <- seq.int(-hist_n + 1L, 0L)
    y_center <- center[j, ]
    y_lb     <- if (is.null(lb)) rep(NA_real_, h) else lb[j, ]
    y_ub     <- if (is.null(ub)) rep(NA_real_, h) else ub[j, ]
    
    yrng <- range(c(y_hist_j, y_center, y_lb, y_ub), na.rm = TRUE)
    plot(x_hist, y_hist_j, type = "l", col = hist_col, lwd = 1,
         ylim = yrng, xlab = "Time (0 = forecast start)",
         ylab = if (isTRUE(show_median) && !is.null(fc$median)) "Median" else "Mean",
         main = paste0("y", j), ...)
    
    if (!all(is.na(y_lb)) && !all(is.na(y_ub))) {
      polygon(c(Hx, rev(Hx)), c(y_lb, rev(y_ub)),
              col = grDevices::adjustcolor(ci_col, alpha.f = 0.4), border = NA)
    }
    lines(Hx, y_center, col = line_col, lwd = 2)
    
    # connector
    segments(0, tail(y_hist_j, 1), 1, y_center[1], col = connect_col, lwd = 2)
    abline(v = 0, lty = vline_lty, col = vline_col)
  }
  
  invisible(NULL)
}