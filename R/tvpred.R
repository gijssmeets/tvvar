#' Forecast from a fitted tvfit model
#'
#' @description
#' Generates multi-step forecasts from a fitted time-varying VAR model.
#' Supports both unpenalized and penalized fits. Allows choice of point
#' forecast type (mean or median) and returns both 68% and 95% intervals.
#'
#' @param object A fitted \code{tvfit} or \code{tvpenfit} object.
#' @param h Integer, forecast horizon.
#' @param B Number of Monte Carlo forecast paths.
#' @param seed Optional random seed.
#' @param point Character, either \code{"mean"} or \code{"median"} for the central forecast.
#' @param use_param_draws Logical; ignored for penalized fits.
#' @return An object of class \code{tvpred} with elements:
#' \itemize{
#'   \item \code{point} — mean or median forecast (N × h)
#'   \item \code{lb68}, \code{ub68} — 68% interval bounds
#'   \item \code{lb95}, \code{ub95} — 95% interval bounds
#'   \item \code{meta} — list with N, h, B, method, and point type.
#' }
#' @export
predict.tvfit <- function(object,
                   h = 8,
                   B = 500,
                   seed = NULL,
                   point = c("mean", "median"),
                   use_param_draws = TRUE) {
  point <- match.arg(point)
  if (!is.null(seed)) set.seed(seed)
  
  fit <- object
  N <- fit$meta$N
  p <- fit$meta$p
  r <- fit$meta$r
  VAR.data <- fit$meta$data
  T.fin <- nrow(VAR.data)
  is_penalized <- !identical(fit$meta$method, "unpenalized")
  
  if (is_penalized) {
    warning("Penalized forecast: parameter uncertainty ignored; only factor & shock variability simulated.")
    Phi_c_est <- object$estimate$Phi_c
    Phi_f_est <- object$estimate$Phi_f
    psi       <- as.numeric(object$estimate$phi_r)
    ev        <- object$eval
  } else {
    zm  <- object$meta$zero.mean
    Phi.f.array <- object$meta$Phi.f.array
    theta <- object$theta
    V     <- object$vcov
    if (is.null(theta) || is.null(V)) {
      theta <- object$optim$par
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
  
  qfun <- function(x, probs) stats::quantile(x, probs = probs, na.rm = TRUE)
  
  lb68 <- apply(Ypaths, c(2, 3), qfun, probs = 0.16)
  ub68 <- apply(Ypaths, c(2, 3), qfun, probs = 0.84)
  lb95 <- apply(Ypaths, c(2, 3), qfun, probs = 0.025)
  ub95 <- apply(Ypaths, c(2, 3), qfun, probs = 0.975)
  med  <- apply(Ypaths, c(2, 3), qfun, probs = 0.50)
  mu   <- apply(Ypaths, c(2, 3), mean, na.rm = TRUE)
  
  point_forecast <- switch(point,
                           mean = t(mu),
                           median = t(med)
  )
  
  fc <- list(
    point = point_forecast,
    lb68  = t(lb68),
    ub68  = t(ub68),
    lb95  = t(lb95),
    ub95  = t(ub95),
    meta  = list(
      N = N, h = h, B = B,
      method = if (is_penalized) "tvpred_pen" else "tvpred",
      point  = point
    )
  )
  
  class(fc) <- c("tvpred", "tvvar_result")
  return(fc)
}



#' Plot forecasts from a tvpred object
#'
#' Visualizes forecasts produced by \code{predict.tvfit()}, showing both 68\% and 95\% prediction intervals
#' (with lighter shading for the 95\% band). The plot displays the last \code{hist_n} observations of the
#' historical sample together with the forecast horizon.
#'
#' If a time index is supplied through \code{hist_time}, it will be used to label the x-axis and to extend
#' the forecast horizon automatically by the same spacing (e.g. monthly or yearly). Otherwise, a simple
#' numeric index (1, 2, …) is used for both history and forecast periods.
#'
#' @param x A `tvpred` object returned by \code{predict.tvfit()}.
#' @param var Integer vector of variables to plot (1..N). Default = all.
#' @param hist_y Optional historical data (T×N matrix or vector) corresponding to the fitted sample.
#' @param hist_time Optional vector of time indices or dates for the x-axis. 
#'   Can be one of:
#'   \itemize{
#'     \item A \code{Date} vector, e.g. \code{seq(as.Date("2000-01-01"), by = "month", length.out = T)}.
#'     \item A \code{POSIXct} vector for higher-frequency data (e.g. daily or intraday).
#'     \item A numeric sequence representing time steps.
#'   }
#'   When provided, the function automatically takes the last \code{hist_n} values and
#'   generates \code{h} additional time points beyond the final observed date using the same spacing.
#' @param hist_n Number of last historical points to display (default = 50).
#' @param main Optional overall plot title.
#' @param ci_col68 Color for the 68\% prediction interval shading (default = "grey75").
#' @param ci_col95 Color for the 95\% prediction interval shading (default = "grey90").
#' @param line_col Color for the forecast line (default = "steelblue").
#' @param hist_col Color for the historical line (default = "black").
#' @param connect_col Color for the connecting line between history and forecast 
#'   (default = same as forecast line).
#' @param vline_col Color for the vertical line at forecast start (default = "grey40").
#' @param vline_lty Linetype for the vertical line (default = 2, dotted).
#' @param ... Additional graphical arguments passed to lower-level plotting functions.
#'
#' @details
#' The time index handling works as follows:
#' \itemize{
#'   \item If \code{hist_time} is provided and of class \code{Date} or \code{POSIXct}, 
#'     the function infers the frequency (e.g. daily, monthly, yearly) from the median spacing 
#'     and extends the forecast horizon accordingly.
#'   \item If no \code{hist_time} is provided, the function defaults to integer indices.
#' }
#'
#' @examples
#' \dontrun{
#' # Example with simulated monthly data
#' simdata$date <- seq(as.Date("2020-01-01"), by = "month", length.out = nrow(simdata$Y))
#' fit <- tvfit(simdata$Y)
#' pred <- predict(fit, h = 12)
#' plot(pred, hist_y = simdata$Y, hist_time = simdata$date)
#' }
#'
#' @import ggplot2
#' @export
plot.tvpred <- function(x,
                        var = NULL,
                        hist_y = NULL,
                        hist_time = NULL,
                        hist_n = 50,
                        main = NULL,
                        ci_col68 = "grey75",
                        ci_col95 = "grey90",
                        line_col = "steelblue",
                        hist_col = "black",
                        connect_col = NULL,
                        vline_col = "grey40",
                        vline_lty = 2,
                        ...) {
  
  `%||%` <- function(a, b) if (is.null(a)) b else a
  connect_col <- connect_col %||% line_col
  stopifnot(is.list(x))
  
  # Extract metadata
  N <- x$meta$N
  h <- x$meta$h
  point_label <- x$meta$point %||% "mean"
  
  # Get forecast components
  point <- as.matrix(x$point)
  lb68  <- as.matrix(x$lb68)
  ub68  <- as.matrix(x$ub68)
  lb95  <- as.matrix(x$lb95)
  ub95  <- as.matrix(x$ub95)
  
  # --- Handle time indexing ---
  if (!is.null(hist_time)) {
    
    # Take only the last hist_n time points
    if (length(hist_time) >= hist_n) {
      hist_time <- tail(hist_time, hist_n)
    } else {
      stop("`hist_time` must be at least as long as the history (hist_n).")
    }
    
    # Generate forecast dates after the last observed time
    if (inherits(hist_time, "Date")) {
      # Determine spacing
      delta <- if (length(hist_time) > 1) median(diff(hist_time)) else 1
      # Extend by same spacing
      fc_time <- seq(from = tail(hist_time, 1) + delta,
                     by = delta, length.out = h)
      fc_time <- as.Date(fc_time, origin = "1970-01-01")
    } else if (inherits(hist_time, "POSIXt")) {
      delta <- if (length(hist_time) > 1) median(diff(hist_time)) else 1
      fc_time <- seq(from = tail(hist_time, 1) + delta,
                     by = delta, length.out = h)
    } else if (is.numeric(hist_time)) {
      delta <- if (length(hist_time) > 1) median(diff(hist_time)) else 1
      fc_time <- seq(from = tail(hist_time, 1) + delta,
                     by = delta, length.out = h)
    } else {
      warning("`hist_time` must be Date, POSIXt, numeric, or convertible; using index instead.")
      hist_time <- seq_len(hist_n)
      fc_time <- seq(max(hist_time) + 1, max(hist_time) + h)
    }
    
  } else {
    # Fallback if no time provided
    hist_time <- seq_len(hist_n)
    fc_time <- seq(max(hist_time) + 1, max(hist_time) + h)
  }
  
  # Combine history and forecast time for completeness
  time_full <- c(tail(hist_time, hist_n), fc_time)
  time_hist <- tail(hist_time, hist_n)
  
  # --- Handle historical data ---
  if (!is.null(hist_y)) {
    hist_y <- as.matrix(hist_y)
    if (ncol(hist_y) == 1 && N > 1)
      hist_y <- cbind(hist_y, matrix(NA, nrow(hist_y), N - 1))
  } else {
    hist_y <- matrix(NA_real_, nrow = hist_n, ncol = N)
  }
  
  if (is.null(var)) var <- seq_len(N)
  
  # --- ggplot2 version (default) ---
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    df_list <- vector("list", length(var))
    
    for (idx in seq_along(var)) {
      j <- var[idx]
      hist_j <- tail(hist_y[, j, drop = TRUE], hist_n)
      
      # Historical portion
      df_hist <- data.frame(
        time   = time_hist,
        value  = hist_j,
        lb68   = NA_real_,
        ub68   = NA_real_,
        lb95   = NA_real_,
        ub95   = NA_real_,
        series = paste0("y", j),
        part   = "history"
      )
      
      # Forecast portion
      df_fc <- data.frame(
        time   = fc_time,
        value  = point[j, ],
        lb68   = lb68[j, ],
        ub68   = ub68[j, ],
        lb95   = lb95[j, ],
        ub95   = ub95[j, ],
        series = paste0("y", j),
        part   = "forecast"
      )
      

      
      df_list[[idx]] <- rbind(df_hist, df_fc)
    }
    
    df <- do.call(rbind, df_list)
    
    gg <- ggplot2::ggplot(df, ggplot2::aes(x = time, y = value)) +
      ggplot2::facet_wrap(~ series, ncol = 1, scales = "free_y") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(
        title = main %||% "Forecasts",
        x = if (!is.null(hist_time)) "Time" else "Index (0 = forecast start)",
        y = paste0("Forecast (", point_label, ")")
      )
    
    # Add ribbons (95% first, then 68%)
    if (!is.null(lb95) && !is.null(ub95)) {
      gg <- gg + ggplot2::geom_ribbon(
        data = subset(df, part == "forecast"),
        ggplot2::aes(ymin = lb95, ymax = ub95),
        fill = ci_col95, alpha = 0.6
      )
    }
    if (!is.null(lb68) && !is.null(ub68)) {
      gg <- gg + ggplot2::geom_ribbon(
        data = subset(df, part == "forecast"),
        ggplot2::aes(ymin = lb68, ymax = ub68),
        fill = ci_col68, alpha = 0.6
      )
    }
    
    # Lines
    gg <- gg +
      ggplot2::geom_line(
        data = subset(df, part == "forecast"),
        ggplot2::aes(y = value),
        color = line_col, linewidth = 1.1
      ) +
      ggplot2::geom_line(
        data = subset(df, part == "history"),
        ggplot2::aes(y = value),
        color = hist_col, linewidth = 0.8
      )
    
    # Connection line (forecast start)
    if (any(df$part == "history")) {
      connectors <- do.call(rbind, lapply(split(df, df$series), function(dsi) {
        dh <- subset(dsi, part == "history")
        dfc <- subset(dsi, part == "forecast")
        if (nrow(dh) == 0 || nrow(dfc) == 0) return(NULL)
        data.frame(
          series = unique(dsi$series),
          x = tail(dh$time, 1),
          xend = dfc$time[1],
          y = tail(dh$value, 1),
          yend = dfc$value[1]
        )
      }))
      gg <- gg + ggplot2::geom_segment(
        data = connectors,
        ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
        linewidth = 1, color = connect_col
      )
    }
    
    if (inherits(hist_time, "Date")) {
      gg <- gg + ggplot2::scale_x_date(date_labels = "%Y", date_breaks = "1 year")
    }
    
    gg <- gg + ggplot2::geom_vline(xintercept = tail(fc_time, h + 1)[1],
                                   color = vline_col, linetype = vline_lty)
    
    return(gg)
  }
  
  # --- Base R fallback ---
  K <- length(var)
  op <- par(mfrow = c(K, 1), mgp = c(2, 0.8, 0), mar = 0.1 + c(3, 3, 2, 1))
  on.exit(par(op), add = TRUE)
  
  for (j in var) {
    hist_j <- tail(hist_y[, j, drop = TRUE], hist_n)
    y_point <- point[j, ]
    y_lb68 <- lb68[j, ]; y_ub68 <- ub68[j, ]
    y_lb95 <- lb95[j, ]; y_ub95 <- ub95[j, ]
    yrng <- range(c(hist_j, y_point, y_lb95, y_ub95), na.rm = TRUE)
    
    plot(seq_along(hist_j), hist_j, type = "l", col = hist_col, ylim = yrng,
         main = paste0("y", j), ylab = paste0("Forecast (", point_label, ")"),
         xlab = "Time", ...)
    polygon(c(seq_len(h), rev(seq_len(h))), c(y_lb95, rev(y_ub95)),
            col = adjustcolor(ci_col95, alpha.f = 0.4), border = NA)
    polygon(c(seq_len(h), rev(seq_len(h))), c(y_lb68, rev(y_ub68)),
            col = adjustcolor(ci_col68, alpha.f = 0.4), border = NA)
    lines(seq_len(h), y_point, col = line_col, lwd = 2)
    abline(v = hist_n, lty = vline_lty, col = vline_col)
  }
  
  invisible(NULL)
}