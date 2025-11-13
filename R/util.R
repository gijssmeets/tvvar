# ---- helper: parameter indexing & names (must match your packing order) ----
.tvvar_param_index <- function(fit) {
  meta <- fit$meta
  N    <- meta$N; p <- meta$p; r <- meta$r
  zm   <- isTRUE(meta$zero.mean)
  arr  <- meta$Phi.f.array    # N x (1+N*p) x r
  
  # mask for free Phi_f entries (same flattening used in packing)
  Phi.f.mat  <- matrix(arr, nrow = N)
  mask       <- (Phi.f.mat != 0)
  n_phi_free <- sum(mask)
  
  Nstar <- N * (N + 1L) / 2L
  needC <- if (zm) N^2 * p else (N^2 * p + N)
  
  # names: A, B, phi_r
  names_head <- c("A", "B")
  if (r > 0) names_head <- c(names_head, paste0("phi_r[", seq_len(r), "]"))
  
  # names: Phi_f free entries
  # column mapping: for col c in 1:( (1+N*p)*r ), j = ((c-1) %% (1+N*p)) + 1, k = floor((c-1)/(1+N*p)) + 1
  build_names_Phi_f <- function() {
    nm <- character(0)
    if (n_phi_free == 0) return(nm)
    idx_cols <- which(colSums(mask) >= 0)  # just to get range
    total_cols <- ncol(mask)
    # loop over all matrix entries; add names only where mask TRUE
    for (c in seq_len(total_cols)) {
      j <- ((c - 1L) %% (1 + N * p)) + 1L  # regressor column 1..(1+N*p)
      k <- (c - 1L) %/% (1 + N * p) + 1L   # factor slice 1..r
      for (i in seq_len(N)) {
        if (mask[i, c]) {
          if (j == 1L) {
            reglab <- "const"
          } else {
            off <- j - 1L
            var <- ((off - 1L) %% N) + 1L
            lag <- (off - 1L) %/% N + 1L
            reglab <- sprintf("y%d_lag%d", var, lag)
          }
          nm <- c(nm, sprintf("Phi_f[%d,%s,f%d]", i, reglab, k))
        }
      }
    }
    nm
  }
  
  names_phi_f <- build_names_Phi_f()
  
  # names: vech(L) entries (Cholesky factor for Omega).  We just label by (row,col).
  L_names <- {
    out <- character(0)
    for (j in seq_len(N)) {
      for (i in j:N) out <- c(out, sprintf("L[%d,%d]", i, j))
    }
    out
  }
  
  # names: Phi_c vector (by row-major unpacking consistent with your packing)
  # Phi_c is N x (1 + N*p) if zero.mean=FALSE, or cbind(0, N x N*p) if zero.mean=TRUE.
  build_names_Phi_c <- function() {
    nm <- character(0)
    if (zm) {
      # no intercept in params; just lagged regressors
      for (i in seq_len(N)) {
        for (lag in seq_len(p)) {
          for (var in seq_len(N)) {
            nm <- c(nm, sprintf("Phi_c[%d,y%d_lag%d]", i, var, lag))
          }
        }
      }
    } else {
      # intercept first, then lagged
      for (i in seq_len(N)) nm <- c(nm, sprintf("Phi_c[%d,const]", i))
      for (i in seq_len(N)) {
        for (lag in seq_len(p)) {
          for (var in seq_len(N)) {
            nm <- c(nm, sprintf("Phi_c[%d,y%d_lag%d]", i, var, lag))
          }
        }
      }
    }
    nm
  }
  names_phi_c <- build_names_Phi_c()
  
  # indices for slicing in theta (match your packing in opti.fct / unpen code)
  n_head      <- 2L + r
  i_head      <- seq_len(n_head)
  i_phi_f     <- if (n_phi_free > 0) (max(i_head) + 1L):(max(i_head) + n_phi_free) else integer(0)
  i_L         <- if (Nstar > 0) (if (length(i_phi_f)) max(i_phi_f) else max(i_head)) + 1L
  i_L         <- if (Nstar > 0) i_L:(i_L + Nstar - 1L) else integer(0)
  i_phi_c     <- if (needC > 0) (if (length(i_L)) max(i_L) else if (length(i_phi_f)) max(i_phi_f) else max(i_head)) + 1L
  i_phi_c     <- if (needC > 0) i_phi_c:(i_phi_c + needC - 1L) else integer(0)
  
  list(
    idx = list(head = i_head, phi_f = i_phi_f, L = i_L, phi_c = i_phi_c),
    names = c(names_head, names_phi_f, L_names, names_phi_c)
  )
}


#' Summarize a fitted tvfit model
#'
#' @description
#' Provides formatted parameter estimates, standard errors, z-values, and p-values
#' for a fitted \code{tvfit} object. Parameters are reported in their original order,
#' ensuring continuous numbering across blocks.
#'
#' @param fit A fitted \code{tvfit} object.
#' @param digits Number of digits to display.
#' @param print Logical; if TRUE, prints the summary to console.
#' @return Invisibly returns a list with:
#'   \itemize{
#'     \item \code{info}: model information and ICs
#'     \item \code{params}: full parameter table
#'     \item \code{blocks}: block indices (for reference)
#'   }
#' @export
summary.tvfit <- function(fit, digits = 3, print = TRUE) {
  stopifnot(is.list(fit), !is.null(fit$meta))
  
  theta <- fit$theta %||% fit$optim$par
  V     <- fit$vcov
  if (is.list(V) && !is.null(V$untransformed)) V <- V$untransformed
  
  map  <- .tvvar_param_index(fit)
  npar <- length(map$names)
  
  # align lengths
  if (!length(theta) || length(theta) != npar) {
    theta <- theta[seq_len(min(length(theta), npar))]
    if (length(theta) < npar)
      theta <- c(theta, rep(NA_real_, npar - length(theta)))
  }
  
  have_V <- is.matrix(V) && all(dim(V) == c(length(theta), length(theta)))
  se_raw <- if (have_V) sqrt(pmax(diag(V), 0)) else rep(NA_real_, length(theta))
  
  est_nat <- theta
  se_nat  <- se_raw
  
  # Transform head scalars to (0,1)
  idx <- map$idx
  if (length(idx$head)) {
    inv_logit <- stats::plogis
    h <- idx$head
    p <- inv_logit(theta[h])
    est_nat[h] <- p
    if (have_V) {
      grad <- p * (1 - p)
      se_nat[h] <- se_raw[h] * grad
    }
  }
  
  # z/p-values
  z   <- if (all(!is.na(se_nat))) est_nat / se_nat else rep(NA_real_, length(theta))
  pvl <- if (all(!is.na(z))) 2 * stats::pnorm(-abs(z)) else rep(NA_real_, length(theta))
  
  # unified parameter table (with continuous numbering)
  params_df <- data.frame(
    parameter = seq_along(map$names),
    name      = map$names,
    estimate  = est_nat,
    std.error = se_nat,
    z.value   = z,
    p.value   = pvl,
    stringsAsFactors = FALSE
  )
  
  # use indices for grouping but keep order from map$names
  blocks <- list(
    scalars = idx$head,
    Phi_c   = idx$phi_c,
    Phi_f   = idx$phi_f,
    L_vech  = idx$L
  )
  
  ic <- fit$ic %||% list()
  meta <- fit$meta %||% list()
  info_df <- data.frame(
    metric = c("Method","N","p","r","zero.mean","nobs","AIC","AICc","BIC","Runtime_sec"),
    value  = c(meta$method %||% NA, meta$N %||% NA, meta$p %||% NA, meta$r %||% NA,
               meta$zero.mean %||% NA, meta$nobs %||% NA,
               ic$AIC %||% NA, ic$AICc %||% NA, ic$BIC %||% NA, meta$time %||% NA),
    row.names = NULL
  )
  
  if (isTRUE(print)) {
    if (requireNamespace("knitr", quietly = TRUE)) {
      cat("## Model info\n")
      print(knitr::kable(info_df, digits = digits, align = "l"))
      cat("\n## Parameter estimates\n")
      print(knitr::kable(params_df, digits = digits))
    } else {
      message("knitr not found; printing with base::print()")
      print(info_df)
      print(params_df)
    }
  }
  
  invisible(list(info = info_df, params = params_df, blocks = blocks))
}




initialize_omega <- function(N, diag_val = 0.3, offdiag_val = 0.2) {
  M <- matrix(offdiag_val, N, N)
  diag(M) <- diag_val
  M
}

# Ensures the phi_f_structure has the correct 3D shape [N x N x r]
.normalize_phi_structure <- function(phi_f_structure, N, r) {
  # Case 1: NULL or scalar — default to ones
  if (is.null(phi_f_structure)) {
    return(array(1, dim = c(N, N, r)))
  }
  
  # Case 2: matrix input (common case, e.g., matrix(1, 2, 2))
  if (is.matrix(phi_f_structure)) {
    # replicate the same matrix across r slices
    return(array(phi_f_structure, dim = c(N, N, r)))
  }
  
  # Case 3: already a 3D array
  if (length(dim(phi_f_structure)) == 3L) {
    dims <- dim(phi_f_structure)
    if (!all(dims[1:2] == N))
      stop(sprintf("phi_f_structure must have dims [%d,%d,*], got [%d,%d,%d]",
                   N, N, dims[1], dims[2], dims[3]))
    if (dims[3] != r)
      stop(sprintf("phi_f_structure has %d factors but r=%d.", dims[3], r))
    return(phi_f_structure)
  }
  
  # Case 4: list of matrices
  if (is.list(phi_f_structure)) {
    if (length(phi_f_structure) != r)
      stop("If phi_f_structure is a list, it must have length r.")
    arr <- array(NA, dim = c(N, N, r))
    for (i in seq_len(r)) arr[,,i] <- phi_f_structure[[i]]
    return(arr)
  }
  
  stop("Invalid phi_f_structure: must be matrix, 3D array, list of matrices, or NULL.")
}

initialize_parameters <- function(VAR.data,
                                  p,
                                  r,
                                  zero.mean = TRUE,
                                  phi_f_structure,
                                  init = c("default", "random", "custom"),
                                  init_list = NULL) {
  init <- match.arg(init)
  VAR.data <- as.matrix(VAR.data)
  N  <- ncol(VAR.data)
  
  # --- Build Φ_f structure and mask ---
  phi_arr      <- .normalize_phi_structure(phi_f_structure, N, r)
  Phi.f.array  <- make.Phi.f(structure = phi_arr, lag.order = p)   # N x (N*p+1) x r
  Phi.f.as.mat <- matrix(Phi.f.array, nrow = N)                     # N x ((N*p+1)*r)
  Phi.f.mask   <- (Phi.f.as.mat != 0)
  n_phi_free   <- sum(Phi.f.mask)
  
  # --- Baseline defaults ---
  A0     <- 0.10
  B0     <- 0.80
  phi_r0 <- if (r > 0) rep(0.95, r) else numeric(0L)
  Omega0 <- matrix(0.2, N, N); diag(Omega0) <- 0.3   # 0.3 diag / 0.2 off-diag
  
  # Φ_c via static VAR (as before)
  Y_named <- as.data.frame(VAR.data); colnames(Y_named) <- paste0("y", seq_len(N))
  B.t <- matrix(0, nrow = N, ncol = (N * p + 1))
  for (j in 1:N) B.t[j, ] <- coefficients(VAR(Y_named, p = p))[[j]][, 1]
  A.t <- B.t[, -(N * p + 1), drop = FALSE]                 # drop intercept
  if (isFALSE(zero.mean)) A.t <- cbind(B.t[, (N * p + 1), drop = FALSE], A.t)
  Phi_c0 <- A.t
  
  # --- Init modes ---
  if (init == "random") {
    set.seed(123)
    # A, B, phi_r
    A0     <- runif(1, 0.05, 0.15)
    B0     <- runif(1, 0.75, 0.85)
    phi_r0 <- if (r > 0) runif(r, 0.90, 0.99) else numeric(0L)
    # Ω
    Omega0 <- diag(runif(N, 0.08, 0.15))
    # Φ_c small noise
    Phi_c0 <- Phi_c0 + 0.01 * matrix(rnorm(length(Phi_c0)), nrow = nrow(Phi_c0))
  } else if (init == "custom") {
    if (!is.null(init_list)) {
      if (!is.null(init_list$A))     A0     <- init_list$A
      if (!is.null(init_list$B))     B0     <- init_list$B
      if (!is.null(init_list$phi_r)) phi_r0 <- init_list$phi_r
      if (!is.null(init_list$Omega)) Omega0 <- init_list$Omega
      if (!is.null(init_list$Phi_c)) Phi_c0 <- init_list$Phi_c
      if (!is.null(init_list$Phi_f)) {
        stopifnot(all(dim(init_list$Phi_f) == dim(Phi.f.array)))
        Phi.f.array  <- init_list$Phi_f
        Phi.f.as.mat <- matrix(Phi.f.array, nrow = N)
        Phi.f.mask   <- (Phi.f.as.mat != 0)
        n_phi_free   <- sum(Phi.f.mask)
      }
    }
  }
  # Φ_f free entries start values
  Phi_f_start <- numeric(n_phi_free)
  if (n_phi_free > 0) {
    if (init == "random") {
      # small random starts for free Φ_f entries
      Phi_f_start <- rnorm(n_phi_free, mean = 0, sd = 0.05)
    } else if (init == "custom" && !is.null(init_list) && !is.null(init_list$Phi_f)) {
      Phi_f_start <- Phi.f.as.mat[Phi.f.mask]   # from provided Φ_f
    } else {
      Phi_f_start[] <- 0.01                     # default small positive
    }
  }
  
  # pack (A,B,phi_r) transformed; Φ_f free; vech(Ω); vec(Φ_c)
  logit <- function(x) log(x/(1 - x))
  head_vec <- if (r > 0) c(logit(A0), logit(B0), atanh(phi_r0)) else c(logit(A0), logit(B0))
  par_init <- c(head_vec, Phi_f_start, vech(Omega0), as.vector(Phi_c0))
  
  list(
    par_init   = par_init,
    Phi.f.array = Phi.f.array,
    Phi.f.array.mat.structure = Phi.f.as.mat,
    Phi_c0     = Phi_c0,
    Omega0     = Omega0,
    A0 = A0, B0 = B0, phi_r0 = phi_r0,
    n_phi_free = n_phi_free,
    Phi_f_mask = Phi.f.mask
  )
}



calc_indices_of_lambda_diagonals <- function(dim.VAR, lag.order, number.factors){
  with_zeros <- rep(0,number.factors)
  without_zeros <- rep(0, number.factors)
  
  count_with <- 0
  count_without <- 0

  for(i in (1:number.factors)){
    with_zeros[i] <- count_with + 1
    without_zeros[i] <- count_without + 1
    
    count_with <- count_with + dim.VAR^2*lag.order + 1
    count_without <- count_without + dim.VAR^2*lag.order + 1 - i
  }
  
  return(list(
    with_zeros = with_zeros,
    without_zeros = without_zeros
  ))
}

params_to_list <- function(params, Phi.f.array, cfg){
  
  dim.VAR = cfg$dim.VAR
  number.factors = cfg$number.factors
  lag.order = cfg$lag.order
  
  A = params[1]
  B = params[2]
  pphi = params[3:(3+number.factors-1)]
  
  #Phi.f.as.matrix <- matrix(Phi.f.array, nrow=dim.VAR)
  Phi.f.as.matrix <- Phi.f.array
  count.params <- 2+number.factors+length(which(Phi.f.as.matrix != 0))

  
  Phi.f <- params[(2+number.factors+1):count.params]
  
  
  Nstar <- dim.VAR*(dim.VAR+1)/2
  momega <- params[(count.params+1):(count.params+Nstar)]

  count.params = count.params+Nstar
  
  Phi.c <- params[(count.params+1):(count.params+dim.VAR^2*lag.order)]
  
  #if(zero.mean == TRUE){
  #  Phi.c <- cbind(rep(0,dim.VAR),matrix(params[(count.params+1):(count.params+dim.VAR^2*lag.order)], nrow=dim.VAR))
  #}else{
  #  Phi.c <- params[(count.params+1):(count.params+dim.VAR^2*lag.order)]
  #}
  
  parlist <- list(A = A, B = B, pphi = pphi, Phi.f = Phi.f, momega = momega, Phi.c = Phi.c)
  
  return(parlist)
}

merge_conditional_params <- function(params, Phi.f.array, cfg, pphi = NaN, Phi.c = NaN){
  parlist <- params_to_list(params, Phi.f.array, cfg)
  
  if(!all(is.nan(pphi))){
    parlist$pphi <- pphi
  }
  if (!all(is.nan(Phi.c))) {
    parlist$Phi.c <- Phi.c
  }
  
  par <- list_to_params(parlist, Phi.f.array)
  return(par)
}

## construct Q-matrix (selection matrix)
Q.matrix.fun <- function(r,N,p){
  aid.vec <- matrix(1:r, ncol=1)
  aid.mat <- matrix(aid.vec %x% diag(N*p+1))
  Qmat <- matrix(0, nrow=r*(N*p+1)^2, ncol=r)
  for(i in 1:r){
    Qmat[(which(aid.mat == i)),i]=1
  }
  Qmat 
}

make.Phi.f <- function(structure, lag.order){
  dim.VAR <- dim(structure)[1]
  number.factors <- dim(structure)[3]
  Phi.f <- array(0, c(dim.VAR,(dim.VAR*lag.order+1), number.factors))
  for(j in 1:number.factors){
    Phi.f[,,j] = cbind(rep(0,dim.VAR), matrix(rep(structure[,,j], lag.order), nrow=dim.VAR))
  }
  return(Phi.f)
}

build.statespace.form <- function(VAR.data, lag.order, number.factors, Phi.f, Phi.c) {
  dim.VAR <- ncol(VAR.data)
  
  ## regressors: (Np x T)-matrix
  T.fin <- nrow(VAR.data)
  Y <- t(VAR.data)
  
  ## calls c++ function
  Y.minus1 =create_Y_minus1(Y, lag.order, T.fin, dim.VAR)
  
  ## Adjust Y to remove initial lags
  Y <- Y[, -(1:lag.order)]
  
  ## Dependent variable for state space form
  ytilde <- Y - Phi.c %*% Y.minus1
  ## Build array of Z-matrices
  Q <- Q.matrix.fun(r = number.factors, N = dim.VAR, p = lag.order)
  Z <- create_Z(Y.minus1, Phi.f, Q, T.fin, dim.VAR, number.factors, lag.order) 
  ## Return the result as a list
  list(Z = Z, ytilde = ytilde, Y.minus1 = Y.minus1, Q = Q)
}

list_to_params <- function(parlist, Phi.f.array){
  A = parlist$A
  B = parlist$B
  pphi = parlist$pphi
  Phi.f = parlist$Phi.f
  momega = parlist$momega
  Phi.c = parlist$Phi.c
  
  
  params = 
    c(A,
      B,
      pphi, 
      Phi.f, 
      momega,
      Phi.c
    )
  
  
  return(params)
}

divide_params <- function(params, pphi, cfg, Phi.c, Phi.f, Phi.f.array){
  
  result <- new.env()
  fixed_par <- list()
  if(!all(is.nan(pphi))){
    fixed_par$pphi <- pphi
  }
  
  if(!all(is.nan(Phi.c))){
    fixed_par$Phi.c <- Phi.c
  }
  
  if(!all(is.nan(Phi.f))){
    fixed_par$Phi.f <- Phi.f
  }
  
  fixed_names <- names(fixed_par)
  parlist = params_to_list(params, Phi.f.array, cfg)
  
  # List all names that are not fixed
  free_names <- names(parlist)[!names(parlist) %in% fixed_names]
  
  # Create a list of free parameters
  free_par <- list()
  for (name in free_names) {
    free_par[[name]] <- parlist[[name]]
  }
  
  result$par_free <- free_par
  result$par_fixed <- fixed_par
  return(result)
}

logistic <- function(x) {
  return(exp(x) / (1 + exp(x)))
}

rebuild_params <- function(par_free, par_fixed) {
  parlist <- list()
  
  # Always rebuild A and B
  parlist$A <- par_free[["A"]]
  parlist$B <- par_free[["B"]]
  
  # Always rebuild momega
  momega_idx <- grep("^momega", names(par_free))
  parlist$momega <- unname(unlist(par_free[momega_idx]))
  
  # Conditionally rebuild Phi.c
  if (!"Phi.c" %in% names(par_fixed)) {
    phi_c_idx <- grep("^Phi\\.c", names(par_free))
    parlist$Phi.c <- unname(unlist(par_free[phi_c_idx]))
  }
  
  # Conditionally rebuild Phi.f
  if (!"Phi.f" %in% names(par_fixed)) {
    phi_f_idx <- grep("^Phi\\.f", names(par_free))
    parlist$Phi.f <- unname(unlist(par_free[phi_f_idx]))
  }
  
  # Conditionally rebuild pphi
  if (!"pphi" %in% names(par_fixed)) {
    pphi_idx <- grep("^pphi", names(par_free))
    parlist$pphi <- unname(unlist(par_free[pphi_idx]))
  }
  
  all_par <- c(parlist, par_fixed)
  params <- list_to_params(all_par, Phi.f.array)
  
  return(params)
}


transform_parameters <- function(opt_sim, cfg_sim, sim){
  opt_sim$par[1:2] <- logistic(opt_sim$par[1:2])
  
  opt_sim$par[3:(2+cfg_sim$number.factors)] <- tanh(opt_sim$par[3:(2+cfg_sim$number.factors)])
  
  count.params <- 2+cfg_sim$number.factors+length(which(matrix(sim$Phi.f.array, nrow=cfg_sim$dim.VAR) != 0))
  Nstar <- cfg_sim$dim.VAR*(cfg_sim$dim.VAR+1)/2
  m1 <- matrix(D.matrix(cfg_sim$dim.VAR) %*% opt_sim$par[(count.params+1):(count.params+Nstar)], ncol = cfg_sim$dim.VAR)
  m1[upper.tri(m1)] <- 0
  momega <- m1 %*% t(m1)
  opt_sim$par[(count.params+1):(count.params+Nstar)] <- momega[lower.tri(momega, diag = TRUE)]
  
  return(opt_sim)
}


optim_Phi_c <- function(opti.eval.em, data, cfg, lambda, weights){
  f <- opti.eval.em$smoothed.state
  last.obs <- length(f[1,]) 
  Y <- t(data)
  Y <- Y[, -(1:cfg$lag.order)] # removes first lag.order observations (as we need Y minus 1 to be this)
  Y <- Y[,-last.obs, drop=FALSE] # removes the last observations as we don't have predictions
  
  Y.minus1 <- opti.eval.em$Y.minus1
  Y.minus1 <- Y.minus1[-1,]
  Y.minus1 <- Y.minus1[,-last.obs, drop=FALSE]

  Z <- opti.eval.em$Z

  T.loop <- length(Y[1,])-1
  
  y <- matrix(0, nrow = nrow(Y), ncol = ncol(Y)) 
  for (t in seq_len(T.loop)) {
    
    if (cfg$number.factors == 1){
      Z.t <- Z[, , t, drop=FALSE]
    } else{
      Z.t <- Z[, , t]
    }
    f.t <- f[, t, drop = FALSE]
    y[, t] <- Y[, t] - Z.t %*% f.t
  }
  
  X <- t(Y.minus1)
  y <- t(y)
  
  # OLS solution
  Phi.c.OLS <- solve(t(X) %*% X) %*% t(X) %*% y
  Phi.c.OLS <- as.vector(t(Phi.c.OLS))
  

  return(list(Phi.c.OLS = Phi.c.OLS
              ))
  
}

optim_pphi_c <- function(opti.eval.em){
  if (length(opti.eval.em$pphi_ols) == 1) {
    pphi <- matrix(atanh(opti.eval.em$pphi_ols), nrow = 1, ncol = 1)
  } else {
    pphi <- atanh(opti.eval.em$pphi_ols)

  }
  return(pphi)
}


optim_Phi_f_full_likelihood <- function(opti.eval.em, params, data, Phi.f.array, cfg, lambda, type = 'regular', weights) {
  smoothed_factors <- opti.eval.em$smoothed.state
  V_array <- opti.eval.em$smoothed.state.variance
  
  current_par_list <- params_to_list(params, Phi.f.array, cfg)
  Phi.c_current <- cbind(rep(0,cfg$dim.VAR),matrix(current_par_list$Phi.c, nrow=cfg$dim.VAR))
  
  loglik_wrapper <- function(phi_f_vec) {
      
      # Get the current full list of parameters
      parlist <- params_to_list(params, Phi.f.array, cfg)
    
      # Substitute the trial Phi.f vector
      parlist$Phi.f <- phi_f_vec
      
      # Rebuild the parameter vector that opti.fct expects
      par_free_for_opti <- parlist
      par_fixed_for_opti <- NaN
      
      result <- opti.fct(
        par_free = par_free_for_opti,
        par_fixed = par_fixed_for_opti,
        VAR.data = data,
        Phi.f.array = Phi.f.array,
        cfg = cfg,
        Smooth = FALSE,
        purpose = "optim" 
      )
      
      return(result)
    }
    
  grad_for_apg <- function(phi_vec, opts) {
    return(numDeriv::grad(func = loglik_wrapper, x = phi_vec, method = "Richardson"))
  }
  
  prox_l1_weighted <- function(x, t, opts) {
    lambda <- opts$lambda # single value
    weights <- opts$weights # vector of weights

    thres <- t * lambda * weights 
    idx.1 <- which(x < -thres)
    idx.2 <- which(x > thres)
    res <- rep(0,length(x))
    if ( length(idx.1)>0 ) res[idx.1] <- x[idx.1] + thres[idx.1]
    if ( length(idx.2)>0 ) res[idx.2]<- x[idx.2] - thres[idx.2]
    return(res)
  }
  
  Phi.f.init <- current_par_list$Phi.f
  
  # Regular LASSO
  if(type == 'regular'){
    opts <- list(X_INIT = Phi.f.init, lambda = lambda, MAX_ITERS = 100, EPS = 1e-3)
    prox_grad_result <- apg::apg(grad_for_apg, apg::prox.l1, length(Phi.f.init), opts)
    phi_est_penalized <- prox_grad_result$x
  }
  
  # Adaptive LASSO
  if(type == 'adaptive'){
    step_size <- 0.1 
    opts_weighted <- list(X_INIT = Phi.f.init, lambda = lambda, MAX_ITERS = 100, EPS = 1e-3, weights = weights, FIXED_STEP_SIZE = TRUE, STEP_SIZE = step_size)
    prox_grad_result_weighted <- apg::apg(grad_for_apg, prox_l1_weighted, length(Phi.f.init), opts_weighted)
    phi_est_penalized <- prox_grad_result_weighted$x
  }
  

  active_idx <- which(abs(phi_est_penalized) != 0)
  
  if (length(active_idx) == 0) {
    return(list(phi_est_penalized = phi_est_penalized))
  }
  
  # Refit
  refit_objective <- function(active_params) {
    phi_full <- numeric(length(phi_est_penalized))
    phi_full[active_idx] <- active_params
    return(loglik_wrapper(phi_full))
  }
  
  phi_est_refitted <- phi_est_penalized 
  refit <- optim(par = phi_est_penalized[active_idx], 
                            fn = refit_objective, 
                            method = "BFGS")
  
  phi_refitted <- numeric(length(phi_est_penalized))
  phi_est_refitted[active_idx] <- refit$par
  
  phi_est_penalized <- phi_est_refitted
  
  return(list(
    phi_est_penalized = phi_est_penalized))
}

initialize_identified_phi <- function(dim.VAR, number.factors, lag.order) {
  
  # Get key dimensions
  r <- number.factors
  N <- dim.VAR
  p <- lag.order
  
  Delta_r <- matrix(0, nrow = r, ncol = r)
  Delta_r[lower.tri(Delta_r, diag = TRUE)] <- 1
  
  full_delta_rows <- N * (N * p)
  Delta_full <- matrix(1, nrow = full_delta_rows, ncol = r)
  Delta_full[1:r, ] <- Delta_r
  
  Phi_f_3D_no_intercept <- array(0, dim = c(N, N * p, r))
  for (j in 1:r) {
    Phi_f_3D_no_intercept[, , j] <- matrix(Delta_full[, j], nrow = N, ncol = N * p)
  }
  
  final_cols <- N * p + 1
  Phi_f_array <- array(0, dim = c(N, final_cols, r))
  
  for (j in 1:r) {
    Phi_f_array[, , j] <- cbind(rep(0, N), Phi_f_3D_no_intercept[, , j])
  }
  
  
  Phi_f_vector <- as.vector(Phi_f_array)
  Phi_f_structure_matrix <- matrix(Phi_f_vector, nrow = dim.VAR)
  
  return(
    list(
      Phi_f_array = Phi_f_array,
      Phi_f_structure_matrix = Phi_f_structure_matrix,
      Phi_f_vector = Phi_f_vector,
      Delta_r = Delta_r
    )
  )
}

optim_pphi_bfgs <- function(opti.eval.em, cfg) {

  V0_star  <- opti.eval.em$V_0_star
  V1_star  <- opti.eval.em$V_minus1_star
  V01_star <- opti.eval.em$V_0_minus1_star
  
  T_minus_1 <- cfg$T.here - 1
  
  q_function_pphi <- function(pphi_untransformed) {
    
    if (cfg$number.factors > 1){
      pphi <- tanh(diag(pphi_untransformed))
      if (any(abs(pphi) >= 1.0)) return(Inf) 
      Sigma_eta <- diag(cfg$number.factors) - pphi %*% t(pphi)
      log_lik_part1 <- - (T_minus_1 / 2) * log(det(Sigma_eta))
      quadratic_part <- V0_star - V01_star %*% t(pphi) - pphi %*% t(V01_star) + pphi %*% V1_star %*% t(pphi)
      log_lik_part2 <- - 0.5 * sum(diag(solve(Sigma_eta) %*% quadratic_part))
      
      }
    else{
      pphi <- tanh(pphi_untransformed)
      if (abs(pphi) >= 1.0) return(Inf)
      log_lik_part1 <- - (T_minus_1 / 2) * log(1 - pphi^2)
      quadratic_part <- V0_star - 2 * pphi * V01_star + pphi^2 * V1_star
      log_lik_part2 <- - (0.5 / (1 - pphi^2)) * quadratic_part
    }
    return(-(log_lik_part1 + log_lik_part2))
  }

  initial_pphi_ols <- opti.eval.em$pphi_ols
  initial_pphi_ols_clipped <- pmin(pmax(initial_pphi_ols, -0.999), 0.999)
  start_val <- atanh(initial_pphi_ols_clipped)
  
  # Use optim with BFGS to find the maximum of the Q-function
  optim_result <- optim(
    par = start_val,
    fn = q_function_pphi,
    method = "BFGS"
  )
  
  if(cfg$number.factors == 1){
    return(matrix(optim_result$par, nrow=1, ncol=1))
  }
  else if (cfg$number.factors > 1) {
    return(optim_result$par)
  }
}


#' @keywords internal
#' @noRd
check_identification <- function(factor.structure, dim.VAR, number.factors, lag.order) {
  N <- dim.VAR
  r <- number.factors
  p <- lag.order
  
  # --- Step 1: Convert to 3D array form [N x (N*p) x r] ---
  if (length(dim(factor.structure)) == 2L) {
    # e.g., a flattened N x (N*p*r) matrix
    if (ncol(factor.structure) != N * p * r)
      stop("If 2D, factor.structure must have N * (N*p*r) columns.")
    factor.structure <- array(factor.structure, dim = c(N, N * p, r))
  } else if (length(dim(factor.structure)) != 3L) {
    stop("`factor.structure` must be either a 3D array [N x (N*p) x r] or equivalent 2D form.")
  }
  
  # --- Step 2: Build Λ = [vec(Phi_1^f), ..., vec(Phi_r^f)] ---
  Lambda <- matrix(NA, nrow = N * N * p, ncol = r)
  for (j in seq_len(r)) {
    Lambda[, j] <- as.vector(factor.structure[, , j])
  }
  
  # --- Step 3: Extract the identification block Λ₀ (first r rows) ---
  Lambda_0 <- Lambda[1:r, , drop = FALSE]
  
  # --- Step 4: Check identification condition up to permutation ---
  perms <- gtools::permutations(r, r)
  identified <- FALSE
  successful_perm <- NULL
  
  for (i in seq_len(nrow(perms))) {
    permuted <- Lambda_0[, perms[i, ]]
    
    # Must be lower triangular (zeros strictly above diag)
    lower_ok <- all(permuted[upper.tri(permuted)] == 0)
    
    # Diagonal must be strictly positive
    diag_ok <- all(diag(permuted) > 0)
    
    if (lower_ok && diag_ok) {
      identified <- TRUE
      successful_perm <- perms[i, ]
      break
    }
  }
  
  # --- Step 5: Output and messages ---
  if (identified) {
    message("✅ Identification check passed: Λ₀ can be permuted into lower-triangular form ",
            "(permutation: ", paste(successful_perm, collapse = ","), ").")
    return(TRUE)
  } else {
    warning("The supplied factor.structure cannot be permuted into an identifiable form.\n",
            "Ensure there are enough zeros to construct a lower-triangular Λ₀ with positive diagonal.")
    return(FALSE)
  }
}

#' @keywords internal
#' @noRd
init.factor.structure <- function(dim.VAR, number.factors, lag.order) {
  
  # Get key dimensions
  r <- number.factors
  N <- dim.VAR
  p <- lag.order
  
  # --- Step 1: Create the target Delta_r matrix with 1s and 0s ---
  # This enforces the lower-triangular identification scheme.
  Delta_r <- matrix(0, nrow = r, ncol = r)
  Delta_r[lower.tri(Delta_r, diag = TRUE)] <- 1
  
  
  # --- Step 2: Build the full Delta matrix for the lag coefficients ---
  # The dimensions here are only for the lag coefficients, not the intercept.
  full_delta_rows <- N * (N * p)
  Delta_full <- matrix(1, nrow = full_delta_rows, ncol = r)
  Delta_full[1:r, ] <- Delta_r
  
  
  # --- Step 3: Reshape Delta back into the 3D Phi array (without intercept) ---
  Phi_f_3D_no_intercept <- array(0, dim = c(N, N * p, r))
  for (j in 1:r) {
    Phi_f_3D_no_intercept[, , j] <- matrix(Delta_full[, j], nrow = N, ncol = N * p)
  }
  
  
  # --- Step 4: Create the final 3D array including the zero intercept column ---
  # This mimics the formatting behavior of your `make.Phi.f` function.
  final_cols <- N * p + 1
  Phi_f_array <- array(0, dim = c(N, final_cols, r))
  
  for (j in 1:r) {
    # The first column is the restricted intercept (all zeros).
    # The remaining columns are filled with the lag coefficients from the structure.
    Phi_f_array[, , j] <- cbind(rep(0, N), Phi_f_3D_no_intercept[, , j])
  }
  
  
  # --- Step 5: Generate the other required output formats ---
  
  # Flatten the 3D array into a single vector
  Phi_f_vector <- as.vector(Phi_f_array)
  
  # Create the 2D matrix version, which matches your 'Phi.f.array.mat.structure'
  Phi_f_structure_matrix <- matrix(Phi_f_vector, nrow = dim.VAR)
  
  
  # --- Step 6: Return all formatted structures in a list ---
  return(
    list(
      Phi_f_array = Phi_f_array,
      Phi_f_structure_matrix = Phi_f_structure_matrix,
      Phi_f_vector = Phi_f_vector,
      Delta_r = Delta_r
    )
  )
}
