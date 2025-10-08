#' Unpenalized estimation (ML or EM) for TV-VAR
#' @export
unpenalized_estimate <- function(data, p, r, zero.mean = TRUE,
                                 phi_f_structure,
                                 method = c("ML","EM"),
                                 control = list(maxit = 2000, trace = 0),
                                 em_control = list(max_iter = 200, tol = 1e-3, trace = TRUE)) {
  method <- match.arg(method)
  
  #Y <- .ensure_TxN(as.matrix(data), fun = "unpenalized_estimate")
  
  
  if (identical(method, "ML")) {
    return(.unpen_ml(data, p, r, zero.mean, phi_f_structure, control))
  } else {
    return(.unpen_em(data, p, r, zero.mean, phi_f_structure, em_control))
  }
}




#' @keywords internal
.unpen_em <- function(data, p, r, zero.mean, phi_f_structure, em_control) {
  Y <- as.matrix(data)
  n_var <- ncol(Y); T.fin <- nrow(Y)
  
  # structure → Phi.f array (N x (N*p+1) x r)
  phi_arr <- .normalize_phi_structure(phi_f_structure, n_var, r)
  Phi.f.array <- make.Phi.f(structure = phi_arr, lag.order = p)
  
  # cfg for opti.fct
  cfg <- list(
    dim.VAR = n_var,
    lag.order = p,
    number.factors = r,
    zero.mean = isTRUE(zero.mean)
  )
  
  # --- initialization (same as ML) ---
  logit <- function(x) log(x/(1 - x))
  A0 <- 0.10; B0 <- 0.80
  phi_r0 <- if (r > 0) rep(0.95, r) else numeric(0L)
  Omega0 <- diag(0.1, n_var)
  
  # OLS VAR for Phi_c
  Ydep <- Y[(p + 1):T.fin, , drop = FALSE]
  Xlag <- do.call(cbind, lapply(1:p, function(L) Y[(p + 1 - L):(T.fin - L), , drop = FALSE]))
  X <- cbind(1, Xlag)
  Bhat <- vapply(seq_len(n_var), function(j) lm.fit(X, Ydep[, j])$coefficients, numeric(ncol(X)))
  Phi_c0 <- t(Bhat); if (isTRUE(zero.mean)) Phi_c0[, 1] <- 0
  
  # number of free Phi.f params (mask)
  Phi.f.as.matrix <- matrix(Phi.f.array, nrow = n_var)
  Phi.f.mask  <- (Phi.f.as.matrix != 0)
  n_phi_free  <- sum(Phi.f.mask)
  
  # packed params
  params <- c(
    if (length(phi_r0)) c(logit(A0), logit(B0), logit(phi_r0)) else c(logit(A0), logit(B0)),
    rep(0.01, n_phi_free),
    vech(Omega0),
    as.vector(Phi_c0)
  )
  
  # bounds only for (A,B,phi_r)
  lower_bounds <- rep(-Inf, length(params)); upper_bounds <- rep( Inf, length(params))
  n_head <- 2 + length(phi_r0)
  if (n_head > 0) { lower_bounds[1:n_head] <- -10; upper_bounds[1:n_head] <- 10 }
  
  tol <- em_control$tol %||% 1e-3
  max_iter <- em_control$max_iter %||% 200
  do_trace <- isTRUE(em_control$trace)
  
  obj_prev <- Inf
  iter <- 1L
  repeat {
    # --- E-step (smoother expectations kept inside opti.fct state) ---
    invisible(opti.fct(
      par_free   = params,
      par_fixed  = NaN,
      VAR.data   = Y,
      Phi.f.array = Phi.f.array,
      cfg = cfg,
      Smooth = TRUE,
      purpose = "eval"
    ))
    
    # --- M-step (maximize expected complete-data loglik) ---
    opt <- optim(
      par = params,
      fn  = opti.fct,
      par_fixed  = NaN,
      VAR.data   = Y,
      Phi.f.array = Phi.f.array,
      cfg = cfg,
      Smooth = FALSE,
      purpose = "optim",
      method = "L-BFGS-B",
      lower = lower_bounds,
      upper = upper_bounds,
      control = list(maxit = 500, trace = 0)
    )
    params <- opt$par
    obj <- opt$value
    
    rel <- 2 * abs(obj - obj_prev) / (abs(obj) + abs(obj_prev))
    if (do_trace) message(sprintf("[EM %03d] obj=%.6f  rel=%.3e", iter, obj, rel))
    if (!is.finite(rel) || rel < tol || iter >= max_iter) break
    obj_prev <- obj; iter <- iter + 1L
  }
  
  # final evaluation
  ev <- opti.fct(
    par_free   = params,
    par_fixed  = NaN,
    VAR.data   = Y,
    Phi.f.array = Phi.f.array,
    cfg = cfg,
    Smooth = FALSE,
    purpose = "eval"
  )
  
  # ---- unpack params to estimates (mirror slicing in opti.fct) ----
  inv_logit <- stats::plogis
  idx <- 1L
  A  <- inv_logit(params[idx]); idx <- idx + 1L
  B  <- inv_logit(params[idx]); idx <- idx + 1L
  phi_r <- if (r > 0) { x <- inv_logit(params[idx:(idx + r - 1L)]); idx <- idx + r; x } else numeric(0L)
  
  # ensure mask & count for unpacking Phi.f
  Phi.f.as.matrix <- matrix(Phi.f.array, nrow = n_var)
  Phi.f.mask  <- (Phi.f.as.matrix != 0)
  n_phi_free  <- sum(Phi.f.mask)
  
  Phi.f.est <- matrix(0, nrow = nrow(Phi.f.as.matrix), ncol = ncol(Phi.f.as.matrix))
  if (n_phi_free > 0) {
    Phi.f.est[Phi.f.mask] <- params[idx:(idx + n_phi_free - 1L)]
    idx <- idx + n_phi_free
  }
  Phi.f.est <- array(Phi.f.est, dim = dim(Phi.f.array))
  
  Nstar <- n_var * (n_var + 1L) / 2L
  L_vec <- params[idx:(idx + Nstar - 1L)]; idx <- idx + Nstar
  Lmat  <- matrix(D.matrix(n_var) %*% L_vec, ncol = n_var)
  Lmat[upper.tri(Lmat)] <- 0
  Omega <- Lmat %*% t(Lmat)
  
  need <- if (isTRUE(zero.mean)) n_var^2 * p else (n_var^2 * p + n_var)
  Phi_c_vec <- params[idx:(idx + need - 1L)]
  Phi_c <- if (isTRUE(zero.mean)) cbind(0, matrix(Phi_c_vec, nrow = n_var)) else matrix(Phi_c_vec, nrow = n_var)
  
  # ICs
  K <- length(c(A, B, phi_r)) + n_phi_free + Nstar + length(Phi_c_vec)
  n_eff <- (T.fin - p) * n_var
  AIC  <-  2 * K + 2 * ev$average.L * (T.fin - p)
  AICc <-  AIC + (2*K^2 + 2*K) / (n_eff - K - 1)
  BIC  <-  log(n_eff) * K + 2 * ev$average.L * (T.fin - p)
  
  out <- list(
    estimate = list(A = A, B = B, phi_r = phi_r,
                    Phi_f = Phi.f.est, Phi_c = Phi_c, Omega = Omega),
    lik = list(avg = ev$average.L, sum = ev$full.L),
    ic  = list(AIC = AIC, AICc = AICc, BIC = BIC),
    eval = ev,
    optim = list(par = params, value = obj, convergence = NA,
                 counts = NA, message = sprintf("EM finished in %d iters", iter)),
    meta = list(N = n_var, p = p, r = r, zero.mean = isTRUE(zero.mean),
                nobs = T.fin, method = "EM")
  )
  class(out) <- "tvvar_fit"
  out
}




#' Unpenalized estimation (ML) for TV-VAR — minimal wrapper
#' @param data numeric matrix (T x N), already transformed (demeaned if zero.mean=TRUE).
#' @param p integer, VAR lag order.
#' @param r integer, number of factors.
#' @param zero.mean logical; TRUE fixes intercept to 0.
#' @param phi_f_structure 3D array [N,N,r], or list of r (N x N) matrices, or block matrix N x (N*r).
#' @param control list for optim (e.g., list(maxit=2000, trace=0)).
#' @return list with estimates, likelihood, ICs, eval, optim info.
#' @export
#' @keywords internal
.unpen_ml <- function(data, p, r, zero.mean, phi_f_structure, control) {
  Y <- as.matrix(data)
  n_var <- ncol(Y); T.fin <- nrow(Y)
  
  # structure → Phi.f array (N x (N*p+1) x r)
  Phi.f.str <- array(1, dim = c(n_var,n_var,r))
  Phi.f.array <- make.Phi.f(structure = Phi.f.str, lag.order = p)
  
  # cfg for opti.fct
  cfg <- list(
    dim.VAR = n_var,
    lag.order = p,
    number.factors = r,
    zero.mean = isTRUE(zero.mean)
  )
  
  # --- initialization (same as EM) ---
  logit <- function(x) log(x/(1 - x))
  A0 <- 0.10; B0 <- 0.80
  phi_r0 <- if (r > 0) rep(0.95, r) else numeric(0L)
  Omega0 <- matrix(0.2, n_var, n_var)   # fill with 0.2 off-diagonal
  diag(Omega0) <- 0.3   
  
  ## --- Estimate Phi^c as coefficient matrix in static VAR (exact Gorgi-style) ---
  Y_named <- as.data.frame(Y)
  colnames(Y_named) <- paste0("y", seq_len(ncol(Y_named)))
  
  B.t <- matrix(0, nrow = n_var, ncol = (n_var * p + 1))
  for (j in 1:n_var) {
    B.t[j, ] <- coefficients(VAR(Y_named, p = p))[[j]][, 1]
  }
  
  A.t <- B.t[, -(n_var * p + 1), drop = FALSE]   # drop intercept column
  
  ## case of nonzero intercept
  if (isFALSE(zero.mean)) {
    A.t <- cbind(B.t[, (n_var * p + 1), drop = FALSE], A.t)
  }
  
  Phi_c0 <- A.t
  
  # number of free Phi.f params (mask)
  Phi.f.as.matrix <- matrix(Phi.f.array, nrow = n_var)
  Phi.f.mask  <- (Phi.f.as.matrix != 0)
  n_phi_free  <- sum(Phi.f.mask)
  
  # packed params
  par_init <- c(
    if (length(phi_r0)) c(logit(A0), logit(B0), logit(phi_r0)) else c(logit(A0), logit(B0)),
    rep(0.01, n_phi_free),
    vech(Omega0),
    as.vector(Phi_c0)
  )
  
  # ML optimize
  opt <- optim(
    par = par_init,
    fn  = opti.fct,
    par_fixed = NaN,
    VAR.data  = Y,
    Phi.f.array = Phi.f.array,
    cfg = cfg,
    Smooth = FALSE,
    purpose = "optim",
    method = "BFGS",
    control = control
  )
  
  # evaluates at ML optimum
  ev <- opti.fct(
    par_free    = opt$par,
    par_fixed   = NaN,
    VAR.data    = Y,
    Phi.f.array = Phi.f.array,
    cfg         = cfg,
    Smooth      = TRUE,
    purpose     = "eval"
  )
  
  p_out <- opt$par
  T.fin <- nrow(Y)  
  my.hess      <- optimHess(p_out, hessian.fct.untr,
                            VAR.data = Y, Phi.f.array = Phi.f.array,
                            zero.mean = isTRUE(zero.mean), cfg = cfg)
  
  cov.matrix   <- solve(my.hess * T.fin)

  
  
  
  # ---- unpack params to estimates (mirror slicing in opti.fct) ----
  inv_logit <- stats::plogis
  idx <- 1L
  A  <- inv_logit(opt$par[idx]); idx <- idx + 1L
  B  <- inv_logit(opt$par[idx]); idx <- idx + 1L
  phi_r <- if (r > 0) inv_logit(opt$par[idx:(idx + r - 1L)]) else numeric(0L)
  idx <- idx + length(phi_r)
  
  # ensure mask & count for unpacking Phi.f
  Phi.f.as.matrix <- matrix(Phi.f.array, nrow = n_var)
  Phi.f.mask  <- (Phi.f.as.matrix != 0)
  n_phi_free  <- sum(Phi.f.mask)
  
  Phi.f.est <- matrix(0, nrow = nrow(Phi.f.as.matrix), ncol = ncol(Phi.f.as.matrix))
  if (n_phi_free > 0) {
    Phi.f.est[Phi.f.mask] <- opt$par[idx:(idx + n_phi_free - 1L)]
    idx <- idx + n_phi_free
  }
  Phi.f.est <- array(Phi.f.est, dim = dim(Phi.f.array))
  
  Nstar <- n_var * (n_var + 1L) / 2L
  L_vec <- opt$par[idx:(idx + Nstar - 1L)]; idx <- idx + Nstar
  Lmat  <- matrix(D.matrix(n_var) %*% L_vec, ncol = n_var)
  Lmat[upper.tri(Lmat)] <- 0
  Omega <- Lmat %*% t(Lmat)
  
  need <- if (isTRUE(zero.mean)) n_var^2 * p else (n_var^2 * p + n_var)
  Phi_c_vec <- opt$par[idx:(idx + need - 1L)]
  Phi_c <- if (isTRUE(zero.mean)) cbind(0, matrix(Phi_c_vec, nrow = n_var)) else matrix(Phi_c_vec, nrow = n_var)
  
  # ICs
  K <- length(c(A, B, phi_r)) + n_phi_free + Nstar + length(Phi_c_vec)
  n_eff <- (T.fin - p) * n_var
  AIC  <-  2 * K + 2 * ev$average.L * (T.fin - p)
  AICc <-  AIC + (2*K^2 + 2*K) / (n_eff - K - 1)
  BIC  <-  log(n_eff) * K + 2 * ev$average.L * (T.fin - p)
  
  out <- list(
    estimate = list(A = A, B = B, phi_r = phi_r,
                    Phi_f = Phi.f.est, Phi_c = Phi_c, Omega = Omega),
    lik = list(avg = ev$average.L, sum = ev$full.L),
    ic  = list(AIC = AIC, AICc = AICc, BIC = BIC),
    optim = list(par = opt$par, value = opt$value, convergence = opt$convergence,
                 counts = opt$counts, message = opt$message),
    vcov = cov.matrix,
    theta = p_out,
    meta = list(N = n_var, p = p, r = r, zero.mean = isTRUE(zero.mean),
                nobs = T.fin, method = "ML", Phi.f.array = Phi.f.array)
  )
  class(out) <- "tvvar_fit"
  out
}

# internal: accept array/list/block matrix; enforce r if provided
#' @keywords internal
#' @noRd
.normalize_phi_structure <- function(phi, n_var, r = NULL) {
  if (is.array(phi) && length(dim(phi)) == 3L) {
    stopifnot(dim(phi)[1L] == n_var, dim(phi)[2L] == n_var)
    if (!is.null(r)) stopifnot(dim(phi)[3L] == r)
    return(phi)
  }
  if (is.list(phi)) {
    r0 <- length(phi); if (!is.null(r)) stopifnot(r0 == r) else r <- r0
    stopifnot(all(vapply(phi, function(m) all(dim(m) == c(n_var, n_var)), logical(1))))
    arr <- array(0, c(n_var, n_var, r))
    for (j in seq_len(r)) arr[, , j] <- phi[[j]]
    return(arr)
  }
  if (is.matrix(phi)) {
    stopifnot(nrow(phi) == n_var, ncol(phi) %% n_var == 0L)
    r0 <- ncol(phi) %/% n_var; if (!is.null(r)) stopifnot(r0 == r) else r <- r0
    arr <- array(0, c(n_var, n_var, r))
    for (j in seq_len(r)) {
      cols <- ((j - 1L) * n_var + 1L):(j * n_var)
      arr[, , j] <- phi[, cols, drop = FALSE]
    }
    return(arr)
  }
  stop("`phi_f_structure` must be a 3D array [N,N,r], a list of N x N matrices, or an N x (N*r) block matrix.")
}






#' Core optimizer for TV-VAR likelihood / evaluation
#'
#' Runs the Kalman filter (and optional smoother) and returns the
#' negative average log-likelihood (for optimization) or detailed
#' state outputs (for evaluation).
#'
#' @details Parameters are unconstrained and mapped internally
#'   (e.g., via `plogis`). Nonzero entries of `Phi.f` are filled
#'   using `Phi.f.mask`.
#'
#' @param params Numeric vector of unconstrained parameters (see Details).
#' @param VAR.data Matrix (T x N).
#' @param Phi.f.mask Matrix or array used to indicate nonzero entries of Phi.f
#'   (same shape as Phi.f when reshaped); values != 0 are filled from `params`.
#' @param zero.mean Logical; if TRUE, intercept is fixed at 0.
#' @param Smooth Logical; if TRUE, run Kalman smoother.
#' @param purpose "optim" (return negative avg loglik) or "eval" (return list).
#' @param r Integer, number of factors.
#' @param p Integer, VAR lag order.
#' @param n_var Integer, number of VAR variables (ncol(VAR.data)).
#'
#' @return If purpose="optim": scalar negative average loglik.
#'   If purpose="eval": list with likelihoods, states, and components.
#' @export
#' @importFrom stats plogis

#' @keywords internal
#' @noRd
opti.fct <- function(par_free, par_fixed, VAR.data, Phi.f.array, cfg, Smooth, purpose){
  
  # Unpack
  zero.mean = cfg$zero.mean
  dim.VAR = cfg$dim.VAR
  lag.order = cfg$lag.order
  number.factors = cfg$number.factors
  
  
  # Combine [vec] par free and [list] par fixed into [vec] params
  if (is.numeric(par_free)) {
    params <- par_free
  } else {
    params <- rebuild_params(par_free, par_fixed)
  }  
  # transform back to (0,1) space
  A <- exp(params[1])/(1+exp(params[1]))
  B <- exp(params[2])/(1+exp(params[2]))
  
  
  if(number.factors>0){
    
    # [FIX] zero mean option
    # Unpacks (conditional) parameter 'pphi'
    
    if(number.factors==1){
      pphi <- matrix(tanh(params[3:(3+number.factors-1)]), nrow=1, ncol=1)} else{
        pphi <- diag(tanh(params[3:(3+number.factors-1)]))
      }
    
    
    # sigma_eta is I - phi * phi. 
    Q <- diag(number.factors)-pphi %*% t(pphi)
    
    # a1 for the factors, is 0 as stationary
    starta <- numeric(number.factors)
    
    # P1 is unconditional variance of the variance, which is Identity due to orthogonal variance of eta.
    startP <- diag(number.factors)
    
  } else if (number.factors==0){
    pphi <- matrix(0, nrow=1, ncol=1)
    Q <- diag(1)-pphi %*% t(pphi)
    starta <- numeric(1)
    startP <- diag(1)
  }
  
  Phi.f.as.matrix <- matrix(Phi.f.array, nrow=dim.VAR)
  count.params <- 2+number.factors+length(which(Phi.f.as.matrix != 0))
  Phi.f <- matrix(0, nrow(Phi.f.as.matrix), ncol(Phi.f.as.matrix))
  Phi.f[which(Phi.f.as.matrix != 0)]=params[(2+number.factors+1):count.params]

  #length(which(Phi.f != 0))
  #count.params <- 3+length(as.vector(Phi.f))
  
  ######## initializations for state, state variance and covariance matrix
  
  # Obtain the base covariance matrix omega
  Nstar <- dim.VAR*(dim.VAR+1)/2
  m1 <- matrix(D.matrix(dim.VAR) %*% params[(count.params+1):(count.params+Nstar)], ncol=dim.VAR)
  m1[upper.tri(m1)]=0
  momega <- m1%*%t(m1)
  
  # Init the covariance matrix H_t (using unconditional variance)
  firstH <- momega/(1-A-B)
  
  # Flag for updating covariance (1, 2 both using the BEKK model for H_t). Can change this!
  covi=1
  
  count.params = count.params+Nstar
  
  
  # [FIX] zero.mean option
  # Unpacks (conditional) parameter 'Phi.c'
  if(zero.mean==TRUE){
    Phi.c <- cbind(rep(0,dim.VAR),matrix(params[(count.params+1):(count.params+dim.VAR^2*lag.order)], nrow=dim.VAR))
  } else{
    #Phi.c <- cbind(rep(0,dim.VAR),matrix(params[(count.params+1):(count.params+dim.VAR^2*lag.order)], nrow=dim.VAR))
    Phi.c <- matrix(params[(count.params+1):(count.params+dim.VAR^2*lag.order+dim.VAR)], nrow=dim.VAR)
  }
  
  
  # Obtain y_tilde, Z from the VAR data by building state space form.
  if(number.factors==0){
    data.ss <- build.statespace.form(VAR.data=VAR.data, lag.order=lag.order, number.factors=1, Phi.f = Phi.f, Phi.c=Phi.c)} else{
      data.ss <- build.statespace.form(VAR.data=VAR.data, lag.order=lag.order, number.factors=number.factors, Phi.f = Phi.f, Phi.c=Phi.c)}
  
  
  ytilde <- data.ss$ytilde
  Z <- data.ss$Z
  T.here <- length(ytilde[1,])-1
  
  
  ### call C++ loop function for the KF
  ### Output of the Kalman filter gives (sumlik, avlik, pstate, pstatevariance, fstate, fstatevariance, F, L, vmat, aH, T)
  
  output <- my_loop_main(ytilde, Z, startP, covi, firstH, momega, pphi, A, B, Q)
  
  last.obs <- length(output$pstate[1,])
  
  if(number.factors==0){} else{
    ### Kalman smoother 
    predicted.a <- matrix(output$pstate[,-last.obs, drop=FALSE],ncol=T.here)
    predicted.P <- output$pstatevariance[,,-last.obs, drop=FALSE]
    filtered.a <- output$fstate[,-last.obs, drop=FALSE]
    filtered.P <- output$fstatevariance[,,-last.obs, drop=FALSE]
    
    L.array <- output$L[,,-last.obs, drop=FALSE]
    F.array <- output$`F`[,,-last.obs, drop=FALSE]
    v.mat <- output$vmat[,-last.obs,drop=FALSE]
    pred.err.var <- output$`F`[,,-last.obs,drop=FALSE]
    
    ### more things to collect
    r.mat <- matrix(0, nrow=number.factors,ncol=T.here)
    N.array <- array(0, dim=c(number.factors,number.factors,T.here))
    alphahat.mat <- matrix(0, nrow=nrow(predicted.a), ncol=ncol(predicted.a))
    V.array <- array(0, dim=c(number.factors,number.factors,T.here))
    smooth.error <- matrix(0, nrow=nrow(ytilde), ncol=T.here)
    
    # Added init of autocovariance matrix (Γ = similar(P, length(a) - 1))
    gamma.array <- array(0, dim=c(number.factors,number.factors,T.here-1))
  }
  
  if(Smooth == FALSE){} else{
    ytilde.short = ytilde[,-last.obs, drop=FALSE]
    
    if(number.factors>0){
      
      for(t in T.here:1){
        r.t <- r.mat[,t]
        N.t <- N.array[,,t]
        
        F.t.inv <- solve(F.array[,,t])
        L.t <- L.array[,,t]
        v.t <- v.mat[,t]
        a.t <- predicted.a[,t]
        P.t <- predicted.P[,,t]
        Z.t = Z[,,t, drop=FALSE]
        Z.t <- matrix(Z.t, nrow = dim(Z)[1], ncol = dim(Z)[2])  # force 2D
        
        r.tminus1 <- t(Z.t) %*% F.t.inv %*% v.t + t(L.t) %*% r.t
        N.tminus1 <- t(Z.t) %*% F.t.inv %*% Z.t + t(L.t) %*% N.t %*% L.t
        
        alphahat.mat[,t] <- a.t + P.t %*% r.tminus1
        V.array[,,t] <- P.t - P.t %*% N.tminus1 %*% P.t
        
        # Extended KS recursions 
        if (t > 1){
          gamma.array[,, t-1] <- diag(number.factors) - P.t %*% N.tminus1}
        if (t < T.here) {
          J_t <- L.t %*% P.t
          gamma.array[,, t] <- gamma.array[,, t] %*% J_t
        }
        
        
        r.mat[,(t-1)] = r.tminus1
        N.array[,,(t-1)] = N.tminus1
        
        smooth.error[,t] <- ytilde.short[,t]-Z.t %*% alphahat.mat[,t]
        
      }
    } else{}
  }
  
  # Add penalisation
  #
  #penalty_phi_c <- lambda/length(Phi.c) * mean(abs(Phi.c))
  #penalty_phi_f <- lambda/length(Phi.f) * mean(abs(Phi.f))
  
  if(purpose=="optim"){
    L = -output$avlik
  } else if(purpose=="eval"){
    L = -output$avlik
    
    
    if(number.factors==0){
      result <- new.env()
      result$average.L <- -output$avlik
      result$full.L <- -output$sumlik
      result$Phi.c <- Phi.c
      result$array.filtered.H <- output$aH
      result$momega =vech(momega)
      result$Nstar =Nstar
      as.list(result)} else{
        
        
        result <- new.env()
        result$average.L <- -output$avlik
        result$full.L <- -output$sumlik
        result$filtered.state <- filtered.a
        result$filtered.state.variance <- filtered.P
        result$predicted.state <- predicted.a
        result$predicted.state.variance <- predicted.P
        result$Phi.c <- Phi.c
        result$array.filtered.H <- output$aH
        result$momega =vech(momega)
        result$Nstar =Nstar
        result$Z=Z
        result$v.mat= v.mat
        result$pred.err.var=pred.err.var
        result$Y.minus1 <- data.ss$Y.minus1
        
        if(Smooth == TRUE){
          result$smoothed.state <- alphahat.mat
          result$smoothed.state.variance <- V.array
          result$smooth.error=smooth.error
          result$ytilde <- ytilde
          result$smoothed.state.autocov <- gamma.array
          
          # Initialize V* matrices
          V_0_minus1_star <- array(0, dim=c(number.factors,number.factors))
          V_0_star <- array(0, dim=c(number.factors,number.factors))
          V_minus1_star <- array(0, dim=c(number.factors,number.factors))
          
          
          # Fill V* matrices
          for (t in 2:T.here){
            V_0_minus1_star <- V_0_minus1_star + gamma.array[,,t-1] + result$smoothed.state[,t] %*% t(result$smoothed.state[,t-1])
            V_0_star <- V_0_star + V.array[,,t] + result$smoothed.state[,t] %*% t(result$smoothed.state[,t])
            V_minus1_star <- V_minus1_star + V.array[,,t-1] + result$smoothed.state[,t-1] %*% t(result$smoothed.state[,t-1])
          }
          
          # Store V* matrices
          result$V_0_minus1_star <- V_0_minus1_star
          result$V_0_star <- V_0_star
          result$V_minus1_star <- V_minus1_star
          
          # pphi  
          pphi_full <- V_0_minus1_star %*% solve(V_minus1_star)
          pphi_full <- pmin(pmax(pphi_full, -0.999), 0.999)
          
          
          pphi_trans <- atanh(pphi_full)
          if (any(is.nan(pphi_trans))) {
            print("pphi_full:")
            print(pphi_full)
            print("Smoothed state:")
            print(result$smoothed.state)
            print("V_0_minus1_star:")
            print(V_0_minus1_star)
            print("V_minus1_star:")
            print(V_minus1_star)
            stop("NaN detected in atanh(pphi_full)")
          }
          
          if (nrow(pphi_full) == 1) {
            result$pphi_ols <- (pphi_full[1, 1])
          } else {
            result$pphi_ols <- (diag(pphi_full))
          }
          
          # Phi.f gradient
          grad <- matrix(0, nrow=dim.VAR, ncol=dim.VAR)
          for (t in 1:T.here){
            Z.t = Z[,,t, drop=FALSE]
            Z.t <- matrix(Z.t, nrow = dim(Z)[1], ncol = dim(Z)[2]) 
            
            grad <- grad + (smooth.error[,t] %*% t(smooth.error[,t]) + Z.t %*% V.array[,,t] %*% t(Z.t))
          }
          grad.Phi.f <- -(0.5/T.here)*grad
          result$grad.Phi.f <- grad.Phi.f
          
        }
        
        as.list(result)
      }
  }}