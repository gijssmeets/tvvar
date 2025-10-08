initialize_omega <- function(N, diag_val = 0.3, offdiag_val = 0.2) {
  M <- matrix(offdiag_val, N, N)
  diag(M) <- diag_val
  M
}



initialize_parameters <- function(VAR.data, Phi.f.array.mat.structure, cfg){
  print('[Initializing parameters..]')
  
  # Unpack config values
  dim.VAR = cfg$dim.VAR
  lag.order = cfg$lag.order
  number.factors = cfg$number.factors
  zero.mean = cfg$zero.mean
  colnames(VAR.data) <- paste0("y", 1:ncol(VAR.data))
  
  # For simulation, initialize as true values to support optimization.
  if(number.factors == 1){
    p1 = c(0.3, 0.6, 0.95)}
  else{
    p1 = c(0.3, 0.6, diag(0.95))
  }
  
  omega <- matrix(0.2, nrow = dim.VAR, ncol = dim.VAR)
  diag(omega) <- 0.3
  start.H <- omega
  print(start.H)

  
  B.t <- matrix(0, nrow=dim.VAR, ncol=(dim.VAR*lag.order+1))
  for(j in 1:dim.VAR){
    B.t[j,] = coefficients(VAR(VAR.data, p=lag.order))[[j]][,1]
  }
  A.t <- B.t[,-(dim.VAR*lag.order+1)]
  
  if(zero.mean==FALSE){
    A.t <- cbind(B.t[,(dim.VAR*lag.order+1)], A.t)
  }
  Phi.c.ini = A.t
  
  if(zero.mean==TRUE){
    Phi.f.ini = Phi.f.array.mat.structure[,-1]
  }
  else{
    Phi.f.ini = Phi.f.array.mat.structure
  }
  
  # Transform for numerical stability
  if(number.factors>0){
    params1 = 
      c(log(p1[1]/(1-p1[1])),
        log(p1[2]/(1-p1[2])),
        atanh(p1[3:(3+number.factors-1)]), 
        Phi.f.array.mat.structure[which(Phi.f.array.mat.structure !=0)]*0.01, 
        vech(start.H),
        as.vector(Phi.c.ini)
      )
  } else if(number.factors==0){
    params1 = 
      c(log(p1[1]/(1-p1[1])),
        log(p1[2]/(1-p1[2])),
        Phi.f.array.mat.structure[which(Phi.f.array.mat.structure !=0)]*0.01, 
        vech(start.H),
        as.vector(Phi.c.ini)
      )}
  
  print('[Parameters initialized]')
  return(params1)
  
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
    print('regular!')
    opts <- list(X_INIT = Phi.f.init, lambda = lambda, MAX_ITERS = 100, EPS = 1e-3)
    prox_grad_result <- apg::apg(grad_for_apg, apg::prox.l1, length(Phi.f.init), opts)
    phi_est_penalized <- prox_grad_result$x
  }
  
  # Adaptive LASSO
  if(type == 'adaptive'){
    print('adaptive!')
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
