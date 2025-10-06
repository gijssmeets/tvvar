em_algorithm <- function(params, data, Phi.f.array, cfg, conditional = FALSE, lambda = c(0.015, 0.015), weights.c = 1, weights.f = 1) {
  
  # Initialize loop variables
  iter <- 1
  obj <- -1e8
  converged <- FALSE
  T.here <- ncol(data) - 1
  
  # Initialize EM variables
  tol <- 1e-4
  max_iter <- 100
  
  # Set bounds for alpha, beta (in BEKK model) to prevent division by 0 in optimization
  lower_bounds <- rep(-Inf, length(params))
  upper_bounds <- rep( Inf, length(params))
  lower_bounds[1:2] <- -10  
  upper_bounds[1:2] <- 10    

  
  while (!converged && iter < max_iter) {
    print(paste0("Iteration: ", iter))
    
    # Format parameters to a vector
    par_free <- params_to_list(params, Phi.f.array, cfg)

    ### Expectation step
    opti.eval.em <- opti.fct(par_free = par_free, par_fixed = NaN,
                             VAR.data = data, Phi.f.array = Phi.f.array,
                             cfg = cfg, Smooth = TRUE, purpose = "eval")

    
    ### Conditional step
    if (conditional) {
      
      # Initialize placeholder variables to store the results of the C-steps
      pphi_updated <- NaN
      Phi.c_updated <- NaN
      Phi.f_updated <- NaN
      
      # [pphi] Optimize pphi
      pphi_updated <- optim_pphi_bfgs(opti.eval.em, cfg)
      
      # [pphi] Merge new pphi estimate with parameters of previous iteration
      params <- merge_conditional_params(params, Phi.f.array, cfg, pphi = pphi_updated, Phi.c = NaN)
  
      # [Phi.c] Optimize Phi.c (does not use penalization, but can be made effective in future research)
      Phi.c.total <- optim_Phi_c(opti.eval.em, data, cfg, lambda[1], weights = weights.c)
      Phi.c_updated <- Phi.c.total$Phi.c.OLS
      
      # [Phi.c] Merge new Phi.c estimate with parameters of previous iteration
      params <- merge_conditional_params(params, Phi.f.array, cfg, pphi = NaN, Phi.c = Phi.c_updated)
      
      # [Phi.f] Optimize Phi.f with penalization and adaptive weights
      Phi.f.results <- optim_Phi_f_full_likelihood(
          opti.eval.em = opti.eval.em,
          params       = params, 
          cfg          = cfg,
          data         = data,
          Phi.f.array  = Phi.f.array,
          lambda       = lambda[2],
          type         = cfg$penalty_type,
          weights      = weights.f)
      Phi.f_updated <- Phi.f.results$phi_est_penalized

      # Divide conditionally estimated parameters (fixed) and free parameters (free)
      par_div <- divide_params(
        params = params, 
        pphi = pphi_updated, 
        cfg = cfg, 
        Phi.c = Phi.c_updated, 
        Phi.f = Phi.f_updated, 
        Phi.f.array = Phi.f.array
      )
      
    } else { 
      
      # Prepare parameters if we want to estimate unconditionally (regular EM)
      par_div <- divide_params(params, pphi = NaN, cfg = cfg, Phi.c = NaN, Phi.f = NaN, Phi.f.array = Phi.f.array)
    }
      
    # Parameters that go into the M-step.
    par_free <- par_div$par_free
    
    # Parameters that are estimated in the C-step
    par_fixed <- par_div$par_fixed
    
    
    ### Maximization step
    opti <- optim(par = unlist(par_free),
                  fn = opti.fct,
                  par_fixed = par_fixed,
                  purpose = "optim",
                  VAR.data = data,
                  Phi.f.array = Phi.f.array,
                  cfg = cfg,
                  Smooth = FALSE,
                  method = "L-BFGS-B",
                  lower = lower_bounds,
                  upper = upper_bounds,
                  hessian = FALSE,
                  control = list(trace = 0, maxit = 1000))
    par_free <- opti$par
    
    # Merge (free) params and (fixed) params into one vector
    params_new <- rebuild_params(par_free, par_fixed)
    
    # Convergence check
    delta_obj <- 2 * abs(opti$value - obj) / (abs(opti$value) + abs(obj))
    converged <- delta_obj < tol
    
    # Update for next iteration
    params <- params_new
    obj <- opti$value
    iter <- iter + 1
  }
  
  # Sign identification of Phi.f.
  converged_par <- params_to_list(params, Phi.f.array, cfg)
  converged_phi_f <- converged_par$Phi.f
  
  if(converged_phi_f[1] < 0) {
    converged_par$Phi.f <- -converged_phi_f
  }
  
  
  params <- list_to_params(converged_par, cfg)
  
  
  return(list(par = params, result = opti.eval.em, obj = obj))
}  