#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
double check_lyapunov(const arma::mat& Phi_c, const arma::cube& Phi_f, const arma::vec& phi_r, int seed, int D = 250, int m = 100) {
  int dim_VAR = Phi_c.n_rows;
  int lag_order = Phi_c.n_cols / dim_VAR;
  
  int r = phi_r.n_elem;
  
  // Create phi.mat (diagonal matrix from phi_r)
  arma::mat phi_mat = diagmat(phi_r);
  
  // Create fac_cov (Identity matrix minus diag(phi_r^2))
  arma::mat fac_cov = eye<arma::mat>(r, r) - diagmat(phi_r % phi_r);
  
  arma::vec sim_lyap(D, fill::zeros);  // store Lyapunov exponents
  
  for (int k = 0; k < D; ++k) {
    
    // Generate errors 
    arma::mat fac_errors = randn(m, r);  
    
    // Simulated factors matrix
    arma::mat sim_factors(m, r, fill::zeros);
    
    // Simulate factors
    for (int i = 1; i < m; ++i) {
      sim_factors.row(i) = (phi_mat * sim_factors.row(i - 1).t()).t() + (sqrt(fac_cov) * fac_errors.row(i).t()).t();
    }
    
    // Initialize Phi_t_array
    arma::cube Phi_t_array(dim_VAR * lag_order, dim_VAR * lag_order, m, fill::zeros);
    
    for (int j = 0; j < m; ++j) {
      arma::vec factors_j = sim_factors.row(j).t();
      
      // First block of Phi_t_array
      arma::mat adjusted_Phi_f = Phi_c;  // Start with Phi_c (no intercept assumption removed)
      
      for (int slice_idx = 0; slice_idx < r; ++slice_idx) {
        adjusted_Phi_f += Phi_f.slice(slice_idx) * factors_j(slice_idx);  // Multiply each slice by factors_j
      }
      
      Phi_t_array.slice(j).submat(0, 0, dim_VAR - 1, adjusted_Phi_f.n_cols - 1) = adjusted_Phi_f;
      
      // Identity block in the second part
      Phi_t_array.slice(j).submat(dim_VAR, 0, dim_VAR * lag_order - 1, dim_VAR - 1) = eye<arma::mat>(dim_VAR, dim_VAR);
    }
    
    arma::mat prod_Phi_t = Phi_t_array.slice(0);
    for (int i = 1; i < m; ++i) {
      prod_Phi_t *= Phi_t_array.slice(i);
    }
    
    sim_lyap(k) = log(norm(prod_Phi_t, 2));  
  }
  
  double exp_lyap = mean(sim_lyap);
  
  return exp_lyap;
}