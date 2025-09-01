#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;



// [[Rcpp::export]]
Rcpp::List irf_generator_loop(arma::cube& PHI, arma::mat& irf_out, arma::mat& PI_out, 
                              const arma::cube& Phi, const arma::mat& psi, const arma::vec& aT, 
                              const arma::mat& FF1, const arma::cube& FF, const arma::mat& CT, 
                              const arma::vec& ek, int lags, int B, int N, int p, int r) {
  
  // IRF matrix and intermediate variables
  arma::mat irf(N, lags, arma::fill::zeros);
  arma::mat PI(N * p, N * p, arma::fill::zeros);
  
  // Loop over B replications
  for (int j = 0; j < B; ++j) {
    // Initialize f matrix
    arma::mat f(r, lags, arma::fill::zeros);
    f.col(0) = aT + FF1.row(j).t();  // Ensuring FF1 row is transposed to match dimensions
    
    // Initialize Phit (ensure it's a matrix with correct dimensions)
    arma::mat Phit = Phi.slice(0);  // Phi.slice(0) should be a matrix of size (N, N * p)
    
    for (int i = 0; i < r; ++i) {
      Phit += Phi.slice(i + 1) * f(i, 0);  // Ensure Phi.slice(i+1) is of compatible dimensions
    }
    
    // Store PHI for the first time step
    // Ensure submatrices and subcubes match dimensions exactly
    PHI.slice(0).submat(0, 0, N - 1, N * p - 1) = Phit;
    
    // Handle the block structure in PHI (if p > 1)
    if (p > 1) {
      PHI.slice(0).submat(N, 0, N * p - 1, N * (p - 1) - 1) = arma::eye(N * (p - 1), N * (p - 1));
    }
    
    // Iterate through lags
    for (int t = 1; t < lags; ++t) {
      // Update factor process
      f.col(t) = psi * f.col(t - 1) + FF.slice(j).col(t - 1) * (1 - std::pow(psi, 2));
      
      // Update Phit for each lag
      Phit = Phi.slice(0);
      for (int i = 0; i < r; ++i) {
        Phit += Phi.slice(i + 1) * f(i, t);
      }
      
      // Store PHI for the current time step
      PHI.slice(t).submat(0, 0, N - 1, N * p - 1) = Phit;
      
      if (p > 1) {
        PHI.slice(t).submat(N, 0, N * p - 1, N * (p - 1) - 1) = arma::eye(N * (p - 1), N * (p - 1));
      }
    }
    
    // Initialize IRF calculation
    irf.col(0) = CT * ek;
    PI = PHI.slice(0);  // Make sure PI gets a matrix from PHI.slice
    
    // Compute IRF for each lag
    for (int l = 1; l < lags; ++l) {
      irf.col(l) = PI.submat(0, 0, N - 1, N - 1) * CT * ek;
      PI = PHI.slice(l) * PI;  // Ensure the matrix multiplication is valid
    }
    
    // Check dimensions before adding to avoid mismatches
    if (irf_out.n_rows == irf.n_rows && irf_out.n_cols == irf.n_cols) {
      irf_out += irf / B;
    } else {
      Rcpp::Rcerr << "Error: Dimension mismatch in irf_out." << std::endl;
      return Rcpp::List::create(Rcpp::Named("error") = "Dimension mismatch in irf_out");
    }
    
    if (PI_out.n_rows == PI.n_rows && PI_out.n_cols == PI.n_cols) {
      PI_out += PI / B;
    } else {
      Rcpp::Rcerr << "Error: Dimension mismatch in PI_out." << std::endl;
      return Rcpp::List::create(Rcpp::Named("error") = "Dimension mismatch in PI_out");
    }
  }
  
  // Return the result as a list
  return Rcpp::List::create(Rcpp::Named("irf") = irf_out, Rcpp::Named("PI") = PI_out);
}