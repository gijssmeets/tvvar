// irf_generator.cpp
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List irf_generator_cpp(arma::cube Phi,
                       const arma::vec& psi,
                       const arma::vec& aT,
                       const arma::mat& PT,
                       const arma::mat& HT,
                       const arma::vec& ek,
                       const int lags,
                       const bool fixed = false,
                       const int B = 500,
                       const int seed_f = 1234)
{
  // drop intercept column
  if (Phi.n_cols < 2) stop("Phi must have intercept + coefficients.");
  Phi.shed_col(0);  // N x (N*p) x (r+1)
  
  const uword N = Phi.n_rows;
  if (N == 0) stop("N == 0.");
  if (lags <= 0) stop("`lags` must be positive.");
  if (ek.n_elem != N) stop("Length of ek must be N.");
  if (HT.n_rows != N || HT.n_cols != N) stop("HT must be N x N.");
  if (Phi.n_cols % N != 0) stop("After dropping intercept, Phi.n_cols must be multiple of N.");
  const uword p = Phi.n_cols / N;
  
  const uword r = aT.n_elem;
  if (Phi.n_slices != r + 1) stop("Phi.n_slices must equal r + 1.");
  if (PT.n_rows != r || PT.n_cols != r) stop("PT must be r x r.");
  
  // psi (allow scalar or length r)
  vec psi_vec(r, fill::ones);
  if (psi.n_elem == 1) psi_vec.fill(psi[0]);
  else if (psi.n_elem == r) psi_vec = psi;
  else stop("psi must be length 1 or r.");
  
  cube PHI(N * p, N * p, lags, fill::zeros);
  mat  irf(N, lags, fill::zeros);
  mat  irf_out(N, lags, fill::zeros);
  mat  PI_out(N * p, N * p, fill::zeros);
  mat  CT = chol(HT, "lower");  // t(chol(HT)) in R
  
  auto build_Phit = [&](const vec& f_col)->mat {
    mat Phit = Phi.slice(0);
    for (uword i = 0; i < r; ++i) Phit += Phi.slice(i + 1) * f_col(i);
    return Phit;
  };
  
  if (fixed) {
    // --- fixed = TRUE (unchanged) ---
    mat fmat(r, lags, fill::none); fmat.fill(datum::nan);
    fmat.col(0) = aT;
    
    PHI.slice(0).submat(0, 0, N - 1, N * p - 1) = build_Phit(fmat.col(0));
    if (p > 1) {
      const uword npp1 = N * (p - 1);
      PHI.slice(0).submat(N, 0,      N * p - 1, npp1 - 1) = eye(npp1, npp1);
      PHI.slice(0).submat(N, npp1,   N * p - 1, N * p - 1).zeros();
    }
    for (uword t = 1; t < (uword)lags; ++t) {
      fmat.col(t) = fmat.col(t - 1);
      PHI.slice(t).submat(0, 0, N - 1, N * p - 1) = build_Phit(fmat.col(t));
      if (p > 1) {
        const uword npp1 = N * (p - 1);
        PHI.slice(t).submat(N, 0,      N * p - 1, npp1 - 1) = eye(npp1, npp1);
        PHI.slice(t).submat(N, npp1,   N * p - 1, N * p - 1).zeros();
      }
    }
    vec shock = CT * ek;
    irf.col(0) = shock;
    mat PI = PHI.slice(0);
    for (uword l = 1; l < (uword)lags; ++l) {
      irf.col(l) = PI.submat(0, 0, N - 1, N - 1) * shock;
      PI = PHI.slice(l) * PI;
    }
    irf_out = irf; PI_out = PI;
    return List::create(_["irf"] = irf_out, _["PI"] = PI_out);
  }
  
  // --- fixed = FALSE: use R's RNG + MASS::mvrnorm for FF1 (exact parity) ---
  {
    RNGScope scope;                               // use R’s RNG state
    Function set_seed("set.seed");
    set_seed(seed_f);
    
    // 1) FF1 <- MASS::mvrnorm(n=B, mu=rep(0,r), Sigma=PT)
    Environment MASS = Environment::namespace_env("MASS");
    Function mvrnorm = MASS["mvrnorm"];
    NumericVector mu(r);                          // zeros
    NumericMatrix FF1_R = mvrnorm(_["n"] = B, _["mu"] = mu, _["Sigma"] = PT, _["empirical"] = false);
    mat FF1 = as<mat>(FF1_R);                     // B x r
    
    // 2) FF <- array(rnorm(r*(lags-1)*B), c(r, lags-1, B))  (same order as R)
    cube FF(r, lags - 1, B, fill::none);
    {
      NumericVector z = Rcpp::rnorm((int)r * (int)(lags - 1) * (int)B);
      std::copy(z.begin(), z.end(), FF.memptr());
    }
    
    const double w = 1.0 / (double)B;
    vec one_minus_psi2 = 1.0 - square(psi_vec);
    
    for (uword j = 0; j < (uword)B; ++j) {
      mat fmat(r, lags, fill::none); fmat.fill(datum::nan);
      fmat.col(0) = aT + FF1.row(j).t();
      
      PHI.slice(0).zeros();
      PHI.slice(0).submat(0, 0, N - 1, N * p - 1) = build_Phit(fmat.col(0));
      if (p > 1) {
        const uword npp1 = N * (p - 1);
        PHI.slice(0).submat(N, 0,      N * p - 1, npp1 - 1) = eye(npp1, npp1);
        PHI.slice(0).submat(N, npp1,   N * p - 1, N * p - 1).zeros();
      }
      
      for (uword t = 1; t < (uword)lags; ++t) {
        fmat.col(t) = psi_vec % fmat.col(t - 1) + FF.slice(j).col(t - 1) % one_minus_psi2;
        PHI.slice(t).zeros();
        PHI.slice(t).submat(0, 0, N - 1, N * p - 1) = build_Phit(fmat.col(t));
        if (p > 1) {
          const uword npp1 = N * (p - 1);
          PHI.slice(t).submat(N, 0,      N * p - 1, npp1 - 1) = eye(npp1, npp1);
          PHI.slice(t).submat(N, npp1,   N * p - 1, N * p - 1).zeros();
        }
      }
      
      vec shock = CT * ek;
      irf.col(0) = shock;
      mat PI = PHI.slice(0);
      for (uword l = 1; l < (uword)lags; ++l) {
        irf.col(l) = PI.submat(0, 0, N - 1, N - 1) * shock;
        PI = PHI.slice(l) * PI;
      }
      
      irf_out += irf * w;
      PI_out  += PI  * w;
    }
  }
  
  return List::create(_["irf"] = irf_out, _["PI"] = PI_out);
}