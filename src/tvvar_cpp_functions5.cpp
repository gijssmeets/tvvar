#include <RcppArmadillo.h>
#include <math.h>



using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List my_loop_main(arma::mat ytilde, arma::cube Z, arma::mat startP, int cov,  
arma::mat startH, arma::mat momega, arma::mat pphi, double A, double B, arma::mat Q){
	
	// dimensions
	double N = ytilde.n_rows ;
	double T = ytilde.n_cols ;
	double r = startP.n_cols ;
	 
	// declare arrays to collect things
	arma::vec vlik=zeros(T-1), starta=zeros(r);
	arma::mat pstate=zeros(r,T), fstate=zeros(r,T), vmat=zeros(N,T);
	arma::cube aH=zeros(N,N,T), pstatevariance=zeros(r,r,T), fstatevariance=zeros(r,r,T), L=zeros(r,r,T), F=zeros(N,N,T);
	
	// declare the rest
	double sumlik, avlik; 
	arma::vec yt, vt, at;
	arma::mat Ht,Ft,Pt, Kt, Lt,Zt;
	 
	// initialization
	pstate.col(0) = starta;
	pstatevariance.slice(0) = startP;
	aH.slice(0) = momega + B*startH;
	
	for(double t = 0; t < T-1; ++t){
		yt = ytilde.col(t);
		Zt = Z.slice(t);

		at = pstate.col(t);
		Pt = pstatevariance.slice(t);

		Ht = aH.slice(t);

		vt = yt - Zt*at;
		Ft = Zt*Pt*trans(Zt) + Ht;
		Kt = pphi*Pt*trans(Zt) * inv(Ft);
		Lt = pphi - Kt*Zt;

		if(cov == 1){
			aH.slice(t+1) = momega +B*Ht + A* (vt*trans(vt));
		} else if(cov==2){
			aH.slice(t+1) = momega +B*Ht + A* (vt*trans(vt));
		}
		
		F.slice(t) = Ft;
		L.slice(t) = Lt;
		vmat.col(t) =vt;
		
		fstate.col(t) = at+Pt*trans(Zt)*inv(Ft)*vt;
		fstatevariance.slice(t) = Pt-Pt*trans(Zt)*inv(Ft)*Zt*Pt;

		pstate.col(t+1) = pphi*at + Kt * vt;
		pstatevariance.slice(t+1) = pphi*Pt*trans(Lt) + Q;

		vlik(t) = -0.5* N*log(2*M_PI) - 0.5*log(det(Ft)) - 0.5*dot(trans(vt)*inv(Ft),vt);
	}
	
	sumlik =sum(vlik);
	avlik = sum(vlik)/vlik.n_elem;
	
	if(avlik != avlik){avlik=-100;}
	if(avlik < -100){avlik= -100;}
	
	return List::create(	
	_["sumlik"] = sumlik,
	_["avlik"] = avlik,
	_["pstate"] = pstate,
	_["pstatevariance"] = pstatevariance,	
	_["fstate"] = fstate,
	_["fstatevariance"] = fstatevariance,	
	_["F"] =F,
	_["L"] =L,
	_["vmat"] = vmat,
	_["aH"] = aH,	
	_["T"] = T );
}





  
// [[Rcpp::export]]  
arma::mat create_Y_minus1(const arma::mat& Y, int lag_order, int T_fin, int dim_VAR) {

  // Preallocate Y_minus1 with appropriate size
  arma::mat Y_minus1(dim_VAR * lag_order + 1, T_fin - lag_order, arma::fill::zeros);
  
  // Populate Y_minus1
  for (int t = lag_order; t < T_fin; ++t) {
    int col_index = t - lag_order; // Adjust for 0-based index
      
    // Set the first element of the column to 1.0
    Y_minus1(0, col_index) = 1.0;
      
    // Fill the rest with the lagged values from Y
    int pos = 1; // Start position in Y_minus1
    for (int lag = 0; lag < lag_order; ++lag) {
      // Copy the values from Y(:, t-lag)
      Y_minus1.submat(pos, col_index, pos + dim_VAR - 1, col_index) = Y.col(t - lag - 1);
      pos += dim_VAR;
    }
  }
    
  return Y_minus1;
}  
  
  
  

  
// Function to create and return the Z cube
// [[Rcpp::export]] 
arma::cube create_Z(const arma::mat& Y_minus1, const arma::mat& Phi_f, const arma::mat& Q_here, int T_fin, int dim_VAR, int number_factors, int lag_order) {
  // Preallocate the Z cube with the correct dimensions
  arma::cube Z(dim_VAR, number_factors, T_fin-lag_order);
  
  // Populate the Z cube
  for (int t = 0; t < T_fin-lag_order; ++t) {
    // Extract the column from Y_minus1 and reshape it into a row vector
    arma::rowvec Y_minus1_t = Y_minus1.col(t).t();
      
    // Compute the Kronecker product of Y_minus1_t and Phi_f
    arma::mat kron_prod = arma::kron(Y_minus1_t, Phi_f);
    
    // Multiply by Q_here to get the matrix for Z(:,:,t)
    Z.slice(t) = kron_prod * Q_here;
  }
    
  return Z;
}  
