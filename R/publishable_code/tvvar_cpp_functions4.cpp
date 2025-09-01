#include <RcppArmadillo.h>
#include <math.h>



using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List my_loop(arma::mat ytilde, arma::cube Z, arma::mat startP, int cov,  
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







