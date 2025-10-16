rm(list=ls())
library(vars)
library(Rcpp)
library(RcppArmadillo)
library(numDeriv)
library(matrixcalc) 
library(dplyr)
library(MASS)


#setwd("C:/Users/Julia/Dropbox/tvvar/R/publishable_code")



load("/Users/gijssmeets/Documents/2. werk/3. vu/tvvar/R/publishable_code/model_output.RData")   # final model from Gorgi et al (2024)


## call c++ function script
sourceCpp('/Users/gijssmeets/Documents/2. werk/3. vu/tvvar/R/publishable_code/tvvar_cpp_functions5.cpp')
sourceCpp('/Users/gijssmeets/Documents/2. werk/3. vu/tvvar/R/publishable_code/lyapunov_cond.cpp')
sourceCpp('/Users/gijssmeets/Documents/2. werk/3. vu/tvvar/R/publishable_code/irf_generator.cpp')

source('/Users/gijssmeets/Documents/2. werk/3. vu/tvvar/R/publishable_code/tvvar_r_functions.R')


#source("irf_functions_cpp.r")  ## R code that produces the same output, just slower


########### Time-varying impulse response functions

T.max = 10
#T.max = 2


B=500 #number of Monte Carlo repetitions
lags=10 # horizon for IRFs
N=dim.VAR


IRF_lb=array(NA,c(N,lags, T.max))
IRF_ub=array(NA,c(N,lags, T.max))
IRF_med=array(NA,c(N,lags, T.max))


## loop through time points
for(t in 1:T.max){

#### settings for IRF generation
T.index=t #time point
theta=p_out #untransformed estimated coefficients
V=cov.matrix_untr # covariance matrix of untransformed estimated coeffs
s.h=1 # see index
ek= c(0,0,1) # shock position (here: third variable)


irf = array(0,c(N,lags,B))
Pi = array(0, c(N*p,N*p,B))
ev <- numeric()
lyapunov <- numeric()

## simulate from asymptotic distribution
for(i in 1:B){
  
  s.h = s.h+1
  index=s.h*t
  set.seed(index)
  par_val = mvrnorm(n = 1, mu = theta, Sigma = V)
  
  #### obtain Phi,psi,aT,PT,HT evaluated at parameter value par_val  
  opti.eval <- opti.fct(params=par_val,VAR.data = VAR.data, Phi.f.array=Phi.f.array, zero.mean=zero.mean, Smooth=FALSE, purpose="eval")
  
  phi_est <- par_val[3:(3+number.factors-1)]
  
  Phi_f_est <- Phi.f.array.mat.structure
  Phi_f_est[which(Phi.f.array.mat.structure !=0)]=par_val[(3+number.factors):(3+number.factors+number.nonzero.phis-1)]
  
  Phi_c_est <- opti.eval$Phi.c
  
  phi_est_tr <- exp(par_val[3:(3+number.factors-1)])/(exp(par_val[3:(3+number.factors-1)])+1)
  Phi_c_noint <- Phi_c_est[,-1]
  Phi_f_noint <- array(Phi_f_est, dim=c(dim.VAR, (dim.VAR*lag.order+1), r))[,-1,]
  seed_ly = 12
  
  lyapunov[i] = check_lyapunov(Phi_c_noint, Phi_f_noint, phi_est_tr, seed_ly)
  
  Phi = array(cbind(Phi_c_est, Phi_f_est), c(dim.VAR,dim.VAR*lag.order+1, (number.factors+1)))
  psi=exp(phi.est)/(1+exp(phi.est))
  #psi=phi.est
  
  aT <- opti.eval$filtered.state[,T.index]
  PT <- opti.eval$filtered.state.variance[,,T.index]
  HT <- opti.eval$array.filtered.H[,,T.index]
  
  irf.gen <- irf_generator_cpp(Phi, psi, aT, PT, HT, ek, lags, fixed = TRUE, B = B, seed_f = 1234)
  
  irf[,,i] = irf.gen$irf
  Pi[,,i] = irf.gen$PI
  ev[i] <- eigen(Pi[,,i])$values[1]
}

unstable <- which(lyapunov>0)
irf[,,unstable]=NA

#print(unstable)

for(i in 1:N){
  for(j in 1:lags){
    IRF_lb[i,j,t] = quantile(irf[i,j,], probs=0.16,na.rm=TRUE)
    IRF_ub[i,j,t] = quantile(irf[i,j,], probs=0.84,na.rm=TRUE)
    IRF_med[i,j,t] = quantile(irf[i,j,], probs=0.5,na.rm=TRUE)
  }
}


print(t)

}





## plot example

t.here=10
#good.draws=good.draws[t.here]
j=1

windows(14,5)
op <- par(mfrow = c(1,2), mgp = c(2,.8,0), mar = .1+c(3,3,3,1),cex.main=1.3)

plot(0:(lags-1),IRF_med[j,,t.here],type="l",  xlab="", ylab="", lwd=2, col=2,
     ylim=c(min(c(IRF_lb[j,,t.here])), max(c(IRF_ub[j,,t.here]))),  
     main =paste("bond spread -> IP growth, t=",t.here, sep=""))
lines(0:(lags-1),IRF_ub[j,,t.here],lty="dashed", lwd=2, col=3)
lines(0:(lags-1),IRF_lb[j,,t.here],lty="dashed", lwd=2, col=3)
abline(h=0, col=1)






