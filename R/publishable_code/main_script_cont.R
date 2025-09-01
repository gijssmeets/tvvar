rm(list=ls())
library(vars)
library(Rcpp)
library(RcppArmadillo)
library(numDeriv)
library(matrixcalc)
library(dplyr)



setwd("C:/Users/Julia/Dropbox/DynamicVAR_privat/Rpackage/publishable_code")

load("C:/Users/Julia/Dropbox/DynamicVAR_privat/Rpackage/publishable_code/model_output.RData")


## call c++ function script
sourceCpp("tvvar_cpp_functions5.cpp")
sourceCpp("lyapunov_cond.cpp")
#sourceCpp("irf_gen.cpp")


source("tvvar_r_functions.r")
source("irf_functions.r")


########### Time-varying impulse response functions

T.max = 20
#T.max = 2


B=500
lags=12
N=dim.VAR

good.draws <-c(NA, T.max)

IRF_out_lb2=array(NA,c(N,lags, T.max))
IRF_out_ub2=array(NA,c(N,lags, T.max))
IRF_out_med2=array(NA,c(N,lags, T.max))



## loop through time points

#for(t in 1:T.max){


t=1
T.index=t
theta=p_out
V=cov.matrix_untr
s.h=1 
ek= c(0,0,1)


irf <- IR.one.t(T.index=t, theta=p_out, V=cov.matrix_untr, lags=lags, N=dim.VAR, B=B, s.h=1, ek= c(0,0,1), fixed=FALSE)



for(i in 1:N){
  for(j in 1:lags){
    IRF_out_lb[i,j,t] = quantile(irf[i,j,], probs=0.16,na.rm=TRUE)
    IRF_out_ub[i,j,t] = quantile(irf[i,j,], probs=0.84,na.rm=TRUE)
    IRF_out_med[i,j,t] = quantile(irf[i,j,], probs=0.5,na.rm=TRUE)
  }
}

print(t)

#}

IRF_out_lb2[,,1]
IRF_out_ub2[,,1]
IRF_out_med2[,,1]









