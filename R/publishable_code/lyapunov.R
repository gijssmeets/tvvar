rm(list=ls())
library(vars)
library(Rcpp)
library(RcppArmadillo)
library(numDeriv)
library(matrixcalc)


setwd("C:/Users/Julia/Dropbox/DynamicVAR")
#setwd("C:/Users/jsg490/Dropbox/DynamicVAR")

setting.number= 4
lag.order=2
zero.mean= TRUE



load(paste("Code/Results/empirics/empirical_setting_",setting.number, "_lags", lag.order,"_zeromean",zero.mean,"_withSE.RData", sep=""))

sourceCpp("C:/Users/Julia/Dropbox/DynamicVAR_privat/Rpackage/publishable_code/lyapunov_cond.cpp")


phi_r = phi.est
r = length(phi_r)

Phi.c.noint <- Phi.c.est[,-1]
Phi.f.noint <- array(Phi.f.est, dim=c(dim.VAR, (dim.VAR*lag.order+1), r))[,-1,]

Phi_c=Phi.c.noint
Phi_f=Phi.f.noint


seed=14
result <- check_lyapunov(Phi_c, Phi_f, phi_r)






D <- 250
sim.lyap <- numeric(D) 
m=100
seed=147


check.lyapunov <- function(Phi_c, Phi_f, phi.r, D=250, m=100, seed=147){

for(k in 1:D){

r = length(phi.r)



phi.mat <- diag(phi.r)
fac.cov <- diag(r)- diag(phi.r^2)
fac.errors <- matrix(rnorm(r*m),nrow=m, ncol=r)
  
set.seed(seed+k)
sim.factors<- matrix(0,nrow=m, ncol=r)

for(i in 2:m){
  sim.factors[i,] = phi.mat %*% sim.factors[(i-1),] +  sqrt(fac.cov) %*% fac.errors[i,]
}



Phi.t.array <- array(0, dim=c(dim.VAR*lag.order, dim.VAR*lag.order, m))

for(j in 1:m){
  factors.j <- as.matrix(sim.factors[j,], ncol=1)
  Phi.t.array[1:dim.VAR,,j] = Phi.c.noint + apply(Phi.f.noint, c(1,2), function(x) x %*% factors.j)
  Phi.t.array[(dim.VAR+1):(dim.VAR*lag.order),1:dim.VAR,j] = diag(dim.VAR)
}


prod.Phi.t <-Reduce("%*%", lapply(1:100, function(i) Phi.t.array[,,i]))

sim.lyap[k] <- log(norm(prod.Phi.t, type = c("2")))

}

exp.lyap <- mean(sim.lyap)
  
return(exp.lyap)
}


check.lyapunov(Phi.c.est, Phi.f.est, phi.est)


