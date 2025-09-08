rm(list=ls())
library(vars)
library(Rcpp)
library(RcppArmadillo)
library(numDeriv)
library(matrixcalc)
library(dplyr)



setwd("C:/Users/Julia/Dropbox/DynamicVAR_privat/Rpackage/publishable_code")
#setwd("C:/Users/jsg490/Dropbox/DynamicVAR")

## call c++ function script
sourceCpp("tvvar_cpp_functions5.cpp")
sourceCpp("lyapunov_cond.cpp")
#sourceCpp("irf_gen.cpp")


source("tvvar_r_functions.r")



lag.order = 2 #VAR lag order
zero.mean = TRUE #include intercept (FALSE) or work with demeaned data (TRUE)
number.factors =4 #number of factors (for tv. coefficient matrix)





## read in data
all.data <- read.table("full_data_updated.csv", header=TRUE, dec=",", sep=";")
size.data <- length(all.data[,1])

### set start date such that sample starts in January 1970
start.date <- as.Date("12-01-69", format="%m-%d-%y")

levels.data <- as.matrix(all.data[c(which(as.Date(all.data[,1])==start.date):(size.data)),2:4])
recessions <- as.matrix(all.data[c(which(as.Date(all.data[,1])==start.date):(size.data)),5])

transformed.data1 <- as.matrix(cbind(
  diff(log(levels.data[,1]))*100, diff(log(levels.data[,2]))*100, levels.data[-1,3]))


if(zero.mean==TRUE){
  transformed.data <- scale(transformed.data1, center=TRUE, scale= FALSE)} else{
    transformed.data <- scale(transformed.data1, center=FALSE, scale= FALSE)
  }
colnames(transformed.data) = c("ip growth", "cpi inflation", "bond spread")


VAR.data =transformed.data
dim.VAR <- ncol(VAR.data)

## discard pre-sample values if not needed (for comparison among 2, 3 and 4 lags)
if(lag.order==2){
  VAR.data =VAR.data[-c(1:2),]} else if(
    lag.order==3){
    VAR.data =VAR.data[-1,]}





## define structure of PHi^f (loadings structure)
Phi.f.str <- array(0, c(dim.VAR,dim.VAR, number.factors))
# baa -> ip and infl: financial -> macro
Phi.f.str[1,3,1] = 1
Phi.f.str[2,3,1] = 1
# ip -> ip and baa -> baa: ip growth and spread persistence
Phi.f.str[1,1,2] = 1
Phi.f.str[3,3,2] = 1
# infl -> infl: inflation persistence
Phi.f.str[2,2,3] = 1
# ip -> infl and infl -> ip: macro-macro spillovers
Phi.f.str[1,2,4] = 1
Phi.f.str[2,1,4] = 1




Phi.f.array <- Phi.f(structure= Phi.f.str, lag.order=lag.order)

Phi.f.array.mat <- matrix(Phi.f.array, nrow=dim.VAR)
Phi.f.array.mat.structure <- matrix(Phi.f.array, nrow=dim.VAR)

number.nonzero.phis <- length(which(Phi.f.array.mat.structure !=0))
Phi.f.array.mat.structure[which(Phi.f.array.mat.structure !=0)]=1


### starting values for estimation (e.g. true parameters)
### ordering: A,B,phi, nonzero elements in Phi^f, Omega, Phi^c
## transform for numerical stability

p1 = c(0.1, 0.8, rep(0.95, number.factors))

start.H <- diag(0.1, dim.VAR)


## estimate Phi^c as coefficient matrix in static VAR
B.t <- matrix(0, nrow=dim.VAR, ncol=(dim.VAR*lag.order+1))
for(j in 1:dim.VAR){
  B.t[j,] = coefficients(VAR(VAR.data, p=lag.order))[[j]][,1]
}
A.t <- B.t[,-(dim.VAR*lag.order+1)]

## case of nonzero intercept
if(zero.mean==FALSE){
  A.t <- cbind(B.t[,(dim.VAR*lag.order+1)], A.t)
}


Phi.c.ini = A.t




if(number.factors>0){
  params1 = 
    c(log(p1[1]/(1-p1[1])),
      log(p1[2]/(1-p1[2])),
      log(p1[3:(3+number.factors-1)]/(1-p1[3:(3+number.factors-1)])), 
      Phi.f.array.mat.structure[which(Phi.f.array.mat.structure !=0)]*0.01, 
      vech(start.H),
      as.vector(Phi.c.ini)
      #c(0.1, 0.3, 2)
    )
} else if(number.factors==0){
  params1 = 
    c(log(p1[1]/(1-p1[1])),
      log(p1[2]/(1-p1[2])),
      Phi.f.array.mat.structure[which(Phi.f.array.mat.structure !=0)]*0.01, 
      vech(start.H),
      as.vector(Phi.c.ini)
      #c(0.1, 0.3, 2)
    )}








########## Optimization

opti <-
  optim(par=params1,
        fn=opti.fct, 
        purpose="optim", 
        VAR.data = VAR.data, 
        Phi.f.array=Phi.f.array,
        zero.mean= zero.mean,
        Smooth = FALSE,
        method='BFGS', 
        hessian = FALSE, 
        control=list(trace=6, maxit=5000 ))









########## Evaluation of estimates


p_out= opti$par
opti.eval <- opti.fct(params=p_out,VAR.data = VAR.data, Phi.f.array=Phi.f.array,zero.mean= zero.mean, Smooth = FALSE, purpose="eval")


### final estimates (transformed back)
pars=c(
  exp(p_out[1])/(exp(p_out[1])+1),
  exp(p_out[2])/(exp(p_out[2])+1),
  exp(p_out[3:(3+number.factors-1)])/(exp(p_out[3:(3+number.factors-1)])+1),
  p_out[(2+number.factors+1):length(p_out)])



###goodness of fit

T.fin <- nrow(VAR.data)
K <- length(pars)
AIC <- 2*K+ 2*opti.eval$average.L*(T.fin-lag.order)
AIC.c <- AIC + (2*K^2+2*K)/((T.fin-lag.order)*dim.VAR-K-1)
BIC <- log((T.fin-lag.order)*dim.VAR)*K + 2*opti.eval$average.L*(T.fin-lag.order)






## Standard errors for final parameter estimates (transformed)
my.hess <- optimHess(pars, hessian.fct.tr, VAR.data = VAR.data, Phi.f.array=Phi.f.array, zero.mean=zero.mean)
my.hess.inv <- solve(my.hess*T.fin)
cov.matrix <- my.hess.inv
se.vector <- sqrt(diag(cov.matrix))



## Covariance matrix of untransformed parameter estimates (for IRFs)
par_out_hess <- optimHess(p_out, hessian.fct.untr, VAR.data = VAR.data, Phi.f.array=Phi.f.array, zero.mean=zero.mean)
cov.matrix_untr <- solve(par_out_hess*T.fin)




## Display estimation results for final specification
result.list <- cbind(round(pars, 3), round(se.vector,3), (abs(pars/se.vector)>1.96))

### ordering: A,B,phi, nonzero elements in Phi^f, Omega, Phi^c
names.vector <- 
  c("A", "B", "phi1", "phi2", "phi3", "phi4",# "phi5","phi6","phi7",
    
    #"Phi_f_int1",  "Phi_f_int2", "Phi_f_int3",
    
    "Phi_f1_13 L1", "Phi_f1_23 L1", 
    "Phi_f1_13 L2", "Phi_f1_23 L2", 
    
    "Phi_f2_11 L1", "Phi_f2_33 L1",
    "Phi_f2_11 L2", "Phi_f2_33 L2",
    
    "Phi_f3_22 L1",
    "Phi_f3_22 L2",
    
    "Phi_f4_21 L1", "Phi_f4_12 L1",  
    "Phi_f4_21 L2", "Phi_f4_12 L2",  
    
    paste("omega_", 1:6, sep=""),
    paste("Phi_c_", 1:3, "1 L1", sep=""),  
    paste("Phi_c_", 1:3, "2 L1", sep=""),  
    paste("Phi_c_", 1:3, "3 L1", sep=""),
    paste("Phi_c_", 1:3, "1 L2", sep=""),  
    paste("Phi_c_", 1:3, "2 L2", sep=""),  
    paste("Phi_c_", 1:3, "3 L2", sep="")
  )

rownames(result.list)=names.vector












