library(MASS)


## construct Q-matrix (selection matrix)
Q.matrix.fun <- function(r,N,p){
aid.vec <- matrix(1:r, ncol=1)
aid.mat <- matrix(aid.vec %x% diag(N*p+1))
Qmat <- matrix(0, nrow=r*(N*p+1)^2, ncol=r)
	for(i in 1:r){
		Qmat[(which(aid.mat == i)),i]=1
	}
Qmat 
}






Phi.f <- function(structure, lag.order){
  dim.VAR <- dim(structure)[1]
  number.factors <- dim(structure)[3]
  Phi.f <- array(0, c(dim.VAR,(dim.VAR*lag.order+1), number.factors))
  for(j in 1:number.factors){
    Phi.f[,,j] = cbind(rep(0,dim.VAR), matrix(rep(structure[,,j], lag.order), nrow=dim.VAR))
  }
  return(Phi.f)
}



###### data transformation for state space form
### inputs: VAR.data, lag.order, number.factors



build.statespace.form <- function(VAR.data, lag.order, number.factors, Phi.f, Phi.c) {
  dim.VAR <- ncol(VAR.data)
  
  ## regressors: (Np x T)-matrix
  T.fin <- nrow(VAR.data)
  Y <- t(VAR.data)
  
  ## calls c++ function
  Y.minus1 =create_Y_minus1(Y, lag.order, T.fin, dim.VAR)
  
  ## Adjust Y to remove initial lags
  Y <- Y[, -(1:lag.order)]
  
  ## Dependent variable for state space form
  ytilde <- Y - Phi.c %*% Y.minus1
  
  ## Build array of Z-matrices
  Q <- Q.matrix.fun(r = number.factors, N = dim.VAR, p = lag.order)
  
  Z <- create_Z(Y.minus1, Phi.f, Q, T.fin, dim.VAR, number.factors, lag.order) 

  ## Return the result as a list
  list(Z = Z, ytilde = ytilde)
}






################################################## 
############ Estimation function


opti.fct <- function(params, VAR.data, Phi.f.array, zero.mean, Smooth, purpose){
	
A <- exp(params[1])/(1+exp(params[1]))
B <- exp(params[2])/(1+exp(params[2]))


if(number.factors>0){
if(number.factors==1){
pphi <- matrix(exp(params[3:(3+number.factors-1)])/(1+exp(params[3:(3+number.factors-1)])), nrow=1, ncol=1)} else{
pphi <- diag(exp(params[3:(3+number.factors-1)])/(1+exp(params[3:(3+number.factors-1)])))}

Q <- diag(number.factors)-pphi %*% t(pphi)
starta <- numeric(number.factors)
startP <- diag(number.factors)

} else if (number.factors==0){
pphi <- matrix(0, nrow=1, ncol=1)
Q <- diag(1)-pphi %*% t(pphi)
starta <- numeric(1)
startP <- diag(1)
}

Phi.f.as.matrix <- matrix(Phi.f.array, nrow=dim.VAR)
count.params <- 2+number.factors+length(which(Phi.f.as.matrix != 0))


Phi.f <- matrix(0, nrow(Phi.f.as.matrix), ncol(Phi.f.as.matrix))
Phi.f[which(Phi.f.as.matrix != 0)]=params[(2+number.factors+1):count.params]

#length(which(Phi.f != 0))
#count.params <- 3+length(as.vector(Phi.f))

######## initializations for state, state variance and covariance matrix

Nstar <- dim.VAR*(dim.VAR+1)/2
m1 <- matrix(D.matrix(dim.VAR) %*% params[(count.params+1):(count.params+Nstar)], ncol=dim.VAR)
m1[upper.tri(m1)]=0
momega <- m1%*%t(m1)
firstH <- momega/(1-A-B)
covi=1

count.params = count.params+Nstar


if(zero.mean==TRUE){
  Phi.c <- cbind(rep(0,dim.VAR),matrix(params[(count.params+1):(count.params+dim.VAR^2*lag.order)], nrow=dim.VAR))
} else{
#Phi.c <- cbind(rep(0,dim.VAR),matrix(params[(count.params+1):(count.params+dim.VAR^2*lag.order)], nrow=dim.VAR))
Phi.c <- matrix(params[(count.params+1):(count.params+dim.VAR^2*lag.order+dim.VAR)], nrow=dim.VAR)
}




if(number.factors==0){
data.ss <- build.statespace.form(VAR.data=VAR.data, lag.order=lag.order, number.factors=1, Phi.f = Phi.f, Phi.c=Phi.c)} else{
data.ss <- build.statespace.form(VAR.data=VAR.data, lag.order=lag.order, number.factors=number.factors, Phi.f = Phi.f, Phi.c=Phi.c)}


ytilde <- data.ss$ytilde
Z <- data.ss$Z
T.here <- length(ytilde[1,])-1

### call C++ loop function
output <- my_loop_main(ytilde, Z, startP, covi, firstH, momega, pphi, A, B, Q)


last.obs <- length(output$pstate[1,])




if(number.factors==0){} else{
### Kalman smoother 
predicted.a <- matrix(output$pstate[,-last.obs],ncol=T.here)
predicted.P <- output$pstatevariance[,,-last.obs]
filtered.a <- output$fstate[,-last.obs]
filtered.P <- output$fstatevariance[,,-last.obs]


L.array <- output$L[,,-last.obs]
F.array <- output$F[,,-last.obs]
v.mat <- output$vmat[,-last.obs]
pred.err.var <- output$F[,,-last.obs]

### more things to collect
r.mat <- matrix(0, nrow=number.factors,ncol=T.here)
N.array <- array(0, dim=c(number.factors,number.factors,T.here))
alphahat.mat <- matrix(0, nrow=nrow(predicted.a), ncol=ncol(predicted.a))
V.array <- array(0, dim=c(number.factors,number.factors,T.here))
smooth.error <- matrix(0, nrow=nrow(ytilde), ncol=T.here)
}



if(Smooth == FALSE){} else{
ytilde.short = ytilde[,-last.obs]

if(number.factors>0){
  
for(t in T.here:1){
r.t <- r.mat[,t]
N.t <- N.array[,,t]
F.t.inv <- solve(F.array[,,t])
L.t <- L.array[,,t]
v.t <- v.mat[,t]
a.t <- predicted.a[,t]
P.t <- predicted.P[,,t]
Z.t = Z[,,t]

r.tminus1 <- t(Z.t) %*% F.t.inv %*% v.t + t(L.t) %*% r.t
N.tminus1 <- t(Z.t) %*% F.t.inv %*% Z.t + t(L.t) %*% N.t %*% L.t

alphahat.mat[,t] <- a.t + P.t %*% r.tminus1
V.array[,,t] <- P.t - P.t %*% N.tminus1 %*% P.t

r.mat[,(t-1)] = r.tminus1
N.array[,,(t-1)] = N.tminus1

smooth.error[,t] <- ytilde.short[,t]-Z.t %*% alphahat.mat[,t]

}
} else{}
}



if(purpose=="optim"){
L = -output$avlik
} else if(purpose=="eval"){

L = -output$avlik


if(number.factors==0){
result <- new.env()
result$average.L <- -output$avlik
result$full.L <- -output$sumlik
result$Phi.c <- Phi.c
result$array.filtered.H <- output$aH
result$momega =vech(momega)
result$Nstar =Nstar
as.list(result)} else{


result <- new.env()
result$average.L <- -output$avlik
result$full.L <- -output$sumlik
result$filtered.state <- filtered.a
result$filtered.state.variance <- filtered.P
result$predicted.state <- predicted.a
result$predicted.state.variance <- predicted.P
result$Phi.c <- Phi.c
result$array.filtered.H <- output$aH
result$momega =vech(momega)
result$Nstar =Nstar
result$Z=Z
result$v.mat= v.mat
result$pred.err.var=pred.err.var

if(Smooth == TRUE){
result$smoothed.state <- alphahat.mat
result$smoothed.state.variance <- V.array
result$smooth.error=smooth.error
}
as.list(result)
}
}
}








hessian.fct.tr <- function(params, VAR.data, Phi.f.array, zero.mean){
  
  A <- params[1]
  B <- params[2]
  
  if(number.factors==1){
    pphi <- matrix(params[3:(3+number.factors-1)], nrow=1, ncol=1)} else{
      pphi <- diag(params[3:(3+number.factors-1)])}

  
  Q <- diag(number.factors)-pphi %*% t(pphi)
  
  
  
  Phi.f.as.matrix <- matrix(Phi.f.array, nrow=dim.VAR)
  count.params <- 2+number.factors+length(which(Phi.f.as.matrix != 0))
  
  
  Phi.f <- matrix(0, nrow(Phi.f.as.matrix), ncol(Phi.f.as.matrix))
  Phi.f[which(Phi.f.as.matrix != 0)]=params[(2+number.factors+1):count.params]
  
  #length(which(Phi.f != 0))
  
  #count.params <- 3+length(as.vector(Phi.f))
  
  ######## initializations for state, state variance and covariance matrix
  starta <- numeric(number.factors)
  startP <- diag(number.factors)
  
  
  Nstar <- dim.VAR*(dim.VAR+1)/2
  m1 <- matrix(D.matrix(dim.VAR) %*% params[(count.params+1):(count.params+Nstar)], ncol=dim.VAR)
  m1[upper.tri(m1)]=0
  momega <- m1%*%t(m1)
  firstH <- momega/(1-A-B)
  covi=1
  
  count.params = count.params+Nstar
  
  if(zero.mean==TRUE){
    Phi.c <- cbind(rep(0,dim.VAR),matrix(params[(count.params+1):(count.params+dim.VAR^2*lag.order)], nrow=dim.VAR))
  } else{
    Phi.c <- matrix(params[(count.params+1):(count.params+dim.VAR^2*lag.order+dim.VAR)], nrow=dim.VAR)
  }
  
  if(number.factors==0){
    data.ss <- build.statespace.form(VAR.data=VAR.data, lag.order=lag.order, number.factors=1, Phi.f = Phi.f, Phi.c=Phi.c)} else{
    data.ss <- build.statespace.form(VAR.data=VAR.data, lag.order=lag.order, number.factors=number.factors, Phi.f = Phi.f, Phi.c=Phi.c)}
  
  
  
  ytilde <- data.ss$ytilde
  Z <- data.ss$Z
  T.here <- length(ytilde[1,])-1
  
  ### call C++ loop function
  output <- my_loop_main(ytilde, Z, startP, covi, firstH, momega, pphi, A, B, Q)
  
  
  L = -output$avlik
  
}









hessian.fct.untr <- function(params, VAR.data, Phi.f.array, zero.mean){
  
  A <- exp(params[1])/(1+exp(params[1]))
  B <- exp(params[2])/(1+exp(params[2]))
  
  
  
  
  if(number.factors>0){
    if(number.factors==1){
      pphi <- matrix(exp(params[3:(3+number.factors-1)])/(1+exp(params[3:(3+number.factors-1)])), nrow=1, ncol=1)} else{
        pphi <- diag(exp(params[3:(3+number.factors-1)])/(1+exp(params[3:(3+number.factors-1)])))}
    
    Q <- diag(number.factors)-pphi %*% t(pphi)
    starta <- numeric(number.factors)
    startP <- diag(number.factors)
    
  } else if (number.factors==0){
    pphi <- matrix(0, nrow=1, ncol=1)
    Q <- diag(1)-pphi %*% t(pphi)
    starta <- numeric(1)
    startP <- diag(1)
  }
  
  
  
  Phi.f.as.matrix <- matrix(Phi.f.array, nrow=dim.VAR)
  count.params <- 2+number.factors+length(which(Phi.f.as.matrix != 0))
  
  
  Phi.f <- matrix(0, nrow(Phi.f.as.matrix), ncol(Phi.f.as.matrix))
  Phi.f[which(Phi.f.as.matrix != 0)]=params[(2+number.factors+1):count.params]
  
  #length(which(Phi.f != 0))
  
  #count.params <- 3+length(as.vector(Phi.f))
  
  ######## initializations for state, state variance and covariance matrix
  starta <- numeric(number.factors)
  startP <- diag(number.factors)
  
  
  Nstar <- dim.VAR*(dim.VAR+1)/2
  m1 <- matrix(D.matrix(dim.VAR) %*% params[(count.params+1):(count.params+Nstar)], ncol=dim.VAR)
  m1[upper.tri(m1)]=0
  momega <- m1%*%t(m1)
  firstH <- momega/(1-A-B)
  covi=1
  
  count.params = count.params+Nstar
  
  if(zero.mean==TRUE){
    Phi.c <- cbind(rep(0,dim.VAR),matrix(params[(count.params+1):(count.params+dim.VAR^2*lag.order)], nrow=dim.VAR))
  } else{
    Phi.c <- matrix(params[(count.params+1):(count.params+dim.VAR^2*lag.order+dim.VAR)], nrow=dim.VAR)
  }
  
  if(number.factors==0){
    data.ss <- build.statespace.form(VAR.data=VAR.data, lag.order=lag.order, number.factors=1, Phi.f = Phi.f, Phi.c=Phi.c)} else{
    data.ss <- build.statespace.form(VAR.data=VAR.data, lag.order=lag.order, number.factors=number.factors, Phi.f = Phi.f, Phi.c=Phi.c)}
  
  
  
  ytilde <- data.ss$ytilde
  Z <- data.ss$Z
  T.here <- length(ytilde[1,])-1
  
  ### call C++ loop function
  output <- my_loop_main(ytilde, Z, startP, covi, firstH, momega, pphi, A, B, Q)
  
  
  L = -output$avlik
  
}

















####### impulse responses

irf.generator = function(Phi1,psi,aT,PT,HT,ek,lags,B=500,seed_f=1234){
  
  Phi=Phi1[,-1,]
  r = length(aT)
  N = dim(Phi)[1]
  p = (dim(Phi)[2])/dim(Phi)[1]
  PHI = array(0,c((N*p),(N*p),lags))
  PI_out = matrix(0,(N*p),(N*p))
  f = matrix(NA,r,lags)
  irf = matrix(0,N,lags)
  irf_out = matrix(0,N,lags)
  CT = t(chol(HT))

  
  ##################################
  ########### Expectation
  
    FF1 = t(CT%*% matrix(rnorm(B*r, 0, 1), nrow=r))
    FF = array(rnorm((r*(lags-1)*B)),c(r,(lags-1),B))
    

    for(j in 1:B){
      
      f = matrix(NA,r,lags)
      f[,1] = aT + FF1[j,]
      
      Phit = Phi[,,1]
      for(i in 1:r){
        Phit = Phit+Phi[,,(i+1)]*f[i,1] 
      }
      
      PHI[1:N,,1] = Phit
      if(p>1){
        PHI[(N+1):(N*p),1:(N*(p-1)),1] = diag(1,(N*(p-1)))
        PHI[(N+1):(N*p),(N*(p-1)+1):(N*p),1] = matrix(0,(N*(p-1)),N)
      }
      
      for(t in 2:lags){
        
        f[,t] = psi*f[,(t-1)] + FF[,(t-1),j]*(1-psi^2)
        
        Phit = Phi[,,1]
        for(i in 1:r){
          Phit = Phit+Phi[,,(i+1)]*f[i,t] 
        }
        
        PHI[1:N,,t] = Phit
        
        if(p>1){
          PHI[(N+1):(N*p),1:(N*(p-1)),t] = diag(1,(N*(p-1)))
          PHI[(N+1):(N*p),(N*(p-1)+1):(N*p),t] = matrix(0,(N*(p-1)),N)
        }
      }
      
      irf[,1] = CT%*%ek
      PI = PHI[,,1]
      
      for(l in 2:lags){
        irf[,l] = PI[1:N,1:N]%*%CT%*%ek
        PI = PHI[,,l]%*%PI
      }
      
      irf_out = irf_out + irf/B
      PI_out = PI_out + PI/B
    }
    
    output = list(irf=irf_out, PI=PI_out)
    return(output)
  
  
}







irf_generator = function(Phi1, psi, aT, PT, HT, ek, lags, B = 500, seed_f = 1234) {
  Phi = Phi1[, -1, ]
  r = length(aT)
  N = dim(Phi)[1]
  p = (dim(Phi)[2]) / dim(Phi)[1]
  PHI = array(NA, c((N * p), (N * p), lags))
  PI_out = matrix(0, (N * p), (N * p))
  irf_out = matrix(0, N, lags)
  CT = t(chol(HT))
  
  ##################################
  ########### Expectation
  
  set.seed(seed_f)
  FF1 = t(PT %*% matrix(rnorm(B * r, 0, 1), nrow = r))
  FF = array(rnorm((r * (lags - 1) * B)), c(r, (lags - 1), B))
  
  # Call C++ function to run the main loop
  result = irf_generator_loop(PHI, irf_out, PI_out, Phi, psi, aT, FF1, FF, CT, ek, lags, B, N, p)
  
  # Return the result
  return(list(irf = result$irf, PI = result$PI))
}

























lyaponov <- function(phi.est, Phi.c.est, Phi.f.est, seedy=15){

phi_r = phi.est
r = length(phi_r)

Phi.c.noint <- Phi.c.est[,-1]
Phi.f.noint <- array(Phi.f.est, dim=c(dim.VAR, (dim.VAR*lag.order+1), r))[,-1,]

Phi_c=Phi.c.noint
Phi_f=Phi.f.noint


seed=seedy
result <- check_lyapunov(Phi_c, Phi_f, phi_r)

return(result)
}









IR.one.t <- function(T.index, theta, V, N, fixed, s.h, B, lags, ek){

  irf = array(NA,c(N,lags,B))
  
  for(i in 1:B){
    
    
    #### re-draw static coefficients from the asymptotic distribution
    s.h = s.h+1
    set.seed(s.h)
    par_val = mvrnorm(n = 1, mu = theta, Sigma = V)

    #### obtain Phi,psi,aT,PT,HT evaluated at parameter value par_val  
    opti.eval <- opti.fct(params=par_val,VAR.data = VAR.data, Phi.f.array=Phi.f.array, zero.mean=zero.mean, Smooth=FALSE, purpose="eval")
    
    phi.est <- par_val[3:(3+number.factors-1)]
    
    Phi.f.est <- Phi.f.array.mat.structure
    Phi.f.est[which(Phi.f.array.mat.structure !=0)]=par_val[(3+number.factors):(3+number.factors+number.nonzero.phis-1)]
    
    Phi.c.est <- opti.eval$Phi.c
    
    Phi1 = array(cbind(Phi.c.est, Phi.f.est), c(dim.VAR,dim.VAR*lag.order+1, (number.factors+1)))
    psi=exp(phi.est)/(1+exp(phi.est))
    
    aT <- opti.eval$filtered.state[,T.index]
    PT <- opti.eval$filtered.state.variance[,,T.index]
    HT <- opti.eval$array.filtered.H[,,T.index]
    
    lyaponov.coeff <- lyaponov(psi, Phi.c.est, Phi.f.est)
    
    if(lyaponov.coeff<0){
      psi.mat <- (diag(psi))
    irf.gen <- irf.generator(Phi1,psi.mat,aT,PT,HT,ek,lags=lags)
    
    
    irf[,,i] = irf.gen$irf
    } else{}
    
    
  }
  
  
  
  irf_generator
  
  
  return(irf)
}



