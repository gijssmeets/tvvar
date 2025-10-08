library(MASS)

opti.fct <- function(par_free, par_fixed, VAR.data, Phi.f.array, cfg, Smooth, purpose){

# Unpack
zero.mean = cfg$zero.mean
dim.VAR = cfg$dim.VAR
lag.order = cfg$lag.order
number.factors = cfg$number.factors
  
# Combine [vec] par_free (free to optimise params) and [list] par_fixed (conditionally maximized params) into [vec] params
params <- rebuild_params(par_free, par_fixed)

# Transform back to (0,1) space
A <- exp(params[1])/(1+exp(params[1]))
B <- exp(params[2])/(1+exp(params[2]))


if(number.factors>0){
if(number.factors==1){
    pphi <- matrix(tanh(params[3:(3+number.factors-1)]), nrow=1, ncol=1)} else{
    pphi <- diag(tanh(params[3:(3+number.factors-1)]))
    }

# Factor unit variance assumption
Q <- diag(number.factors)- pphi %*% t(pphi)

# Initialise state mean as 0 (assumed stationary)
starta <- numeric(number.factors)

# Initialise state variance as unconditional variance (unit variance)
startP <- diag(number.factors)

} else if (number.factors==0){
pphi <- matrix(0, nrow=1, ncol=1)
Q <- diag(1)-pphi %*% t(pphi)
starta <- numeric(1)
startP <- diag(1)
}

# Unpack Phi.f matrix
Phi.f.as.matrix <- matrix(Phi.f.array, nrow=dim.VAR)
count.params <- 2+number.factors+length(which(Phi.f.as.matrix != 0))
Phi.f <- matrix(0, nrow(Phi.f.as.matrix), ncol(Phi.f.as.matrix))
Phi.f[which(Phi.f.as.matrix != 0)]=params[(2+number.factors+1):count.params]

# Unpack Omega matrix (for BEKK equation)
Nstar <- dim.VAR*(dim.VAR+1)/2
m1 <- matrix(D.matrix(dim.VAR) %*% params[(count.params+1):(count.params+Nstar)], ncol=dim.VAR)
m1[upper.tri(m1)]=0
momega <- m1%*%t(m1)

# Initialize the covariance matrix H_t (using unconditional variance)
firstH <- momega/(1-A-B)

# Flag for updating covariance (1, 2 both using the BEKK model for H_t). Can change this!
covi=1

# Update count.params
count.params = count.params+Nstar

# Unpacks Phi.c
if(zero.mean==TRUE){
    Phi.c <- cbind(rep(0,dim.VAR),matrix(params[(count.params+1):(count.params+dim.VAR^2*lag.order)], nrow=dim.VAR))
  } else{
    Phi.c <- matrix(params[(count.params+1):(count.params+dim.VAR^2*lag.order+dim.VAR)], nrow=dim.VAR)
  }


# Build statespace form
if(number.factors==0){
data.ss <- build.statespace.form(VAR.data=VAR.data, lag.order=lag.order, number.factors=1, Phi.f = Phi.f, Phi.c=Phi.c)} else{
data.ss <- build.statespace.form(VAR.data=VAR.data, lag.order=lag.order, number.factors=number.factors, Phi.f = Phi.f, Phi.c=Phi.c)}

# Store y_tilde, Z
ytilde <- data.ss$ytilde
Z <- data.ss$Z
T.here <- length(ytilde[1,])-1

# Run the Kalman filter
output <- my_loop_main(ytilde, Z, startP, covi, firstH, momega, pphi, A, B, Q)
last.obs <- length(output$pstate[1,])

if(number.factors==0){} else{
# Obtain Kalman filter elements
predicted.a <- matrix(output$pstate[,-last.obs, drop=FALSE],ncol=T.here)
predicted.P <- output$pstatevariance[,,-last.obs, drop=FALSE]
filtered.a <- output$fstate[,-last.obs, drop=FALSE]
filtered.P <- output$fstatevariance[,,-last.obs, drop=FALSE]
L.array <- output$L[,,-last.obs, drop=FALSE]
F.array <- output$`F`[,,-last.obs, drop=FALSE]
v.mat <- output$vmat[,-last.obs,drop=FALSE]
pred.err.var <- output$`F`[,,-last.obs,drop=FALSE]

# Initializations for the Kalman smoother
r.mat <- matrix(0, nrow=number.factors,ncol=T.here)
N.array <- array(0, dim=c(number.factors,number.factors,T.here))
alphahat.mat <- matrix(0, nrow=nrow(predicted.a), ncol=ncol(predicted.a))
V.array <- array(0, dim=c(number.factors,number.factors,T.here))
smooth.error <- matrix(0, nrow=nrow(ytilde), ncol=T.here)
gamma.array <- array(0, dim=c(number.factors,number.factors,T.here-1))
}

# Run the Kalman smoother
if(Smooth == FALSE){} else{
ytilde.short = ytilde[,-last.obs, drop=FALSE]

if(number.factors>0){
  
for(t in T.here:1){
r.t <- r.mat[,t]
N.t <- N.array[,,t]

F.t.inv <- solve(F.array[,,t])
L.t <- L.array[,,t]
v.t <- v.mat[,t]
a.t <- predicted.a[,t]
P.t <- predicted.P[,,t]
Z.t = Z[,,t, drop=FALSE]
Z.t <- matrix(Z.t, nrow = dim(Z)[1], ncol = dim(Z)[2])  # force 2D

r.tminus1 <- t(Z.t) %*% F.t.inv %*% v.t + t(L.t) %*% r.t
N.tminus1 <- t(Z.t) %*% F.t.inv %*% Z.t + t(L.t) %*% N.t %*% L.t

alphahat.mat[,t] <- a.t + P.t %*% r.tminus1
V.array[,,t] <- P.t - P.t %*% N.tminus1 %*% P.t

# Extended KS recursions for autocovariance matrix
if (t > 1){
gamma.array[,, t-1] <- diag(number.factors) - P.t %*% N.tminus1}
if (t < T.here) {
  J_t <- L.t %*% P.t
  gamma.array[,, t] <- gamma.array[,, t] %*% J_t
}

r.mat[,(t-1)] = r.tminus1
N.array[,,(t-1)] = N.tminus1
smooth.error[,t] <- ytilde.short[,t]-Z.t %*% alphahat.mat[,t]

}
} else{}
}


# Store log likelihood
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
result$Y.minus1 <- data.ss$Y.minus1

if(Smooth == TRUE){
result$smoothed.state <- alphahat.mat
result$smoothed.state.variance <- V.array
result$smooth.error=smooth.error
result$ytilde <- ytilde
result$smoothed.state.autocov <- gamma.array

# Initialize V* matrices
V_0_minus1_star <- array(0, dim=c(number.factors,number.factors))
V_0_star <- array(0, dim=c(number.factors,number.factors))
V_minus1_star <- array(0, dim=c(number.factors,number.factors))


# Fill V* matrices
for (t in 2:T.here){
  V_0_minus1_star <- V_0_minus1_star + gamma.array[,,t-1] + result$smoothed.state[,t] %*% t(result$smoothed.state[,t-1])
  V_0_star <- V_0_star + V.array[,,t] + result$smoothed.state[,t] %*% t(result$smoothed.state[,t])
  V_minus1_star <- V_minus1_star + V.array[,,t-1] + result$smoothed.state[,t-1] %*% t(result$smoothed.state[,t-1])
}

# Store V* matrices
result$V_0_minus1_star <- V_0_minus1_star
result$V_0_star <- V_0_star
result$V_minus1_star <- V_minus1_star

# Solve 'naive' pphi OLS, not taking into account the dependency of Sigma_eta.
# This is a good initial guess for the optimization of pphi in the conditional maximization step.
pphi_full <- V_0_minus1_star %*% solve(V_minus1_star)
pphi_full <- pmin(pmax(pphi_full, -0.999), 0.999)
if (nrow(pphi_full) == 1) {
  result$pphi_ols <- (pphi_full[1, 1])
} else {
  result$pphi_ols <- (diag(pphi_full))
}


}

as.list(result)
}
}
}








hessian.fct.tr <- function(params, VAR.data, Phi.f.array, zero.mean, cfg){
  number.factors <- cfg$number.factors
  dim.VAR <- cfg$dim.VAR
  lag.order <- cfg$lag.order
  
  
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


hessian.fct.untr <- function(params, VAR.data, Phi.f.array, zero.mean, cfg){
  

  
  number.factors <- cfg$number.factors
  dim.VAR <- cfg$dim.VAR
  lag.order <- cfg$lag.order
  
  
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





