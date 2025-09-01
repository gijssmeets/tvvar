library(MASS)

irf_generator = function(Phi,psi,aT,PT,HT,ek,lags,fixed=FALSE,B=500,seed_f=1234){
  
  Phi=Phi[,-1,]
  r = length(aT)
  N = dim(Phi)[1]
  p = (dim(Phi)[2])/dim(Phi)[1]
  PHI = array(NA,c((N*p),(N*p),lags))
  PI_out = matrix(0,(N*p),(N*p))
  f = matrix(NA,r,lags)
  irf = matrix(0,N,lags)
  irf_out = matrix(0,N,lags)
  CT = t(chol(HT))
  
  ##################################
  ########### Fixed
  
  if(fixed==TRUE){
    f[,1] = aT
    
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
      
      f[,t] = f[,(t-1)]
      
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
    
    irf_out = irf
    PI_out = PI
    
    output = list(irf=irf_out, PI=PI_out)
    return(output)
  }
  
  ##################################
  ########### Expectation
  
  if(fixed==FALSE){
    
    set.seed(seed_f)
    
    FF1 = mvrnorm(n=B, rep(0,r), PT)
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
  
}



#inputs: theta, MLE, V, estimated asymptotic variance of MLE

IR.one.t <- function(T.index, theta, V, N, fixed, s.h, B, lags, ek){

  irf = array(0,c(N,lags,B))
  
  for(i in 1:B){
    
    s.h = s.h+1
    set.seed(s.h)
    par_val = mvrnorm(n = 1, mu = theta, Sigma = V)

    #### obtain Phi,psi,aT,PT,HT evaluated at parameter value par_val  
    opti.eval <- opti.fct(params=par_val,VAR.data = VAR.data, Phi.f.array=Phi.f.array, zero.mean=zero.mean, purpose="eval")
    
    phi.est <- par_val[3:(3+number.factors-1)]
    
    Phi.f.est <- Phi.f.array.mat.structure
    Phi.f.est[which(Phi.f.array.mat.structure !=0)]=par_val[(3+number.factors):(3+number.factors+number.nonzero.phis-1)]
    
    Phi.c.est <- opti.eval$Phi.c
    
    Phi = array(cbind(Phi.c.est, Phi.f.est), c(dim.VAR,dim.VAR*lag.order+1, (number.factors+1)))
    psi=exp(phi.est)/(1+exp(phi.est))
    #psi=phi.est
    
    aT <- opti.eval$filtered.state[,T.index]
    PT <- opti.eval$filtered.state.variance[,,T.index]
    HT <- opti.eval$array.filtered.H[,,T.index]
    
    irf.gen <- irf_generator(Phi,psi,aT,PT,HT,ek,lags=lags,fixed=fixed)
    
    irf[,,i] = irf.gen$irf
    
  }
  
  list(irf)
}
