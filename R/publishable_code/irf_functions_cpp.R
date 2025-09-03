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