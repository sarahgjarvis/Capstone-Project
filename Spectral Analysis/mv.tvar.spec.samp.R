mv.tvar.spec.samp = function(fit,lim,wo.wn){
    library(mvtnorm)
    
    p = dim(fit$m)[1]
    T = dim(fit$m)[2]
    n = floor(T/2)
    j = 1:n
    w = 2*pi*j/T
    i = complex(imaginary = 1)
  
  if(fit$m[1,(p+1)] == fit$m[1,floor(T/2)]){
      
      m = fit$m[,1]
      n = fit$n[1]
      C = fit$C[,,1]
      s = fit$s[1]
      n.samp = 500
      samp = rmvt(n.samp, df = n, delta = m, type="shifted", sigma = C)
      
      spec = matrix(NA,n.samp,length(w))
      for(k in 1:n.samp){
        for(j in 1:length(w)){
          spec[k,j] = s/(2*pi*Mod(sum(c(1,samp[k,])*c(1,rep(-1,p))*exp(-i*seq(0,p)*w[j])))^2)
        }
      }
      
      if(wo.wn){ spec = 2*pi*spec/s }
      
      int = matrix(NA,length(w),3)
      for(j in 1:length(w)){
        int[j,] = c(quantile(spec[,j],0.025),mean(spec[,j]),quantile(spec[,j],0.975))
      }
  }else{
    
    n.samp = 200
    spec = array(NA,c(T-p,n.samp,length(w)))
    int = array(NA,c(T-p,length(w),3))
    for(t in 1:(T-p)){
        m = fit$m[,t+p]
        n = fit$n[t+p]
        C = fit$C[,,t+p]
        s = fit$s[t+p]
        samp = rmvt(n.samp, df = n, delta = m, type="shifted", sigma = C)
        
        for(k in 1:n.samp){
          for(j in 1:length(w)){
            spec[t,k,j] = s/(2*pi*Mod(sum(c(1,samp[k,])*c(1,rep(-1,p))*exp(-i*seq(0,p)*w[j])))^2)
          }
          if(wo.wn){ spec[t,k,] = 2*pi*spec[t,k,]/s }
        }
        
        for(j in 1:length(w)){
          int[t,j,] = c(quantile(spec[t,,j],0.025),mean(spec[t,,j]),quantile(spec[t,,j],0.975))
        }
    }
    
  }

  return(list(w = w, spec.samps = spec, int=int))
}
