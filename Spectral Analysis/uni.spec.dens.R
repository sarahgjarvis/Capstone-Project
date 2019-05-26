uni.spec.dens = function(fit,lim,wo.wn,freq){
  
  if(freq){
  m = fit$m
  s = fit$s
  p = dim(m)[1]
  T = dim(m)[2]
  n = floor(T/2)
  j = 1:n
  w = 2*pi*j/T
  i = complex(imaginary = 1)
  
  spec = rep(NA,length(w))
  spec.wo.wn = rep(NA,length(w))
    for(j in 1:length(w)){
      spec[j] = s[50]/(2*pi*Mod(sum(c(1,m[,50])*c(1,rep(-1,p))*exp(-i*seq(0,p)*w[j])))^2)
      spec.wo.wn[j] = 1/(Mod(sum(c(1,m[,50])*c(1,rep(-1,p))*exp(-i*seq(0,p)*w[j])))^2)
    }
  
  if(wo.wn){ spec = spec.wo.wn }
  return(list(w=w,spec=spec))
  }else{
    m = fit$m
    s = fit$s
    p = dim(m)[1]
    T = dim(m)[2]
    lam=seq(2,lim,length.out=150)
    w = 2*pi/lam
    i = complex(imaginary = 1)
    
    spec = rep(NA,length(w))
    spec.wo.wn = rep(NA,length(w))
    for(j in 1:length(w)){
      spec[j] = s[50]/(2*pi*Mod(sum(c(1,m[,50])*c(1,rep(-1,p))*exp(-i*seq(0,p)*w[j])))^2)
      spec.wo.wn[j] = 1/(Mod(sum(c(1,m[,50])*c(1,rep(-1,p))*exp(-i*seq(0,p)*w[j])))^2)
    }
    
    if(wo.wn){ spec = spec.wo.wn }
    return(list(lam=lam,spec=spec))
  }
}
