uni.spec.dens.intervals = function(fit,color,add,ttl,wo.wn,lims, linty){
  library(ggplot2)
  
  #calulate 95% intervals for phis
  library(mvtnorm)
  m = as.vector(fit$m[,1])
  p = length(m)
  T = ncol(fit$m)
  n = fit$n[1]
  s = fit$s[1]
  C = as.matrix(fit$C[,,1])
  n.samp = 500
  samp = rmvt(n.samp, df = n, delta = m, type="shifted", sigma = C)
  
  n = floor(T/2)
  j = 1:n
  w = 2*pi*j/T
  i = complex(imaginary = 1)
  spec.wo.wn = matrix(NA,n.samp,length(w))
  for(k in 1:n.samp){
    for(j in 1:length(w)){
      spec.wo.wn[k,j] = 1/(Mod(sum(c(1,samp[k,])*c(1,rep(-1,p))*exp(-i*seq(0,p)*w[j])))^2)
    }
  }
  
  int = matrix(NA,length(w),3)
  for(i in 1:length(w)){
    int[i,] = log10(c(quantile(spec.wo.wn[,i],0.025),quantile(spec.wo.wn[,i],0.5),quantile(spec.wo.wn[,i],0.975)))
  }
  
  if(!wo.wn){ int = s*int/(2*pi) }
  
  #plot spectrum
  if(!add){
    plot(w,int[,2],ylim=c(min(c(lims[1],int[,1])),max(c(lims[2],int[,3]))),type="l",xlab="period (years)",ylab="log10 spectra",xaxt = "n")
    title(ttl,col.main=color)
  }else{
    #lines(w,int[,2],col=color)
    title(ttl,col.main=color)
  }
  polygon(c(w, rev(w)), c(int[,3], rev(int[,1])),
          col=alpha(color,0.2), border = NA)
  lines(w,int[,2],col=color,lwd=2, lty = linty)
}
