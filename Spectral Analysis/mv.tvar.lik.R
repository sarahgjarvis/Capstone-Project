mv.tvar.lik = function(x,plims,dfmins,step,s0){
  
  dn = seq(dfmins[1],1,step)
  bn = seq(dfmins[2],1,step)
  pn = seq(plims[1],plims[2],1)
  
  x = as.matrix(x)
  T = nrow(x)
  R = ncol(x)
  
  llik = array(NA,c(length(pn),length(dn),length(bn)) )
  aic = array(NA,c(length(pn),length(dn),length(bn)) )
  bic = array(NA,c(length(pn),length(dn),length(bn)) )
  # se = rep(NA,T)
  
  maxll = -1e3000
  
  for(i in 1:length(pn)){
    p = pn[i]
    # m0,C0,s0,n0
    m0 = rep(0,p) 
    C0 = diag(1,p) 
    s0 = s0 
    n0 = 1 
    for(j in 1:length(dn)){
      d = dn[j]
      for(k in 1:length(bn)){
        b = bn[k]
        
        temp.llik = 0
        mt = m0
        Ct = C0
        st = s0
        nt = n0
        for(t in (plims[2]+1):T){
          F = matrix(x[seq((t-1),(t-p),-1),],p,R)
          A = matrix(Ct%*%F/d,p,R)
          Q = t(F)%*%A + st*diag(1,R)
          inv.Q = solve(Q)
          A = A%*%inv.Q
          e = x[t,] - t(F)%*%mt
          nt = b*nt 
          
          temp.llik = temp.llik + lgamma((nt+R)/2)-lgamma(nt/2)-R*log(pi*nt)/2-log(det(as.matrix(Q)))/2-(nt+R)*log(1+t(e)%*%inv.Q%*%e/nt)/2
          # se[t] = t(e)%*%e #squared error for alternate AIC calculation
          
          r = nt + t(e)%*%inv.Q%*%e
          mt = mt + A%*%e
          r = as.vector(r/(nt+1))
          st = st*r
          nt = nt + 1
          Ct = r*(Ct/d - A%*%Q%*%t(A))
          Ct = (Ct + t(Ct))/2
          
        }
        
        llik[i,j,k] = temp.llik
        aic[i,j,k] =   2*p - 2*temp.llik  #2*p + (T-plims[2])*log(sum(se))
        bic[i,j,k] = log(T-plims[2])*p - 2*temp.llik
        if(temp.llik>maxll){  
          maxll = temp.llik
          optim.p = p
          optim.d = d
          optim.b = b
        }
      }
    }
  }
  
  #image.plot(pn,bn,llik[,1,])
  return(list(llik=llik,aic=aic,bic=bic,p=optim.p,d=optim.d,b=optim.b))
  
}
