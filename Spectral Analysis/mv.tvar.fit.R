mv.tvar.fit = function(x,p,df,m0,C0,s0,n0){
  T = nrow(x)
  R = ncol(x)
  d = df[1]
  b = df[2]
  
  m = matrix(NA,p,T)
  C = array(NA,c(p,p,T))
  s = rep(NA,T)
  n = rep(NA,T)
  one.step.resids = matrix(NA,R,T)
  
  mt = m0
  Ct = C0
  st = s0
  nt = n0
  
  # filtering
  for(t in (p+1):T){
    F = x[seq((t-1),(t-p),-1),]
    A = Ct%*%F/d
    Q = t(F)%*%A + st*diag(1,R)
    inv.Q = solve(Q)
    A = A%*%inv.Q
    e = x[t,] - t(F)%*%mt
    mt = mt + A%*%e
    r = b*nt + t(e)%*%inv.Q%*%e
    nt = b*nt + 1
    r = as.vector(r/nt)
    st = st*r
    Ct = r*(Ct/d - A%*%Q%*%t(A))
    Ct = (Ct + t(Ct))/2
    
    m[,t] = mt
    C[,,t] = Ct
    n[t] = nt
    s[t] = st
    one.step.resids[,t] = e 
  }
  
  m[,1:p] = m[,(p+1)]
  C[,,1:p] = C[,,(p+1)]
  n[1:p] = n[(p+1)]
  s[1:p] = s[(p+1)]
  
  # smoothing
  e = matrix(NA,R,T)
  for(t in seq((T-1),(p+1),-1)){
    m[,t] = (1-d)*m[,t] + d*m[,(t+1)]
    e[,t] = x[t,] - t(x[seq((t-1),(t-p),-1),])%*%m[,t]
    n[t] = (1-b)*n[t] + b*n[(t+1)]
    st = s[t]
    s[t] = 1/((1-b)/st + b/s[(t+1)])
    C[,,t] = s[t]*((1-d)*C[,,t]/st + d*d*C[,,(t+1)]/s[(t+1)])
  }
  
  m[,1:p] = m[,(p+1)]
  C[,,1:p] = C[,,(p+1)]
  n[1:p] = n[(p+1)]
  s[1:p] = s[(p+1)]
  
  return(list(m=m,C=C,n=n,s=s,smooth.resids=e,one.step.resids=one.step.resids))
}
